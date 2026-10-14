from pathlib import Path
import os
W=Path('/Users/dudu/Documents/Codex/2026-10-09/ggn-target-ggn-target-opt-anaconda3');os.environ['NUMBA_CACHE_DIR']=str(W/'work/numba_cache');os.environ['MPLCONFIGDIR']=str(W/'work/mplconfig')
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
import h5py,anndata as ad,numpy as np,pandas as pd,json,datetime,hashlib
from scipy import sparse
from scanpy.tools._score_genes import _score_genes_bins
B=Path('/Users/dudu/全部数据库/GGN_target');R=B/'results/50_screen';O=R/'step5/HLCA';F=B/'raw/HLCA/HLCA_core_688185ad-11c2-4172-a53a-f4f1f4076860.h5ad';T=W/'work/hlca';now=lambda:datetime.datetime.now().astimezone().isoformat()
m=json.loads((O/'HLCA_methods_before_analysis.json').read_text());cfg=json.loads((O/'HLCA_annotation_scope.json').read_text());fr=json.loads((R/'step4/freeze_manifest.json').read_text());frozen=pd.read_csv(fr['frozen_files'][0]['path'],sep='\t',low_memory=False);pool=frozen[frozen.eligible_topN_pool].copy();assert len(pool)==82
print(now(),'START expression analysis',flush=True)
with h5py.File(F) as h:
 obs=ad.io.read_elem(h['obs']);var=ad.io.read_elem(h['raw/var']);g=h['raw/X'];assert g.attrs['encoding-type']=='csr_matrix';shape=tuple(g.attrs['shape']);assert shape[0]==len(obs)
 symbols=pd.Index(var[cfg['gene_symbol_column']].astype(str));names=pd.Index(var.index.astype(str));assert names.is_unique
 keep=obs.index.isin(pd.read_csv(O/'included_cell_ids.tsv.gz',sep='\t').cell_id);assert keep.sum()==cfg['analyzed_cells'];z=obs.loc[keep].copy();z['cell_type_analysis']=z[cfg['cell_type_column']].astype(str);z['donor_analysis']=z[cfg['donor_column']].astype(str)
 symbol_counts=symbols.value_counts();lookup={s:i for i,s in enumerate(symbols) if symbol_counts[s]==1};ids={str(s).split('.')[0]:i for i,s in enumerate(var.index)}
 maps=[]
 for r in pool.itertuples():
  j=lookup.get(r.gene);basis='exact_symbol'
  if j is None:
   possible=[ids[e] for e in str(r.ensembl_ids).split(';') if e in ids]
   if len(set(possible))==1:j=possible[0];basis='frozen_Ensembl_ID'
  maps.append({'gene':r.gene,'lineage':r.lineage,'index':j,'HLCA_symbol':symbols[j] if j is not None else '', 'basis':basis if j is not None else 'unmapped'})
 mp=pd.DataFrame(maps);mp.to_csv(O/'pool82_identifier_mapping.tsv',sep='\t',index=False);assert mp['index'].notna().all(),mp[mp['index'].isna()]
 programs={};coverage=[]
 for program in ['KAC','CAS']:
  indices=[]
  for symbol in m[program+'_genes']:
   j=lookup.get(symbol);basis='exact_symbol'
   if j is None:
    matches=frozen[frozen.gene.eq(symbol)].ensembl_ids.dropna().unique();poss={ids[e] for value in matches for e in str(value).split(';') if e in ids}
    if len(poss)==1:j=next(iter(poss));basis='frozen_Ensembl_ID'
   coverage.append({'program':program,'gene':symbol,'present':j is not None,'Ensembl_ID':names[j] if j is not None else '', 'HLCA_symbol':symbols[j] if j is not None else '', 'basis':basis if j is not None else 'unmapped'})
   if j is not None:indices.append(j)
  programs[program]=indices;assert len(set(indices))==len(indices)
 pd.DataFrame(coverage).to_csv(O/'program_gene_coverage.tsv',sep='\t',index=False)

 ptr=g['indptr'][:];n=len(z);ng=shape[1];sums=np.zeros(ng,np.float64);expr=np.lib.format.open_memmap(T/'pool82_log1pCP10k.npy',mode='w+',dtype='float32',shape=(n,82));lib=np.zeros(n);pmean={p:np.zeros(n) for p in programs};progress=0;qc={'noninteger_values':0,'negative_values':0,'nonfinite_values':0,'zero_library_cells':0}
 def chunks():
  outrow=0
  for a in range(0,shape[0],4000):
   b=min(a+4000,shape[0]);sel=keep[a:b]
   if not sel.any():continue
   lo,hi=int(ptr[a]),int(ptr[b]);x=sparse.csr_matrix((g['data'][lo:hi],g['indices'][lo:hi],ptr[a:b+1]-lo),shape=(b-a,ng))[sel].astype(np.float32)
   yield outrow,x
   outrow+=int(sel.sum())
 for start,x in chunks():
  end=start+x.shape[0];qc['noninteger_values']+=int(np.count_nonzero(x.data!=np.rint(x.data)));qc['negative_values']+=int(np.count_nonzero(x.data<0));qc['nonfinite_values']+=int(np.count_nonzero(~np.isfinite(x.data)));counts=np.asarray(x.sum(axis=1,dtype=np.float64)).ravel();qc['zero_library_cells']+=int((counts==0).sum());assert (counts>0).all();lib[start:end]=counts;x.data*=np.repeat((10000/counts).astype('float32'),np.diff(x.indptr));np.log1p(x.data,out=x.data);sums+=np.asarray(x.sum(axis=0,dtype=np.float64)).ravel();expr[start:end]=x[:,mp['index'].astype(int)].toarray()
  for p,indices in programs.items():pmean[p][start:end]=np.asarray(x[:,indices].mean(axis=1,dtype=np.float64)).ravel()
  if end-progress>=50000:print(now(),'PASS1',end,'/',n,flush=True);progress=end
 expr.flush();assert not any(qc.values()),qc
 avg=pd.Series(sums/n,index=names);avg.to_csv(O/'score_control_pool_gene_means.tsv',sep='\t',header=['mean_log1pCP10k'],index_label='gene');ctrl={}
 for p in programs:
  np.random.seed(0);gs=pd.Index([names[j] for j in programs[p]],dtype='string');control=pd.Index([],dtype='string')
  for c in _score_genes_bins(gs,names,ctrl_as_ref=True,ctrl_size=50,n_bins=25,get_subset=lambda genes:avg.loc[genes].to_numpy()[None,:]):control=control.union(c)
  ctrl[p]=[ids[s] for s in control];assert len(control)>0;pd.DataFrame({'gene':control,'feature_name':[symbols[ids[c]] for c in control]}).to_csv(O/f'{p}_score_control_genes.tsv',sep='\t',index=False)
 scores={p:pmean[p].copy() for p in programs};progress=0
 for start,x in chunks():
  end=start+x.shape[0];x.data*=np.repeat((10000/lib[start:end]).astype('float32'),np.diff(x.indptr));np.log1p(x.data,out=x.data)
  for p in programs:scores[p][start:end]-=np.asarray(x[:,ctrl[p]].mean(axis=1,dtype=np.float64)).ravel()
  if end-progress>=100000:print(now(),'PASS2',end,'/',n,flush=True);progress=end
z['library_counts']=lib;thresholds={}
for p in programs:
 s=scores[p];assert np.isfinite(s).all();q1,q3=np.quantile(s,[.25,.75],method='linear');thr=q3+1.5*(q3-q1);z[p+'_score']=s;z[p+'_high']=s>thr;thresholds[p]={'q1':float(q1),'q3':float(q3),'threshold':float(thr),'genes_available':len(programs[p]),'genes_requested':len(m[p+'_genes']),'control_genes':len(ctrl[p]),'high_n':int((s>thr).sum()),'high_fraction':float((s>thr).mean())}
(O/'HLCA_thresholds_descriptive.json').write_text(json.dumps(thresholds,ensure_ascii=False,indent=2));cols=['cell_type_analysis','donor_analysis',cfg['tissue_column'],'library_counts','KAC_score','CAS_score','KAC_high','CAS_high'];z[cols].to_csv(O/'HLCA_cell_scores.tsv.gz',sep='\t',index_label='cell_id')
# Pure summaries; no inferential functions or p values.
srows=[]
for scope,cols in [('all',[]),('cell_type',['cell_type_analysis']),('donor',['donor_analysis']),('donor_cell_type',['donor_analysis','cell_type_analysis'])]:
 groups=[('all',z.index)] if not cols else list(z.groupby(cols,observed=True).groups.items())
 for key,ix in groups:
  key=key if isinstance(key,tuple) else (key,);sub=z.loc[ix]
  for p in programs:
   vals=sub[p+'_score'].to_numpy();quant=np.quantile(vals,[.25,.5,.75,.9,.95,.99]);rec={'scope':scope,'donor':'','cell_type':'','program':p,'n_cells':len(sub),'n_donors':sub.donor_analysis.nunique(),'small_stratum_under50':len(sub)<50,'mean_score':float(vals.mean()),'Q1':quant[0],'median':quant[1],'Q3':quant[2],'p90':quant[3],'p95':quant[4],'p99':quant[5],'threshold':thresholds[p]['threshold'],'high_cells':int(sub[p+'_high'].sum()),'high_fraction':float(sub[p+'_high'].mean())}
   if scope=='donor':rec['donor']=key[0]
   if scope=='cell_type':rec['cell_type']=key[0]
   if scope=='donor_cell_type':rec.update(donor=key[0],cell_type=key[1])
   srows.append(rec)
pd.DataFrame(srows).to_csv(O/'HLCA_scores_descriptive_summary.tsv',sep='\t',index=False)
erows=[];position=pd.Series(np.arange(n),index=z.index)
for scope,cols in [('cell_type',['cell_type_analysis']),('donor_cell_type',['donor_analysis','cell_type_analysis'])]:
 for key,ix in z.groupby(cols,observed=True).groups.items():
  key=key if isinstance(key,tuple) else (key,);idx=position.loc[ix].to_numpy();a=expr[idx];quant=np.quantile(a,[.5,.9],axis=0);means=a.mean(axis=0,dtype=np.float64);det=(a>0).mean(axis=0)
  for j,r in enumerate(pool.itertuples()):erows.append({'scope':scope,'cell_type':key[-1],'donor':key[0] if scope=='donor_cell_type' else '', 'gene':r.gene,'lineage':r.lineage,'n_cells':len(idx),'n_donors':z.loc[ix,'donor_analysis'].nunique(),'small_stratum_under50':len(idx)<50,'detection_fraction':det[j],'mean_log1pCP10k':means[j],'median_log1pCP10k':quant[0,j],'p90_log1pCP10k':quant[1,j]})
e=pd.DataFrame(erows);e[e.scope.eq('cell_type')].to_csv(O/'pool82_expression_by_cell_type.tsv',sep='\t',index=False);e[e.scope.eq('donor_cell_type')].to_csv(O/'pool82_expression_by_donor_cell_type.tsv.gz',sep='\t',index=False)
# Donor-equal distributions, not standard errors or significance tests.
d=e[e.scope.eq('donor_cell_type')];dd=d.groupby(['gene','lineage','cell_type'],observed=True).agg(donors_observed=('donor','nunique'),donor_median_detection=('detection_fraction','median'),donor_min_detection=('detection_fraction','min'),donor_max_detection=('detection_fraction','max'),donor_median_mean_log1pCP10k=('mean_log1pCP10k','median')).reset_index();dd.to_csv(O/'pool82_expression_donor_summary.tsv',sep='\t',index=False)
ss=pd.DataFrame(srows);ss=ss[ss.scope.eq('donor_cell_type')];ss.groupby(['cell_type','program'],observed=True).agg(donors_observed=('donor','nunique'),donor_median_high_fraction=('high_fraction','median'),donor_min_high_fraction=('high_fraction','min'),donor_max_high_fraction=('high_fraction','max'),donor_median_score=('median','median')).reset_index().to_csv(O/'HLCA_score_donor_summary.tsv',sep='\t',index=False)
qa={'completed_at':now(),'analyzed_cells':n,'donors':z.donor_analysis.nunique(),'cell_types':z.cell_type_analysis.nunique(),'raw_shape':list(map(int,shape)),'raw_count_QC':qc,'program_thresholds':thresholds,'expression_cells_by82_shape':[n,82],'inferential_tests_performed':False,'cross_dataset_comparisons_performed':False,'step6_started':False,'source_download_receipt':json.loads((B/'raw/HLCA/download_receipt.json').read_text())};(O/'HLCA_analysis_QA.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2));print(now(),'COMPLETE',json.dumps(qa,ensure_ascii=False),flush=True)
