#!/opt/anaconda3/envs/ggn_target/bin/python
"""Frozen direction validation only. No mutation of ranks, weights, gates or pools."""
from pathlib import Path
import os
W=Path('/Users/dudu/Documents/Codex/2026-10-09/ggn-target-ggn-target-opt-anaconda3')
for k,v in {'MPLCONFIGDIR':str(W/'work/mplconfig'),'NUMBA_CACHE_DIR':str(W/'work/numba_cache'),'XDG_CACHE_HOME':str(W/'work/cache'),'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1','MKL_NUM_THREADS':'1','VECLIB_MAXIMUM_THREADS':'1'}.items():os.environ[k]=v
import argparse,json,hashlib,datetime,re,gzip,math,collections,gc
import numpy as np,pandas as pd,scipy.sparse as sp
from scipy.io import mmread
from scipy.stats import mannwhitneyu,wilcoxon,spearmanr
from statsmodels.stats.multitest import multipletests
import h5py
B=Path('/Users/dudu/全部数据库/GGN_target');R=B/'results/50_screen';O=R/'step5'
def sha(p):
 with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def now():return datetime.datetime.now().astimezone().isoformat()
def log(x):
 t=now()+' '+str(x);print(t,flush=True)
 with (O/'execution.log').open('a') as f:f.write(t+'\n')
def save(x,name):x.to_csv(O/name,sep='\t',index=False,float_format='%.17g')
def bh(x,method='fdr_bh'):
 a=pd.Series(x).astype(float);b=pd.Series(np.nan,index=a.index);ok=a.notna()&np.isfinite(a)
 if ok.any():b.loc[ok]=multipletests(a[ok],method=method)[1]
 return b
M=json.loads((O/'validation_methods_before_analysis.json').read_text());F=M['freeze_manifest']
def guard():
 assert all(sha(v['path'])==v['sha256'] for v in F['frozen_files']),'FROZEN INPUT CHANGED'
 reg=json.loads((O/'spatial_registration.json').read_text())
 assert sha(reg['hypothesis_path'])==reg['hypothesis_sha256']
 assert sha(O/'validation_methods_before_analysis.json')==reg['methods_sha256']
 return pd.read_csv(F['frozen_files'][0]['path'],sep='\t',low_memory=False)
def symbol(g):return {'GBA':'GBA1'}.get(str(g),str(g))
def verdict(effect,qpositive,q2):
 return np.where(np.isfinite(effect)&(effect>0)&(qpositive<.05),'方向一致',np.where(np.isfinite(effect)&(effect<0)&(q2<.05),'方向相反','无法判定'))
def bulk(reuse_models=False):
 frozen=guard();genes=sorted(frozen.gene.unique());assert len(genes)==3405
 log('GSE282617 bulk start; freeze hashes verified; no paired inference')
 raw=pd.read_csv(B/'raw/bulk_GSE282617_AIS_MIA_IAC/GSE282617_processed_data.csv.gz')
 raw['symbol']=[symbol(g) if pd.notna(g) and str(g).strip() else str(i) for g,i in zip(raw.gene_name,raw.gene_id)]
 mapping=raw[['gene_id','gene_name','symbol']];save(mapping,'GSE282617_gene_symbol_mapping.tsv')
 # Preserve frozen symbols and resolve nomenclature through unambiguous stable IDs.
 stable=raw.gene_id.str.replace(r'\.\d+$','',regex=True);lookup_rows=[]
 for fr in frozen.drop_duplicates('gene').itertuples():
  ids=re.findall(r'ENSG\d+',str(fr.ensembl_ids));names=sorted(raw.loc[stable.isin(ids),'symbol'].unique())
  if len(names)==1:key=names[0];basis='frozen Ensembl ID matches raw gene_id'
  elif fr.gene in set(raw.symbol):key=fr.gene;basis='exact symbol; nonunique IDs not pooled'
  else:key='';basis='unresolved'
  lookup_rows.append(dict(gene=fr.gene,bulk_model_symbol=key,identifier_basis=basis))
 gene_lookup=pd.DataFrame(lookup_rows).set_index('gene');model_genes=gene_lookup.reindex(genes).bulk_model_symbol.tolist()
 assert all(model_genes) and len(set(model_genes))==len(model_genes)
 save(gene_lookup.reset_index(),'GSE282617_resolved_candidate_lookup.tsv')
 meta=pd.read_csv(R/'GSE282617_sample_mapping.tsv',sep='\t').set_index('matrix_sample')
 meta['included']=meta.mapping_status.eq('exact')&~meta.stage.eq('LUSC')&~meta.index.isin(['WJR18','WYA12'])
 meta['exclusion_reason']=np.where(meta.stage.eq('LUSC'),'FLA肺鳞癌',np.where(~meta.mapping_status.eq('exact'),'身份未解决，不猜测更名',''))
 save(meta.reset_index(),'GSE282617_sample_eligibility.tsv');meta=meta[meta.included].copy()
 assert meta.stage.value_counts().to_dict()=={'IAC':23,'AIS':19,'MIA':16,'Normal':5}
 counts=raw.groupby('symbol',sort=True)[['count.'+z for z in meta.index]].sum().T;counts.index=meta.index
 fpkm=raw.groupby('symbol',sort=True)[['FPKM.'+z for z in meta.index]].sum().T;fpkm.index=meta.index
 assert np.isfinite(counts.values).all() and np.all(counts.values>=0) and np.allclose(counts,np.rint(counts))
 assert np.isfinite(fpkm.values).all() and np.all(fpkm.values>=0)
 counts=counts.astype(np.int64);save(pd.DataFrame({'sample':counts.index,'library_counts':counts.sum(axis=1).values,'stage':meta.stage}),'GSE282617_libraries.tsv')
 from pydeseq2.dds import DeseqDataSet
 from pydeseq2.ds import DeseqStats
 allrows=[]
 for step,target,reference in [('I',['AIS','MIA'],['Normal']),('II',['IAC'],['AIS','MIA'])]:
  mm=meta[meta.stage.isin(target+reference)].copy();mm['group']=np.where(mm.stage.isin(target),'target','reference');mm.group=pd.Categorical(mm.group,categories=['reference','target']);c=counts.loc[mm.index]
  n1=int((mm.group=='target').sum());n0=int((mm.group=='reference').sum());assert min(n1,n0)>=5
  save(mm.reset_index(),'GSE282617_'+step+'_units.tsv')
  keep=((c>=10).sum(axis=0)>=math.ceil(len(c)/3))&(c.sum(axis=0)>=50)
  save(pd.DataFrame({'gene':c.columns,'retained_for_model':keep.values,'total':c.sum(axis=0).values}),'GSE282617_'+step+'_model_gene_eligibility.tsv')
  log(f'GSE282617 {step}:target{n1}/reference{n0};fit{int(keep.sum())} genes')
  if reuse_models:
   res=pd.read_csv(O/f'GSE282617_{step}_all_genes_Wald.tsv.gz',sep='\t',index_col=0)
   assert set(res.index)==set(c.columns[keep]);log('Reuse unchanged full-genome Wald fit; correct identifier extraction/BH only')
  else:
   dds=DeseqDataSet(counts=c.loc[:,keep],metadata=mm[['group']],design='~ group',refit_cooks=True,n_cpus=4,quiet=True)
   dds.deseq2();st=DeseqStats(dds,contrast=['group','target','reference'],alpha=.05,n_cpus=4,quiet=True);st.summary();res=st.results_df
   res.to_csv(O/f'GSE282617_{step}_all_genes_Wald.tsv.gz',sep='\t')
   normalized=pd.DataFrame(dds.layers['normed_counts'],index=dds.obs_names,columns=dds.var_names);normalized.to_csv(O/f'GSE282617_{step}_normalized_counts.tsv.gz',sep='\t')
   fitqc=pd.DataFrame(index=dds.var_names)
   for key in ['_LFC_converged','_beta_converged','_genewise_converged','_MAP_converged','non_zero','replaced','refitted']:
    if key in dds.var:fitqc[key]=dds.var[key]
    elif key in dds.varm and np.asarray(dds.varm[key]).size==len(dds.var_names):fitqc[key]=np.asarray(dds.varm[key]).ravel()
   fitqc.to_csv(O/f'GSE282617_{step}_fit_diagnostics.tsv',sep='\t')
  rr=res.reindex(model_genes).copy();rr.index=genes;rr.index.name='gene';rr=rr.reset_index();rr['bulk_model_symbol']=model_genes;rr['identifier_basis']=gene_lookup.reindex(genes).identifier_basis.values
  rr['p_greater']=np.where(rr.log2FoldChange>0,rr.pvalue/2,1-rr.pvalue/2);rr['q_greater_all_unique']=bh(rr.p_greater);rr['q_two_sided_all_unique']=bh(rr.pvalue)
  rr['ci95_low']=rr.log2FoldChange-1.96*rr.lfcSE;rr['ci95_high']=rr.log2FoldChange+1.96*rr.lfcSE
  ff=np.log2(fpkm.loc[mm.index].reindex(columns=model_genes)+1);a=ff.loc[mm.group=='target'].to_numpy();b=ff.loc[mm.group=='reference'].to_numpy()
  rr['FPKM_log2_median_target']=np.nanmedian(a,axis=0);rr['FPKM_log2_median_reference']=np.nanmedian(b,axis=0);rr['FPKM_log2_median_difference']=rr.FPKM_log2_median_target-rr.FPKM_log2_median_reference
  pg=np.full(len(genes),np.nan);p2=pg.copy();valid=np.isfinite(a).all(axis=0)&np.isfinite(b).all(axis=0)
  pg[valid]=mannwhitneyu(a[:,valid],b[:,valid],axis=0,alternative='greater',method='asymptotic',use_continuity=True).pvalue
  p2[valid]=mannwhitneyu(a[:,valid],b[:,valid],axis=0,alternative='two-sided',method='asymptotic',use_continuity=True).pvalue
  rr['MW_p_greater']=pg;rr['MW_p_two_sided']=p2;rr['MW_q_greater_all_unique']=bh(pg);rr['MW_q_two_sided_all_unique']=bh(p2)
  rr['gene_in_original_matrix']=rr.bulk_model_symbol.isin(counts.columns);rr['model_retained']=rr.bulk_model_symbol.isin(c.columns[keep]);rr['step']=step;rr['target_n']=n1;rr['reference_n']=n0;rr['patient_ID_available']=False
  save(rr,'GSE282617_'+step+'_unique_gene_results.tsv')
  mapped=frozen[['gene','lineage','rank_all_display','rank_regression_display','eligible_topN_pool','I_upstream_eligible','II_upstream_eligible']].merge(rr,on='gene',how='left',validate='many_to_one')
  for lineage,idx in mapped.groupby('lineage').groups.items():
   for src,out in [('p_greater','q_greater_lineage'),('pvalue','q_two_sided_lineage'),('MW_p_greater','MW_q_greater_lineage'),('MW_p_two_sided','MW_q_two_sided_lineage')]:mapped.loc[idx,out]=bh(mapped.loc[idx,src]).values
  mapped['direction_verdict']=verdict(mapped.log2FoldChange,mapped.q_greater_lineage,mapped.q_two_sided_lineage)
  mapped['MW_direction_verdict']=verdict(mapped.FPKM_log2_median_difference,mapped.MW_q_greater_lineage,mapped.MW_q_two_sided_lineage)
  mapped['raw_effect_direction']=np.where(mapped.log2FoldChange>0,'正向',np.where(mapped.log2FoldChange<0,'反向','无有效方向'))
  mapped['method_direction_conflict']=(mapped.log2FoldChange*mapped.FPKM_log2_median_difference)<0
  mapped['interpretation']='bulk混合组织方向；不是该谱系细胞内复现；仅并列报告，不计D3或改变排序'
  mapped['limitations']='无患者ID不能排除伪重复；第I步正常仅5个仅大效应可检出；非配对受患者间混杂；WJR18/WYA12身份未解决'
  save(mapped,'GSE282617_'+step+'_candidate_direction.tsv');allrows.append(mapped)
  log('GSE282617 '+step+' verdict '+str(mapped.direction_verdict.value_counts().to_dict()))
  if not reuse_models:del dds,st
  gc.collect()
 save(pd.concat(allrows,ignore_index=True),'GSE282617_all_candidate_directions.tsv');guard();log('bulk complete; immutable hashes unchanged')

def read_matrix(matrix,features,barcodes):
 gf=pd.read_csv(features,sep='\t',header=None);gg=np.array([symbol(x) for x in gf[1]])
 bf=pd.read_csv(barcodes,sep='\t',header=None)[0].astype(str).to_numpy()
 with gzip.open(matrix,'rb') as f:X=mmread(f).T.tocsr()
 if gf.shape[1]>=3:
  ok=gf[2].eq('Gene Expression').to_numpy();X=X[:,ok];gg=gg[ok]
 return collapse(X,gg),bf

def collapse(X,genes):
 assert X.shape[1]==len(genes)
 assert np.isfinite(X.data).all() and (X.data>=0).all() and np.allclose(X.data,np.rint(X.data))
 unique,inv=np.unique(genes,return_inverse=True)
 if len(unique)<len(genes):X=X@sp.csr_matrix((np.ones(len(inv),dtype=np.int64),(np.arange(len(inv)),inv)),shape=(len(inv),len(unique)))
 else:X=X[:,np.argsort(genes)]
 return X.tocsr(),unique

def read_h5(p):
 with h5py.File(p,'r') as f:
  q=f['matrix'];shape=tuple(q['shape'][:]);X=sp.csc_matrix((q['data'][:],q['indices'][:],q['indptr'][:]),shape=shape).T.tocsr();gg=np.array([x.decode() for x in q['features/name'][:]])
  if 'feature_type' in q['features']:
   ok=np.array([x.decode()=='Gene Expression' for x in q['features/feature_type'][:]]);X=X[:,ok];gg=gg[ok]
  bc=np.array([x.decode() for x in q['barcodes'][:]])
 return collapse(X,gg),bc

def spatial_inventory():
 metadata=pd.read_csv(O/'GEO_sample_metadata.tsv',sep='\t');rows=[]
 for p in sorted((B/'raw/ST_GSE307534/extracted').glob('*/*/filtered_feature_bc_matrix/matrix.mtx.gz')):
  gsm=re.search(r'GSM\d+',str(p))[0];name=p.parent.parent.name;r=metadata[(metadata.dataset=='GSE307534')&(metadata.gsm==gsm)].iloc[0].to_dict();r.update(matrix=str(p),features=str(p.parent/'features.tsv.gz'),barcodes=str(p.parent/'barcodes.tsv.gz'),positions=str(p.parent.parent/'spatial/tissue_positions_list.csv'),sample_name=name)
  if not Path(r['positions']).exists():r['positions']=str(p.parent.parent/'spatial/tissue_positions.csv')
  assert re.search(r'P\d+',name)[0]==r['patient'];rows.append(r)
 for p in sorted((B/'raw/ST_GSE189487').glob('*_matrix.mtx.gz')):
  gsm=re.search(r'GSM\d+',p.name)[0];prefix=str(p).removesuffix('_matrix.mtx.gz');r=metadata[(metadata.dataset=='GSE189487')&(metadata.gsm==gsm)].iloc[0].to_dict();r.update(matrix=str(p),features=prefix+'_features.tsv.gz',barcodes=prefix+'_barcodes.tsv.gz',positions=prefix+'_tissue_positions_list.csv.gz',sample_name=r['patient']);rows.append(r)
 t=pd.DataFrame(rows);assert len(t)==62 and not t.gsm.duplicated().any();return t

def positions(p):
 t=pd.read_csv(p,header=None)
 if str(t.iloc[0,0])=='barcode':t=t.iloc[1:].copy()
 assert t.shape[1]==6
 t.columns=['barcode','in_tissue','array_row','array_col','pixel_row','pixel_col'];t.barcode=t.barcode.astype(str)
 for c in t.columns[1:]:t[c]=pd.to_numeric(t[c],errors='coerce')
 assert not t.barcode.duplicated().any();return t.set_index('barcode')

def correlate(x,y):
 if len(x)<20 or not np.isfinite(y).all() or np.std(x)==0 or np.std(y)==0:return np.nan
 return float(spearmanr(x,y).statistic)

def spatial():
 frozen=guard();inv=spatial_inventory();save(inv,'spatial_sample_inventory.tsv');targets=pd.read_csv(O/'spatial_targets_frozen_top20.tsv',sep='\t');markers=M['spatial']['markers'];rows=[];qc=[];score_qc=[]
 (O/'spatial_arrays').mkdir(exist_ok=True)
 log('spatial analysis start AFTER registered hypothesis SHA; 62 samples; no spot-level P')
 for r in inv.itertuples(index=False):
  (X,gg),bc=read_matrix(r.matrix,r.features,r.barcodes);xy=positions(r.positions).reindex(bc);lib=np.asarray(X.sum(axis=1)).ravel();valid=xy.in_tissue.eq(1).to_numpy()&(lib>0)&xy[['pixel_row','pixel_col']].notna().all(axis=1).to_numpy();raw_n=len(bc);X=X[valid];xy=xy.iloc[np.where(valid)[0]];lib=lib[valid];bc=bc[valid];n=len(bc)
  norm=X.astype(np.float64).multiply((1e4/lib)[:,None]).tocsr();norm.data=np.log1p(norm.data);lookup={g:i for i,g in enumerate(gg)};score={};used={}
  for name,genes in markers.items():
   ii=[lookup[g] for g in genes if g in lookup and X[:,lookup[g]].sum()>0];used[name]=[gg[i] for i in ii];score[name]=np.asarray(norm[:,ii].mean(axis=1)).ravel() if len(ii)>=3 else np.full(n,np.nan)
   score_qc.append(dict(dataset=r.dataset,gsm=r.gsm,program=name,n_markers=len(ii),markers=';'.join(used[name]),valid=len(ii)>=3))
  mat=np.column_stack([norm[:,lookup[g]].toarray().ravel() if g in lookup else np.full(n,np.nan) for g in targets.gene])
  for j,t in enumerate(targets.itertuples(index=False)):
   val=mat[:,j];own=score[t.lineage];good=n>=20 and np.isfinite(own).all() and np.isfinite(val).all();rho=correlate(val,own) if good else np.nan;kac=correlate(val,score['KAC']) if np.isfinite(val).all() else np.nan;cdc=correlate(val,score['cDC_descriptive']) if t.lineage=='Myeloid' and np.isfinite(val).all() else np.nan
   nh=nl=0;diff=np.nan
   if good:
    low,high=np.quantile(own,[.25,.75]);lo=own<=low;hi=own>=high
    if low<high:
     nh=int(hi.sum());nl=int(lo.sum())
     if min(nh,nl)>=5:diff=float(val[hi].mean()-val[lo].mean())
   rows.append(dict(dataset=r.dataset,gsm=r.gsm,patient=r.patient,stage=r.stage,gene=t.gene,lineage=t.lineage,n_spots=n,gene_present=t.gene in lookup,gene_detected_spots=int(np.count_nonzero(val>0)),detection_fraction=float(np.mean(val>0)) if np.isfinite(val).all() else np.nan,mean_log1pCP10k=float(np.mean(val)) if np.isfinite(val).all() else np.nan,localization_spearman=rho,localization_high_minus_low=diff,high_spots=nh,low_spots=nl,KAC_colocalization_spearman=kac,cDC_descriptive_spearman=cdc,inference_unit='patient/sample; no spot P'))
  np.savez_compressed(O/'spatial_arrays'/f'{r.dataset}_{r.gsm}.npz',barcodes=bc.astype(str),genes=targets.gene.to_numpy(str),expression=mat.astype(np.float32),pixel_row=xy.pixel_row.to_numpy(),pixel_col=xy.pixel_col.to_numpy(),library_counts=lib,**{k:np.asarray(v,dtype=np.float32) for k,v in score.items()})
  qc.append(dict(dataset=r.dataset,gsm=r.gsm,patient=r.patient,stage=r.stage,raw_filtered_spots=raw_n,tissue_positive_library_coordinate_spots=n,removed=raw_n-n,total_counts=int(lib.sum()),detected_genes=int(np.asarray(X.sum(axis=0)>0).sum()),sample_usable=n>=20))
  log(f'{r.dataset} {r.gsm} {r.patient} {r.stage}:spots{raw_n}->{n}')
  save(pd.DataFrame(rows),'spatial_sample_effects.tsv');save(pd.DataFrame(qc),'spatial_sample_QC.tsv');save(pd.DataFrame(score_qc),'spatial_marker_coverage.tsv')
  del X,norm,mat;gc.collect()
 summarize_spatial();guard();log('spatial complete; freeze hashes unchanged')

def bootmedian(x,seed):
 if len(x)<2:return (np.nan,np.nan)
 rng=np.random.default_rng(seed);v=np.median(rng.choice(x,(2000,len(x)),replace=True),axis=1);return tuple(np.quantile(v,[.025,.975]))

def summarize_spatial():
 t=pd.read_csv(O/'spatial_sample_effects.tsv',sep='\t');metrics=['localization_spearman','localization_high_minus_low','KAC_colocalization_spearman'];out=[];patient_rows=[]
 for ds in ['GSE307534','GSE189487']:
  d=t[t.dataset.eq(ds)&~t.stage.eq('Normal')].copy()
  for (gene,lineage),z in d.groupby(['gene','lineage']):
   for metric in metrics:
    a=z.groupby('patient')[metric].mean();a=a[np.isfinite(a)];v=a.values
    for patient,eff in a.items():patient_rows.append(dict(dataset=ds,gene=gene,lineage=lineage,metric=metric,patient=patient,effect=eff))
    p=pt=np.nan
    if ds=='GSE307534' and len(v)>=8:
     p=1. if np.all(v==0) else float(wilcoxon(v,alternative='greater',zero_method='wilcox',method='auto').pvalue)
     pt=1. if np.all(v==0) else float(wilcoxon(v,alternative='two-sided',zero_method='wilcox',method='auto').pvalue)
    ci=bootmedian(v,2026100905) if len(v)>=8 else (np.nan,np.nan)
    out.append(dict(dataset=ds,gene=gene,lineage=lineage,metric=metric,n_patients=len(v),n_samples=z[metric].notna().sum(),median_patient_effect=float(np.median(v)) if len(v) else np.nan,positive_patients=int((v>0).sum()),negative_patients=int((v<0).sum()),ci95_low=ci[0],ci95_high=ci[1],p_greater=p,p_two_sided=pt,interpretation='共享发现患者，非独立复现' if ds=='GSE307534' else '独立队列仅6样本，只报方向、不作显著性推断'))
 out=pd.DataFrame(out)
 for (ds,metric),idx in out.groupby(['dataset','metric']).groups.items():
  out.loc[idx,'p_holm_greater']=bh(out.loc[idx,'p_greater'],'holm').values;out.loc[idx,'p_holm_two_sided']=bh(out.loc[idx,'p_two_sided'],'holm').values
 out['direction']=np.where(out.median_patient_effect>0,'正向',np.where(out.median_patient_effect<0,'反向','无法估计'))
 out['inference_verdict']=np.where(out.dataset=='GSE189487','仅描述，不作判定',verdict(out.median_patient_effect,out.p_holm_greater,out.p_holm_two_sided))
 save(out,'spatial_patient_summary.tsv');save(pd.DataFrame(patient_rows),'spatial_patient_effects.tsv')
 stage=t.groupby(['dataset','stage','gene','lineage']).agg(n_samples=('gsm','nunique'),median_localization=('localization_spearman','median'),median_KAC_colocalization=('KAC_colocalization_spearman','median'),median_detection_fraction=('detection_fraction','median')).reset_index();save(stage,'spatial_stage_descriptive.tsv')

def descriptive():
 frozen=guard();genes=sorted(frozen.gene.unique());meta=pd.read_csv(O/'GEO_sample_metadata.tsv',sep='\t');rows=[];qc=[]
 for ds,root,pat in [('GSE300827','ST_GSE300827_AIS_MIA','*.h5'),('GSE277206','ST_GSE277206','extracted/*.h5')]:
  for p in sorted((B/'raw'/root).glob(pat)):
   gsm=re.search(r'GSM\d+',p.name)[0];r=meta[(meta.dataset==ds)&(meta.gsm==gsm)].iloc[0];(X,gg),bc=read_h5(p);v=np.asarray(X.sum(axis=0)).ravel();total=float(v.sum());n=X.shape[0];nz=np.asarray((X>0).sum(axis=0)).ravel();idx={g:i for i,g in enumerate(gg)}
   for gene in genes:
    j=idx.get(gene);cnt=float(v[j]) if j is not None else np.nan
    rows.append(dict(dataset=ds,gsm=gsm,title=r.title,patient='unknown',stage=r.stage,gene=gene,total_sample_gene_counts=cnt,log2CPM=float(np.log2(cnt/total*1e6+1)) if j is not None and total>0 else np.nan,spot_detection_fraction=float(nz[j]/n) if j is not None else np.nan,has_coordinates=False,verdict='仅描述，不作判定',P_value=np.nan,reason='无Normal；AIS仅2未达每组5' if ds=='GSE300827' else '只有MIA两样本，无第二组'))
   qc.append(dict(dataset=ds,gsm=gsm,title=r.title,stage=r.stage,filtered_spots=n,total_counts=int(total),n_genes=len(gg),patient='unknown',coordinates_present=False));log(f'{ds} {gsm} descriptive sample pseudobulk only')
 t=pd.DataFrame(rows);save(t,'nonspatial_reference_sample_pseudobulk.tsv.gz');save(pd.DataFrame(qc),'nonspatial_reference_sample_QC.tsv')
 summary=t.groupby(['dataset','stage','gene']).agg(n_samples=('gsm','nunique'),median_log2CPM=('log2CPM','median'),min_log2CPM=('log2CPM','min'),max_log2CPM=('log2CPM','max'),median_spot_detection_fraction=('spot_detection_fraction','median')).reset_index();summary['interpretation']='仅描述，不作判定；无跨阶段显著性检验';save(summary,'nonspatial_reference_gene_summary.tsv');guard()

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('mode',choices=['bulk','bulk_remap','spatial','descriptive','spatial_summary']);a=ap.parse_args();guard();log('SCRIPT SHA256='+sha(__file__)+' mode='+a.mode)
 {'bulk':bulk,'bulk_remap':lambda:bulk(reuse_models=True),'spatial':spatial,'descriptive':descriptive,'spatial_summary':summarize_spatial}[a.mode]()
