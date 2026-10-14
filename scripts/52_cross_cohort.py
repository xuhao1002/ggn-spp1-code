#!/opt/anaconda3/envs/ggn_target/bin/python
"""Independent study step 2. Inputs only raw/ and results/50_screen/. No legacy result imports.
Commands: audit (cached GEO); preprocess DATASET; bulk; sc DATASET; summarize. No step 3 operations.
"""
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']: os.environ[k]='1'
from pathlib import Path
CACHE=Path('/Users/dudu/Documents/Codex/2026-10-09/ggn-target-ggn-target-opt-anaconda3/work')
os.environ['NUMBA_CACHE_DIR']=str(CACHE/'numba_cache');os.environ['MPLCONFIGDIR']=str(CACHE/'mpl');os.environ['XDG_CACHE_HOME']=str(CACHE/'cache')
import sys,re,gzip,json,time,hashlib,gc
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.io import mmread
from scipy.stats import wilcoxon,binom
from statsmodels.stats.multitest import multipletests
B=Path('/Users/dudu/全部数据库/GGN_target');O=B/'results/50_screen/step2'
VL={'AT2L':['AT2','Club/secretory'],'Myeloid':['Macrophage','Monocyte','cDC','pDC','Neutrophil','Proliferating myeloid'],'Macrophage':['Macrophage'],'Fibroblast':['Fibroblast'],'T_NK':['T cell','NK'],'Endothelial':['Endothelial']}
EPI=['AT2','AT1','Club/secretory','Ciliated','Basal']
def log(x):
 s=time.strftime('%Y-%m-%d %H:%M:%S')+' '+str(x); print(s,flush=True)
 with (O/'execution.log').open('a') as f:f.write(s+'\n')
def table(x,p):x.to_csv(p,sep='\t',index=False)

MARKERS={'AT2': ['SFTPC', 'SFTPB', 'SFTPA1', 'NAPSA', 'LAMP3', 'ABCA3'], 'AT1': ['AGER', 'RTKN2', 'CLIC5', 'EMP2'], 'Club/secretory': ['SCGB1A1', 'SCGB3A2', 'SCGB3A1', 'BPIFB1'], 'Ciliated': ['FOXJ1', 'CAPS', 'TPPP3', 'RSPH1'], 'Basal': ['KRT5', 'TP63', 'KRT15'], 'Endothelial': ['PECAM1', 'VWF', 'CLDN5', 'CDH5'], 'Lymphatic EC': ['PROX1', 'CCL21', 'MMRN1'], 'Fibroblast': ['COL1A1', 'COL1A2', 'DCN', 'LUM'], 'SMC/pericyte': ['ACTA2', 'MYH11', 'RGS5', 'PDGFRB'], 'Mesothelial': ['MSLN', 'UPK3B', 'WT1'], 'T cell': ['CD3E', 'CD3D', 'TRAC', 'CD2'], 'NK': ['NKG7', 'GNLY', 'KLRD1', 'KLRF1'], 'B cell': ['MS4A1', 'CD79A', 'BANK1'], 'Plasma': ['JCHAIN', 'MZB1', 'IGHG1'], 'Macrophage': ['CD68', 'C1QA', 'C1QB', 'MARCO', 'APOE'], 'Monocyte': ['CD14', 'VCAN', 'S100A8'], 'cDC': ['CLEC9A', 'CD1C', 'FCER1A', 'CLEC10A'], 'pDC': ['LILRA4', 'IL3RA', 'CLEC4C'], 'Mast': ['CPA3', 'KIT', 'MS4A2', 'HDC'], 'Neutrophil': ['CSF3R', 'FCGR3B', 'CXCR2'], 'Platelet': ['PPBP', 'PF4', 'ITGA2B', 'TUBB1']}

def preprocess(ds):
 import scanpy as sc, anndata as ad
 d=O/ds;d.mkdir(exist_ok=True)
 audit=pd.read_csv(O/'sample_audit.tsv',sep='\t').fillna('')
 s=audit[(audit.dataset==ds)&audit.include].copy()
 if ds=='GSE288479':s=s[s.gsm.str[3:].astype(int)<=8768575]
 parts=[];qc=[]
 for _,r in s.iterrows():
  if ds=='GSE164789':
   root=B/'raw/scRNA_GSE164789_precursor';pre=root/(r.gsm+'_'+r.title)
   gp=str(pre)+'.genes.tsv.gz';bp=str(pre)+'.barcodes.tsv.gz';mp=str(pre)+'.matrix.mtx.gz'
  elif ds=='GSE288479':
   root=B/'raw/multi_GSE288479_partsolid';pre=root/(r.gsm+'_'+r.patient+'_'+r.site+'_scRNA')
   gp=str(pre)+'_features.tsv.gz';bp=str(pre)+'_barcodes.tsv.gz';mp=str(pre)+'_matrix.mtx.gz'
  else:
   root=B/'raw/scRNA_GSE203360_mGGO';mp=root/(r.gsm+'_'+r.patient+'-'+r.site+'.expression_matrix.txt.gz')
  if ds!='GSE203360':
   gf=pd.read_csv(gp,sep='\t',header=None);genes=gf[1 if gf.shape[1]>1 else 0].astype(str).values
   bf=pd.read_csv(bp,sep='\t',header=None)
   barcodes=bf[0].astype(str).values if bf.shape[1]==1 else np.array(['matrix_column_'+str(i+1) for i in range(len(bf))])
   with gzip.open(mp,'rb') as f:X=mmread(f).T.tocsr().astype(np.float32)
  else:
   # Chunk dense text by gene rows to bound peak memory.
   blocks=[];gg=[];barcodes=None
   for ch in pd.read_csv(mp,sep='\t',index_col=0,chunksize=500):
    if barcodes is None:barcodes=ch.columns.to_numpy()
    gg.extend(ch.index.astype(str));blocks.append(sp.csr_matrix(ch.values,dtype=np.float32))
   genes=np.array(gg);X=sp.vstack(blocks).T.tocsr();del blocks
  if ds=='GSE288479':
   assert bf.shape[1]==16 and np.allclose(np.asarray(X.sum(axis=1)).ravel(),bf[1])
   assert np.array_equal(np.asarray((X>0).sum(axis=1)).ravel(),bf[2])
  assert X.shape==(len(barcodes),len(genes))
  assert np.all(X.data>=0) and np.all(X.data==np.floor(X.data)),'Counts must be nonnegative integers'
  # Duplicate symbols are summed, never renamed as fictitious genes.
  uniq,inv=np.unique(genes,return_inverse=True)
  X=X@sp.csr_matrix((np.ones(len(inv)),(np.arange(len(inv)),inv)),shape=(len(inv),len(uniq)))
  a=ad.AnnData(X.astype(np.float32),obs=pd.DataFrame(index=[r.gsm+':'+str(x) for x in barcodes]),var=pd.DataFrame(index=uniq))
  for c in ['patient','site','histology','gsm']:a.obs[c]=r[c]
  a.obs['sample']=r.gsm;a.var['mt']=a.var_names.str.startswith('MT-')
  sc.pp.calculate_qc_metrics(a,qc_vars=['mt'],percent_top=None,log1p=False,inplace=True)
  n0=a.n_obs;keep=(a.obs.n_genes_by_counts>=200)&(a.obs.total_counts>=500)&(a.obs.pct_counts_mt<20)
  a=a[keep].copy();n1=a.n_obs
  sc.pp.scrublet(a,random_state=0)
  a=a[~a.obs.predicted_doublet].copy()
  qc.append(dict(sample=r.gsm,patient=r.patient,site=r.site,raw_cells=n0,after_qc=n1,after_doublet=a.n_obs,mt_genes=int(a.var.mt.sum())))
  parts.append(a);log(f'{ds} {r.gsm} cells {n0}->{n1}->{a.n_obs}')
  table(pd.DataFrame(qc),d/'sample_qc.tsv')
 a=ad.concat(parts,join='outer',fill_value=0);del parts;gc.collect()
 a.X=sp.csr_matrix(a.X);sc.pp.filter_genes(a,min_cells=10);a.layers['counts']=a.X.copy()
 sc.pp.highly_variable_genes(a,n_top_genes=3000,flavor='seurat_v3',batch_key='sample',span=.3)
 sc.pp.normalize_total(a,target_sum=1e4);sc.pp.log1p(a)
 hv=a[:,a.var.highly_variable].copy();sc.pp.scale(hv,max_value=10);sc.tl.pca(hv,n_comps=50,random_state=0)
 a.obsm['X_pca']=hv.obsm['X_pca'];del hv;gc.collect()
 import harmonypy
 ho=harmonypy.run_harmony(a.obsm['X_pca'],a.obs,'sample',random_state=0,max_iter_harmony=20)
 z=ho.Z_corr;z=z.detach().cpu().numpy() if hasattr(z,'detach') else np.asarray(z)
 a.obsm['X_pca_harmony']=z if z.shape[0]==a.n_obs else z.T
 sc.pp.neighbors(a,use_rep='X_pca_harmony',n_neighbors=30,random_state=0)
 sc.tl.leiden(a,resolution=1.,key_added='cluster',flavor='igraph',n_iterations=2,directed=False,random_state=0)
 for label,genes in MARKERS.items():sc.tl.score_genes(a,[g for g in genes if g in a.var_names],score_name=label,random_state=0)
 cs=a.obs.groupby('cluster',observed=True)[list(MARKERS)].mean();cs.to_csv(d/'cluster_marker_scores.tsv',sep='\t')
 ann=pd.DataFrame({'label':cs.idxmax(axis=1),'n_cells':a.obs.groupby('cluster',observed=True).size()})
 ann['second_label']=cs.apply(lambda r:r.nlargest(2).index[1],axis=1)
 ann['score_margin']=np.sort(cs.values,axis=1)[:,-1]-np.sort(cs.values,axis=1)[:,-2]
 check=list(dict.fromkeys(sum(MARKERS.values(),[])+['EPCAM','PTPRC','MKI67','TOP2A','KRT8','KRT19']))
 check=[g for g in check if g in a.var_names]
 cm=pd.DataFrame(a[:,check].X.toarray(),index=a.obs_names,columns=check).groupby(a.obs.cluster,observed=True).mean()
 cm.to_csv(d/'cluster_marker_expression.tsv',sep='\t');ann.to_csv(d/'cluster_annotation.tsv',sep='\t')
 a.obs['label']=a.obs.cluster.map(ann.label).astype(str)
 a.obs.to_csv(d/'cell_annotations.tsv.gz',sep='\t')
 a.write_h5ad(d/'rebuilt.h5ad',compression='lzf');log(f'{ds} PREPROCESS COMPLETE {a.shape}; annotation requires review before inference')


def bh(p):
 p=np.asarray(p,dtype=float);q=np.full(p.shape,np.nan);ok=np.isfinite(p)
 if ok.any():q[ok]=multipletests(p[ok],method='fdr_bh')[1]
 return q

def median_ci(x):
 x=np.sort(x);n=len(x)
 if not n:return np.nan,np.nan
 k=[k for k in range(1,n//2+1) if 1-2*binom.cdf(k-1,n,.5)>=.95]
 return (x[max(k)-1],x[n-max(k)]) if k else (-np.inf,np.inf)

def paired_summary(delta):
 rows=[]
 for g,rr in delta.iterrows():
  x=rr.dropna().values;lo,hi=median_ci(x)
  rows.append(dict(gene=g,n_pairs=len(x),median_patient_effect=np.median(x) if len(x) else np.nan,mean_patient_effect=np.mean(x) if len(x) else np.nan,positive_patients=int((x>0).sum()),negative_patients=int((x<0).sum()),median_ci95_low=lo,median_ci95_high=hi))
 return pd.DataFrame(rows).set_index('gene')

def finish(r,ds,axis,lineage,n,minpairs,weight):
 r=r.copy();r['dataset']=ds;r['axis']=axis;r['lineage']=lineage;r['n_pairs']=n;r['weight']=weight
 r['q_greater']=bh(r['p_greater']);r['q_two_sided']=bh(r['p_two_sided'])
 r['status']=np.where(r.p_greater.isna(),'无检测值',np.where((r.effect>0)&(r.q_greater<.05),'支持',np.where((r.effect<0)&(r.q_two_sided<.05),'方向相反','未通过')))
 if n<minpairs:r['status']='无法验证：配对数不足';r[['p_greater','p_two_sided','q_greater','q_two_sided']]=np.nan
 r['support']=r.status.eq('支持');r.index.name='gene'
 name=f'{ds}_{axis}_{lineage}';r.to_csv(O/(name+'.tsv'),sep='\t');return r

def wilcox_test(a,b,genes,ds,axis):
 assert list(a.columns)==list(b.columns)
 delta=a.reindex(genes)-b.reindex(genes);r=paired_summary(delta)
 r['effect']=r.median_patient_effect;r['p_greater']=np.nan;r['p_two_sided']=np.nan
 n=len(a.columns)
 if n>=7:
  for g,x in delta.iterrows():
   if x.isna().any():continue
   if (x==0).all():pg=p2=1.
   else:pg=wilcoxon(x.values,alternative='greater').pvalue;p2=wilcoxon(x.values,alternative='two-sided').pvalue
   r.loc[g,['p_greater','p_two_sided']]=[pg,p2]
 delta.to_csv(O/f'{ds}_{axis}_patient_effects.tsv.gz',sep='\t')
 return finish(r,ds,axis,'bulk',n,7,.5)

def deseq(counts,meta,genes,ds,axis,lineage,minpairs,weight):
 from pydeseq2.dds import DeseqDataSet
 from pydeseq2.ds import DeseqStats
 n=meta.patient.nunique();name=f'{ds}_{axis}_{lineage}'
 meta.to_csv(O/(name+'_units.tsv'),sep='\t');counts.to_csv(O/(name+'_counts.tsv.gz'),sep='\t')
 # Descriptive paired log2 CPM used even when below inferential minimum.
 cpm=np.log2(counts.div(counts.sum(axis=1),axis=0)*1e6+1)
 a=cpm.loc[meta.condition=='target'].copy();a.index=meta.loc[a.index,'patient']
 b=cpm.loc[meta.condition=='reference'].copy();b.index=meta.loc[b.index,'patient']
 delta=(a-b.reindex(a.index)).T.reindex(genes)
 r=paired_summary(delta);r['effect']=r.median_patient_effect;r['p_greater']=np.nan;r['p_two_sided']=np.nan
 delta.to_csv(O/(name+'_patient_log2CPM_differences.tsv.gz'),sep='\t')
 if n>=minpairs:
  keep=((counts>=10).sum(axis=0)>=np.ceil(len(counts)/3))&(counts.sum(axis=0)>=50)
  c=counts.loc[:,keep].astype(int);m=meta[['patient','condition']].copy()
  m.condition=pd.Categorical(m.condition,categories=['reference','target'])
  dds=DeseqDataSet(counts=c,metadata=m,design='~ patient + condition',refit_cooks=True,n_cpus=4,quiet=True)
  dds.deseq2();st=DeseqStats(dds,contrast=['condition','target','reference'],n_cpus=4,quiet=True);st.summary()
  rr=st.results_df;rr.to_csv(O/(name+'_all_genes_Wald.tsv.gz'),sep='\t')
  z=rr.reindex(genes);r['effect']=z.log2FoldChange;r['lfcSE']=z.lfcSE;r['effect_ci95_low']=z.log2FoldChange-1.96*z.lfcSE;r['effect_ci95_high']=z.log2FoldChange+1.96*z.lfcSE
  r['p_two_sided']=z.pvalue;r['p_greater']=np.where(z.log2FoldChange>0,z.pvalue/2,1-z.pvalue/2)
  norm=pd.DataFrame(dds.layers['normed_counts'],index=c.index,columns=c.columns)
  aa=norm.loc[meta.condition=='target'].copy();aa.index=meta.loc[aa.index,'patient']
  bb=norm.loc[meta.condition=='reference'].copy();bb.index=meta.loc[bb.index,'patient']
  nd=(np.log2(aa+.5)-np.log2(bb.reindex(aa.index)+.5)).T.reindex(genes)
  ps=paired_summary(nd);r[ps.columns]=ps
  nd.to_csv(O/(name+'_patient_normalized_log2FC.tsv.gz'),sep='\t')
 return finish(r,ds,axis,lineage,n,minpairs,weight)

def bulk():
 cand=pd.read_csv(O.parent/'step1/candidates_step1.tsv',sep='\t');genes=sorted(cand.gene.unique())
 audit=pd.read_csv(O/'sample_audit.tsv',sep='\t').fillna('')
 # GSE102511 file is already log2-normalized/batch-corrected per current GEO metadata.
 x=pd.read_csv(B/'raw/bulk/GSE102511_Smruthy_etal_allsamples_TPM.txt.gz',sep='\t',index_col=0).groupby(level=0).mean()
 ps=sorted(p for p in {c.split('-')[0] for c in x.columns} if p+'-AAH' in x and p+'-NL' in x)
 a=x[[p+'-AAH' for p in ps]].copy();b=x[[p+'-NL' for p in ps]].copy();a.columns=b.columns=ps
 wilcox_test(a,b,genes,'GSE102511','I_AAH_vs_NL');log(f'GSE102511 completed {len(ps)} pairs')
 x=pd.read_csv(B/'raw/bulk_GSE319666_HTAN_PML/GSE319666_counts.tsv.gz',sep='\t',index_col=0)
 # Library denominators include all input counts, including genes without a symbol.
 lib=x.sum(axis=0)
 mp={}
 with gzip.open(B/'raw/annotation/gencode.v44.basic.annotation.gtf.gz','rt') as f:
  for l in f:
   if l.startswith('#'):continue
   v=l.rstrip().split('\t')
   if v[2]=='gene':mp[re.search(r'gene_id "([^".]+)',v[8])[1]]=re.search(r'gene_name "([^"]+)',v[8])[1]
 x.index=[mp.get(i.split('.')[0],None) for i in x.index];x=x.loc[x.index.notna()].groupby(level=0).sum()
 lc=np.log2(x.div(lib,axis=1)*1e6+1)
 sm=audit[audit.dataset=='GSE319666'].copy();assert sm.matrix_column.is_unique
 sm['in_matrix']=sm.matrix_column.isin(lc.columns);table(sm,O/'GSE319666_column_match.tsv')
 assert sm.in_matrix.all() and set(sm.matrix_column)==set(lc.columns)
 def means(h):
  ss=sm[sm.histology==h];return lc[ss.matrix_column].T.groupby(ss.patient.values).mean().T
 nn=means('N')
 for hist,axis in [('AAH','I_AAH_vs_N'),('AIS/MIA','I_AISMIA_vs_N')]:
  aa=means(hist);ps=sorted(set(aa.columns)&set(nn.columns));wilcox_test(aa[ps],nn[ps],genes,'GSE319666',axis);log(f'GSE319666 {axis} {len(ps)} pairs')
 x=pd.read_csv(B/'raw/multi_GSE288479_partsolid/GSE288479_raw_counts_All_samples.txt.gz',index_col=0)
 x=x.drop(columns=['gene_id','gene_type']).groupby('gene_name').sum()
 sm=audit[(audit.dataset=='GSE288479')&(audit.gsm.str[3:].astype(int)>=8768576)&(audit.site!='Visium')]
 ps=sorted(sm.patient.unique());assert len(ps)==7
 for tar,ref,axis in [('GG','DN','I_GG_vs_DN'),('S','GG','regional_S_vs_GG')]:
  cols=[p+'_'+s+'_counts' for p in ps for s in [ref,tar]]
  meta=pd.DataFrame([dict(patient=p,condition='reference' if s==ref else 'target') for p in ps for s in [ref,tar]],index=cols)
  deseq(x[cols].T,meta,genes,'GSE288479',axis,'bulk',7,.5);log(f'GSE288479 {axis} {len(ps)} pairs')


def sc_inference(ds):
 import anndata as ad
 d=O/ds
 assert (d/'annotation_review.json').exists(),'Current annotation must be reviewed before any contrasts'
 a=ad.read_h5ad(d/'rebuilt.h5ad');review=json.loads((d/'annotation_review.json').read_text())
 ann=pd.read_csv(d/'cluster_annotation.tsv',sep='\t',index_col=0);ann.index=ann.index.astype(str)
 for cl,label in review.get('overrides',{}).items():ann.loc[str(cl),'label']=label
 a.obs['label']=a.obs.cluster.astype(str).map(ann.label)
 ep=np.asarray(a[:, 'EPCAM'].layers['counts'].toarray()).ravel() if 'EPCAM' in a.var_names else np.zeros(a.n_obs)
 pt=np.asarray(a[:, 'PTPRC'].layers['counts'].toarray()).ravel() if 'PTPRC' in a.var_names else np.zeros(a.n_obs)
 isepi=a.obs.label.isin(EPI).values
 # Same technical contamination cleanup as source-code validation definition.
 clean=~((isepi&(pt>=2))|(~isepi&(ep>=2)))
 a.obs['technical_clean']=clean
 a.obs[['patient','site','sample','cluster','label','technical_clean']].to_csv(d/'reviewed_cell_annotations.tsv.gz',sep='\t')
 cand=pd.read_csv(O.parent/'step1/candidates_step1.tsv',sep='\t')
 axes=[('T','A','I_precursor_vs_adjacent')] if ds=='GSE164789' else [('GG','DN','I_GG_vs_DN'),('S','GG','regional_S_vs_GG')] if ds=='GSE288479' else [('S','GG','regional_S_vs_GG')]
 cells=[]
 for lin,labs in VL.items():
  lm=a.obs.label.isin(labs).values&clean;genes=sorted(cand.loc[cand.lineage==lin,'gene'].unique())
  for tar,ref,axis in axes:
   counts=[];meta=[];ex=[]
   for patient in sorted(a.obs.patient.unique()):
    ix={s:np.where(lm&(a.obs.patient.values==patient)&(a.obs.site.values==s))[0] for s in [ref,tar]}
    valid=all(len(z)>=50 for z in ix.values())
    for site,indices in ix.items():cells.append(dict(dataset=ds,lineage=lin,axis=axis,patient=patient,site=site,n_cells=len(indices),pair_eligible=valid,reason='' if valid else 'one/both sides <50 cells'))
    if not valid:continue
    for site,indices in ix.items():
     counts.append(np.asarray(a.layers['counts'][indices].sum(axis=0)).ravel())
     meta.append(dict(unit=patient+'|'+site,patient=patient,condition='target' if site==tar else 'reference',n_cells=len(indices)))
   if not meta:
    r=pd.DataFrame(index=genes,columns=['effect','p_greater','p_two_sided'],dtype=float);finish(r,ds,axis,lin,0,8,1.)
   else:
    m=pd.DataFrame(meta).set_index('unit');c=pd.DataFrame(np.vstack(counts),index=m.index,columns=a.var_names)
    deseq(c,m,genes,ds,axis,lin,8,1.)
   log(f'{ds} {axis} {lin}: {len(meta)//2} pairs')
 table(pd.DataFrame(cells),d/'patient_cell_counts_and_exclusions.tsv')


def summarize():
 cand=pd.read_csv(O.parent/'step1/candidates_step1.tsv',sep='\t')
 evid=[]
 for f in sorted(O.glob('GSE*.tsv')):
  header=f.open().readline().strip().split('\t')
  if not {'dataset','axis','support','q_greater','status','lineage'}.issubset(header):continue
  r=pd.read_csv(f,sep='\t');evid.append(r)
 # Explicit unavailable results; no inferential substitute allowed for unpaired GSE189357.
 for lin,g in cand.groupby('lineage'):
  evid.append(pd.DataFrame({'gene':g.gene,'dataset':'GSE189357','axis':'II_IAC_vs_AISMIA','lineage':lin,'n_pairs':0,'weight':1.,'status':'无法验证：无可核实配对','support':False,'effect':np.nan,'p_greater':np.nan,'p_two_sided':np.nan,'q_greater':np.nan,'q_two_sided':np.nan}))
 ev=pd.concat(evid,ignore_index=True);table(ev,O/'validation_evidence_long.tsv')
 rows=[]
 by_gene={g:t for g,t in ev.groupby('gene',sort=False)}
 for flag in ['I_upstream_eligible','II_upstream_eligible']:cand[flag]=cand[flag].eq(True)
 for _,c in cand.iterrows():
  eg=by_gene[c.gene];e=eg[eg.lineage.isin(['bulk',c.lineage])]
  # Accession-level union: each independent cohort contributes at most its specified weight.
  sp=e[e.support.fillna(False)]
  score=sp.groupby('dataset').weight.max().sum()
  si=sp[sp.axis.str.startswith('I_')].groupby('dataset').weight.max().sum()
  sii=sp[sp.axis.str.startswith('II_')].groupby('dataset').weight.max().sum()
  sr=sp[sp.axis.str.startswith('regional_')].groupby('dataset').weight.max().sum()
  opp=e[e.status.eq('方向相反')]
  r={'gene':c.gene,'lineage':c.lineage,'I_upstream_eligible':c.I_upstream_eligible,'II_upstream_eligible':c.II_upstream_eligible,'D3':score,'D3_I':si,'D3_II_histology':sii,'D3_regional':sr,'D3_gate_pass':bool(score>=1.),'D3_I_ge1':bool(si>=1.),'D3_II_histology_ge1':bool(sii>=1.),'support_datasets':';'.join(sorted(sp.dataset.unique())),'support_axes':';'.join(sp.dataset+':'+sp.axis),'opposite_axes':';'.join(opp.dataset+':'+opp.axis),'has_opposite_significant':len(opp)>0,'evidence_scope':'I precursor comparisons and regional evidence separately; S>GG is not histological II','intervention_direction':'表达升高不推断抑制有益','D7_status':'未执行第3步；不沿用旧安全否决'}
  rows.append(r)
 out=pd.DataFrame(rows);table(out,O/'D3_candidates.tsv')
 full=cand.drop(columns=['D3_status']).merge(out.drop(columns=['I_upstream_eligible','II_upstream_eligible','intervention_direction']),on=['gene','lineage'],validate='one_to_one')
 full['D3_status']=np.where(full.D3_gate_pass,'D3数值门槛通过；轴别与反向证据须另读','D3数值门槛未通过')
 full['all_hard_gates_pass']=np.where(full.D3_gate_pass,'尚不能判定；D7未运行','D3未通过；D7未运行')
 table(full,O/'candidates_step2.tsv')
 summaries=[]
 for step in ['I','II']:
  for lin,g in out[out[step+'_upstream_eligible']].groupby('lineage'):
   summaries.append(dict(step=step,lineage=lin,n_upstream=len(g),D3_gate_pass=int(g.D3_gate_pass.sum()),I_axis_D3_ge1=int(g.D3_I_ge1.sum()),II_histology_D3_ge1=int(g.D3_II_histology_ge1.sum()),regional_any=int((g.D3_regional>0).sum()),significant_opposite_any=int(g.has_opposite_significant.sum())))
 table(pd.DataFrame(summaries),O/'validation_by_discovery_step.tsv')
 table(out[out.gene.isin(['SPP1','TREM2'])],O/'SPP1_TREM2_D3.tsv')
 # Assertions encode the scientific gates, not just file shape.
 assert len(out)==5207 and not out[['gene','lineage']].duplicated().any()
 assert not ev.loc[ev.n_pairs<ev.lineage.eq('bulk').map({True:7,False:8}),'support'].any()
 assert not ev.loc[ev.dataset.isin(['GSE189357','GSE203360']),'support'].any()
 assert not ev.loc[ev.status.eq('方向相反'),'support'].any()
 assert (out.D3<=2.5).all() # 1 scRNA + 3 independent bulk cohorts.
 assert (out.D3_II_histology==0).all()
 log('Step2 summaries complete; no ranking/no step3.')


def audit_cached_geo():
 rows=[]
 for ds in ['GSE164789','GSE189357','GSE288479','GSE102511','GSE319666','GSE203360']:
  p=O/(ds+'.soft'); text=p.read_text()
  for ch in text.split('^SAMPLE = ')[1:]:
   gsm=ch.splitlines()[0]; title=re.search(r'!Sample_title = (.*)',ch)[1]
   chars=dict(x.split(': ',1) for x in re.findall(r'!Sample_characteristics_ch1 = (.*)',ch) if ': ' in x)
   row=dict(dataset=ds,gsm=gsm,title=title,patient='',site='',histology='',matrix_column='',include=False,reason='')
   if ds=='GSE164789':
    if 'TCR' in title:continue
    row.update(patient=title.split('-')[0],site='T' if re.search(r'!Sample_source_name_ch1 = (.*)',ch)[1].endswith('tumor') else 'A',histology=chars['histology'],include=chars['histology'] in ['AIS','MIA'])
    row['reason']='I: AIS/MIA vs adjacent; multiple same-patient samples pooled' if row['include'] else 'IAC excluded from precursor comparison'
   elif ds=='GSE189357': row.update(patient='',site='T',histology=chars['histolgical type'],reason='No verifiable paired patient IDs; cannot validate; no D3')
   elif ds=='GSE288479':
    patient=chars['individuals']; site=re.search(r'\((DN|AN|GG|S)\)',title)
    row.update(patient=patient,site=site[1] if site else 'Visium',include=bool(site),reason='scRNA 4 patients: direction only' if int(gsm[3:])<=8768575 else 'bulk 7 patients' if site else 'Visium not used')
   elif ds=='GSE102511':row.update(patient=chars['donor'],site=title.split('-')[1],matrix_column=title,include=title.endswith(('-AAH','-NL')),reason='AAH vs NL; input already log2 transformed and batch corrected')
   elif ds=='GSE319666':row.update(patient=re.search('patient (.*?), derived',title)[1],histology=chars['sample histology'],matrix_column=re.search(r'\[(.*?)\]',title)[1],include=chars['sample histology'] in ['N','AAH','AIS/MIA'],reason='Mean log2 CPM per patient per histology')
   else:row.update(patient=title.split('-')[0],site=title.split('-')[1].split(',')[0],include=True,reason='3 patients: direction only, no D3')
   rows.append(row)
 a=pd.DataFrame(rows);a.to_csv(O/'sample_audit.tsv',sep='\t',index=False)

if __name__=='__main__':
 if sys.argv[1]=='preprocess':preprocess(sys.argv[2])
 elif sys.argv[1]=='bulk':bulk()
 elif sys.argv[1]=='sc':sc_inference(sys.argv[2])
 elif sys.argv[1]=='summarize':summarize()
 elif sys.argv[1]=='audit':audit_cached_geo()
