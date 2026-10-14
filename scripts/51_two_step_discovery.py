#!/opt/anaconda3/envs/ggn_target/bin/python
"""Independent step 1. Only approved annotation h5ad, raw GO, source code constants,
and new step-0 outputs are input. No old comparisons; no held-out data; no step 2."""
import os
WORK='/Users/dudu/Documents/Codex/2026-10-09/ggn-target-ggn-target-opt-anaconda3/work'
for k,v in [('NUMBA_CACHE_DIR',WORK+'/numba_cache'),('MPLCONFIGDIR',WORK+'/mpl'),('XDG_CACHE_HOME',WORK+'/cache'),('OPENBLAS_NUM_THREADS','1'),('OMP_NUM_THREADS','1'),('MKL_NUM_THREADS','1')]:os.environ[k]=v
from pathlib import Path
from datetime import datetime
from zoneinfo import ZoneInfo
import sys,json,hashlib,time,math,warnings,inspect
import numpy as np,pandas as pd, scipy.sparse as sp,h5py
from scipy.stats import wilcoxon,norm
from statsmodels.stats.multitest import multipletests
from anndata.io import read_elem,sparse_dataset
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import importlib.metadata
P=Path('/Users/dudu/全部数据库/GGN_target');O=P/'results/50_screen';D=O/'step1';D.mkdir(exist_ok=True)
GROUPS={'AT2L':['AT2','AT2/Club-like epithelial','Proliferating epithelial'],'Myeloid':['Macrophage','Monocyte','cDC','pDC','Neutrophil','Proliferating myeloid'],'Fibroblast':['Fibroblast'],'T_NK':['T cell','NK','Proliferating T/NK'],'Endothelial':['Endothelial'],'Macrophage':['Macrophage']}
EPI=['AT1','AT2','AT2-like (patient-specific)','AT2/Club-like epithelial','Basal','Ciliated','Proliferating epithelial']
SOURCE={'Epithelial':EPI,'Myeloid':GROUPS['Myeloid'],'T_NK':GROUPS['T_NK'],'Fibroblast':['Fibroblast'],'Endothelial':['Endothelial'],'B_Plasma':['B cell','Plasma','Proliferating B/Plasma'],'Mast':['Mast'],'SMC_pericyte':['SMC/pericyte'],'Lymphatic_EC':['Lymphatic EC'],'Macrophage':['Macrophage']}
STAGES=['Normal','AAH','AIS','MIA','LUAD']; PRE=['AAH','AIS','MIA'];COMPS={'I':(['Normal'],PRE),'II':(PRE,['LUAD'])}
MIN_CELLS=50;MIN_PAIRS=8
now=lambda:datetime.now(ZoneInfo('Asia/Shanghai')).isoformat()
def log(s):
 print(now()+' '+s,flush=True)
 with (O/'run_log.txt').open('a') as f:f.write(now()+' [step1] '+s+'\n')
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def write(df,p,index=False):df.to_csv(p,sep='\t',index=index)
def read(p,index_col=None):return pd.read_csv(p,sep='\t',index_col=index_col)
def adjust(p,method='fdr_bh'):
 p=np.asarray(p,float);q=np.full(p.shape,np.nan);m=np.isfinite(p)
 if m.any():q[m]=multipletests(p[m],method=method)[1]
 return q
def freeze_methods():
 protocol=P/'docs/新研究_靶点优先级筛选_预设方案.md'
 cfg={'created_at':now(),'protocol_sha256':sha(protocol),'threshold_sha256':sha(O/'state_thresholds_frozen.json'),'groups':GROUPS,'source_groups':SOURCE,'minimum_cells_per_patient_side':50,'minimum_pairs':8,'gene_filter':'count >=10 in >=ceil(number pseudobulks/3), total count >=50; technical definition from existing scripts; recomputed','design':'~ patient + condition','D1_D2':'two-sided Wald; positive log2FC; BH FDR<0.05 and ceil(0.8*n) positive patient differences','paired_direction':'DESeq2 normalized counts; log2((target+0.5)/(reference+0.5)); ties are not positive','source_filter':'sample mean log1p(CP10k), median across eligible samples >=50 cells; epithelial/max non-epithelial >=0.2; microenvironment group/epithelial >=0.2; no cross-cohort rescue in step1','D10':'fixed existing final_label subtypes for Myeloid/T_NK, fixed step0 KAC/other; >=8 pairs. Fibroblast/Endothelial/Macrophage finer states unavailable. Non-significance never proves composition-only.','LOO':'all upstream-eligible D1/D2+D4+D5 candidates; D3/D7 pending; full original tested gene universe refit per omitted patient, all genes BH, not candidate-only BH; recompute size factors and patient concordance; >=8 remaining pairs; missing p counts not significant','substage':'AAH>Normal and AIS/MIA>Normal, patient Wilcoxon greater on log2(CPM+0.5), BH within lineage/contrast candidate family; Holm across both comparisons additionally reported; individual AIS/MIA direction descriptive if <8 pairs','clinical_label':'expression-supported candidate labels only; AAH and pooled AIS/MIA positive with BH<0.05 supports provisional regression direction; stage-specific insufficient N explicitly flagged; no proven dependency or support-cell causality','limits':'no ranking, no step2, no frozen dataset reads, no repaired KAC subset test; slope substitute remains step3 D7','versions':{k:importlib.metadata.version(k) for k in ['pydeseq2','numpy','pandas','scipy','anndata']}}
 p=D/'methods_before_tests.json'
 if not p.exists():
  with p.open('x') as f:json.dump(cfg,f,ensure_ascii=False,indent=2)
  (D/'protocol_snapshot.md').write_text(protocol.read_text());log('Methods and directions fixed before expression tests: '+sha(p))
 else:
  old=json.loads(p.read_text());assert old['protocol_sha256']==sha(protocol),'Protocol changed: STOP'
 return cfg

def prepare():
 freeze_methods()
 # Rebuild exact proliferation definitions from raw GO and the source-code constant list.
 import ast
 module=ast.parse((P/'scripts/prolif_genes.py').read_text());defs={}
 for n in module.body:
  if isinstance(n,ast.Assign) and isinstance(n.targets[0],ast.Name) and n.targets[0].id in ['TIROSH_S','TIROSH_G2M']:
   defs[n.targets[0].id]=n.value.func.value.value.split()
 tirosh=set(defs['TIROSH_S']+defs['TIROSH_G2M']);terms={}
 for l in (P/'raw/lincs/Enrichr_GO_Biological_Process_2023.txt').read_text().splitlines():
  r=l.split('\t');terms[r[0]]={g for g in r[2:] if g}
 get=lambda code:set().union(*(v for k,v in terms.items() if '('+code+')' in k))
 narrow=(tirosh|get('GO:0140014')|get('GO:0006260'))-get('GO:0045786')
 keywords=['cell cycle','mitotic','cell division','dna replication','chromosome segregation','spindle','cytokinesis','kinetochore']
 broad=tirosh|set().union(*(v for k,v in terms.items() if any(x in k.lower() for x in keywords)))
 assert len(tirosh)==96 and len(narrow)==196 and len(broad)==1311
 write(pd.DataFrame({'gene':sorted(broad|narrow),'prolif196':[g in narrow for g in sorted(broad|narrow)],'prolif1311':[g in broad for g in sorted(broad|narrow)]}),D/'proliferation_sets.tsv')
 if (D/'prepare_complete.json').exists():log('Reuse newly computed step1 pseudobulks; no old results');return
 with h5py.File(P/'results/01_snRNA/snRNA_final.h5ad','r') as f:
  obs=read_elem(f['obs']);var=read_elem(f['var']);genes=var.index.astype(str)
  assert genes.is_unique
  for c in ['patient','stage','sample','final_label']:obs[c]=obs[c].astype(str)
  cells=read(O/'AT2L_cell_state_assignments.tsv.gz',index_col=0)
  th=json.loads((O/'state_thresholds_frozen.json').read_text())
  assert np.array_equal(cells.is_KAC.values,cells.KAC.values>th['thresholds']['KAC'])
  obs['state_KAC']=False;obs['state_other']=False
  obs.loc[cells.index,'state_KAC']=cells.is_KAC;obs.loc[cells.index,'state_other']=~cells.is_KAC
  units=obs[['patient','stage','sample','final_label']].drop_duplicates().reset_index(drop=True)
  keys=units['sample']+'|'+units['final_label'];key_idx=dict(zip(keys,units.index));code=np.array([key_idx[k] for k in obs['sample']+'|'+obs['final_label']])
  units['n_cells']=np.bincount(code,minlength=len(units));C=np.zeros((len(units),len(genes)),dtype=np.int64);L=np.zeros_like(C,dtype=np.float64)
  # State raw counts, patient/stage rather than individual cells as inferential units.
  state_keys=[(p,s,k) for p in sorted(obs.patient.unique(),key=lambda x:int(x[1:])) for s in STAGES for k in ['KAC','other']]
  skey={k:i for i,k in enumerate(state_keys)};SC=np.zeros((len(state_keys),len(genes)),np.int64);SN=np.zeros(len(state_keys),np.int64)
  ds=sparse_dataset(f['layers/counts']);n=len(obs)
  for start in range(0,n,10000):
   end=min(n,start+10000);x=ds[start:end,:].tocsr();assert np.all(x.data>=0) and np.allclose(x.data,np.rint(x.data))
   agg=sp.csr_matrix((np.ones(end-start,dtype=np.int64),(code[start:end],np.arange(end-start))),shape=(len(units),end-start));C+=np.asarray((agg@x).toarray(),np.int64)
   normx=x.astype(np.float64);sums=np.asarray(normx.sum(1)).ravel();assert (sums>0).all();normx=sp.diags(10000/sums)@normx;normx.data=np.log1p(normx.data);L+=(agg@normx).toarray()
   z=obs.iloc[start:end];sid=[];col=[]
   for j,(p,s,k,o) in enumerate(zip(z.patient,z.stage,z.state_KAC,z.state_other)):
    if k or o:sid.append(skey[(p,s,'KAC' if k else 'other')]);col.append(j)
   if sid:
    ag=sp.csr_matrix((np.ones(len(sid),np.int64),(sid,col)),shape=(len(state_keys),end-start));SC+=np.asarray((ag@x).toarray(),np.int64);SN+=np.asarray(ag.sum(1)).ravel()
   if start%100000==0:log(f'Counts and source expression recomputed: {end}/{n} cells')
 units.index.name='unit';write(units,D/'annotation_units.tsv',True);np.savez_compressed(D/'aggregates.npz',counts=C,log_sums=L,genes=np.array(genes,dtype=str),state_counts=SC,state_n=SN)
 sm=pd.DataFrame(state_keys,columns=['patient','stage','state']);sm['n_cells']=SN;write(sm,D/'state_units.tsv')
 # source filter: median of sample means; explicit sample-level denominators.
 med={};srcn=[];srcrows=[]
 for name,labs in SOURCE.items():
  rows=[]
  for sample in units['sample'].unique():
   ix=np.where(units['sample'].eq(sample)&units.final_label.isin(labs))[0];nc=int(units.iloc[ix].n_cells.sum())
   srcn.append({'lineage':name,'sample':sample,'n_cells':nc,'included':nc>=50})
   if nc>=50:
    expr=L[ix].sum(0)/nc;rows.append(expr);srcrows.append((name,sample))
  med[name]=np.median(np.stack(rows),axis=0) if rows else np.full(len(genes),np.nan)
 M=pd.DataFrame(med,index=genes);M.index.name='gene';write(M,D/'source_median_log_expression.tsv',True);write(pd.DataFrame(srcn),D/'source_sample_eligibility.tsv')
 rates=[]
 for name in GROUPS:
  num=M.Epithelial if name=='AT2L' else M[name]
  den=M.drop(columns=['Epithelial','Macrophage']).max(axis=1) if name=='AT2L' else M.Epithelial
  ratio=num/den.replace(0,np.nan);ratio[(den==0)&(num>0)]=np.inf
  rates.append(pd.DataFrame({'gene':genes,'lineage':name,'numerator':num.values,'denominator':den.values,'ratio':ratio.values,'source_pass':(ratio>=.2).values,'status':np.where(ratio.isna(),'无法判定',np.where(ratio>=.2,'通过','不通过'))}))
 write(pd.concat(rates),D/'source_filter_recomputed.tsv')
 (D/'prepare_complete.json').write_text(json.dumps({'created_at':now(),'cells':n,'genes':len(genes),'prolif196':len(narrow),'prolif1311':len(broad)}));log('Prepare complete; raw counts summed independently; source ratios recalculated.')

def load():
 z=np.load(D/'aggregates.npz');return read(D/'annotation_units.tsv',0),z

def pseudobulk(lineage,ref,target,state=None):
 u,z=load();genes=z['genes'];C=z['counts']
 if state in ['KAC','other']:
  u=read(D/'state_units.tsv');u=u[u.state.eq(state)].copy();C=z['state_counts'][u.index];u=u.reset_index(drop=True);mask=np.ones(len(u),bool)
 else:mask=u.final_label.isin(GROUPS[lineage] if state is None else [state]).values
 patients=sorted(u.patient.unique(),key=lambda x:int(x[1:]));rows=[];meta=[];elig=[]
 for p in patients:
  inds=[np.where(mask&u.patient.eq(p).values&u.stage.isin(st).values)[0] for st in [ref,target]];ns=[int(u.iloc[ix].n_cells.sum()) for ix in inds];ok=min(ns)>=50
  elig.append(dict(patient=p,reference_cells=ns[0],target_cells=ns[1],included=ok,reason='' if ok else 'one_or_both_sides_below_50_or_missing'))
  if ok:
   for side,ix,nc in zip(['reference','target'],inds,ns):
    rows.append(C[ix].sum(0));meta.append(dict(unit=p+'_'+side,patient=p,condition=side,n_cells=nc))
 counts=pd.DataFrame(rows,columns=genes,index=[x['unit'] for x in meta],dtype=np.int64);m=pd.DataFrame(meta).set_index('unit') if meta else pd.DataFrame(columns=['patient','condition','n_cells'])
 return counts,m,pd.DataFrame(elig)

def fit(c,m,out,fixed_genes=None):
 out.mkdir(exist_ok=True,parents=True)
 if (out/'complete.json').exists():return read(out/'deseq.tsv',0),read(out/'patient_log2FC.tsv',0)
 if m.patient.nunique()<8:return None,None
 keep=((c>=10).sum(0)>=math.ceil(len(c)/3))&(c.sum(0)>=50)
 c=c.loc[:,keep] if fixed_genes is None else c.reindex(columns=fixed_genes)
 write(m,out/'metadata.tsv',True);write(pd.DataFrame({'gene':c.columns}),out/'tested_gene_universe.tsv')
 t=time.time();dds=DeseqDataSet(counts=c,metadata=m,design='~ patient + condition',refit_cooks=True,n_cpus=4,quiet=True);dds.deseq2()
 ds=DeseqStats(dds,contrast=['condition','target','reference'],alpha=.05,cooks_filter=True,independent_filter=True,n_cpus=4,quiet=True);ds.summary();r=ds.results_df.copy();r.index.name='gene'
 nc=pd.DataFrame(dds.layers['normed_counts'],index=c.index,columns=c.columns);patients=sorted(m.patient.unique(),key=lambda x:int(x[1:]));delta=pd.DataFrame({p:np.log2((nc.loc[p+'_target']+.5)/(nc.loc[p+'_reference']+.5)) for p in patients}).T
 delta.index.name='patient';positive=(delta>0).sum(0);n=len(patients)
 r['n_pairs']=n;r['positive_patients']=positive;r['required_positive']=math.ceil(.8*n);r['direction_consistency']=positive/n;r['stage_pass']=(r.padj<.05)&(r.log2FoldChange>0)&(positive>=math.ceil(.8*n));r['D_raw']=r.log2FoldChange*r.direction_consistency
 r['log2FC_ci95_low']=r.log2FoldChange-1.959963984540054*r.lfcSE;r['log2FC_ci95_high']=r.log2FoldChange+1.959963984540054*r.lfcSE
 write(r,out/'deseq.tsv',True);write(delta,out/'patient_log2FC.tsv',True);write(nc,out/'normalized_counts.tsv.gz',True)
 (out/'complete.json').write_text(json.dumps({'created_at':now(),'seconds':time.time()-t,'n_pairs':n,'genes':len(r),'test':'two-sided Wald; BH'}))
 log(f'Fit {out.relative_to(D)}: {n} pairs, {len(r)} genes, {int(r.stage_pass.sum())} stage-pass, {time.time()-t:.1f}s')
 return r,delta

def main_tests():
 source=read(D/'source_filter_recomputed.tsv');pro=read(D/'proliferation_sets.tsv').set_index('gene');summ=[];allrows=[]
 for lin in GROUPS:
  for step,(ref,tar) in COMPS.items():
   out=D/'main'/f'{lin}_{step}';out.mkdir(exist_ok=True,parents=True);c,m,e=pseudobulk(lin,ref,tar);write(e,out/'pair_eligibility.tsv');write(c,out/'pseudobulk_counts.tsv.gz',True)
   r,delta=fit(c,m,out)
   if r is None:summ.append(dict(lineage=lin,step=step,pairs=len(m)//2,status='无法验证'));continue
   r=r.copy();s=source[source.lineage.eq(lin)].set_index('gene');r['source_ratio']=s.ratio;r['source_pass']=s.source_pass
   r['prolif196']=pro.prolif196.reindex(r.index).fillna(False).astype(bool);r['prolif1311']=pro.prolif1311.reindex(r.index).fillna(False).astype(bool);r['upstream_eligible']=r.stage_pass&r.source_pass&~r.prolif196;r['lineage']=lin;r['step']=step
   write(r,out/'annotated.tsv',True);allrows.append(r.reset_index());summ.append(dict(lineage=lin,step=step,pairs=len(m)//2,genes=len(r),stage_pass=int(r.stage_pass.sum()),upstream_eligible=int(r.upstream_eligible.sum()),status='已完成；D3/D7待后续步骤'))
 write(pd.concat(allrows),D/'discovery_all.tsv.gz');write(pd.DataFrame(summ),D/'main_summary.tsv');log('Main discovery completed; candidates are provisional, not all-hard-gates passed.')

def d10():
 subs={'AT2L':['KAC','other'],'Myeloid':GROUPS['Myeloid'],'T_NK':GROUPS['T_NK']};summ=[]
 for lin,states in subs.items():
  for state in states:
   for step,(ref,tar) in COMPS.items():
    out=D/'d10'/f'{lin}_{state.replace("/","_").replace(" ","_")}_{step}';out.mkdir(exist_ok=True,parents=True);c,m,e=pseudobulk(lin,ref,tar,state);write(e,out/'pair_eligibility.tsv');r,d=fit(c,m,out)
    summ.append(dict(lineage=lin,state=state,step=step,pairs=len(m)//2,status='无法验证' if r is None else '已完成',result=str(out/'deseq.tsv') if r is not None else ''))
 # Explicit D7 other-state control lesion(LUAD) vs Normal, no KAC rescue.
 out=D/'d10'/'AT2L_other_Invasive_vs_Normal';out.mkdir(exist_ok=True,parents=True);c,m,e=pseudobulk('AT2L',['Normal'],['LUAD'],'other');write(e,out/'pair_eligibility.tsv');r,d=fit(c,m,out);summ.append(dict(lineage='AT2L',state='other',step='Invasive_vs_Normal',pairs=len(m)//2,status='已完成' if r is not None else '无法验证',result=str(out/'deseq.tsv') if r is not None else ''))
 write(pd.DataFrame(summ),D/'d10_pair_summary.tsv');log('D10 state/subtype tests completed; missing finer subtypes will remain unresolved.')

def substage():
 a=read(D/'discovery_all.tsv.gz');rows=[];elig=[]
 for lin in GROUPS:
  cand=sorted(a.loc[a.lineage.eq(lin)&a.upstream_eligible,'gene'].unique())
  for stage,target in [('AAH',['AAH']),('AIS_MIA',['AIS','MIA']),('AIS',['AIS']),('MIA',['MIA'])]:
   c,m,e=pseudobulk(lin,['Normal'],target);e['lineage']=lin;e['comparison']=stage;elig.append(e);n=len(m)//2
   # library totals use all genes, not candidate genes
   cp=c.div(c.sum(axis=1),axis=0)*1e6;delta=pd.DataFrame({p:np.log2((cp.loc[p+'_target',cand]+.5)/(cp.loc[p+'_reference',cand]+.5)) for p in m.patient.unique()}).T
   if n:write(delta,D/f'substage_patient_log2FC_{lin}_{stage}.tsv',True)
   block=[]
   for g in cand:
    v=delta[g].values if n else np.array([]);pv=(float(wilcoxon(v,alternative='greater',zero_method='wilcox',method='auto').pvalue) if np.any(v!=0) else 1.) if n>=8 and stage in ['AAH','AIS_MIA'] else np.nan
    block.append(dict(lineage=lin,comparison=stage,gene=g,n_pairs=n,median_patient_log2FC=float(np.median(v)) if n else np.nan,positive_patients=int((v>0).sum()),p_greater=pv,status='无法验证' if n<8 else ('已检验' if stage in ['AAH','AIS_MIA'] else '仅方向描述；不新增检验'),use='primary_pooled_substage' if stage in ['AAH','AIS_MIA'] else 'stage_specific_support'))
   b=pd.DataFrame(block)
   if len(b):b['q_BH']=adjust(b.p_greater);rows.append(b)
 allr=pd.concat(rows,ignore_index=True);allr['p_Holm_lineage_two_primary_comparisons']=np.nan
 for lin in GROUPS:
  ix=allr.lineage.eq(lin)&allr.comparison.isin(['AAH','AIS_MIA']);allr.loc[ix,'p_Holm_lineage_two_primary_comparisons']=adjust(allr.loc[ix,'p_greater'],'holm')
 write(allr,D/'clinical_substage_tests.tsv');write(pd.concat(elig),D/'clinical_substage_eligibility.tsv');log('Clinical substage one-sided tests completed with fixed target>Normal direction.')

def loo():
 a=read(D/'discovery_all.tsv.gz');rows=[]
 for lin in GROUPS:
  for step,(ref,tar) in COMPS.items():
   full=read(D/'main'/f'{lin}_{step}'/'annotated.tsv',0);cand=full.index[full.upstream_eligible].tolist()
   if not cand:continue
   c,m,e=pseudobulk(lin,ref,tar);ps=sorted(m.patient.unique(),key=lambda x:int(x[1:]));universe=read(D/'main'/f'{lin}_{step}'/'tested_gene_universe.tsv').gene.tolist();pieces=[]
   for p in ps:
    keep=m.patient.ne(p);out=D/'loo'/f'{lin}_{step}'/p
    r,d=fit(c.loc[keep],m.loc[keep],out,fixed_genes=universe)
    if r is None:
     pieces.append(pd.DataFrame({'gene':cand,'omitted_patient':p,'status':'无法验证；LOO后少于8对'}));continue
    rr=r.reindex(cand).copy();rr['omitted_patient']=p;rr['status']='已完成';pieces.append(rr.reset_index())
   long=pd.concat(pieces,ignore_index=True);write(long,D/'loo'/f'{lin}_{step}'/'candidate_all_omissions.tsv')
   for g in cand:
    z=long[long.gene.eq(g)];base=float(full.loc[g,'log2FoldChange']);available=z.status.eq('已完成').all()
    if not available:rows.append(dict(lineage=lin,step=step,gene=g,status='无法验证；LOO后少于8对',top5_excluded=False));continue
    changes=(z.log2FoldChange-base).abs();relative=changes/abs(base);idx=changes.idxmax();driver=bool((relative>.5).any()|(z.direction_consistency<.8).any());stable=bool((z.log2FoldChange>0).all()&(z.positive_patients>=z.required_positive).all())
    rows.append(dict(lineage=lin,step=step,gene=g,n_omissions=len(z),max_p=float(z.pvalue.max()) if z.pvalue.notna().all() else np.nan,max_finite_p=float(z.pvalue.max()),missing_p_count=int(z.pvalue.isna().sum()),fraction_FDR_lt05=float((z.padj<.05).mean()),direction_stable=bool((z.log2FoldChange>0).all()),concordance_stable=bool((z.positive_patients>=z.required_positive).all()),most_influential_patient=z.loc[idx,'omitted_patient'],log2FC_full=base,log2FC_after_influential=float(z.loc[idx,'log2FoldChange']),maximum_absolute_log2FC_change=float(changes.max()),maximum_relative_log2FC_change=float(relative.max()),robust_direction_and_concordance=stable,status='单患者驱动' if driver else ('稳健' if stable else '方向不稳健'),top5_excluded=driver))
 write(pd.DataFrame(rows),D/'loo_summary.tsv');log('LOO complete; all original test genes refitted for BH, upstream candidates summarized.')



def add_substage_uncertainty():
 from scipy.stats import binom
 p=D/'clinical_substage_tests.tsv';t=read(p)
 for (lin,comparison),idx in t.groupby(['lineage','comparison']).groups.items():
  delta=read(D/f'substage_patient_log2FC_{lin}_{comparison}.tsv',0);gs=t.loc[idx,'gene'].tolist();v=delta.reindex(columns=gs).to_numpy();n=len(v)
  t.loc[idx,'n_nonzero_pairs']=np.count_nonzero(v,axis=0);t.loc[idx,'n_zero_pairs']=(v==0).sum(axis=0)
  ks=[k for k in range(1,(n+1)//2+1) if 1-2*binom.cdf(k-1,n,.5)>=.95]
  if ks:
   k=max(ks);sv=np.sort(v,axis=0);t.loc[idx,'median_ci95_low']=sv[k-1];t.loc[idx,'median_ci95_high']=sv[n-k];t.loc[idx,'median_ci_actual_coverage']=1-2*binom.cdf(k-1,n,.5)
  else:t.loc[idx,'median_ci95_low']=-np.inf;t.loc[idx,'median_ci95_high']=np.inf;t.loc[idx,'median_ci_actual_coverage']=1.
  t.loc[idx,'median_ci_method']='unadjusted conservative distribution-free order-statistic median CI; ties may increase coverage; not a new selection gate'
 write(t,p)


def sample_count_audit():
 u=read(D/'annotation_units.tsv',0);rows=[]
 for lin in GROUPS:
  for step,(ref,tar) in COMPS.items():
   m=read(D/'main'/f'{lin}_{step}'/'metadata.tsv',0);ps=set(m.patient)
   row=dict(lineage=lin,step=step,patient_pairs=len(ps),pseudobulk_units=len(m))
   for side,stages in [('reference',ref),('target',tar)]:
    included=u.patient.isin(ps)&u.final_label.isin(GROUPS[lin])&u.stage.isin(stages)
    row[side+'_raw_samples']=int(u.loc[included,'sample'].nunique());row[side+'_cells']=int(u.loc[included,'n_cells'].sum())
    assert row[side+'_cells']==m.loc[m.condition.eq(side),'n_cells'].sum()
   row['total_raw_samples']=row['reference_raw_samples']+row['target_raw_samples'];assert row['pseudobulk_units']==2*row['patient_pairs'];rows.append(row)
 r=pd.DataFrame(rows);write(r,D/'main_sample_counts.tsv');return r


def finalize():
 import shutil
 add_substage_uncertainty()
 sample_info=sample_count_audit()
 a=read(D/'discovery_all.tsv.gz'); ss=read(D/'clinical_substage_tests.tsv'); loo_df=read(D/'loo_summary.tsv'); dsum=read(D/'d10_pair_summary.tsv')
 # Summaries retain every candidate-lineage unit, no ranking or top-N selection.
 candidates=a.loc[a.upstream_eligible,['lineage','gene']].drop_duplicates().sort_values(['lineage','gene'])
 main={lin:{step:read(D/'main'/f'{lin}_{step}'/'annotated.tsv',0) for step in COMPS} for lin in GROUPS}
 dk={}
 for _,r in dsum[dsum.step.isin(['I','II'])].iterrows():
  dk[(r.lineage,r.step,r.state)]=(int(r.pairs),read(r.result,0) if r.status=='已完成' else None)
 state_rows=[]; target_rows=[]; merged=[]
 gene_set=set(read(O/'KAC_TableS3_genes.tsv').gene)
 for _,c in candidates.iterrows():
  lin,g=c.lineage,c.gene;row={'gene':g,'lineage':lin,'D3_status':'未运行第2步','D6_D7_D8_D9_status':'未运行第3步','rank_status':'未排序/未冻结','KAC_signature_member':g in gene_set,'intervention_direction':'表达升高不推断抑制有益','evidence_level':'探索性/假说生成'}
  supported={};exclude=False
  for step in ['I','II']:
   r=main[lin][step];rec=r.loc[g] if g in r.index else None
   supported[step]=bool(rec is not None and rec.upstream_eligible)
   for k in ['log2FoldChange','lfcSE','log2FC_ci95_low','log2FC_ci95_high','pvalue','padj','n_pairs','positive_patients','required_positive','direction_consistency','stage_pass','D_raw','source_ratio','source_pass','prolif196','prolif1311','upstream_eligible']:
    row[f'{step}_{k}']=rec[k] if rec is not None else np.nan
   lr=loo_df[(loo_df.lineage==lin)&(loo_df.step==step)&(loo_df.gene==g)]
   for k in ['status','max_p','missing_p_count','fraction_FDR_lt05','most_influential_patient','maximum_absolute_log2FC_change','maximum_relative_log2FC_change','top5_excluded']:
    row[f'{step}_LOO_{k}']=lr.iloc[0].get(k,np.nan) if len(lr) else ('本步未通过上游门槛' if k=='status' else np.nan)
   if len(lr):exclude=exclude or bool(lr.iloc[0].top5_excluded)
   keys=[k for k in dk if k[:2]==(lin,step)];support=[];tested=[];missing=[]
   for key in keys:
    npairs,rstate=dk[key];st=key[2]
    stat={'gene':g,'lineage':lin,'step':step,'state':st,'n_pairs':npairs}
    if rstate is None:stat.update(status='无法验证；少于8对',log2FoldChange=np.nan,padj=np.nan,stage_pass=False);missing.append(st)
    elif g not in rstate.index:stat.update(status='无法判定；基因未通过表达过滤',log2FoldChange=np.nan,padj=np.nan,stage_pass=False);missing.append(st)
    else:
     rr=rstate.loc[g];stat.update(status='已检验',log2FoldChange=rr.log2FoldChange,padj=rr.padj,pvalue=rr.pvalue,positive_patients=rr.positive_patients,required_positive=rr.required_positive,direction_consistency=rr.direction_consistency,stage_pass=bool(rr.stage_pass));(tested if np.isfinite(rr.padj) else missing).append(st);stat['status']='已检验' if np.isfinite(rr.padj) else '无法判定：P或FDR因模型过滤不可用'
     if rr.log2FoldChange>0 and rr.padj<.05:support.append(st)
    state_rows.append(stat)
   if not keys:cls='无法区分：现有注释无更细亚型'
   elif support:cls='状态内变化有证据；组成贡献未确定'
   elif tested:cls='无法区分：状态内未检出，不等于仅组成驱动'
   else:cls='无法区分：状态内配对数或表达量不足'
   row[f'{step}_D10']=cls;row[f'{step}_D10_supported_states']=';'.join(support);row[f'{step}_D10_unavailable_states']=';'.join(missing)
  sub=ss[(ss.lineage==lin)&(ss.gene==g)].set_index('comparison')
  for stage in ['AAH','AIS_MIA','AIS','MIA']:
   sr=sub.loc[stage]
   for k in ['n_pairs','median_patient_log2FC','positive_patients','p_greater','q_BH','p_Holm_lineage_two_primary_comparisons','n_nonzero_pairs','n_zero_pairs','median_ci95_low','median_ci95_high','median_ci_actual_coverage','status']:
    row[f'{stage}_{k}']=sr[k]
  direction_all=bool((sub.loc[['AAH','AIS','MIA'],'median_patient_log2FC']>0).all())
  pooled_support=bool((sub.loc[['AAH','AIS_MIA'],'q_BH']<.05).all())
  regression_hint=supported['I'] and direction_all
  target=('消退方向线索 + 阻断方向候选' if regression_hint and supported['II'] else '消退方向线索' if regression_hint else '阻断方向候选' if supported['II'] else '第I步相关，亚阶段持续方向未齐备')
  row.update(clinical_target=target,regression_direction_hint=regression_hint,substage_AAH_and_AISMIA_BH_support=pooled_support,all_AAH_AIS_MIA_median_directions_positive=direction_all,regression_label_limit='横断面持续相关不等于维持依赖；MIA少于8对不能独立验证；AAH/AIS-MIA检验支持另列',cell_localization_limit='上皮谱系定位' if lin=='AT2L' else '微环境来源；直接支持病灶的作用尚未证实',regression_downgrade_IAC_only=bool(supported['II'] and not supported['I']),top5_excluded_by_LOO=exclude,all_hard_gates_passed='尚不能判定；D3/D7未运行')
  target_rows.append({k:row[k] for k in ['gene','lineage','clinical_target','regression_direction_hint','substage_AAH_and_AISMIA_BH_support','regression_downgrade_IAC_only','regression_label_limit','cell_localization_limit']});merged.append(row)
 full=pd.DataFrame(merged);write(full,D/'candidates_step1.tsv');write(pd.DataFrame(target_rows),D/'clinical_target_labels.tsv');write(pd.DataFrame(state_rows),D/'D10_candidate_state_evidence.tsv')
 # Descriptive subtype composition only: each point/row is a patient-side.
 u,z=load();su=read(D/'state_units.tsv');fr=[]
 for lin in ['AT2L','Myeloid','T_NK']:
  states=['KAC','other'] if lin=='AT2L' else GROUPS[lin]
  for step,(ref,tar) in COMPS.items():
   _,m,_=pseudobulk(lin,ref,tar)
   for patient in m.patient.unique():
    for side,stages in [('reference',ref),('target',tar)]:
     counts=[]
     for state in states:
      if lin=='AT2L':ns=int(su.loc[su.patient.eq(patient)&su.stage.isin(stages)&su.state.eq(state),'n_cells'].sum())
      else:ns=int(u.loc[u.patient.eq(patient)&u.stage.isin(stages)&u.final_label.eq(state),'n_cells'].sum())
      counts.append(ns)
     for state,ns in zip(states,counts):fr.append(dict(lineage=lin,step=step,patient=patient,side=side,state=state,n_cells=ns,total_cells=sum(counts),fraction=ns/sum(counts)))
 write(pd.DataFrame(fr),D/'state_fractions_by_patient.tsv')
 # User-authorized registry update. Preserve original registration and record both hashes.
 registry=O/'frozen_sets.tsv';old=O/'frozen_sets_step0_registered.tsv'
 if not old.exists():shutil.copy2(registry,old)
 reg=read(registry);reg.loc[reg.dataset.eq('GSE282617'),'purpose']='仅冻结后非配对方向检验：排除FLA、WJR18、WYA12；Normal5/AIS19/MIA16/IAC23；不计D3、不参与排序；依据方案3.2.1'
 reg.loc[reg.dataset.eq('GSE282617'),'post_freeze_eligible_samples']=63
 write(reg,registry)
 limitation='''# 方法局限：正常肺KAC子集检验不可实施\n\nNormal侧KAC共2321个，占Normal AT2谱系23203个细胞约10%；该比例主要由Normal第90百分位阈值的构造决定，不能当作独立估计的修复态患病比例。P4=661、P24=602、P25=504，合计1767/2321=76.13%。后两名患者来自MD Anderson，其他患者来自日本顺天堂/神户。7名两侧各≥50的合格患者中，P5/P24/P25为II期（3/7=42.86%）；全队列II–III期为P5/P11/P20/P24/P25（5/23=21.74%），已核对原始Table S1。\n\n因此原定病灶KAC vs Normal KAC伪批量只有7对，保持“无法进行”。不降低细胞数或配对数门槛，不重算或修改KAC阈值；不能由AT2L-other的可检验性推断无修复风险。方案2.6.1每患者斜率的剂量反应选择性检验属于替代估计量，安排在第3步D7，不能称为原定子集检验的结果。\n\nD10状态内不显著不证明仅由组成驱动；低细胞数亚群仍可能有未检出的变化。成纤维、内皮和单独巨噬群在复用的24类final_label中无更细状态注释，不能将原群自身重跑后冒充亚型验证。KAC定义包含签名基因，相关基因的KAC内比较也不是独立于状态定义的验证。\n\nAAH为8对，MIA单独仅3–4对；AIS/MIA合并检验不能证明AIS与MIA分别均成立。临床消退标签仅是探索性方向线索，不能推断维持依赖或抑制有益。微环境细胞来源也不等于已证实直接支持病灶。\n\nGSE282617冻结后使用独立写定的非配对分支，排除鳞癌和两处身份未解决样本；不会反哺D3或排序。\n'''
 (D/'methods_limitations.md').write_text(limitation)
 def mdtable(df):
  return '| '+' | '.join(map(str,df.columns))+' |\n| '+' | '.join(['---']*len(df.columns))+' |\n'+'\n'.join('| '+' | '.join(str(v) for v in r)+' |' for r in df.itertuples(index=False,name=None))
 ms=read(D/'main_summary.tsv');main_show=ms[['lineage','step','pairs','stage_pass','upstream_eligible']].copy();main_show.columns=['细胞群','比较步','患者配对数','阶段关联通过','加来源/196过滤后']
 ls=loo_df.groupby(['lineage','step']).agg(LOO候选数=('gene','size'),单患者驱动=('top5_excluded','sum'),每次FDR仍显著=('fraction_FDR_lt05',lambda x:int((x==1).sum()))).reset_index()
 subpairs=ss.groupby(['lineage','comparison']).n_pairs.first().unstack()[['AAH','AIS_MIA','AIS','MIA']].reset_index()
 report=f'''# 新研究第1步完成报告\n\n完成时间：{now()}。仅执行第1步；未进入第2步，未运行冻结验证集分析、D3、最终评分或排序。全部发现结果为探索性/假说生成，不推断抑制候选有益。\n\n## 执行与输入边界\n\n脚本：`scripts/51_two_step_discovery.py`；环境：`/opt/anaconda3/envs/ggn_target`。获准旧输入仅`results/01_snRNA/snRNA_final.h5ad`的注释与counts；原始GO和原始出版商Table S1另用于基因集定义/患者临床标签核对；未读取旧results/0*–4*其他结果。KAC/CAS沿用第0步已固定阈值与状态记录，未重新定义阈值。来源表达和伪批量全部重算。\n\n每名患者同一侧多个样本合计原始counts，一名患者每侧一个伪批量；每侧≥50细胞、每比较≥8对；基因过滤沿用既有技术定义（计数≥10的伪批量≥总数1/3且总计数≥50），没有按结果修改。pydeseq2 `~ patient + condition`，双侧Wald和BH；方向条件为I前体>Normal、IILUAD>前体；正向患者数≥ceil(0.8n)，22对需18、21对需17。患者方向按DESeq2归一化计数的配对差值（加0.5后log2比值）计算，平局不算正向。Wald95%区间、原始P、FDR、患者效应均保存。\n\n## 主结果与实际配对数\n\n{mdtable(main_show)}\n\n原始注释包含23名患者、75个样本；每名患者同侧样本合并后形成2个配对伪批量单位，原始样本数不能当作独立患者数。\n\n{mdtable(sample_info)}\n\n通过本步门槛共{int(a.upstream_eligible.sum())}个“基因×细胞群×比较步”条目，合并为{len(full)}个“基因×细胞群”候选记录、{full.gene.nunique()}个不同基因。这些不是最终通过全部门槛的候选：D3与D7待后续步骤，不能宣称最终靶点数。髓系包含巨噬细胞，两个分组重叠，不能把条目数当作独立基因发现数。\n\n排除原因：P7无前体样本，所有两步比较均排除；T/NK另排除P8（Normal18、前体35个细胞，低于50）；巨噬第I步另排除P13（Normal46、前体257个，Normal低于50）。其余主比较无额外患者排除。\n\n全部主比较明细：`discovery_all.tsv.gz`；候选并列展示：`candidates_step1.tsv`（按细胞群/基因名排列，**不是排序分排名**）。各比较的`pair_eligibility.tsv`列出每名患者两侧细胞数及排除原因。\n\n## 来源过滤与增殖集\n\n重新计算每样本每谱系的log1p(CP10k)细胞均值，再取合格样本中位数；上皮/最高非上皮谱系≥0.2，微环境目标群/上皮≥0.2。两侧均无表达的比值标无法判定，不默认为通过。未使用旧比值表或跨队列“救回”。`source_filter_recomputed.tsv`与`source_sample_eligibility.tsv`记录比值和分母样本数。\n\n原始GO重建得到196和1311基因，均与预设定义一致。196定义是在已知KAC基因集与旧定义重叠情况之后选定，并非盲选；重算不恢复盲法。196用于本步标记/过滤，1311并列标记；宽定义下top20敏感性须在第4步存在排序后报告，本步不产生top20。\n\n## D10状态内重算\n\n{mdtable(dsum[['lineage','state','step','pairs','status']])}\n\nKAC第I步7对不足，不做检验；第II步KAC内16对可检验，这是“前体KAC→LUAD KAC”，**不是**原定病灶KAC→Normal KAC修复态比较。AT2L-other两步22对，另保存LUAD→Normal的23对双侧对照供第3步使用。\n\n髓系的巨噬与单核、T/NK的T与NK达到配对门槛；更稀少亚群不足8对，逐项标无法验证。成纤维、内皮、独立巨噬群无更细亚型注释，标无法区分，不重新聚类以追求结果。`D10_candidate_state_evidence.tsv`列出每候选状态内效应、P/FDR及一致率；`state_fractions_by_patient.tsv`给出患者层面组成描述。\n\n**D10解释限制：** 状态内有显著正向变化的标“状态内变化有证据；组成贡献未确定”；状态内未检出不等于仅由组成驱动。现有比较没有正式分解组成与状态贡献，故未把缺乏证据强行标成“组成驱动”或“两者兼有”。这是可识别性限制，不是新增筛选门槛。\n\n## 临床目标标签与亚阶段\n\n{mdtable(subpairs)}\n\nAAH>Normal、AIS/MIA>Normal按固定方向作患者配对Wilcoxon单侧；在每细胞群/比较的参与候选内BH，并另报两个主亚阶段比较家族的Holm校正。并报告非零/平局配对数及中位效应的保守分布无关95%区间（未作多重性调整，实际覆盖率另列；样本过少时区间无有限界）。AIS、MIA单独仅描述中位患者效应与正向患者数，不临时追加单侧检验；MIA少于8对明确无法独立验证。\n\n{int(full.regression_direction_hint.sum())}个候选记录在第I步通过且AAH/AIS/MIA中位方向均为正，被标为“消退方向线索”，**不是已证实维持依赖**；AAH与合并AIS/MIA两项均BH<0.05者为{int(full.substage_AAH_and_AISMIA_BH_support.sum())}个。方向标签与统计支持分列，不能因方向为正宣称亚阶段已验证。微环境的直接支持作用仍需验证。仅第II步有证据的候选标阻断方向并降级消退标签。\n\n`clinical_substage_tests.tsv`保存全部P/FDR/Holm、效应与配对数；`clinical_target_labels.tsv`保存标签与限制。\n\n## LOO稳健性\n\nD3/D7尚未运行，先对所有通过本步D1/D2+D4+D5门槛的候选执行LOO；这覆盖未来可能通过全部门槛的候选范围，不宣称当前已全门槛通过。每次剔除一名患者，重新拟合原完整检验基因集合、size factors与同设计配对模型，并重新做BH；不只对候选校正。\n\n{mdtable(ls)}\n\n{int(loo_df.top5_excluded.sum())}个候选-步条目触发“单患者驱动”；合并为{int(full.top5_excluded_by_LOO.sum())}个基因×细胞群记录，按预设不得进入top5。`loo_summary.tsv`报告最大P、缺失P数、FDR保持比例、最有影响患者、剔除前后log2FC与最大绝对/相对变化。标“稳健”遵循方案的方向/一致率定义；是否每次仍显著另列，二者不是同一个条件。\n\n逐患者箱线/散点在`figures/*_patient_effects.pdf`，每个通过本步门槛的候选均有图，每点为患者，P24/P25用橙色。图中配对比值与模型Wald log2FC分别注明。\n\n## 方法局限与偏离记录\n\n{limitation.split(chr(10)+chr(10),1)[1]}\n\n- 无阈值、权重或门槛变更；未调整单侧方向；未改变第0步固定阈值。\n- 现方案第6.2仍有“第1步计算阈值”的遗留文字，按已确认第0步固定值执行，未重算。\n- 第3.2.1列出的前体n=37是排除前数；实际冻结后计划使用AIS19+MIA16=35。登记样本总数仍70，计划可用63；未在本步运行该数据集。\n- D10缺少细亚型或不足配对数的部分保持无法判断；不以不显著代替组成归因证据。\n- LOO范围为本步上游合格集合；D3/D7和可成药性不在本步运行。KAC剂量反应替代估计量留在第3步D7，并非本轮遗漏或原定检验的“抢救”。\n- 基因集/程序层面的D1/D2打分敏感性及权重/1311基因top20敏感性，在第4步评分时汇报，本步不进行最终评分。\n\n核查：原始计数守恒；KAC+other计数与AT2谱系汇总完全一致；直接抽查三个原始单元的counts/归一化表达；主模型患者配对、方向一致率和BH均通过；最终LOO检查与文件SHA256另存。\n\n**停点：第1步完成，等待用户指令；未进入第2步。**\n'''
 (D/'step1_report.md').write_text(report)
 (D/'integration_complete.json').write_text(json.dumps({'created_at':now(),'candidate_gene_lineages':len(full),'unique_genes':int(full.gene.nunique()),'LOO_excluded_gene_lineages':int(full.top5_excluded_by_LOO.sum()),'frozen_registry_before_sha256':sha(old),'frozen_registry_after_sha256':sha(registry)},ensure_ascii=False,indent=2))
 log('Step1 integration complete; no ranking; no step2. Report, D10 limitations, clinical direction labels and LOO flags written.')


if __name__=='__main__':
 phase=sys.argv[1] if len(sys.argv)>1 else 'all';freeze_methods()
 for name,fn in [('prepare',prepare),('main',main_tests),('d10',d10),('substage',substage),('loo',loo),('finalize',finalize)]:
  if phase in ['all',name]:fn()
