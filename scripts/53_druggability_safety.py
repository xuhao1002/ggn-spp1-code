#!/opt/anaconda3/envs/ggn_target/bin/python
"""Step3 only: D3 axis-matched pool, D6-D9 raw-source annotations and D7 replacement estimator.
No ranking, no legacy result reuse except approved snRNA annotations/counts, no frozen-set analysis.
"""
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
CACHE=Path('/Users/dudu/Documents/Codex/2026-10-09/ggn-target-ggn-target-opt-anaconda3/work')
os.environ['MPLCONFIGDIR']=str(CACHE/'mpl');os.environ['NUMBA_CACHE_DIR']=str(CACHE/'numba_cache');os.environ['XDG_CACHE_HOME']=str(CACHE/'cache')
import sys,json,gzip,zipfile,re,hashlib,time,warnings
import numpy as np,pandas as pd,scipy.sparse as sp
from scipy.stats import wilcoxon,binom
from statsmodels.stats.multitest import multipletests
B=Path('/Users/dudu/全部数据库/GGN_target');R=B/'results/50_screen';O=R/'step3';S=O/'sources';O.mkdir(exist_ok=True);S.mkdir(exist_ok=True)
def log(x):
 s=time.strftime('%Y-%m-%d %H:%M:%S')+' '+str(x);print(s,flush=True)
 with (O/'execution.log').open('a') as f:f.write(s+'\n')
def write(x,name,index=False):x.to_csv(O/name,sep='\t',index=index)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def pool():
 c=pd.read_csv(R/'step2/candidates_step2.tsv',sep='\t')
 # Old cross-axis scalar retained only in historical step2, not used in any current eligibility.
 c=c.rename(columns={'D3':'D3_cross_axis_step2_archival','D3_gate_pass':'D3_cross_axis_gate_archival'})
 c['I_axis_matched_D3_pass']=c.I_upstream_eligible.eq(True)&c.D3_I.ge(1)
 c['II_axis_matched_D3_pass']=c.II_upstream_eligible.eq(True)&c.D3_II_histology.ge(1)
 c['pool_before_D7']=c.I_axis_matched_D3_pass|c.II_axis_matched_D3_pass
 ev=pd.read_csv(R/'step2/validation_evidence_long.tsv',sep='\t');ev=ev[ev.support & ev.axis.str.startswith('I_')]
 es={g:t for g,t in ev.groupby('gene')};ds=[]
 for r in c.itertuples():
  t=es.get(r.gene,pd.DataFrame(columns=['lineage','dataset']));ss=set(t.loc[t.lineage.isin(['bulk',r.lineage]),'dataset'])
  ds.append(';'.join(sorted(ss)))
 c['I_support_datasets']=ds
 def tier(r):
  if not r.pool_before_D7:return '不在轴匹配候选池'
  d=set(r.I_support_datasets.split(';'))
  if 'GSE164789' in d:return '甲'
  if 'GSE319666' in d:return '乙'
  if len(d)>=2 and d<={'GSE102511','GSE288479','GSE319666'}:return '丙'
  return '丁'
 c['cohort_suitability_tier']=c.apply(tier,axis=1)
 c['top5_excluded_opposite']=c.has_opposite_significant.eq(True)
 c['top5_excluded_I_LOO']=c.I_LOO_status.eq('单患者驱动')
 c['top5_excluded_II_LOO']=c.II_LOO_status.eq('单患者驱动')
 c['top5_excluded_existing']=c.top5_excluded_opposite|c.top5_excluded_by_LOO.eq(True)
 p=c[c.pool_before_D7]
 assert len(c)==5207 and len(p)==199 and c.II_axis_matched_D3_pass.sum()==0
 assert p.cohort_suitability_tier.value_counts().to_dict()=={'甲':135,'丙':59,'乙':5}
 assert p.top5_excluded_opposite.sum()==5 and p.top5_excluded_I_LOO.sum()==6
 write(c,'candidates_axis_matched.tsv');write(p,'pool199_before_D7.tsv')
 log(f'Axis-matched pool {len(p)}, {p.gene.nunique()} genes; tiers {p.cohort_suitability_tier.value_counts().to_dict()}')

def safety():
 c=pd.read_csv(O/'candidates_axis_matched.tsv',sep='\t');genes=sorted(c.gene.unique())
 raw=B/'raw/safety_druggability'
 gt=pd.read_csv(raw/'GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz',sep='\t',skiprows=2)
 core=['Heart - Left Ventricle','Heart - Atrial Appendage','Liver','Kidney - Cortex','Kidney - Medulla'];tissues=core+['Whole Blood']
 ident=pd.read_csv(O/'gene_identity_map.tsv',sep='\t').fillna('')
 gt['stable_id']=gt.Name.str.split('.').str[0]
 gtrows=[]
 for g in genes:
  ids=ident.loc[ident.gene.eq(g),'ensembl_ids'].iloc[0].split(';');sub=gt[gt.Description.eq(g)|gt.stable_id.isin(ids)]
  gtrows.append(sub[tissues].max().rename(g))
 gt=pd.DataFrame(gtrows);gt.index.name='gene';gt.reindex(genes).to_csv(O/'GTEx_candidate_tissue_medians.tsv',sep='\t')
 z=zipfile.ZipFile(raw/'normal_ihc_data.tsv.zip');ih=pd.read_csv(z.open(z.namelist()[0]),sep='\t')
 imap=[]
 for rr in ident.itertuples():
  for ens in rr.ensembl_ids.split(';'):
   if ens:imap.append({'Gene':ens,'candidate_gene':rr.gene})
 ih=ih.merge(pd.DataFrame(imap),on='Gene',how='inner')
 ih=ih.rename(columns={'Gene name':'source_gene_symbol','candidate_gene':'Gene name'})
 write(ih,'HPA_candidate_records.tsv')
 ih=ih[ih.Reliability.isin(['Approved','Supported','Enhanced'])]
 tc={'Heart muscle':['cardiomyocytes'],'Liver':['hepatocytes'],'Kidney':['cells in glomeruli','cells in tubules','proximal tubules (cell body)','proximal tubules (microvilli)','distal tubules','collecting ducts'],'Bone marrow':['hematopoietic cells']}
 ih=ih[ih.apply(lambda r:r['Cell type'] in tc.get(r.Tissue,[]),axis=1)]
 hgroups={g:t for g,t in ih.groupby('Gene name')}
 # Fresh step1 aggregates are allowed; do not load old results/04 ratios.
 med=pd.read_csv(R/'step1/source_median_log_expression.tsv',sep='\t',index_col=0)
 u=pd.read_csv(R/'step1/annotation_units.tsv',sep='\t',index_col=0);a=np.load(R/'step1/aggregates.npz');gn=a['genes'];ix=pd.Index(gn).get_indexer(genes);assert (ix>=0).all()
 L=a['log_sums'][:,ix];nc=u.n_cells.values
 # Macrophage and Myeloid overlap: construct disjoint remainder instead of comparing target to a group containing itself.
 labs=['Monocyte','cDC','pDC','Neutrophil','Proliferating myeloid'];vals=[];elig=[]
 for sample in u['sample'].unique():
  ind=np.where(u['sample'].eq(sample)&u.final_label.isin(labs))[0];n=int(nc[ind].sum());elig.append(dict(sample=sample,n_cells=n,included=n>=50))
  if n>=50:vals.append(L[ind].sum(axis=0)/n)
 med['Other_myeloid']=pd.Series(np.median(np.stack(vals),axis=0),index=genes)
 write(pd.DataFrame(elig),'other_myeloid_sample_eligibility.tsv')
 med.reindex(genes).to_csv(O/'D7_immune_source_medians.tsv',sep='\t')
 rows=[]
 for r in c.itertuples():
  g,lin=r.gene,r.lineage;v=gt.loc[g] if g in gt.index else pd.Series(np.nan,index=tissues)
  h=hgroups.get(g,pd.DataFrame(columns=['Tissue','Cell type','Level','Reliability']))
  hc=h[h.Tissue.ne('Bone marrow')];hb=h[h.Tissue.eq('Bone marrow')]
  corehigh=hc[hc.Level.eq('High')];bmhigh=hb.Level.eq('High').any()
  reasons=[]
  for t in core:
   if v[t]>=50:reasons.append(f'GTEx {t} {v[t]:.6g} TPM>=50')
  if len(corehigh):reasons.append('HPA High: '+';'.join(corehigh.Tissue+' / '+corehigh['Cell type']))
  immune=lin in ['Myeloid','Macrophage','T_NK'];fold=np.nan;others=[]
  if immune:
   others=['T_NK','B_Plasma','Mast'] if lin=='Myeloid' else ['Myeloid','B_Plasma','Mast'] if lin=='T_NK' else ['Other_myeloid','T_NK','B_Plasma','Mast']
   num=float(med.loc[g,lin]);den=float(med.loc[g,others].max())
   fold=num/den if den>0 else np.inf if num>0 and den==0 else np.nan
   if np.isfinite(fold) and fold<3:reasons.append(f'免疫谱系选择性不足：{fold:.6g}<3')
   if bmhigh:reasons.append('HPA 骨髓 High')
  else:
   if v['Whole Blood']>=50:reasons.append(f'GTEx Whole Blood {v["Whole Blood"]:.6g} TPM>=50')
   if bmhigh:reasons.append('HPA 骨髓 High')
  unknown=bool(v.isna().any() or (immune and pd.isna(fold)))
  status='预设安全否决' if reasons else '无法判定：关键数值缺失' if unknown else '未触发预设安全否决'
  row=dict(gene=g,lineage=lin,D7_safety_veto=bool(reasons),D7_status=status,D7_data_incomplete=unknown,D7_reasons='；'.join(reasons),D7_immune_ratio=fold,D7_immune_comparator=';'.join(others),D7_hpa_core_high=bool(len(corehigh)),D7_hpa_bone_marrow_high=bool(bmhigh),D7_hpa_core_records=len(hc),D7_hpa_bone_marrow_records=len(hb),D7_hpa_missing_protein_evidence=h.empty,D7_safety_interpretation='未触发预设否决不等于长期用药安全；HPA缺失不等于低表达')
  row.update({'GTEx_'+t:v[t] for t in tissues});rows.append(row)
 out=pd.DataFrame(rows);write(out,'D7_organ_safety.tsv');p=c[c.pool_before_D7].drop(columns=[k for k in out.columns if k in c.columns and k not in ['gene','lineage']]).merge(out,on=['gene','lineage'])
 write(p,'pool199_organ_safety.tsv');log(f'D7 organs all {out.D7_status.value_counts().to_dict()}; pool {p.D7_status.value_counts().to_dict()}')


def depmap_lincs():
 c=pd.read_csv(O/'candidates_axis_matched.tsv',sep='\t');genes=sorted(c.gene.unique())
 model=pd.read_csv(B/'raw/DepMap/Model.csv');luad=model[model.OncotreeSubtype.eq('Lung Adenocarcinoma')].copy()
 header=pd.read_csv(B/'raw/DepMap/CRISPRGeneEffect.csv',nrows=0).columns
 rawmap={re.sub(r' \(\d+\)$','',h):h for h in header[1:]}
 ident=pd.read_csv(O/'gene_identity_map.tsv',sep='\t').fillna('').set_index('gene');mapping={};source_symbols={}
 for g in genes:
  matches=[g] if g in rawmap else sorted(set(ident.loc[g,'canonical_symbols'].split(';'))&set(rawmap))
  if len(matches)==1:mapping[rawmap[matches[0]]]=g;source_symbols[g]=matches[0]
 cols=list(mapping)
 arr=[];allarr=[]
 for ch in pd.read_csv(B/'raw/DepMap/CRISPRGeneEffect.csv',usecols=[header[0]]+cols,chunksize=100,index_col=0):
  sub=ch[ch.index.isin(luad.ModelID)];arr.append(sub);allarr.append(ch)
 a=pd.concat(arr);whole=pd.concat(allarr);a.columns=[mapping[h] for h in a.columns];whole.columns=a.columns
 # Duplicate symbols would need explicit identifier review, not silent averaging.
 assert a.columns.is_unique
 a.to_csv(O/'DepMap_LUAD_gene_effect_matrix.tsv.gz',sep='\t');write(luad[luad.ModelID.isin(a.index)],'DepMap_LUAD_models.tsv')
 rows=[]
 for g in genes:
  v=a[g].dropna() if g in a else pd.Series(dtype=float);w=whole[g].dropna() if g in whole else pd.Series(dtype=float)
  rows.append(dict(gene=g,D9_DepMap_source_symbol=source_symbols.get(g,''),D9_DepMap_n_LUAD=len(v),D9_DepMap_median=v.median(),D9_DepMap_mean=v.mean(),D9_DepMap_q25=v.quantile(.25),D9_DepMap_q75=v.quantile(.75),D9_DepMap_min=v.min(),D9_DepMap_max=v.max(),D9_DepMap_fraction_negative=(v<0).mean(),D9_DepMap_all_model_n=len(w),D9_DepMap_all_model_median=w.median(),D9_DepMap_status='细胞系敲除生长效应已测量；早期病灶干预方向未知' if len(v) else '未检索到本次LUAD Chronos测量',D9_DepMap_limit='细胞系存活依赖偏向增殖，与早期病灶不对口；Chronos负值不证明既有前驱病变消退；不入排序、不作门槛'))
 write(pd.DataFrame(rows),'D9_DepMap_gene_annotations.tsv')
 # Expression is annotation only: match model IDs, no inferential cell-line test.
 eh=pd.read_csv(B/'raw/DepMap/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv',nrows=0).columns
 erm={re.sub(r' \(\d+\)$','',h):h for h in eh if re.search(r' \(\d+\)$',h)};em={}
 for g in genes:
  matches=[g] if g in erm else sorted(set(ident.loc[g,'canonical_symbols'].split(';'))&set(erm))
  if len(matches)==1:em[erm[matches[0]]]=g
 ec=list(em)
 idcol='ModelID' if 'ModelID' in eh else eh[0]
 ee=[]
 for ch in pd.read_csv(B/'raw/DepMap/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv',usecols=[idcol]+ec,chunksize=100):ee.append(ch[ch[idcol].isin(luad.ModelID)])
 ex=pd.concat(ee).set_index(idcol);ex.columns=[em[h] for h in ex.columns]
 ex.groupby(level=0).mean().median(axis=0).rename('LUAD_median_supplied_log2_TPM_plus1').to_csv(O/'DepMap_LUAD_expression_summary.tsv',sep='\t')
 sig={}
 for line in (B/'raw/lincs/Enrichr_LINCS_L1000_CRISPR_KO_Consensus_Sigs.txt').read_text().splitlines():
  v=line.split('\t');m=re.match(r'^(.*?) (Up|Down)$',v[0])
  if m:sig.setdefault(m[1],{})[m[2]]=set(x for x in v[2:] if x)
 kac=set(pd.read_csv(R/'step1/candidates_step1.tsv',sep='\t').query('KAC_signature_member == True').gene)
 # Full fixed signature file rather than candidate-only subset.
 fixed=R/'state_thresholds_frozen.json';conf=json.loads(fixed.read_text())
 kpath=R/'KAC_genes_used.txt'
 if kpath.exists():kac=set(kpath.read_text().split())
 # Recover available 91 signature symbols from step0 JSON when present.
 for key in ['KAC_genes_used','kac_genes_used']:
  if key in conf:kac=set(conf[key])
 res=[]
 for r in c.itertuples():
  matches=[r.gene] if r.gene in sig else sorted(set(ident.loc[r.gene,'canonical_symbols'].split(';'))&set(sig))
  matched=matches[0] if len(matches)==1 else '';z=sig.get(matched,{});up=z.get('Up',set());dn=z.get('Down',set())
  prog=set(c.loc[c.lineage.eq(r.lineage)&c.I_upstream_eligible.eq(True),'gene'])
  res.append(dict(gene=r.gene,lineage=r.lineage,D9_LINCS_source_symbol=matched,D9_LINCS_KO_signature_present=bool(z),D9_LINCS_up_n=len(up),D9_LINCS_down_n=len(dn),D9_LINCS_up_genes=';'.join(sorted(up)),D9_LINCS_down_genes=';'.join(sorted(dn)),D9_LINCS_I_program_up_overlap=len(up&prog),D9_LINCS_I_program_down_overlap=len(dn&prog),D9_LINCS_limit='KO共识转录特征；集合Up/Down为KO后响应方向，不等于病灶效应；本资源缺逐模型/剂量/时点信息；不据此认定早期干预有益'))
 write(pd.DataFrame(res),'D9_LINCS_candidate_annotations.tsv');log(f'DepMap models {len(a)} LUAD; measured candidates {len(cols)}; LINCS signatures {len(sig)}')


def prepare_kac():
 import h5py
 from anndata.io import read_elem,sparse_dataset
 c=pd.read_csv(O/'candidates_axis_matched.tsv',sep='\t');genes=sorted(c.gene.unique())
 cells=pd.read_csv(R/'AT2L_cell_state_assignments.tsv.gz',sep='\t',index_col=0)
 cfg=json.loads((R/'state_thresholds_frozen.json').read_text());assert (cells.is_KAC==(cells.KAC>cfg['thresholds']['KAC'])).all()
 with h5py.File(B/'results/01_snRNA/snRNA_final.h5ad','r') as f:
  obs=read_elem(f['obs']);var=read_elem(f['var']);gidx=var.index.get_indexer(genes);assert (gidx>=0).all()
  ds=sparse_dataset(f['layers/counts']);blocks=[];names=[]
  for start in range(0,len(obs),10000):
   end=min(start+10000,len(obs));idx=obs.index[start:end];keep=idx.isin(cells.index)
   if not keep.any():continue
   raw=ds[start:end,:].tocsr()[keep,:];lib=np.asarray(raw.sum(axis=1)).ravel();assert (lib>0).all()
   x=sp.diags(1e4/lib)@raw[:,gidx].astype(np.float64);x.data=np.log1p(x.data)
   blocks.append(x.tocsr());names.extend(idx[keep]);
  xx=sp.vstack(blocks).tocsr();assert xx.shape==(len(cells),len(genes));assert set(names)==set(cells.index)
  sp.save_npz(O/'AT2L_candidate_log1pCP10k.npz',xx)
  cells=cells.loc[names].copy();cells['side']=np.where(cells.stage.eq('Normal'),'Normal',np.where(cells.stage.eq('LUAD'),'LUAD','Precursor'))
  cells['matrix_row']=np.arange(len(cells));cells.to_csv(O/'AT2L_fit_cells.tsv.gz',sep='\t');(O/'AT2L_fit_genes.json').write_text(json.dumps(genes))
  um=cells.groupby(['patient','side']).size().rename('n_cells').reset_index();um['eligible']=um.n_cells.ge(50);write(um,'KAC_slope_units_before_fit.tsv')
  log(f'KAC replacement prepared {xx.shape}; units {um.groupby("side").eligible.sum().to_dict()}')

def huber_many(x,Y,maxiter=100,tol=1e-8):
 """Vectorized HuberT IRLS + MAD, intercept; no cell-level inferential P."""
 n,g=Y.shape;xm=x.mean();xs=x.std();z=(x-xm)/xs
 b0=Y.mean(axis=0);b1=z@Y/n
 constant=np.ptp(Y,axis=0)==0;converged=constant.copy();invalid=np.zeros(g,bool);iterations=np.zeros(g,int)
 scale=np.zeros(g);eps=np.finfo(float).eps*100*np.maximum(1,np.max(np.abs(Y),axis=0))
 active=~constant
 for it in range(maxiter):
  ids=np.where(active)[0]
  if not len(ids):break
  yy=Y[:,ids];res=yy-b0[ids]-z[:,None]*b1[ids]
  ss=np.median(np.abs(res),axis=0)/.6744897501960817;scale[ids]=ss
  bad=(ss<=eps[ids]);perfect=np.max(np.abs(res),axis=0)<=eps[ids]
  invalid[ids[bad&~perfect]]=True;converged[ids[bad&perfect]]=True;active[ids[bad]]=False
  good=~bad
  if not good.any():continue
  ids=ids[good];yy=yy[:,good];res=res[:,good];ss=ss[good]
  w=np.minimum(1.,1.345*ss[None,:]/np.maximum(np.abs(res),np.finfo(float).tiny))
  sw=w.sum(axis=0);sx=z@w;sxx=(z*z)@w;sy=(w*yy).sum(axis=0);sxy=(w*yy*z[:,None]).sum(axis=0)
  det=sw*sxx-sx*sx;bad=det<=np.finfo(float).eps
  nb0=np.divide(sy*sxx-sx*sxy,det,out=np.full(len(ids),np.nan),where=~bad)
  nb1=np.divide(sw*sxy-sx*sy,det,out=np.full(len(ids),np.nan),where=~bad)
  done=(np.maximum(np.abs(nb0-b0[ids]),np.abs(nb1-b1[ids]))<=tol*(1+np.maximum(np.abs(b0[ids]),np.abs(b1[ids]))))&~bad
  b0[ids]=nb0;b1[ids]=nb1;converged[ids[done]]=True;invalid[ids[bad]]=True;active[ids[done|bad]]=False;iterations[ids]=it+1
 slope=b1/xs;valid=converged&~invalid
 return slope,valid,iterations,scale,constant

_ROB_X=None;_ROB_C=None;_ROB_G=None
def fit_init():
 global _ROB_X,_ROB_C,_ROB_G
 _ROB_X=sp.load_npz(O/'AT2L_candidate_log1pCP10k.npz');_ROB_C=pd.read_csv(O/'AT2L_fit_cells.tsv.gz',sep='\t',index_col=0);_ROB_G=json.loads((O/'AT2L_fit_genes.json').read_text())
def fit_job(job):
 patient,side,ix=job;d=O/'slopes';d.mkdir(exist_ok=True);fp=d/f'{patient}_{side}.tsv'
 if fp.exists():return patient,side,'cached'
 x=_ROB_C.iloc[ix].KAC.values;X=_ROB_X[ix];rows=[]
 for start in range(0,len(_ROB_G),128):
  end=min(start+128,len(_ROB_G));y=X[:,start:end].toarray();sl,valid,nit,scale,const=huber_many(x,y)
  for j,g in enumerate(_ROB_G[start:end]):
   rows.append(dict(gene=g,patient=patient,side=side,n_cells=len(ix),slope_raw=sl[j] if valid[j] else np.nan,slope=np.round(sl[j],10) if valid[j] else np.nan,fit_valid=bool(valid[j]),constant_response=bool(const[j]),iterations=int(nit[j]),residual_scale=scale[j],mean_expression=y[:,j].mean(),fraction_detected=(y[:,j]>0).mean()))
 pd.DataFrame(rows).to_csv(fp,sep='\t',index=False);return patient,side,int(sum(r['fit_valid'] for r in rows))

def fit_kac():
 from concurrent.futures import ProcessPoolExecutor,as_completed
 import multiprocessing
 cells=pd.read_csv(O/'AT2L_fit_cells.tsv.gz',sep='\t',index_col=0);jobs=[]
 for (p,s),ii in cells.groupby(['patient','side']).indices.items():
  if len(ii)>=50:jobs.append((p,s,ii))
 with ProcessPoolExecutor(max_workers=3,mp_context=multiprocessing.get_context('spawn'),initializer=fit_init) as pool:
  for f in as_completed([pool.submit(fit_job,j) for j in jobs]):log('KAC robust fit '+str(f.result()))


def identities():
 c=pd.read_csv(O/'candidates_axis_matched.tsv',sep='\t');genes=sorted(c.gene.unique())
 z=zipfile.ZipFile(S/'proteinatlas.tsv.zip');h=pd.read_csv(z.open(z.namelist()[0]),sep='\t').fillna('')
 gt=pd.read_csv(B/'raw/safety_druggability/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz',sep='\t',skiprows=2,usecols=['Name','Description'])
 gt['ens']=gt.Name.str.split('.').str[0];gtm=gt.groupby('Description').ens.apply(set).to_dict();gm={}
 with gzip.open(B/'raw/annotation/gencode.v44.basic.annotation.gtf.gz','rt') as fp:
  for l in fp:
   if l.startswith('#'):continue
   v=l.rstrip().split('\t')
   if v[2]=='gene':gm.setdefault(re.search(r'gene_name "([^"]+)',v[8])[1],set()).add(re.search(r'gene_id "([^".]+)',v[8])[1])
 alias={}
 for r in h.itertuples(index=False,name=None):
  d=dict(zip(h.columns,r))
  for a in [d['Gene']]+[t.strip() for t in d['Gene synonym'].split(',') if t.strip()]:alias.setdefault(a,set()).add(d['Ensembl'])
 rows=[]
 for g in genes:
  ids=set(gtm.get(g,set()))|set(gm.get(g,set()));src=[]
  if g in gtm:src.append('GTEx exact symbol to Ensembl')
  if g in gm:src.append('GENCODE44 exact symbol to Ensembl')
  if not ids and len(alias.get(g,set()))==1:ids=alias[g];src.append('HPA unique symbol/synonym')
  hit=h[h.Ensembl.isin(ids)]
  rows.append(dict(gene=g,ensembl_ids=';'.join(sorted(ids)),canonical_symbols=';'.join(sorted(set(hit.Gene))),uniprot_accessions=';'.join(sorted({u.strip() for t in hit.Uniprot for u in t.split(',') if u.strip()})),identity_sources=';'.join(src),identity_ambiguous=len(ids)>1))
 write(pd.DataFrame(rows),'gene_identity_map.tsv');log('Gene identities resolved; old symbols mapped by stable Ensembl IDs.')


def adjust(p,method='fdr_bh'):
 p=np.asarray(p,float);q=np.full(len(p),np.nan);ok=np.isfinite(p)
 if ok.any():q[ok]=multipletests(p[ok],method=method)[1]
 return q

def median_ci(x):
 x=np.sort(x);n=len(x);k=[i for i in range(1,n//2+1) if 1-2*binom.cdf(i-1,n,.5)>=.95]
 return (x[max(k)-1],x[n-max(k)]) if k else (-np.inf,np.inf)

def kac_tests():
 ss=pd.concat([pd.read_csv(f,sep='\t') for f in sorted((O/'slopes').glob('*.tsv'))],ignore_index=True)
 assert not ss[['gene','patient','side']].duplicated().any()
 write(ss,'KAC_patient_slopes_all.tsv.gz');rows=[];deltas=[]
 for g,t in ss.groupby('gene'):
  wide=t.pivot(index='patient',columns='side',values='slope');wide_raw=t.pivot(index='patient',columns='side',values='slope_raw')
  for side in ['Precursor','LUAD']:
   dd=(wide[side]-wide.Normal).dropna();rr=(wide_raw[side]-wide_raw.Normal).reindex(dd.index);n=len(dd)
   p=1. if n>=8 and (dd==0).all() else wilcoxon(dd.values,alternative='greater').pvalue if n>=8 else np.nan
   praw=1. if n>=8 and (rr==0).all() else wilcoxon(rr.values,alternative='greater').pvalue if n>=8 else np.nan
   lo,hi=median_ci(dd.values)
   rows.append(dict(gene=g,comparison=side+'_vs_Normal',n_valid_pairs=n,median_slope_difference=dd.median(),mean_slope_difference=dd.mean(),positive_patients=int((dd>0).sum()),negative_patients=int((dd<0).sum()),zero_patients=int((dd==0).sum()),p_greater=p,p_greater_unrounded_numeric_sensitivity=praw,median_ci95_low=lo,median_ci95_high=hi,normal_median_slope=wide.loc[dd.index,'Normal'].median(),lesion_median_slope=wide.loc[dd.index,side].median(),status='已检验（替代估计量）' if n>=8 else '无法检验：有效患者对不足8',interpretation='病灶>正常斜率为替代选择性关联，不等于原定KAC子集检验；不显著不等于修复态共有'))
   for patient,v in dd.items():deltas.append(dict(gene=g,comparison=side+'_vs_Normal',patient=patient,slope_difference=v,unrounded_difference=rr.loc[patient]))
 out=pd.DataFrame(rows);out['q_BH_within_comparison']=out.groupby('comparison').p_greater.transform(lambda x:adjust(x))
 out['p_Holm_all_genes_two_comparisons']=adjust(out.p_greater,'holm')
 out['BH_direction_supported']=out.q_BH_within_comparison.lt(.05)&out.median_slope_difference.gt(0)
 out['Holm_direction_supported']=out.p_Holm_all_genes_two_comparisons.lt(.05)&out.median_slope_difference.gt(0)
 out['numerically_tiny_median_effect']=out.median_slope_difference.abs().lt(1e-8)
 write(out,'D7_KAC_replacement_tests.tsv');write(pd.DataFrame(deltas),'KAC_patient_slope_differences.tsv.gz')
 qa=ss.groupby(['patient','side']).agg(n_cells=('n_cells','first'),genes=('gene','size'),valid_fits=('fit_valid','sum'),constant_responses=('constant_response','sum'),max_iterations=('iterations','max')).reset_index();write(qa,'KAC_fit_quality_by_patient.tsv')
 log('KAC replacement tests '+str(out.groupby('comparison').agg(tested=('p_greater','count'),BH=('BH_direction_supported','sum'),Holm=('Holm_direction_supported','sum')).to_dict()))

def other_controls():
 from pydeseq2.dds import DeseqDataSet
 from pydeseq2.ds import DeseqStats
 z=np.load(R/'step1/aggregates.npz');su=pd.read_csv(R/'step1/state_units.tsv',sep='\t');genes=z['genes'];C=z['state_counts'];N=z['state_n']
 cg=sorted(pd.read_csv(O/'candidates_axis_matched.tsv',sep='\t').gene.unique());results=[]
 for tar,stages in [('Precursor',['AAH','AIS','MIA']),('LUAD',['LUAD'])]:
  rows=[];meta=[]
  for p in sorted(su.patient.unique()):
   indices={side:np.where(su.patient.eq(p)&su.state.eq('other')&su.stage.isin(st))[0] for side,st in [('reference',['Normal']),('target',stages)]}
   if not all(N[ii].sum()>=50 for ii in indices.values()):continue
   for side,ii in indices.items():rows.append(C[ii].sum(axis=0));meta.append(dict(unit=p+'|'+side,patient=p,condition=side,n_cells=int(N[ii].sum())))
  m=pd.DataFrame(meta).set_index('unit');cnt=pd.DataFrame(np.vstack(rows),index=m.index,columns=genes);npairs=len(m)//2
  assert npairs==({'Precursor':22,'LUAD':23}[tar]);m.to_csv(O/f'AT2L_other_{tar}_units.tsv',sep='\t')
  keep=(cnt.ge(10).sum(axis=0)>=np.ceil(len(cnt)/3))&cnt.sum(axis=0).ge(50);cc=cnt.loc[:,keep]
  mm=m[['patient','condition']].copy();mm.condition=pd.Categorical(mm.condition,categories=['reference','target'])
  dds=DeseqDataSet(counts=cc.astype(int),metadata=mm,design='~ patient + condition',refit_cooks=True,n_cpus=4,quiet=True);dds.deseq2()
  st=DeseqStats(dds,contrast=['condition','target','reference'],n_cpus=4,quiet=True);st.summary();r=st.results_df;r.to_csv(O/f'AT2L_other_{tar}_all_genes.tsv.gz',sep='\t')
  r=r.reindex(cg).copy();r['gene']=r.index;r['comparison']=tar+'_vs_Normal';r['n_pairs']=npairs;r['effect_ci95_low']=r.log2FoldChange-1.96*r.lfcSE;r['effect_ci95_high']=r.log2FoldChange+1.96*r.lfcSE;r['status']=np.where(r.pvalue.notna(),'AT2L-other双侧对照已检验','计数过滤/模型无检测值');results.append(r.reset_index(drop=True));log(f'AT2L-other {tar} control {npairs} pairs complete')
 write(pd.concat(results),'D7_AT2L_other_controls.tsv')


def druggability():
 c=pd.read_csv(O/'candidates_axis_matched.tsv',sep='\t');genes=set(c.gene)
 ident=pd.read_csv(O/'gene_identity_map.tsv',sep='\t').fillna('');sym_to_candidate={};uni_to_candidate={}
 for r in ident.itertuples():
  for a in {r.gene}|set(r.canonical_symbols.split(';')):
   if a:sym_to_candidate.setdefault(a,set()).add(r.gene)
  for a in r.uniprot_accessions.split(';'):
   if a:uni_to_candidate.setdefault(a,set()).add(r.gene)
 mech=sum([json.loads(f.read_text())['mechanisms'] for f in S.glob('chembl_mechanism_page*.json')],[])
 tar=sum([json.loads(f.read_text())['targets'] for f in S.glob('chembl_target_batch*.json')],[])
 mol=sum([json.loads(f.read_text())['molecules'] for f in S.glob('chembl_molecule_batch*.json')],[])
 assert len(mech)==json.loads((S/'chembl_mechanism_page0.json').read_text())['page_meta']['total_count']
 tm={t['target_chembl_id']:t for t in tar};mm={m['molecule_chembl_id']:m for m in mol}
 assert {x['target_chembl_id'] for x in mech if x['target_chembl_id']}<=set(tm)
 assert {x['molecule_chembl_id'] for x in mech}<=set(mm)
 inhibitory={'INHIBITOR','ANTAGONIST','BLOCKER','NEGATIVE ALLOSTERIC MODULATOR','NEGATIVE MODULATOR','INVERSE AGONIST','ANTISENSE INHIBITOR','RNAI INHIBITOR','DEGRADER','GENE EDITING NEGATIVE MODULATOR','ALLOSTERIC ANTAGONIST'}
 ambiguous={'','BINDING AGENT','MODULATOR','OTHER','DISRUPTING AGENT','SEQUESTERING AGENT','CROSS-LINKING AGENT'}
 rows=[]
 for m in mech:
  t=tm.get(m['target_chembl_id']);drug=mm[m['molecule_chembl_id']]
  if t is None or t.get('organism')!='Homo sapiens':continue
  hit=set();how=[]
  for co in t.get('target_components',[]):
   u=co.get('accession','')
   if u in uni_to_candidate:hit|=uni_to_candidate[u];how.append('UniProt '+u)
   for sy in co.get('target_component_synonyms',[]):
    if sy['syn_type']=='GENE_SYMBOL' and sy['component_synonym'] in sym_to_candidate:hit|=sym_to_candidate[sy['component_synonym']];how.append('gene symbol '+sy['component_synonym'])
  act=m.get('action_type') or ''
  for g in hit:
   rows.append(dict(gene=g,source='ChEMBL37 direct mechanism API',record_id=m['mec_id'],source_version='ChEMBL_37 (2026-05-01)',drug_id=m['molecule_chembl_id'],drug_name=drug.get('pref_name') or m['molecule_chembl_id'],drug_type=drug.get('molecule_type'),target_id=m['target_chembl_id'],target_type=t.get('target_type'),target_name=t.get('pref_name'),action_type=act,explicit_inhibitory=act in inhibitory,type_pending=act in ambiguous,clinical_max_phase=m.get('max_phase'),compound_max_phase=drug.get('max_phase'),withdrawn_flag=drug.get('withdrawn_flag'),direct_interaction=m.get('direct_interaction'),mechanism=m.get('mechanism_of_action'),mapping=';'.join(sorted(set(how))),mechanism_url='https://www.ebi.ac.uk/chembl/explore/compound/'+m['molecule_chembl_id'],primary_refs=';'.join(r.get('ref_url','') or '' for r in m.get('mechanism_refs',[])),interpretation='复合物/家族组分关联，不能推断单基因选择性抑制' if t.get('target_type')!='SINGLE PROTEIN' else '单蛋白机制记录；不是早期肺病灶疗效证据'))
 dg=pd.read_csv(B/'raw/safety_druggability/dgidb_interactions.tsv',sep='\t').fillna('')
 # TTD is retained with provenance as a secondary transfer, not falsely called a freshly verified direct TTD record.
 inh={'inhibitor','antagonist','blocker','negative modulator','inverse agonist','suppressor','antisense oligonucleotide'}
 unk={'','modulator','binder','other/unknown','antibody'}
 for r in dg[dg.interaction_source_db_name.eq('TTD')].itertuples():
  hits=sym_to_candidate.get(r.gene_name,set());act=str(r.interaction_type).lower()
  for g in hits:
   rows.append(dict(gene=g,source='TTD via local DGIdb raw snapshot',record_id='',source_version=r.interaction_source_db_version,drug_id=r.drug_concept_id,drug_name=r.drug_name or r.drug_claim_name,drug_type='待核',target_id='',target_type='待核',target_name=r.gene_name,action_type=act,explicit_inhibitory=act in inh,type_pending=act in unk,clinical_max_phase=np.nan,compound_max_phase=np.nan,withdrawn_flag=np.nan,direct_interaction=np.nan,mechanism='TTD转录相互作用；原始TTD逐条核验待完成',mapping='DGIdb gene identity',mechanism_url='https://dgidb.org/genes/'+r.gene_name,primary_refs='',approved_in_snapshot=str(r.approved).lower()=='true',interpretation='approved布尔值不能替代临床最高阶段；抗体/binder不自动解释为抑制'))
 ev=pd.DataFrame(rows).drop_duplicates();write(ev,'D6_mechanism_records.tsv')
 out=[]
 for g in sorted(genes):
  a=ev[ev.gene.eq(g)];h=a[a.explicit_inhibitory];direct=h[h.source.str.startswith('ChEMBL')];single=direct[direct.target_type.eq('SINGLE PROTEIN')&direct.direct_interaction.eq(1)]
  phase=pd.to_numeric(h.clinical_max_phase,errors='coerce');pending=a[a.type_pending]
  out.append(dict(gene=g,D6_inhibitory_record_present=not h.empty,D6_ChEMBL_inhibitory_present=not direct.empty,D6_ChEMBL_direct_single_protein=not single.empty,D6_TTD_transferred_inhibitory_present=bool((h.source.str.startswith('TTD')).any()),D6_clinical_highest_phase=phase.max(),D6_phase_missing_with_record=bool(not h.empty and phase.isna().all()),D6_type_pending=not pending.empty,D6_inhibitory_drug_names=';'.join(sorted(h.drug_name.astype(str).unique())),D6_pending_drug_names=';'.join(sorted(pending.drug_name.astype(str).unique())),D6_withdrawn_record_present=bool(h.withdrawn_flag.eq(True).any()),D6_status='有抑制类机制记录' if not h.empty else '仅有类型待核记录' if not pending.empty else '仅检出其他作用方向' if not a.empty else '所检来源未检出记录',D6_limit='药物最高阶段非本适应证阶段；机制存在不证明抑制候选对既有前驱病变有益；未检出不等于不可成药'))
 write(pd.DataFrame(out),'D6_gene_annotations.tsv');log(f'D6 {len(ev)} source records; {sum(r["D6_inhibitory_record_present"] for r in out)} genes with inhibitory records')


def feasibility():
 c=pd.read_csv(O/'candidates_axis_matched.tsv',sep='\t');ids=pd.read_csv(O/'gene_identity_map.tsv',sep='\t').fillna('').set_index('gene')
 z=zipfile.ZipFile(S/'proteinatlas.tsv.zip');hp=pd.read_csv(z.open(z.namelist()[0]),sep='\t').fillna('')
 hom=pd.read_csv(S/'ensembl_mouse_homology.tsv',sep='\t').fillna('');hom.columns=['human_ens','human_symbol','mouse_ens','mouse_symbol','type','confidence']
 drug=pd.read_csv(O/'D6_gene_annotations.tsv',sep='\t').set_index('gene');dep=pd.read_csv(O/'D9_DepMap_gene_annotations.tsv',sep='\t').set_index('gene')
 cells=pd.read_csv(O/'AT2L_fit_cells.tsv.gz',sep='\t',index_col=0);xx=sp.load_npz(O/'AT2L_candidate_log1pCP10k.npz');gn=json.loads((O/'AT2L_fit_genes.json').read_text())
 means=[];state_counts=[]
 for state,mask in [('AT2L',np.ones(len(cells),bool)),('KAC',cells.is_KAC.values)]:
  for (p,side),ii in cells.groupby(['patient','side']).indices.items():
   jj=ii[mask[ii]];state_counts.append(dict(state=state,patient=p,side=side,n_cells=len(jj),included=len(jj)>=50))
   if len(jj)>=50:
    vals=np.asarray(xx[jj].mean(axis=0)).ravel();means.append(pd.DataFrame({'gene':gn,'patient':p,'side':side,'state':state,'mean_log_expression':vals}))
 ms=pd.concat(means,ignore_index=True);write(pd.DataFrame(state_counts),'D8_AT2_KAC_expression_units.tsv')
 sm=ms.groupby(['gene','side','state']).agg(median_patient_mean_log_expression=('mean_log_expression','median'),n_patients=('patient','size')).reset_index();write(sm,'D8_AT2_KAC_expression.tsv')
 ex=sm.pivot(index='gene',columns=['side','state'],values='median_patient_mean_log_expression')
 rows=[]
 for r in c.itertuples():
  g=r.gene;ii=set(ids.loc[g,'ensembl_ids'].split(';'));h=hp[hp.Ensembl.isin(ii)]
  antibodies=sorted({a.strip() for v in h.Antibody for a in v.split(',') if a.strip()});ih=set(h['Reliability (IH)']);reliable=bool(ih&{'Approved','Supported','Enhanced'})
  mo=hom[hom.human_ens.isin(ii)&hom.mouse_ens.ne('')];d=drug.loc[g];modeln=int(dep.loc[g,'D9_DepMap_n_LUAD']);imm=r.lineage in ['Myeloid','Macrophage','T_NK']
  at=ex.loc[g,('Precursor','AT2L')] if ('Precursor','AT2L') in ex else np.nan;ka=ex.loc[g,('Precursor','KAC')] if ('Precursor','KAC') in ex else np.nan
  row=dict(gene=g,lineage=r.lineage,D8_IHC_antibody_ids=';'.join(antibodies),D8_IHC_gene_reliability=';'.join(sorted(ih-{''})),D8_IHC_resource_record=bool(antibodies),D8_IHC_feasibility='HPA有抗体及基因层面IHC支持；具体clone、商品在售与mIF兼容性待核' if antibodies and reliable else 'HPA有抗体记录但IHC支持待核' if antibodies else '所检HPA无抗体记录；不等于市场无抗体',D8_cell_line_KO_models=modeln,D8_cell_line_inhibitory_record=bool(d.D6_inhibitory_record_present),D8_cell_line_feasibility=('应在免疫/间质模型或共培养扰动；LUAD肿瘤细胞系不代表来源群' if r.lineage!='AT2L' else 'LUAD模型有CRISPR测量，可作为功能路线线索；早期适配待核')+'；具体sgRNA/shRNA序列、独立试剂和恢复实验须验证',D8_organoid_AT2L_precursor_expression=at,D8_organoid_KAC_precursor_expression=ka,D8_organoid_inhibitory_record=bool(d.D6_inhibitory_record_present),D8_organoid_feasibility='需来源细胞与肺泡类器官共培养，不能用上皮单培养替代' if r.lineage!='AT2L' else '有AT2/KAC表达量可查；遗传扰动条件可行，培养中表达、递送及工具选择性待核',D8_mouse_ortholog_record=not mo.empty,D8_mouse_symbols=';'.join(sorted(set(mo.mouse_symbol))),D8_mouse_ensembl_ids=';'.join(sorted(set(mo.mouse_ens))),D8_mouse_orthology_types=';'.join(sorted(set(mo.type))),D8_mouse_orthology_confidence=';'.join(sorted(set(mo.confidence.astype(str)))),D8_mouse_knockout_resource='原始论文存在Trem2缺失模型；品系获取与条件敲除适配待核' if g=='TREM2' else '原始论文存在Spp1缺失模型；品系获取与条件敲除适配待核' if g=='SPP1' else '现成敲除品系未逐条核实',D8_mouse_feasibility='有小鼠同源基因；药物交叉反应、品系与建立病灶后干预设计待核' if not mo.empty else '所检Ensembl未匹配小鼠同源记录；需核查同源性/人源化路线',D8_mouse_resource_source='https://pmc.ncbi.nlm.nih.gov/articles/PMC11088947/' if g=='TREM2' else 'https://pubmed.ncbi.nlm.nih.gov/31740923/' if g=='SPP1' else 'https://www.ensembl.org/biomart/martview',D8_is_hard_gate=False,D8_assessment_type='资源记录+条件可行性推断；不是试剂/模型验证完成')
  rows.append(row)
 write(pd.DataFrame(rows),'D8_four_path_annotations.tsv');log('D8 four paths annotated for all '+str(len(rows))+' records')


def literature():
 # Screened corpus is an index, not an exhaustive per-gene systematic review.
 papers={}
 for f in S.glob('literature_targeted_page*.json'):
  for p in json.loads(f.read_text()).get('resultList',{}).get('result',[]):
   if p.get('pmid'):papers[p['pmid']]=p
 # Tuple: gene, PMID, intervention, finding, direction, model. Curated from primary abstracts/full text.
 specs=[
 ('SPP1','31740923','Spp1缺失/上皮表达操纵','缺失减少化学致癌肺腺癌数量和大小；支持肿瘤发生作用，未证明建立前驱病灶后消退','抑制','化学致癌及KRAS驱动肺肿瘤'),
 ('TREM2','37081148','Trem2缺失、阻断抗体','缺失恢复NK细胞并抑制肺肿瘤；抗体于肿瘤细胞接种后14天开始；非AAH/AIS自然病灶','抑制','肺肿瘤细胞接种小鼠；免疫微环境'),
 ('CEACAM6','31474761','pHLIP递送siRNA','沉默抑制A549及异种移植瘤生长','抑制','A549/异种移植瘤'),
 ('CEACAM6','40263546','敲低/再表达','敲低可促EGFR突变细胞凋亡；TKI耐药细胞再表达可恢复药物敏感性；方向取决于药物背景','情境依赖','EGFR突变及TKI耐药细胞'),
 ('BCAT1','PMC9483610','shRNA/过表达','BCAT1敲低减少LUAD细胞增殖与克隆，过表达增强相关表型；下游NF-kB药物不等于BCAT1特异工具','抑制','A549/移植瘤'),
 ('ABHD17C','40582697','敲低','降低肿瘤细胞生长并增加凋亡','抑制','LUAD肿瘤细胞'),
 ('CHMP4C','40050481','敲低','降低LUAD细胞增殖、迁移和侵袭；分子对接药物不是已验证抑制剂','抑制','LUAD细胞'),
 ('KDELR3','40411680','敲低/恢复','促进ER应激凋亡并抑制侵袭；动物中抑制生长/转移','抑制','LUAD细胞/动物肿瘤'),
 ('NUAK2','39068449','敲低及ferrostatin救援','减少A549增殖/干性/移植瘤生长，部分效应可被铁死亡抑制剂逆转','抑制','A549/皮下移植瘤'),
 ('PAFAH1B3','35534807','沉默','抑制增殖、侵袭及移植瘤生长','抑制','LUAD细胞/异种移植瘤'),
 ('PYCR1','32133692','敲低','降低LUAD细胞增殖、迁移和侵袭；不以摘要的下游机制作为确定结论','抑制','LUAD细胞'),
 ('ETV4','31253395','敲低及MSI2恢复','敲低抑制增殖和侵袭，MSI2异位表达可救援','抑制','LUAD细胞'),
 ('XPR1','39928150','敲低/过表达救援','直接敲低降低LUAD细胞增殖侵袭；动物分支主要操纵RBM15，不误写为直接XPR1动物敲除','抑制','A549/PC-9细胞'),
 ('SRD5A3','42401760','敲低','降低A549/H1299恶性表型并促进铁死亡；不等同于现有5α还原酶药物具有同一选择性','抑制','A549/H1299细胞'),
 ('PLA2G10','41271053','抑制剂联合anti-PD-1','联合治疗改善肺癌模型，机制指向CD4 T细胞；摘要不足以判定药物单药效应/特异性','抑制（联合治疗）','肺癌免疫微环境模型；摘要级'),
 ('HOXD1','39450540','过表达','增强HOXD1抑制LUAD细胞及小鼠肿瘤；不支持见到病灶上调即抑制','增强','LUAD细胞/小鼠肿瘤'),
 ('SPRED1','36629984','过表达','增强SPRED1抑制LUAD细胞恶性表型及体内生长','增强','LUAD细胞/小鼠肿瘤'),
 ('PHLDA3','34006890','敲低/过表达','报道PHLDA3促进Wnt相关增殖侵袭；敲低后效应可被GSK3β抑制剂恢复','抑制','LUAD细胞'),
 ('PHLDA3','41886456','UFMylation及突变操纵','UFL1-PHLDA3轴抑制AKT和肿瘤生长；与2021促癌结论存在模型/修饰背景差异','增强/保留功能','LUAD细胞及患者来源移植瘤'),
 ('ZMAT3','33082333','Zmat3胚系缺失','KrasG12D肺腺癌模型中缺失对发生速度/严重度影响小；无益处不能写成抑制有效','未证明有益','KrasG12D肺腺癌发生模型')]
 fallback={'31740923':('Osteopontin drives KRAS-mutant lung adenocarcinoma.',''), 'PMC9483610':('BCAT1 promotes lung adenocarcinoma progression through enhanced mitochondrial function and NF-κB pathway activation','10.1631/jzus.B2100985')}
 rows=[]
 for g,k,iv,result,direction,model in specs:
  p=papers.get(k,{});title=p.get('title',fallback.get(k,('',))[0]);doi=p.get('doi',fallback.get(k,('',''))[1]);url='https://pmc.ncbi.nlm.nih.gov/articles/'+k+'/' if k.startswith('PMC') else 'https://pubmed.ncbi.nlm.nih.gov/'+k+'/'
  rows.append(dict(gene=g,source_id=k,title=title,doi=doi,url=url,intervention=iv,paraphrased_result=result,direction=direction,model=model,verified_level='原始全文及方法' if g in ['TREM2','BCAT1'] else '原始论文摘要',established_precursor_regression_proven=False,limit='原始报告未在本项目复现；不等于既有AAH/AIS消退或前体至浸润阻断证据'))
 ev=pd.DataFrame(rows);write(ev,'D9_literature_curated_evidence.tsv')
 lead=pd.read_csv(O/'D9_literature_screening_leads.tsv',sep='\t');lc=lead.groupby('gene').size();genes=pd.read_csv(O/'gene_identity_map.tsv',sep='\t').gene;out=[]
 for g in genes:
  a=ev[ev.gene.eq(g)];dirs=set(a.direction)
  state='方向有证据（仅所列模型）' if dirs and dirs!={'未证明有益'} else '方向未知'
  if g in ['CEACAM6','PHLDA3']:decision='情境依赖/方向冲突，不能统一推荐抑制'
  elif '增强' in dirs:decision='增强有抑瘤证据；不支持直接推定应抑制'
  elif dirs and dirs!={'未证明有益'}:decision='抑制有模型证据；早期既有病灶方向未确定'
  elif dirs:decision='所核扰动未证明有益；不能据此宣布基因无作用'
  else:decision='早期病灶方向未知；未逐篇核实的线索不能记为方向证据'
  out.append(dict(gene=g,D9_direction_class=state,D9_direction_annotation=decision,D9_curated_paper_n=len(a),D9_curated_directions=';'.join(sorted(dirs)),D9_literature_urls=';'.join(a.url),D9_literature_lead_n=int(lc.get(g,0)),D9_model_relevance='包含肺肿瘤发生模型；不是建立病灶后干预' if g in ['SPP1','ZMAT3'] else '仅肿瘤细胞系/移植或接种瘤模型（方案降级：仅晚期模型）' if len(a) else '未知/尚未完成逐篇方向核实',D9_established_precursor_benefit='未获本次核实证据',D9_literature_coverage='全候选标题摘要符号匹配；少量原始论文方向核实；非全基因穷尽综述，0不等于无文献',D9_in_ranking=False,D9_is_gate=False))
 write(pd.DataFrame(out),'D9_direction_gene_annotations.tsv');log(f'D9 curated {len(ev)} primary papers/records for {ev.gene.nunique()} genes; remaining direction unknown, not absence of evidence')

def integrate():
 c=pd.read_csv(O/'candidates_axis_matched.tsv',sep='\t')
 # Replace stale status/eligibility placeholders while retaining prior evidence columns explicitly as archival.
 for col in ['all_hard_gates_passed','all_hard_gates_pass','D3_status','intervention_direction','evidence_level']:
  if col in c:c=c.rename(columns={col:col+'_pre_step3_archival'})
 for name,keys in [('gene_identity_map.tsv',['gene']),('D6_gene_annotations.tsv',['gene']),('D7_organ_safety.tsv',['gene','lineage']),('D8_four_path_annotations.tsv',['gene','lineage']),('D9_DepMap_gene_annotations.tsv',['gene']),('D9_LINCS_candidate_annotations.tsv',['gene','lineage']),('D9_direction_gene_annotations.tsv',['gene'])]:
  d=pd.read_csv(O/name,sep='\t');assert not d.duplicated(keys).any();overlap=[x for x in d if x in c and x not in keys];c=c.drop(columns=overlap).merge(d,on=keys,how='left',validate='many_to_one')
 for name,prefix in [('D7_KAC_replacement_tests.tsv','D7_KAC'),('D7_AT2L_other_controls.tsv','D7_other')]:
  d=pd.read_csv(O/name,sep='\t')
  for side in ['Precursor','LUAD']:
   sub=d[d.comparison.eq(side+'_vs_Normal')].drop(columns='comparison');sub=sub.rename(columns={k:prefix+'_'+side+'_'+k for k in sub if k!='gene'});c=c.merge(sub,on='gene',validate='many_to_one')
 c['pool_after_D7']=c.pool_before_D7&c.D7_safety_veto.eq(False)&~c.D7_data_incomplete.eq(True)
 c['top5_clear_of_existing_exclusions']=c.pool_after_D7&~c.top5_excluded_existing
 c['D7_original_KAC_subset_test']='无法进行：前体/浸润各7对，低于8对；未放宽阈值'
 c['D7_repair_risk_status']='尚不能确定：原定子集检验无法进行；替代检验及AT2L-other对照已完成，健康肺参照待完成'
 c['D7_selectivity_interpretation']='正常侧表达/斜率为描述；替代显著不证明长期安全，不显著不证明修复态共有；非AT2来源候选的AT2检验仅反映上皮旁系风险'
 c['D7_KAC_self_correlation_caution']=c.KAC_signature_member.eq(True)
 c['D8_mouse_inhibitory_resource_record']=c.D6_inhibitory_record_present
 c['D8_mouse_tool_cross_reactivity']='待核；人体机制记录不等于小鼠可用工具'
 c['D6_D7_D8_D9_status']='第3步来源与计算注释完成；待核项显式保留'
 c['rank_status']='尚未排序；未进入第4步'
 c['all_hard_gates_passed']=c.pool_after_D7
 c['intervention_direction']=c.D9_direction_annotation
 assert len(c)==5207 and not c.duplicated(['gene','lineage']).any() and c.pool_after_D7.sum()==82
 assert c.pool_before_D7.sum()==199 and c.II_axis_matched_D3_pass.sum()==0
 assert not any(col.endswith('_x') or col.endswith('_y') for col in c)
 write(c,'candidates_step3.tsv');write(c[c.pool_before_D7],'pool199_annotated.tsv');write(c[c.pool_after_D7],'pool82_after_D7_unranked.tsv')
 compact=['gene','lineage','I_upstream_eligible','II_upstream_eligible','D3_I','D3_II_histology','D3_regional','I_support_datasets','cohort_suitability_tier','D7_status','D7_reasons','D6_status','D6_clinical_highest_phase','D6_inhibitory_drug_names','D6_withdrawn_record_present','D9_direction_class','D9_direction_annotation','D9_literature_urls','D7_repair_risk_status','KAC_signature_member','top5_excluded_opposite','top5_excluded_I_LOO','top5_excluded_II_LOO','pool_after_D7','rank_status']
 write(c[c.pool_before_D7][compact],'pool199_review_summary.tsv')
 log('Integrated 5207 candidate records; 199 axis-matched, 82 organ screen survivors; unranked, no step4')

if __name__=='__main__':
 # Protocol2.6.2 supersedes these original estimator/current-table writers.
 # Audit snapshot remains runnable only by explicitly using the archived script.
 if len(sys.argv)>1 and sys.argv[1] in {'pool','safety','druggability','feasibility','kac_tests','fit_kac','integrate'}:
  raise SystemExit('Deprecated by protocol2.6.2/decisions25-28: do not overwrite reviewed current annotations or restore biased KAC statistics. Use step4/refresh_annotations.py for annotations and step3_correction for corrected estimator. Original script is audit-only.')
 if sys.argv[1]=='pool':pool()
 elif sys.argv[1]=='identities':identities()
 elif sys.argv[1]=='safety':safety()
 elif sys.argv[1]=='depmap':depmap_lincs()
 elif sys.argv[1]=='prepare_kac':prepare_kac()
 elif sys.argv[1]=='fit_kac':fit_kac()
 elif sys.argv[1]=='kac_tests':kac_tests()
 elif sys.argv[1]=='other_controls':other_controls()
 elif sys.argv[1]=='druggability':druggability()
 elif sys.argv[1]=='feasibility':feasibility()

 elif sys.argv[1]=='literature':literature()
 elif sys.argv[1]=='integrate':integrate()
