#!/opt/anaconda3/envs/ggn_target/bin/python
"""Step 0 only. No expression testing or ranking; no legacy result inputs except annotated h5ad."""
from pathlib import Path
from datetime import datetime
from zoneinfo import ZoneInfo
import csv,gzip,re,json,hashlib,sys,os
os.environ.setdefault("NUMBA_CACHE_DIR", "/Users/dudu/Documents/Codex/2026-10-09/ggn-target-ggn-target-opt-anaconda3/work/numba_cache")
import numpy as np
import pandas as pd
import h5py,openpyxl,scanpy as sc,anndata as ad
from anndata.io import read_elem,sparse_dataset
P=Path('/Users/dudu/全部数据库/GGN_target'); O=P/'results/50_screen'; O.mkdir(exist_ok=True)
now=lambda:datetime.now(ZoneInfo('Asia/Shanghai')).isoformat()
def log(s):
 print(s,flush=True)
 with (O/'run_log.txt').open('a') as f:f.write(now()+' '+s+'\n')
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for x in iter(lambda:f.read(8*1024*1024),b''):h.update(x)
 return h.hexdigest()
def save(df,name):df.to_csv(O/name,sep='\t',index=False)
log('START step 0; user confirmation required before step 1. No thresholds/weights changed. No inferential tests.')
# GEO metadata and header only: frozen expression values are not analyzed.
s=(O/'GSE282617_GEO_full.soft').read_text()
records=[]
for block in re.split(r'^\^SAMPLE = ',s,flags=re.M)[1:]:
 def field(k):
  m=re.search(r'^!Sample_'+k+r' = (.+)$',block,re.M);return m[1].strip() if m else ''
 group=re.search(r'!Sample_characteristics_ch1 = group: (\w+)',block)[1]
 records.append(dict(gsm=block.splitlines()[0],sample=field('title'),group=group,stage={'FLA':'LUSC','JRX':'IAC','WJR':'MIA','YWA':'AIS','ZCF':'Normal'}[group],tissue=field('source_name_ch1'),patient='unknown',pair_id='unknown',platform=field('platform_id')))
geo=pd.DataFrame(records);save(geo,'GSE282617_GEO_samples.tsv')
with gzip.open(P/'raw/bulk_GSE282617_AIS_MIA_IAC/GSE282617_processed_data.csv.gz','rt') as f:header=next(csv.reader(f))
fpkm=[x[5:] for x in header if x.startswith('FPKM.')];count=[x[6:] for x in header if x.startswith('count.')]
assert len(header)==142 and len(fpkm)==len(count)==70 and set(fpkm)==set(count)
lookup=geo.set_index('sample');rows=[]
for x in fpkm:
 r=dict(matrix_sample=x,fpkm_column='FPKM.'+x,count_column='count.'+x,patient='unknown',pair_id='unknown',mapping_status='exact' if x in lookup.index else 'unresolved')
 if x in lookup.index:r.update(lookup.loc[x].to_dict())
 else:r.update(dict(gsm='',stage='unknown',group='',tissue='lung',note='WYA12 is a possible YWA12 typo; not silently corrected' if x=='WYA12' else 'WJR18 exists in matrix; GEO has WJR14 instead; identity unresolved'))
 rows.append(r)
save(pd.DataFrame(rows),'GSE282617_sample_mapping.tsv')
log('GSE282617: 70 samples, group prefixes are NOT patients; GEO stages '+str(geo.stage.value_counts().to_dict())+'; missing patient pairing; 68 exact names, WYA12 and WJR18 unresolved. Both steps unable to validate under paired bulk >=7 rule.')
# Registry uses filenames/metadata only, no held-out expression analysis.
registry=[];manifest=[]
specs=[('GSE282617','bulk_GSE282617_AIS_MIA_IAC','GSE282617_processed_data.csv.gz','bulk RNA-seq; Illumina NovaSeq 6000; GPL24676','I/II pathological groups present; paired validation unavailable',70),('GSE307534','ST_GSE307534','extracted/*/*/filtered_feature_bc_matrix/matrix.mtx.gz','10x Visium (protocol); local spatial files present','post-freeze spatial localization/co-localization; patient/sample unit',None),('GSE189487','ST_GSE189487','*_matrix.mtx.gz','10x Visium (protocol); coordinates present, images absent','post-freeze independent spatial directional reference; no spot-level inference',None),('GSE300827','ST_GSE300827_AIS_MIA','*.h5','10x feature-barcode HDF5; spatial assay per protocol; coordinates absent','post-freeze non-spatial sample pseudobulk reference',None),('GSE277206','ST_GSE277206','extracted/*.h5','10x feature-barcode HDF5; spatial assay per protocol; coordinates absent','post-freeze non-spatial sample pseudobulk reference',None)]
for acc,d,pat,platform,use,n in specs:
 files=sorted((P/'raw'/d).glob(pat));n=len(files) if n is None else n
 registry.append(dict(dataset=acc,path=str(P/'raw'/d),sample_count=n,sample_count_basis='GEO and matrix header' if acc=='GSE282617' else 'local expression file inventory; patient count not inferred',platform=platform,purpose=use,status='registered_frozen',restriction='排序冻结前不得用于任何选择、排序或候选表达检验；第0步只核查结构',patient_pairing='unknown' if acc!='GSE307534' else 'patient IDs in filenames; shared Peng discovery donors, not independent cohort'))
 for f in files:manifest.append(dict(dataset=acc,path=str(f),size_bytes=f.stat().st_size,mtime_ns=f.stat().st_mtime_ns))
registry.append(dict(dataset='HLCA',path='NOT_DOWNLOADED',sample_count='',sample_count_basis='not verified',platform='integrated lung scRNA-seq reference',purpose='post-freeze healthy-lung descriptive reference only; no cross-dataset P values',status='download undecided',restriction='第5步前先报告规模和磁盘需求；排序冻结前不得用于任何选择',patient_pairing='not applicable'))
save(pd.DataFrame(registry),'frozen_sets.tsv');save(pd.DataFrame(manifest),'frozen_file_inventory.tsv')
# KAC gene set is re-extracted from publisher supplement, not legacy results.
supp=Path('/Users/dudu/Documents/Codex/2026-10-08/peng-f-cancer-cell-2026-44/outputs/肺结节文献_12篇/01_Peng_2026_Cancer_Cell/补充材料_Cell官网/Table_S3_mmc4.xlsx')
w=openpyxl.load_workbook(supp,read_only=True,data_only=True)
genes=[r[0].strip() for r in list(w.active.values)[3:] if isinstance(r[0],str) and r[0].strip()]
assert len(genes)==len(set(genes))==102
# Explicit legacy annotation mapping only; no legacy score, comparison, or threshold imported.
labels=['AT2','AT2/Club-like epithelial','Proliferating epithelial']
source=P/'results/01_snRNA/snRNA_final.h5ad'
log('AT2 lineage definition: '+str(labels)+'; patient-specific AT2-like label not included; technical annotation mapping follows existing code, no old result read.')
with h5py.File(source,'r') as f:
 obs=read_elem(f['obs']); var=read_elem(f['var']); mask=obs.final_label.isin(labels).values
 # Sparse counts subset avoids reading old normalized scores or any other result.
 X=sparse_dataset(f['layers/counts'])[mask,:]
a=ad.AnnData(X=X,obs=obs.loc[mask,['patient','sample','stage','final_label','gsm']].copy(),var=var.iloc[:,:0].copy())
assert np.all(a.X.data>=0) and np.allclose(a.X.data,np.rint(a.X.data))
sc.pp.normalize_total(a,target_sum=1e4);sc.pp.log1p(a)
pres=[g for g in genes if g in a.var_names]; assert len(pres)==91
save(pd.DataFrame({'gene':genes,'present':[g in a.var_names for g in genes]}),'KAC_TableS3_genes.tsv')
log(f'Recomputed log1p(CP10k) from counts: {a.n_obs} AT2-lineage cells, {a.obs.patient.nunique()} patients, {a.obs["sample"].nunique()} samples; KAC {len(pres)}/102 genes.')
for name,gs in [('KAC',pres),('CAS',['SCGB3A2','CEACAM6','SPINK1','CEACAM5'])]:
 sc.tl.score_genes(a,gs,score_name=name,random_state=0,ctrl_size=50,n_bins=25,gene_pool=list(a.var_names),use_raw=False,ctrl_as_ref=True)
normal=a.obs.stage.eq('Normal')
thresholds={name:float(np.quantile(a.obs.loc[normal,name].values,.9,method='linear')) for name in ['KAC','CAS']}
freeze=dict(created_at=now(),thresholds=thresholds,quantile=.9,quantile_method='linear',comparison='strictly greater than',normal_cells=int(normal.sum()),normal_patients=int(a.obs.loc[normal,'patient'].nunique()),normal_samples=int(a.obs.loc[normal,'sample'].nunique()),lineage_labels=labels,excluded_label='AT2-like (patient-specific)',score_input='counts; normalize_total target_sum=10000; log1p; all AT2 lineage cells; all 18061 genes control pool',score_parameters=dict(random_state=0,ctrl_size=50,n_bins=25,ctrl_as_ref=True,use_raw=False),scanpy_version=sc.__version__,numpy_version=np.__version__,gene_set_source=str(supp),gene_set_sha256=sha(supp),annotation_path=str(source),annotation_size=source.stat().st_size,annotation_mtime_ns=source.stat().st_mtime_ns,user_override='threshold fixed at step 0, superseding document step 1 timing',statistical_note='cell quantile used only for state definition; comparisons use patients, not cells')
# Exclusive creation: never overwrite a frozen threshold on rerun.
with (O/'state_thresholds_frozen.json').open('x') as f:json.dump(freeze,f,ensure_ascii=False,indent=2)
log('THRESHOLDS FROZEN '+str(thresholds)+' SHA256='+sha(O/'state_thresholds_frozen.json'))
a.obs['is_KAC']=a.obs.KAC>thresholds['KAC'];a.obs['is_other']=~a.obs.is_KAC
save(a.obs.reset_index(names='cell_id'),'AT2L_cell_state_assignments.tsv.gz')
sample=a.obs.groupby(['patient','sample','stage'],observed=True).agg(AT2L=('is_KAC','size'),KAC=('is_KAC','sum'),AT2L_other=('is_other','sum')).reset_index();save(sample,'repair_counts_by_sample.tsv')
rows=[]
for p in sorted(a.obs.patient.unique(),key=lambda v:int(v[1:])):
 sub=a.obs[a.obs.patient==p];row={'patient':p}
 for side,st in [('Normal',['Normal']),('Precursor',['AAH','AIS','MIA']),('Invasive',['LUAD']),('All_lesions',['AAH','AIS','MIA','LUAD'])]:
  z=sub[sub.stage.isin(st)];row.update({side+'_samples':z['sample'].nunique(),side+'_AT2L':len(z),side+'_KAC':int(z.is_KAC.sum()),side+'_other':int(z.is_other.sum())})
 for side in ['Precursor','Invasive','All_lesions']:
  for state in ['KAC','other']:
   row[f'{side}_{state}_paired_ge50']=row['Normal_'+state]>=50 and row[side+'_'+state]>=50
 rows.append(row)
t=pd.DataFrame(rows);save(t,'repair_counts_by_patient.tsv')
summary=[]
for side in ['Precursor','Invasive','All_lesions']:
 for state in ['KAC','other']:
  ok=t[f'{side}_{state}_paired_ge50'];n=int(ok.sum());summary.append(dict(comparison=side+'_vs_Normal',state=state,pairs=n,minimum_pairs=8,status='可进行' if n>=8 else '无法进行',included_patients=','.join(t.loc[ok,'patient']),excluded_patients_and_counts=';'.join(f'{r.patient}:Normal={r["Normal_"+state]},lesion={r[side+"_"+state]}' for _,r in t.loc[~ok].iterrows()),note='All_lesions is descriptive inventory only; no pooling across stages for inference' if side=='All_lesions' else 'same-patient samples aggregated within stated side'))
save(pd.DataFrame(summary),'repair_pair_summary.tsv');log(pd.DataFrame(summary).iloc[:,:5].to_string(index=False))
log('COMPLETE step 0 computations. No differential tests, cross-cohort validation, score/ranking or step 1 run.')
