"""冻结方案 P1 敏感性中此前未跑的项：切片层面（患者聚类）、仅 I 期、剔除吸烟包年缺失。使用已冻结评分（results/70_stas/cptac_scores.tsv）。"""
import numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score
import statsmodels.api as sm
R='/Users/dudu/全部数据库/GGN_target/'; RAW=R+'raw/'
s=pd.read_csv(R+'results/70_stas/cptac_scores.tsv',sep='\t').set_index('patient')
sl=pd.read_excel(RAW+'STAS_labels/STAS_labels_CPTAC_from_PanLR_20261011.xlsx')
sl=sl[sl['患者编号'].isin(s.index)].copy(); sl['score']=s.loc[sl['患者编号'],'score_full'].values
def cluster_boot_auc(d,B=2000,seed=1):
    rng=np.random.default_rng(seed); pts=d['患者编号'].unique(); g={p:d[d['患者编号']==p] for p in pts}; a=[]
    for _ in range(B):
        dd=pd.concat([g[p] for p in rng.choice(pts,len(pts))])
        if dd.STAS.nunique()==2: a.append(roc_auc_score(dd.STAS,dd.score))
    return np.percentile(a,[2.5,97.5])
def slide_level(d,tag):
    auc=roc_auc_score(d.STAS,d.score); lo,hi=cluster_boot_auc(d)
    X=sm.add_constant(((d.score-d.score.mean())/d.score.std()).values); gee=sm.GEE(d.STAS.values,X,groups=d['患者编号'].values,family=sm.families.Binomial(),cov_struct=sm.cov_struct.Exchangeable()).fit()
    print(f'{tag}: 切片{len(d)} 患者{d["患者编号"].nunique()} 阳性切片{int(d.STAS.sum())} 切片AUC={auc:.3f} (患者聚类bootstrap 95%CI {lo:.3f}-{hi:.3f}) GEE OR/SD={np.exp(gee.params[1]):.2f} P={gee.pvalues[1]:.3f}')
slide_level(sl,'切片层面（患者聚类）')
# 患者层面的子集敏感性
cl=pd.read_csv(RAW+'STAS_clinical/CPTAC-3.clinical.tsv.gz',sep='\t',dtype=str).drop_duplicates('submitter_id').set_index('submitter_id')
st=cl['ajcc_pathologic_stage.diagnoses'].reindex(s.index); sI=st.str.contains(r'Stage I[ABC]?\d?$',na=False)&~st.str.contains('II',na=False)
pk=pd.to_numeric(cl['pack_years_smoked.exposures'].reindex(s.index),errors='coerce')
def pat(mask,tag):
    d=s[mask.values]; y=d.stas.values
    if len(set(y))<2: print(tag,'不可评估'); return
    rng=np.random.default_rng(1); bs=[roc_auc_score(y[i],d.score_full.values[i]) for i in (rng.integers(0,len(y),len(y)) for _ in range(2000)) if y[i].min()!=y[i].max()]
    print(f'{tag}: n={len(d)} 阳{int(y.sum())}/阴{int((1-y).sum())} AUC={roc_auc_score(y,d.score_full):.3f} (95%CI {np.percentile(bs,2.5):.3f}-{np.percentile(bs,97.5):.3f})')
pat(sI,'仅 I 期'); pat(pk.notna(),'剔除吸烟包年缺失者')
