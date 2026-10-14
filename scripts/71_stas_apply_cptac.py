"""P1 第2步：把冻结的 TCGA 评分一次性应用到 CPTAC RNA（预设方案 P1）。
样本规则（运行前声明）：每名患者取 GDC 标注为 Primary Tumor 且在 RNA 矩阵中的样本，多个时取编号最小者；CPTAC 内 z 标准化。"""
import hashlib, json, numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score
from sklearn.linear_model import LogisticRegression
R='/Users/dudu/全部数据库/GGN_target/'; O=R+'results/70_stas/'
# 校验冻结文件未被改动
for l in open(O+'frozen_sha256.txt'):
    h,fn=l.split(); assert hashlib.sha256(open(O+fn,'rb').read()).hexdigest()==h, fn
meta=json.load(open(O+'train_meta.json')); top=meta['top5000_genes']
full=pd.read_csv(O+'score_full_coefficients.tsv',sep='\t').set_index('ensembl').coef
c20=pd.read_csv(O+'score_compact20_coefficients.tsv',sep='\t').set_index('ensembl').coef
lab_s=pd.read_excel(R+'raw/STAS_labels/STAS_labels_CPTAC_from_PanLR_20261011.xlsx')
y_pt=lab_s.groupby('患者编号').STAS.max()
nslide=lab_s.groupby('患者编号').size(); incons=lab_s.groupby('患者编号').STAS.nunique()>1
f=R+'raw/STAS_clinical/'
cl=pd.read_csv(f+'CPTAC-3.clinical.tsv.gz',sep='\t',dtype=str)
cols=[x for x in pd.read_csv(f+'CPTAC-3.star_tpm.tsv.gz',sep='\t',nrows=0).columns if x.count('-')==2]
t=cl[(cl['sample_type.samples']=='Primary Tumor')&cl['sample'].isin(cols)&cl['submitter_id'].isin(y_pt.index)]
pick=t.sort_values('sample').groupby('submitter_id')['sample'].first()
pts=sorted(pick.index); print('患者',len(pts),'阳',int(y_pt[pts].sum()),'阴',len(pts)-int(y_pt[pts].sum()))
tpm=pd.read_csv(f+'CPTAC-3.star_tpm.tsv.gz',sep='\t',usecols=['Ensembl_ID']+list(pick[pts]),index_col=0)
tpm.columns=pts; tpm.index=tpm.index.str.split('.').str[0]; tpm=tpm[~tpm.index.duplicated()]
genes=[g for g in top if g in tpm.index]; print('基因匹配',len(genes),'/5000')
X=np.log2(tpm.loc[genes].T+1); Z=(X-X.mean())/X.std(ddof=0).replace(0,1)
y=y_pt[pts].values.astype(int)
def ev(score,yy,tag,B=2000,seed=1):
    a=roc_auc_score(yy,score); rng=np.random.default_rng(seed); bs=[]
    for _ in range(B):
        i=rng.integers(0,len(yy),len(yy))
        if yy[i].min()==yy[i].max(): continue
        bs.append(roc_auc_score(yy[i],score[i]))
    lo,hi=np.percentile(bs,[2.5,97.5])
    m=LogisticRegression(C=1e6,max_iter=1000).fit(((score-score.mean())/score.std()).reshape(-1,1),yy)
    print(f'{tag}: n={len(yy)} 阳{int(yy.sum())}/阴{int((1-yy).sum())} AUC={a:.3f} (95%CI {lo:.3f}-{hi:.3f}) 校准斜率(每SD)={m.coef_[0][0]:.3f}'); return a,lo,hi
s_full=(Z[full.index.intersection(Z.columns)]*full[full.index.intersection(Z.columns)]).sum(1).values
s_20=(Z[c20.index.intersection(Z.columns)]*c20[c20.index.intersection(Z.columns)]).sum(1).values
print('精简版匹配基因',len(c20.index.intersection(Z.columns)),'/20')
res={'full':ev(s_full,y,'P1 主分析 完整评分'),'compact20':ev(s_20,y,'精简20基因')}
ns=nslide.loc[pts].values; ic=incons.loc[pts].values
m1=ns==1; print('--敏感性')
res['single']=ev(s_full[m1],y[m1],'仅单切片患者')
m2=~ic; res['consistent']=ev(s_full[m2],y[m2],'去除切片间不一致患者')
json.dump({k:list(map(float,v)) for k,v in res.items()},open(O+'cptac_p1_result.json','w'))
pd.DataFrame({'patient':pts,'sample':pick[pts].values,'stas':y,'score_full':s_full,'score_c20':s_20}).to_csv(O+'cptac_scores.tsv',sep='\t',index=False)
