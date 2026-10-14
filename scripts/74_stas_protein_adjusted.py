import numpy as np, pandas as pd, re, os
from scipy import stats
from scipy.stats import spearmanr, binomtest
R='/Users/dudu/全部数据库/GGN_target/'; D=R+'raw/CPTAC_LUAD_PDC/'; O=R+'results/74_stas_protein_adj/'; os.makedirs(O,exist_ok=True)
lab=pd.read_excel(R+'raw/STAS_labels/STAS_labels_CPTAC_from_PanLR_20261011.xlsx')
y_pt=lab.groupby('患者编号').STAS.max(); ns=lab.groupby('患者编号').size()
m=pd.read_csv(D+'aliquot_map_labeled151.tsv',sep='\t'); t=m[m.sample_type=='Primary Tumor']
a2c=t.drop_duplicates('case').set_index('aliquot').case.to_dict()
def load(f):
    d=pd.read_csv(D+f,sep='\t',index_col=0); d=d.drop([i for i in ['Mean','Median','StdDev'] if i in d.index])
    cols={}
    for c in d.columns:
        if 'Unshared' in c: continue
        mm=re.search(r'CPT\d{10}',c)
        if mm and mm.group(0) in a2c and a2c[mm.group(0)] not in cols: cols[a2c[mm.group(0)]]=c
    return d[list(cols.values())].set_axis(list(cols.keys()),axis=1)
dis=load([x for x in os.listdir(D) if x.startswith('PDC000153') and x.endswith('tmt10.tsv')][0])
con=load([x for x in os.listdir(D) if x.startswith('PDC000489') and x.endswith('tmt11.tsv')][0])
con=con.drop(columns=[c for c in con.columns if c in dis.columns])
def fit(d,tag):
    pts=list(d.columns); y=y_pt[pts].values.astype(float); k=np.minimum(ns[pts].values,4).astype(float); sr=np.array([p[:3]=='C3N' for p in pts],float)
    rows=[]
    for g,v in d.iterrows():
        x=v.values.astype(float); ok=~np.isnan(x)
        if ok.mean()<0.7 or ok.sum()<20: continue
        X=np.column_stack([np.ones(ok.sum()),y[ok],sr[ok],k[ok]])
        b=np.linalg.lstsq(X,x[ok],rcond=None)[0]; dof=ok.sum()-4; s2=((x[ok]-X@b)**2).sum()/dof
        cov=s2*np.linalg.inv(X.T@X); se=np.sqrt(np.diag(cov))
        rows.append((g,b[1],b[1]/se[1],2*stats.t.sf(abs(b[1]/se[1]),dof),b[2]/se[2],ok.sum()))
    r=pd.DataFrame(rows,columns=['gene','beta','t','p','t_series','n']).set_index('gene'); r['fdr']=stats.false_discovery_control(r.p.values)
    print(tag,'n患者',len(pts),'(C3L',int((1-sr).sum()),'C3N',int(sr.sum()),') 阳',int(y.sum()),'阴',int((1-y).sum()),'蛋白',len(r),'FDR<0.10:',int((r.fdr<0.10).sum()),'P<0.05:',int((r.p<0.05).sum()))
    r.to_csv(O+f'de_adj_{tag}.tsv',sep='\t'); return r
rd=fit(dis,'Discovery'); rc=fit(con,'Confirmatory')
com=rd.index.intersection(rc.index)
print('校正后 STAS t 值 两研究 Spearman',round(spearmanr(rd.loc[com,'t'],rc.loc[com,'t'])[0],3),'P',spearmanr(rd.loc[com,'t'],rc.loc[com,'t'])[1],len(com))
print('系列效应 t 值 两研究 Spearman',round(spearmanr(rd.loc[com,'t_series'],rc.loc[com,'t_series'])[0],3),'P',spearmanr(rd.loc[com,'t_series'],rc.loc[com,'t_series'])[1])
sig=rd[rd.fdr<0.10].index.intersection(rc.index)
if len(sig):
    same=np.sign(rd.loc[sig,'beta'])==np.sign(rc.loc[sig,'beta']); rep=same&(rc.loc[sig,'p']<0.05)
    print('Discovery 发现(可比)',len(sig),'方向一致',int(same.sum()),'二项P',round(binomtest(int(same.sum()),len(sig),0.5).pvalue,4),'复现',int(rep.sum()))
else: print('Discovery 校正后无 FDR<0.10 蛋白')
# 诊断：原 S4 的 STAS t 与系列效应 t 的相关（各研究内）
for tag,rr in [('Discovery',rd),('Confirmatory',rc)]:
    o=pd.read_csv(R+f'results/72_stas_protein/de_{tag}.tsv',sep='\t',index_col=0)
    ix=o.index.intersection(rr.index)
    print(tag,'原S4 STAS t 与 系列效应 t 的 Spearman',round(spearmanr(o.loc[ix,'t'],rr.loc[ix,'t_series'])[0],3),'P',spearmanr(o.loc[ix,'t'],rr.loc[ix,'t_series'])[1])
