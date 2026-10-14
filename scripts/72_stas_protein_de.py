import numpy as np, pandas as pd, re
from scipy import stats
from scipy.stats import spearmanr, binomtest
R='/Users/dudu/全部数据库/GGN_target/'; D=R+'raw/CPTAC_LUAD_PDC/'; O=R+'results/72_stas_protein/'
import os; os.makedirs(O,exist_ok=True)
lab_s=pd.read_excel(R+'raw/STAS_labels/STAS_labels_CPTAC_from_PanLR_20261011.xlsx')
y_pt=lab_s.groupby('患者编号').STAS.max(); ns=lab_s.groupby('患者编号').size()
m=pd.read_csv(D+'aliquot_map_labeled151.tsv',sep='\t'); t=m[m.sample_type=='Primary Tumor']
a2c=t.drop_duplicates('case').set_index('aliquot').case.to_dict()
a2c.update(t.set_index('aliquot').case.to_dict())
def load(f):
    d=pd.read_csv(D+f,sep='\t',index_col=0)
    d=d.drop([i for i in ['Mean','Median','StdDev'] if i in d.index])
    cols={}
    for c in d.columns:
        if c.startswith('Unshared') or 'Unshared' in c: continue
        mm=re.search(r'CPT\d{10}',c)
        if not mm: continue
        a=mm.group(0)
        if a in a2c and a2c[a] not in cols: cols[a2c[a]]=c
    return d[list(cols.values())].set_axis(list(cols.keys()),axis=1)
dis=load([x for x in os.listdir(D) if x.startswith('PDC000153') and x.endswith('tmt10.tsv')][0])
con=load([x for x in os.listdir(D) if x.startswith('PDC000489') and x.endswith('tmt11.tsv')][0])
con=con.drop(columns=[c for c in con.columns if c in dis.columns])
print('Discovery',dis.shape,'Confirmatory(去重后)',con.shape)
def de(d,tag,single=False):
    pts=list(d.columns)
    if single: pts=[p for p in pts if ns[p]==1]
    y=y_pt[pts].values.astype(float); k=np.minimum(ns[pts].values,4).astype(float)
    rows=[]
    for g,v in d[pts].iterrows():
        x=v.values.astype(float); ok=~np.isnan(x)
        if ok.mean()<0.7 or ok.sum()<20: continue
        Xd=np.column_stack([np.ones(ok.sum()),y[ok]] + ([] if single else [k[ok]]))
        b,res,rk,_=np.linalg.lstsq(Xd,x[ok],rcond=None)
        dof=ok.sum()-Xd.shape[1]; rss=((x[ok]-Xd@b)**2).sum(); s2=rss/dof
        cov=s2*np.linalg.inv(Xd.T@Xd); tt=b[1]/np.sqrt(cov[1,1]); p=2*stats.t.sf(abs(tt),dof)
        rows.append((g,b[1],tt,p,ok.sum()))
    r=pd.DataFrame(rows,columns=['gene','beta','t','p','n']).set_index('gene')
    r['fdr']=stats.false_discovery_control(r.p.values,method='bh') if len(r) else []
    print(tag,'n患者',len(pts),'阳',int(y.sum()),'阴',int((1-y).sum()),'检验蛋白',len(r),'FDR<0.10:',int((r.fdr<0.10).sum()),'名义P<0.05:',int((r.p<0.05).sum()))
    r.to_csv(O+f'de_{tag}.tsv',sep='\t'); return r
rd=de(dis,'Discovery'); rc=de(con,'Confirmatory')
com=rd.index.intersection(rc.index)
rho,pr=spearmanr(rd.loc[com,'t'],rc.loc[com,'t']);print('两研究 t 值 Spearman',round(rho,3),'P',pr,'共同蛋白',len(com))
sig=rd[(rd.fdr<0.10)].index.intersection(rc.index)
if len(sig):
    same=(np.sign(rd.loc[sig,'beta'])==np.sign(rc.loc[sig,'beta']))
    rep=same&(rc.loc[sig,'p']<0.05)
    print('Discovery 发现(可比)',len(sig),'方向一致',int(same.sum()),'二项P(vs50%)',round(binomtest(int(same.sum()),len(sig),0.5).pvalue,4),'复现(同向且P<0.05)',int(rep.sum()))
    rd.loc[sig].assign(conf_beta=rc.loc[sig,'beta'],conf_p=rc.loc[sig,'p'],replicated=rep).sort_values('p').to_csv(O+'discovery_hits_with_replication.tsv',sep='\t')
else: print('Discovery 无 FDR<0.10 的蛋白，无法进入复现')
print('--敏感性：仅单切片患者')
sd=de(dis,'Discovery_single',True); sc=de(con,'Confirmatory_single',True)
com=sd.index.intersection(sc.index);rho,pr=spearmanr(sd.loc[com,'t'],sc.loc[com,'t']);print('单切片 两研究 t 值 Spearman',round(rho,3),'P',pr,len(com))
