import numpy as np, pandas as pd, re, os
from scipy import stats
from scipy.stats import spearmanr, binomtest, chi2_contingency
R='/Users/dudu/全部数据库/GGN_target/'; D=R+'raw/CPTAC_LUAD_PDC/'
lab=pd.read_excel(R+'raw/STAS_labels/STAS_labels_CPTAC_from_PanLR_20261011.xlsx')
y_pt=lab.groupby('患者编号').STAS.max(); ns=lab.groupby('患者编号').size()
m=pd.read_csv(D+'aliquot_map_labeled151.tsv',sep='\t'); t=m[m.sample_type=='Primary Tumor']
a2c=t.drop_duplicates('case').set_index('aliquot').case.to_dict()
def plexmap(f):
    mp={}
    for i,line in enumerate(open(D+f)):
        if i==0: continue
        parts=line.rstrip('\n').split('\t')
        for x in parts[2:]:
            if re.fullmatch(r'CPT\d{10}',x): mp[x]=parts[1]
    return mp
def load(f,sf):
    pm=plexmap(sf); d=pd.read_csv(D+f,sep='\t',index_col=0); d=d.drop([i for i in ['Mean','Median','StdDev'] if i in d.index])
    cols={};plex={}
    for c in d.columns:
        if 'Unshared' in c: continue
        mm=re.search(r'CPT\d{10}',c)
        if mm and mm.group(0) in a2c and a2c[mm.group(0)] not in cols:
            cols[a2c[mm.group(0)]]=c; plex[a2c[mm.group(0)]]=pm.get(mm.group(0),'NA')
    return d[list(cols.values())].set_axis(list(cols.keys()),axis=1), pd.Series(plex)
fl=os.listdir(D)
dis,pd_=load([x for x in fl if x.startswith('PDC000153') and x.endswith('tmt10.tsv')][0],[x for x in fl if x.startswith('PDC000153') and x.endswith('sample.txt')][0])
con,pc_=load([x for x in fl if x.startswith('PDC000489') and x.endswith('tmt11.tsv')][0],[x for x in fl if x.startswith('PDC000489') and x.endswith('sample.txt')][0])
dup=[c for c in con.columns if c in dis.columns]; con=con.drop(columns=dup); pc_=pc_.drop(index=dup)
def de(d,plex,tag):
    pts=list(d.columns); g=plex.loc[pts]
    cnt=g.value_counts(); keep=[p for p in pts if cnt[g[p]]>=2 and g[p]!='NA']
    # 仅保留 plex 内 STAS 有差异的 plex 对估计有贡献；其余 plex 的组内 STAS 方差=0 不影响估计但占自由度
    y=y_pt[keep].values.astype(float); k=np.minimum(ns[keep].values,4).astype(float); sr=np.array([p[:3]=='C3N' for p in keep],float); gg=g.loc[keep].values
    ct=pd.crosstab(g.loc[pts],y_pt[pts]); chi=chi2_contingency(ct.values)[1] if ct.shape[0]>1 else np.nan
    var_plex=sorted({p for p in set(gg) if 0<y[gg==p].mean()<1})
    print(tag,'带标签',len(pts),'plex数',g.nunique(),'纳入(plex内≥2)',len(keep),'有STAS差异的plex',len(var_plex),'plex×STAS 卡方P=',round(chi,4))
    def demean(v,gl):
        out=v.copy().astype(float)
        for p in set(gl): out[gl==p]-=out[gl==p].mean()
        return out
    Yd,Kd,Sd=demean(y,gg),demean(k,gg),demean(sr,gg); npl=len(set(gg))
    rows=[]
    for gname,v in d[keep].iterrows():
        x=v.values.astype(float); ok=~np.isnan(x)
        if ok.mean()<0.7 or ok.sum()<20: continue
        # 缺失处理：在 ok 子集上重新组内去均值
        gk=gg[ok]; xo=demean(x[ok],gk); yo=demean(y[ok],gk); ko=demean(k[ok],gk); so=demean(sr[ok],gk)
        X=np.column_stack([yo,so,ko]); b=np.linalg.lstsq(X,xo,rcond=None)[0]; dof=ok.sum()-len(set(gk))-3
        if dof<10: continue
        s2=((xo-X@b)**2).sum()/dof; cov=s2*np.linalg.pinv(X.T@X); tt=b[0]/np.sqrt(cov[0,0]); rows.append((gname,b[0],tt,2*stats.t.sf(abs(tt),dof)))
    r=pd.DataFrame(rows,columns=['gene','beta','t','p']).set_index('gene'); r['fdr']=stats.false_discovery_control(r.p.values)
    print(tag,'蛋白',len(r),'FDR<0.10:',int((r.fdr<0.10).sum()),'P<0.05:',int((r.p<0.05).sum())); return r
rd=de(dis,pd_,'Discovery'); rc=de(con,pc_,'Confirmatory')
com=rd.index.intersection(rc.index); rho,p=spearmanr(rd.loc[com,'t'],rc.loc[com,'t']); print('plex校正后 两研究 STAS t 值 Spearman',round(rho,3),'P',p,len(com))
sig=rd[rd.fdr<0.10].index.intersection(rc.index)
if len(sig):
    same=np.sign(rd.loc[sig,'beta'])==np.sign(rc.loc[sig,'beta']); rep=same&(rc.loc[sig,'p']<0.05)
    print('Discovery 发现(可比)',len(sig),'方向一致',int(same.sum()),'二项P',round(binomtest(int(same.sum()),len(sig),0.5).pvalue,4),'复现',int(rep.sum()))
else: print('Discovery 校正后无 FDR<0.10 蛋白')
