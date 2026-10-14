import numpy as np, pandas as pd
from scipy import stats
from scipy.stats import fisher_exact, mannwhitneyu
import statsmodels.api as sm
R='/Users/dudu/全部数据库/GGN_target/raw/'
lt=pd.read_excel(R+'STAS_labels/STAS_labels_TCGA_from_PanLR_20261011.xlsx'); y=lt.groupby('患者编号').label.agg(['max','size'])
m=pd.read_csv(R+'TCGA_LUAD/TCGA-LUAD.somaticmutation_wxs.tsv.gz',sep='\t')
m=m[m['sample'].str.endswith('-01A')|m['sample'].str.contains(r'-01[A-Z]$')]
m['pt']=m['sample'].str[:12]; m['effect']=m['effect'].fillna('')
syn=m.effect.str.contains('synonymous|intron|UTR|non_coding|upstream|downstream|intergenic',regex=True)&~m.effect.str.contains('missense|stop_gained|frameshift|splice_acceptor|splice_donor|inframe|start_lost|stop_lost')
cod=m[~syn]
pts=sorted(set(y.index)&set(m.pt)); print('有突变数据的带标签患者',len(pts),'阳',int(y.loc[pts,'max'].sum()))
d=pd.DataFrame({'stas':y.loc[pts,'max'],'ns':y.loc[pts,'size']})
d['tmb']=np.log2(cod.groupby('pt').size().reindex(pts).fillna(0)+1)
for g in ['TP53','KRAS','EGFR','STK11','KEAP1']: d[g]=pts and pd.Series({p:int(p in set(cod[cod.gene==g].pt)) for p in pts})
cl=pd.read_csv(R+'TCGA_LUAD/TCGA-LUAD.clinical.tsv.gz',sep='\t'); cl=cl[cl['sample_type.samples']=='Primary Tumor'].drop_duplicates('submitter_id').set_index('submitter_id')
def st(s):
    s=str(s)
    for k,v in [('IV',3),('III',3),('II',2),('I',1)]:
        if f'Stage {k}' in s: return v
    return np.nan
d['stage']=cl['ajcc_pathologic_stage.diagnoses'].reindex(pts).map(st)
def one(df,tag):
    print('==',tag,'n',len(df),'STAS+',int(df.stas.sum()))
    a=df[df.stas==1];b=df[df.stas==0]
    # TP53
    t=pd.crosstab(df.TP53,df.stas).reindex(index=[0,1],columns=[0,1],fill_value=0).values
    orr,p2=fisher_exact(t); p1=fisher_exact(t,alternative='greater')[1]  # 行=TP53,列=STAS; OR>1 => TP53突变在STAS+更多（富集）
    X=sm.add_constant(df[['TP53','stage']].dropna()); yy=df.stas.loc[X.index]
    lr=sm.Logit(yy,X).fit(disp=0); 
    print(f'TP53 突变率 STAS+ {a.TP53.mean():.2f} ({int(a.TP53.sum())}/{len(a)}) vs STAS- {b.TP53.mean():.2f} ({int(b.TP53.sum())}/{len(b)}); OR={orr:.2f} Fisher单侧P={p1:.3f}; 校正分期 OR={np.exp(lr.params["TP53"]):.2f} P(单侧)={stats.norm.sf(lr.params["TP53"]/lr.bse["TP53"]):.3f}')
    u=mannwhitneyu(a.tmb,b.tmb,alternative='greater'); print(f'TMB 中位 STAS+ {2**a.tmb.median()-1:.0f} vs STAS- {2**b.tmb.median()-1:.0f}; 单侧P={u.pvalue:.3f}')
    rows=[]
    for g in ['KRAS','EGFR','STK11','KEAP1']:
        tt=pd.crosstab(df[g],df.stas).reindex(index=[0,1],columns=[0,1],fill_value=0).values; o,p=fisher_exact(tt); rows.append((g,round(a[g].mean(),2),round(b[g].mean(),2),round(o,2),p))
    r=pd.DataFrame(rows,columns=['gene','STAS+率','STAS-率','OR','P双侧']); r['BH']=stats.false_discovery_control(r.P双侧.values); print(r.round(3).to_string(index=False))
one(d,'TCGA 全部')
one(d[d.ns==1],'TCGA 仅单切片')
# CPTAC 描述
cli=pd.read_csv(R+'CPTAC_LUAD_linkedomics/HS_CPTAC_LUAD_cli.tsi',sep='\t',index_col=0).iloc[1:]
lc=pd.read_excel(R+'STAS_labels/STAS_labels_CPTAC_from_PanLR_20261011.xlsx').groupby('患者编号').STAS.max(); lc.index=lc.index.str.replace('-','.')
c=cli.join(lc.rename('stas'),how='inner'); print('== CPTAC LinkedOmics 临床 带标签',len(c),'阳',int(c.stas.sum()))
for g in ['TP53','KRAS','STK11','EGFR','KEAP1']:
    col=f'{g}.mutation.status'; c[col]=pd.to_numeric(c[col],errors='coerce'); k=c.dropna(subset=[col])
    t=pd.crosstab(k[col],k.stas).reindex(index=[0,1],columns=[0,1],fill_value=0).values; o,p=fisher_exact(t)
    print(f'{g}: STAS+ {k[k.stas==1][col].mean():.2f} vs STAS- {k[k.stas==0][col].mean():.2f}; OR={o:.2f} 双侧P={p:.3f} (n={len(k)})')
