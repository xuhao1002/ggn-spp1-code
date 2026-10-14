import sys, re, gzip, numpy as np, pandas as pd
sys.path.insert(0,'/Users/dudu/全部数据库/GGN_target/scripts')
from prolif_genes import prolif_set_stage2, TIROSH
from scipy import stats
from scipy.stats import mannwhitneyu
import statsmodels.api as sm
R='/Users/dudu/全部数据库/GGN_target/'; RAW=R+'raw/'
names={}
with gzip.open(RAW+'annotation/gencode.v44.basic.annotation.gtf.gz','rt') as f:
    for line in f:
        if '\tgene\t' in line:
            m=re.search(r'gene_id "([^.]+)\.[^"]*"; gene_type "[^"]*"; gene_name "([^"]+)"',line)
            if m: names[m.group(1)]=m.group(2)
sets={'prolif196':set(prolif_set_stage2()),'tirosh97':set(TIROSH.split() if isinstance(TIROSH,str) else TIROSH)}
print({k:len(v) for k,v in sets.items()})
def stage3(s):
    s=str(s)
    for k,v in [('IV',3),('III',3),('II',2),('I',1)]:
        if f'Stage {k}' in s: return v
    return np.nan
def scores(tpm):
    X=np.log2(tpm.T+1); X.columns=[names.get(c,c) for c in X.columns]; X=X.loc[:,~X.columns.duplicated()]
    Z=(X-X.mean())/X.std(ddof=0).replace(0,1); out={}
    for k,v in sets.items():
        g=[x for x in v if x in Z.columns]; out[k]=Z[g].mean(1); print('  ',k,'基因匹配',len(g),'/',len(v))
    return pd.DataFrame(out)
# TCGA
lt=pd.read_excel(RAW+'STAS_labels/STAS_labels_TCGA_from_PanLR_20261011.xlsx'); yt=lt.groupby('患者编号').label.agg(['max','size'])
hdr=pd.read_csv(RAW+'TCGA_LUAD/TCGA-LUAD.star_tpm.tsv.gz',sep='\t',nrows=0).columns
pick={}
for c in hdr[1:]:
    if c[:12] in yt.index and c.split('-')[3].startswith('01') and c[:12] not in pick: pick[c[:12]]=c
pts=sorted(pick); tpm=pd.read_csv(RAW+'TCGA_LUAD/TCGA-LUAD.star_tpm.tsv.gz',sep='\t',usecols=['Ensembl_ID']+[pick[p] for p in pts],index_col=0)
tpm.columns=pts; tpm.index=tpm.index.str.split('.').str[0]; tpm=tpm[~tpm.index.duplicated()]
print('TCGA'); st=scores(tpm); st['stas']=yt.loc[pts,'max'].values; st['ns']=yt.loc[pts,'size'].values
cl=pd.read_csv(RAW+'TCGA_LUAD/TCGA-LUAD.clinical.tsv.gz',sep='\t'); cl=cl[cl['sample_type.samples']=='Primary Tumor'].drop_duplicates('submitter_id').set_index('submitter_id')
st['stage']=cl['ajcc_pathologic_stage.diagnoses'].reindex(pts).map(stage3).values
# CPTAC
lc=pd.read_excel(RAW+'STAS_labels/STAS_labels_CPTAC_from_PanLR_20261011.xlsx'); yc=lc.groupby('患者编号').STAS.agg(['max','size'])
f=RAW+'STAS_clinical/'; cc=pd.read_csv(f+'CPTAC-3.clinical.tsv.gz',sep='\t',dtype=str)
cols=[x for x in pd.read_csv(f+'CPTAC-3.star_tpm.tsv.gz',sep='\t',nrows=0).columns if x.count('-')==2]
t=cc[(cc['sample_type.samples']=='Primary Tumor')&cc['sample'].isin(cols)&cc['submitter_id'].isin(yc.index)]
pk=t.sort_values('sample').groupby('submitter_id')['sample'].first(); pc=sorted(pk.index)
tp=pd.read_csv(f+'CPTAC-3.star_tpm.tsv.gz',sep='\t',usecols=['Ensembl_ID']+list(pk[pc]),index_col=0); tp.columns=pc; tp.index=tp.index.str.split('.').str[0]; tp=tp[~tp.index.duplicated()]
print('CPTAC'); sc=scores(tp); sc['stas']=yc.loc[pc,'max'].values; sc['ns']=np.minimum(yc.loc[pc,'size'].values,4)
sc['c3n']=[p[:3]=='C3N' for p in pc]; cd=cc.drop_duplicates('submitter_id').set_index('submitter_id')
sc['grade']=cd['tumor_grade.diagnoses'].reindex(pc).map({'G1':0,'G2':0,'G3':1}).values
sc['stage']=cd['ajcc_pathologic_stage.diagnoses'].reindex(pc).map(stage3).values
res={}
def run(df,key,tag,covs):
    a=df[df.stas==1][key];b=df[df.stas==0][key]; u=mannwhitneyu(a,b,alternative='greater')
    d=df.dropna(subset=[key]+covs); X=sm.add_constant(pd.concat([(d[key]-d[key].mean())/d[key].std()]+[d[c].astype(float) for c in covs],axis=1).astype(float)); X.columns=['const','score']+covs
    m=sm.Logit(d.stas.astype(float),X).fit(disp=0); b1,se=m.params['score'],m.bse['score']
    print(f'{tag} [{key}] n={len(d)} STAS+均值{a.mean():.3f} vs STAS-{b.mean():.3f}; MWU单侧P={u.pvalue:.3f}; 校正后 OR/SD={np.exp(b1):.2f} (95%CI {np.exp(b1-1.96*se):.2f}-{np.exp(b1+1.96*se):.2f}) 单侧P={stats.norm.sf(b1/se):.3f}')
    return b1,se
print('== 管线阳性对照')
g=sc.dropna(subset=['grade']); u=mannwhitneyu(g[g.grade==1].prolif196,g[g.grade==0].prolif196,alternative='greater'); print(f'CPTAC G3 vs G1-2 prolif196: 均值 {g[g.grade==1].prolif196.mean():.2f} vs {g[g.grade==0].prolif196.mean():.2f}, 单侧P={u.pvalue:.4f} (n={len(g)})')
s=st.dropna(subset=['stage']); r=stats.spearmanr(s.stage,s.prolif196); print(f'TCGA 分期(1/2/3) 与 prolif196 Spearman={r[0]:.3f} P={r[1]:.4f}')
print('== 主要检验')
r1=run(st,'prolif196','TCGA',['stage']); r2=run(sc,'prolif196','CPTAC',['stage','ns','c3n'])
w=np.array([1/r1[1]**2,1/r2[1]**2]); bp=(w*np.array([r1[0],r2[0]])).sum()/w.sum(); sp=w.sum()**-.5
print(f'合并 OR/SD={np.exp(bp):.2f} (95%CI {np.exp(bp-1.96*sp):.2f}-{np.exp(bp+1.96*sp):.2f}) 单侧P={stats.norm.sf(bp/sp):.3f}')
print('== 次要/敏感性')
run(st,'tirosh97','TCGA',['stage']); run(sc,'tirosh97','CPTAC',['stage','ns','c3n'])
run(st[st.ns==1],'prolif196','TCGA 仅单切片',['stage']); run(sc[sc.ns==1],'prolif196','CPTAC 仅单切片',['stage'])
run(st,'prolif196','TCGA 未校正',[]); run(sc,'prolif196','CPTAC 未校正',[])
