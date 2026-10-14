"""预设方案_敏感性与蛋白验证_2026-10-12.md 第 B 节：IL-1β 相关阴性结论的检验力（单位定义同 scripts/40）。环境 ggn_target。"""
import os, json, hashlib
import numpy as np, pandas as pd, scanpy as sc
from scipy import stats
P=os.path.expanduser("~/全部数据库/GGN_target"); OUT=f"{P}/results/sens_protein/"; os.makedirs(OUT,exist_ok=True)
sha=hashlib.sha256(open(f"{P}/docs/预设方案_敏感性与蛋白验证_2026-10-12.md","rb").read()).hexdigest(); assert sha.startswith("b15055f8"),sha
LINS={"AT2L":["AT2","AT2/Club-like epithelial","Proliferating epithelial"],"Myeloid":["Macrophage","Monocyte","cDC","pDC","Neutrophil","Proliferating myeloid"]}
STEPS={"AAH":["AAH"],"AISMIA":["AIS","MIA"]}; MIN=50
def summ(d):
    d=np.asarray(d,float); d=d[~np.isnan(d)]; n=len(d); m=d.mean(); sd=d.std(ddof=1); se=sd/np.sqrt(n)
    t975=stats.t.ppf(0.975,n-1); t95=stats.t.ppf(0.95,n-1); t80=stats.t.ppf(0.80,n-1)
    return {"n":n,"mean_diff":float(m),"sd_diff":float(sd),"ci95":[float(m-t975*se),float(m+t975*se)],"upper95_one_sided":float(m+t95*se),"MDE80_one_sided":float((t95+t80)*se)}
res={"plan_sha256":sha}
fr=pd.read_csv(f"{P}/results/40_drivers/stepwise/stepwise_fractions_per_patient.tsv",sep="\t")
for (step,ind),g in fr.groupby(["step","indicator"]): res[f"{step}|{ind}"]=summ(g["diff"])
a=sc.read_h5ad(f"{P}/results/01_snRNA/snRNA_final.h5ad")
o=a.obs; pt=o["patient"].astype(str).values; stg=o["stage"].astype(str).values; fl=o["final_label"].astype(str).values
epcam=np.asarray(a[:,"EPCAM"].layers["counts"].todense()).ravel(); X=a.layers["counts"].tocsr()
patients=sorted(set(pt),key=lambda x:int(x[1:]))
for gene,lin,filt in [("IL1R1","AT2L",False),("IL1B","Myeloid",True)]:
    m_=np.isin(fl,LINS[lin]) & ((epcam<2) if filt else True)
    gi=a.var_names.get_loc(gene)
    for step in STEPS:
        diffs={}
        for p in patients:
            sp=pt==p
            nrm=np.where(m_&sp&(stg=="Normal"))[0]; pre=np.where(m_&sp&np.isin(stg,STEPS[step]))[0]
            if len(nrm)<MIN or len(pre)<MIN: continue
            def l(ix):
                tot=X[ix].sum(); return np.log1p(X[ix][:,gi].sum()/tot*1e6)
            diffs[p]=l(pre)-l(nrm)
        res[f"{step}|{gene}_{lin}_dlog1pCPM"]=dict(summ(list(diffs.values())),patients=list(diffs))
json.dump(res,open(OUT+"B_il1b_power.json","w"),indent=1,ensure_ascii=False)
for k,v in res.items():
    if k!="plan_sha256": print(k,"n=%d"%v["n"],"mean=%.3f"%v["mean_diff"],"CI=[%.3f,%.3f]"%tuple(v["ci95"]),"up95=%.3f"%v["upper95_one_sided"],"MDE80=%.3f"%v["MDE80_one_sided"])
print("DONE")
