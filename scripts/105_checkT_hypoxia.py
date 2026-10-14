"""docs/预设方案_补充检查T_2026-10-12.md T3、T4。环境 ggn_target（需联网取 Enrichr 库）。"""
import matplotlib; matplotlib.use("Agg")
import os, sys, json, hashlib, importlib.util
import numpy as np, pandas as pd
from scipy import stats
from scipy.stats import norm
from statsmodels.duration.hazard_regression import PHReg
P=os.path.expanduser("~/全部数据库/GGN_target"); OUT=f"{P}/results/sens_protein/"
sha=hashlib.sha256(open(f"{P}/docs/预设方案_补充检查T_2026-10-12.md","rb").read()).hexdigest(); assert sha.startswith("b29834f7"),sha
spec=importlib.util.spec_from_file_location("l4",f"{P}/scripts/09_layer4_clinical.py"); l4=importlib.util.module_from_spec(spec); spec.loader.exec_module(l4)
MYE="CD68 C1QA C1QB MARCO APOE CD14 VCAN S100A8 CLEC9A CD1C FCER1A CLEC10A LILRA4 IL3RA CLEC4C CSF3R FCGR3B CXCR2".split()
HIF=["CA9","VEGFA","SLC2A1","PGK1","LDHA","BNIP3","NDRG1","ADM","ANKRD37","P4HA1","ENO1","ALDOA","PDK1","EGLN3"]
try:
    import gseapy as gp
    lib=gp.get_library("MSigDB_Hallmark_2020"); HALL=lib.get("Hypoxia"); assert HALL and len(HALL)>100
except Exception as e:
    HALL=None; print("HALLMARK_HYPOXIA 取不到：",e)
top=pd.read_csv(f"{P}/results/30_spatial/clinical/SPP1_niche_geneset.tsv",sep="\t").gene.tolist()
cohorts={n:(l4.tcga() if n=="TCGA" else l4.cohort_def(n)) for n in ["GSE31210","GSE50081","GSE30219","GSE37745","GSE68465","GSE72094","TCGA"]}
def zr(m): return m.sub(m.mean(1),axis=0).div(m.std(1)+1e-9,axis=0)
def score(g,genes):
    pres=[x for x in genes if x in g.index]; s=zr(g.loc[pres].astype(float)).mean(0); return (s-s.mean())/(s.std()+1e-9),len(pres)
def cox(d,x,extra):
    X0=d[["age","male","stage"]].copy(); X0["age"]=(X0.age-X0.age.mean())/(X0.age.std()+1e-9); X0=X0.loc[:,X0.std()>0]
    for k,v in extra.items(): X0[k]=v.reindex(d.index).values
    r=PHReg(d.time.values,pd.concat([x.reindex(d.index).rename("s"),X0],axis=1).values,status=d.event.values,ties="breslow").fit(disp=False)
    return float(r.params[0]),float(r.bse[0])
def meta(lh,se):
    lh,se=np.array(lh),np.array(se); w=1/se**2; fe=(w*lh).sum()/w.sum(); Q=(w*(lh-fe)**2).sum(); k=len(lh)
    tau2=max(0,(Q-(k-1))/(w.sum()-(w**2).sum()/w.sum())); w2=1/(se**2+tau2); mu=(w2*lh).sum()/w2.sum(); s=np.sqrt(1/w2.sum())
    pi=stats.t.ppf(0.975,k-2)*np.sqrt(tau2+s**2)
    return {"HR":float(np.exp(mu)),"lo":float(np.exp(mu-1.96*s)),"hi":float(np.exp(mu+1.96*s)),"p_one_sided":float(norm.sf(mu/s)),"p_two_sided":float(2*norm.sf(abs(mu/s))),"I2":float(max(0,(Q-(k-1))/Q)) if Q>0 else 0.0,"PI95":[float(np.exp(mu-pi)),float(np.exp(mu+pi))]}
res={"plan_sha256":sha,"HALLMARK_n_genes":len(HALL) if HALL else None}
rho={"HIF14":{},"HALLMARK":{}}; lh=[];se=[]; per={}
for n,(d,g,*_) in cohorts.items():
    s,_=score(g,top); h1,_=score(g,HIF); rho["HIF14"][n]=float(stats.spearmanr(s,h1)[0])
    if HALL:
        h2,np_=score(g,HALL); rho["HALLMARK"][n]=float(stats.spearmanr(s,h2)[0])
        b,e=cox(d,s,{"mye":score(g,MYE)[0],"hyp":h2}); lh.append(b); se.append(e); per[n]=float(np.exp(b))
res["T4_rho"]=rho; res["T4_median_rho"]={k:(float(np.median(list(v.values()))) if v else None) for k,v in rho.items()}
if HALL: res["T3_main_plus_HALLMARK"]=dict(meta(lh,se),per_cohort_HR=per)
else: res["T3_main_plus_HALLMARK"]="无法检验"
json.dump(res,open(OUT+"T_hypoxia.json","w"),indent=1,ensure_ascii=False); print(json.dumps(res,indent=1,ensure_ascii=False)); print("DONE")
