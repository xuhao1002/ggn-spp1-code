"""docs/预设方案_预后稳健性与Xenium串扰校准_2026-10-12.md 第 R 节。环境 ggn_target。"""
import matplotlib; matplotlib.use("Agg")
import os, sys, json, hashlib, importlib.util
import numpy as np, pandas as pd
from scipy import stats
from scipy.stats import norm
from statsmodels.duration.hazard_regression import PHReg
P=os.path.expanduser("~/全部数据库/GGN_target"); OUT=f"{P}/results/sens_protein/"; os.makedirs(OUT,exist_ok=True)
sha=hashlib.sha256(open(f"{P}/docs/预设方案_预后稳健性与Xenium串扰校准_2026-10-12.md","rb").read()).hexdigest(); assert sha.startswith("ec46b0c7"),sha
sys.path.insert(0,f"{P}/scripts")
spec=importlib.util.spec_from_file_location("l4",f"{P}/scripts/09_layer4_clinical.py"); l4=importlib.util.module_from_spec(spec); spec.loader.exec_module(l4)
from prolif_genes import TIROSH
MYE="CD68 C1QA C1QB MARCO APOE CD14 VCAN S100A8 CLEC9A CD1C FCER1A CLEC10A LILRA4 IL3RA CLEC4C CSF3R FCGR3B CXCR2".split()
CAF=["FAP","POSTN","CTHRC1","COL11A1","COL10A1","MMP11","INHBA"]; EPI=["EPCAM","KRT8","KRT18","KRT19","CDH1"]; IMM=["PTPRC"]
HYP=["CA9","VEGFA","SLC2A1","PGK1","LDHA","BNIP3","NDRG1","ADM","ANKRD37","P4HA1","ENO1","ALDOA","PDK1","EGLN3"]
top=pd.read_csv(f"{P}/results/30_spatial/clinical/SPP1_niche_geneset.tsv",sep="\t").gene.tolist(); assert len(top)==30
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
    tk=stats.t.ppf(0.975,k-2); pi=tk*np.sqrt(tau2+s**2)
    return {"HR":float(np.exp(mu)),"lo":float(np.exp(mu-1.96*s)),"hi":float(np.exp(mu+1.96*s)),"p_one_sided":float(norm.sf(mu/s)),"p_two_sided":float(2*norm.sf(abs(mu/s))),
            "I2":float(max(0,(Q-(k-1))/Q)) if Q>0 else 0.0,"tau2":float(tau2),"PI95":[float(np.exp(mu-pi)),float(np.exp(mu+pi))],"k":k}
COV={"mye":MYE,"epi":EPI,"imm":IMM,"hyp":HYP,"prolif":TIROSH,"caf":CAF}
def run(adj):
    lh=[];se=[];per={};cov_n={}
    for n,c in cohorts.items():
        d,g=c[0],c[1]; s,_=score(g,top); ex={}
        for a in adj:
            ex[a],cov_n[f"{n}:{a}"]=score(g,COV[a])
        b,e=cox(d,s,ex); lh.append(b); se.append(e); per[n]=float(np.exp(b))
    out=meta(lh,se); out["per_cohort_HR"]=per; out["genes_present_per_covariate"]=cov_n; return out
res={"plan_sha256":sha}
for name,adj in [("main",["mye"]),("R1_purity",["mye","epi","imm"]),("R2_hypoxia",["mye","hyp"]),("R3_all",["mye","epi","imm","hyp","prolif","caf"])]:
    res[name]=run(adj); v=res[name]
    print(name,"HR %.3f (%.3f-%.3f) p1=%.4f p2=%.4f I2=%.2f PI95=[%.2f,%.2f]"%(v["HR"],v["lo"],v["hi"],v["p_one_sided"],v["p_two_sided"],v["I2"],*v["PI95"]),flush=True)
json.dump(res,open(OUT+"R_prognosis_robustness.json","w"),indent=1,ensure_ascii=False); print("DONE")
