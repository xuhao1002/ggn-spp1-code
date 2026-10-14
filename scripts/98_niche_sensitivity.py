"""预设方案_敏感性与蛋白验证_2026-10-12.md 第 A 节：SPP1+ 生态位预后敏感性（A1 逐一剔除队列、A2 加 CAF 评分、A3 去非髓系来源基因）。环境 ggn_target。"""
import matplotlib; matplotlib.use("Agg")
import os, sys, json, hashlib, importlib.util
import numpy as np, pandas as pd
from scipy.stats import norm
from statsmodels.duration.hazard_regression import PHReg
P=os.path.expanduser("~/全部数据库/GGN_target"); OUT=f"{P}/results/sens_protein/"; os.makedirs(OUT,exist_ok=True)
sha=hashlib.sha256(open(f"{P}/docs/预设方案_敏感性与蛋白验证_2026-10-12.md","rb").read()).hexdigest(); assert sha.startswith("b15055f8"),sha
spec=importlib.util.spec_from_file_location("l4",f"{P}/scripts/09_layer4_clinical.py"); l4=importlib.util.module_from_spec(spec); spec.loader.exec_module(l4)
MYE="CD68 C1QA C1QB MARCO APOE CD14 VCAN S100A8 CLEC9A CD1C FCER1A CLEC10A LILRA4 IL3RA CLEC4C CSF3R FCGR3B CXCR2".split()
CAF=["FAP","POSTN","CTHRC1","COL11A1","COL10A1","MMP11","INHBA"]
top=pd.read_csv(f"{P}/results/30_spatial/clinical/SPP1_niche_geneset.tsv",sep="\t").gene.tolist(); assert len(top)==30
rat=pd.read_csv(f"{P}/results/04_source_filter/ratio_peng.tsv",sep="\t",index_col=0)
cols=["epi"]+[c for c in rat.columns if c.startswith("med_")]
best=rat[cols].idxmax(axis=1)
dropped=[g for g in top if g in best.index and best[g]!="med_Myeloid"]; missing=[g for g in top if g not in best.index]
top_mye=[g for g in top if g not in dropped]
print("A3 dropped (top lineage not myeloid):",{g:best[g] for g in dropped},"not in table:",missing,flush=True)
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
    return {"HR":float(np.exp(mu)),"lo":float(np.exp(mu-1.96*s)),"hi":float(np.exp(mu+1.96*s)),"p_one_sided":float(norm.sf(mu/s)),"I2":float(max(0,(Q-(k-1))/Q)) if Q>0 else 0.0,"k":k}
def run(genes,adj,names):
    lh=[];se=[];per={}
    for n in names:
        d,g=cohorts[n][0],cohorts[n][1]; s,_=score(g,genes); ex={}
        if "mye" in adj: ex["mye"]=score(g,MYE)[0]
        if "caf" in adj: ex["caf"]=score(g,CAF)[0]
        b,e=cox(d,s,ex); lh.append(b); se.append(e); per[n]=float(np.exp(b))
    out=meta(lh,se); out["per_cohort_HR"]=per; return out
allc=list(cohorts)
res={"plan_sha256":sha,"main_reproduced":run(top,["mye"],allc)}
res["A1_leave_one_out"]={f"drop_{n}":run(top,["mye"],[c for c in allc if c!=n]) for n in allc}
res["A2_plus_CAF"]=run(top,["mye","caf"],allc)
res["A3_myeloid_source_genes_only"]=dict(run(top_mye,["mye"],allc),genes_dropped=dropped,n_genes=len(top_mye))
json.dump(res,open(OUT+"A_niche_sensitivity.json","w"),indent=1,ensure_ascii=False)
print("main",res["main_reproduced"]["HR"],res["main_reproduced"]["p_one_sided"])
for k,v in res["A1_leave_one_out"].items(): print(k,round(v["HR"],3),round(v["lo"],3),round(v["hi"],3),"p1=%.4f"%v["p_one_sided"])
for k in ["A2_plus_CAF","A3_myeloid_source_genes_only"]: v=res[k]; print(k,round(v["HR"],3),round(v["lo"],3),round(v["hi"],3),"p1=%.4f"%v["p_one_sided"],"I2=%.2f"%v["I2"])
print("DONE")
