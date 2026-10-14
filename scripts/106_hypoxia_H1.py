"""docs/预设方案_缺氧H1H2_2026-10-12.md H1（snRNA）。环境 ggn_target。"""
import os, json, hashlib
import numpy as np, pandas as pd, scanpy as sc
from scipy import stats
P=os.path.expanduser("~/全部数据库/GGN_target"); OUT=f"{P}/results/hypoxia/"; os.makedirs(OUT,exist_ok=True)
sha=hashlib.sha256(open(f"{P}/docs/预设方案_缺氧H1H2_2026-10-12.md","rb").read()).hexdigest(); assert sha.startswith("d1839e14"),sha
niche=pd.read_csv(f"{P}/results/30_spatial/clinical/SPP1_niche_geneset.tsv",sep="\t").gene.tolist()
HIF=["CA9","VEGFA","SLC2A1","PGK1","LDHA","BNIP3","NDRG1","ADM","ANKRD37","P4HA1","ENO1","ALDOA","PDK1","EGLN3"]
import gseapy as gp
HALL=gp.get_library("MSigDB_Hallmark_2020")["Hypoxia"]
HIFc=[g for g in HIF if g not in niche]; HALLc=[g for g in HALL if g not in niche]
STATE=[g for g in niche if g!="SPP1" and g not in HIF and g not in HALL]
json.dump({"HIF14_used":HIFc,"HALLMARK_used":HALLc,"state_genes":STATE},open(OUT+"gene_sets_used.json","w"),indent=1)
print("HIF used",len(HIFc),"HALL used",len(HALLc),"state genes",len(STATE),STATE,flush=True)
a=sc.read_h5ad(f"{P}/results/01_snRNA/snRNA_final.h5ad")
o=a.obs; epcam=np.asarray(a[:,"EPCAM"].layers["counts"].todense()).ravel()
m=((o["final_label"].astype(str)=="Macrophage").values)&(epcam<2)
b=a[m].copy(); del a
spp1c=np.asarray(b[:,"SPP1"].layers["counts"].todense()).ravel()
for nm,gs in [("HIF",HIFc),("HALL",HALLc),("STATE",STATE)]:
    sc.tl.score_genes(b,[g for g in gs if g in b.var_names],score_name=nm,random_state=0)
ob=b.obs.copy(); ob["patient"]=ob["patient"].astype(str); ob["stage"]=ob["stage"].astype(str); ob["spp1"]=spp1c
ob["grp"]=np.where(ob.stage=="LUAD","LUAD",np.where(ob.stage.isin(["AAH","AIS","MIA"]),"Precursor","Other"))
ob=ob[ob.grp!="Other"]
def onesided(x):
    x=np.asarray(x,float); x=x[~np.isnan(x)]
    return {"n":int(len(x)),"mean":float(x.mean()),"n_pos":int((x>0).sum()),"p_one_sided":float(stats.wilcoxon(x,alternative="greater").pvalue)} if len(x)>=5 else {"n":int(len(x))}
res={"plan_sha256":sha}; units=[]
for (p,g),u in ob.groupby(["patient","grp"]):
    r={"patient":p,"grp":g,"n":len(u)}
    pos=u[u.spp1>0]; neg=u[u.spp1==0]
    for sc_ in ["HIF","HALL"]:
        r[f"H1a_i_{sc_}"]=float(pos[sc_].mean()-neg[sc_].mean()) if (len(pos)>=20 and len(neg)>=20) else np.nan
        q1,q3=u["STATE"].quantile([.25,.75]); hi=u[u.STATE>=q3]; lo=u[u.STATE<=q1]
        r[f"H1a_ii_{sc_}"]=float(hi[sc_].mean()-lo[sc_].mean()) if (len(hi)>=20 and len(lo)>=20) else np.nan
        r[f"mean_{sc_}"]=float(u[sc_].mean()) if len(u)>=50 else np.nan
    units.append(r)
U=pd.DataFrame(units); U.to_csv(OUT+"H1_units.tsv",sep="\t",index=False)
for k in ["H1a_i_HIF","H1a_i_HALL","H1a_ii_HIF","H1a_ii_HALL"]:
    res[k]=onesided(U.groupby("patient")[k].mean().values)
    for g in ["Precursor","LUAD"]: res[f"{k}|{g}"]=onesided(U[U.grp==g][k].values)
for sc_ in ["HIF","HALL"]:
    w=U.pivot(index="patient",columns="grp",values=f"mean_{sc_}").dropna()
    res[f"H1b_{sc_}"]=onesided((w["LUAD"]-w["Precursor"]).values)
json.dump(res,open(OUT+"H1_results.json","w"),indent=1,ensure_ascii=False)
for k,v in res.items():
    if k!="plan_sha256": print(k,v)
print("DONE")
