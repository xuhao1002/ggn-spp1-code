"""docs/预设方案_缺氧H1c_2026-10-12.md。环境 ggn_target。"""
import os, json, hashlib
import numpy as np, pandas as pd, scanpy as sc
from scipy import stats
P=os.path.expanduser("~/全部数据库/GGN_target"); OUT=f"{P}/results/hypoxia/"
sha=hashlib.sha256(open(f"{P}/docs/预设方案_缺氧H1c_2026-10-12.md","rb").read()).hexdigest(); assert sha.startswith(open(f"{P}/results/hypoxia_prereg_sha256.txt").read().split()[-2][:8]),sha
G=json.load(open(OUT+"gene_sets_used.json"))
a=sc.read_h5ad(f"{P}/results/01_snRNA/snRNA_final.h5ad")
epcam=np.asarray(a[:,"EPCAM"].layers["counts"].todense()).ravel()
b=a[((a.obs["final_label"].astype(str)=="Macrophage").values)&(epcam<2)].copy(); del a
spp1=np.asarray(b[:,"SPP1"].layers["counts"].todense()).ravel()
for nm,gs in [("HIF",G["HIF14_used"]),("HALL",G["HALLMARK_used"]),("STATE",G["state_genes"])]:
    sc.tl.score_genes(b,[g for g in gs if g in b.var_names],score_name=nm,random_state=0)
o=b.obs.copy(); o["patient"]=o.patient.astype(str); o["stage"]=o.stage.astype(str); o["spp1"]=spp1
o["grp"]=np.where(o.stage=="LUAD","LUAD",np.where(o.stage.isin(["AAH","AIS","MIA"]),"Precursor","Other")); o=o[o.grp!="Other"]
rows=[]
for (p,g),u in o.groupby(["patient","grp"]):
    r={"patient":p,"grp":g}; pos=u[u.spp1>0]; q3=u.STATE.quantile(.75); hi=u[u.STATE>=q3]
    for s in ["HIF","HALL"]:
        r[f"i_{s}"]=float(pos[s].mean()) if len(pos)>=10 else np.nan
        r[f"ii_{s}"]=float(hi[s].mean()) if len(hi)>=20 else np.nan
    rows.append(r)
U=pd.DataFrame(rows); res={"plan_sha256":sha}
for k in ["i_HIF","i_HALL","ii_HIF","ii_HALL"]:
    w=U.pivot(index="patient",columns="grp",values=k).dropna(); d=(w.LUAD-w.Precursor).values
    res["H1c_"+k]={"n":int(len(d)),"mean":float(d.mean()),"n_pos":int((d>0).sum()),"p_one_sided":float(stats.wilcoxon(d,alternative="greater").pvalue)} if len(d)>=5 else {"n":int(len(d))}
json.dump(res,open(OUT+"H1c_results.json","w"),indent=1); [print(k,v) for k,v in res.items()]; print("DONE")
