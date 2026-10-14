"""docs/预设方案_调控互作RCTD_2026-10-13.md  R 部分（TF 活性）。环境 ggn_target。"""
import os, json, hashlib
import numpy as np, pandas as pd, scanpy as sc, scipy.sparse as sp
from scipy import stats
P=os.path.expanduser("~/全部数据库/GGN_target"); OUT=f"{P}/results/mech/"; os.makedirs(OUT,exist_ok=True)
sha=hashlib.sha256(open(f"{P}/docs/预设方案_调控互作RCTD_2026-10-13.md","rb").read()).hexdigest()
assert sha==open(OUT+"prereg_sha256.txt").read().split()[0], "plan changed after lock: "+sha
print("plan sha ok",sha[:8],flush=True)

# ---- CollecTRI ----
import decoupler as dc
try: net=dc.op.collectri(organism="human")          # decoupler >=2
except AttributeError: net=dc.get_collectri(organism="human",split_complexes=False)  # 1.x
net=net[["source","target","weight"]].copy(); net.to_csv(OUT+"collectri_human.tsv",sep="\t",index=False)
net_sha=hashlib.sha256(open(OUT+"collectri_human.tsv","rb").read()).hexdigest()
print("CollecTRI",net.shape,"TFs",net.source.nunique(),flush=True)

G=json.load(open(f"{P}/results/hypoxia/gene_sets_used.json"))
HIF14=["CA9","VEGFA","SLC2A1","PGK1","LDHA","BNIP3","NDRG1","ADM","ANKRD37","P4HA1","ENO1","ALDOA","PDK1","EGLN3"]
import gseapy as gp
HALL=gp.get_library("MSigDB_Hallmark_2020")["Hypoxia"]

# ---- macrophages (same definition as H1) ----
a=sc.read_h5ad(f"{P}/results/01_snRNA/snRNA_final.h5ad")
epcam=np.asarray(a[:,"EPCAM"].layers["counts"].todense()).ravel()
b=a[((a.obs["final_label"].astype(str)=="Macrophage").values)&(epcam<2)].copy(); del a
C=b.layers["counts"]; C=sp.csr_matrix(C)
spp1=np.asarray(C[:,list(b.var_names).index("SPP1")].todense()).ravel()
sc.tl.score_genes(b,[g for g in G["state_genes"] if g in b.var_names],score_name="STATE",random_state=0)
o=b.obs.copy(); o["patient"]=o.patient.astype(str); o["stage"]=o.stage.astype(str); o["spp1"]=spp1
o["grp"]=np.where(o.stage=="LUAD","LUAD",np.where(o.stage.isin(["AAH","AIS","MIA"]),"Precursor","Other"))
keep=(o.grp!="Other").values
tot=np.asarray(C.sum(axis=0)).ravel(); det=tot>=10
genes=np.array(b.var_names)[det]; Cd=C[:,det]
print("macrophages",keep.sum(),"detected genes",len(genes),flush=True)

# ---- coverage QC ----
gset=set(genes); cov=net.groupby("source").target.agg(lambda t: (len(t), sum(x in gset for x in t)))
covdf=pd.DataFrame([(tf,n,k,k/n) for tf,(n,k) in cov.items()],columns=["TF","n_targets","n_detected","frac"])
covdf["eligible"]=(covdf.frac>=0.5)&(covdf.n_detected>=10); covdf.to_csv(OUT+"R_tf_coverage.tsv",sep="\t",index=False)
elig=covdf[covdf.eligible].TF.tolist(); hifcov=covdf[covdf.TF=="HIF1A"].iloc[0].to_dict()
print("eligible TFs",len(elig),"HIF1A coverage",hifcov,flush=True)

def build_W(netsub):
    ns=netsub[netsub.target.isin(gset)]
    tfs=sorted(ns.source.unique()); gi={g:i for i,g in enumerate(genes)}; ti={t:i for i,t in enumerate(tfs)}
    W=np.zeros((len(genes),len(tfs)))
    for s,t,w in ns.itertuples(index=False): W[gi[t],ti[s]]=w
    return W,tfs
def ulm(Y,W):
    """Y: samples x genes (log1p CPM); W: genes x TFs. t-value of slope in y ~ 1 + w."""
    n=W.shape[0]; Wc=W-W.mean(0); ssw=(Wc**2).sum(0)
    Yc=Y-Y.mean(1,keepdims=True)
    beta=(Yc@Wc)/ssw
    resid_ss=(Yc**2).sum(1,keepdims=True)-beta**2*ssw
    se=np.sqrt(np.maximum(resid_ss,1e-12)/(n-2)/ssw)
    return beta/se
def pb(idx):
    v=np.asarray(Cd[idx].sum(axis=0)).ravel(); return np.log1p(v/v.sum()*1e6)

net_e=net[net.source.isin(elig)]; W,tfs=build_W(net_e)
hifT=net[(net.source=="HIF1A")]; hifT_dc=hifT[~hifT.target.isin(set(HIF14)|set(HALL))]
n_dc=int(hifT_dc.target.isin(gset).sum())
Wdc,_=build_W(hifT_dc) if n_dc>=10 else (None,None)

ob=o[keep]; rows=[]; arr=np.arange(len(o))
for (p,g),u in ob.groupby(["patient","grp"]):
    ii=arr[o.index.get_indexer(u.index)]; sp1=u.spp1.values>0
    q1,q3=u.STATE.quantile([.25,.75]); hi=(u.STATE>=q3).values; lo=(u.STATE<=q1).values
    r={"patient":p,"grp":g,"n":len(u)}
    if len(u)>=50: r["all"]=pb(ii)
    if sp1.sum()>=20 and (~sp1).sum()>=20: r["pos"]=pb(ii[sp1]); r["neg"]=pb(ii[~sp1])
    if hi.sum()>=20 and lo.sum()>=20: r["hi"]=pb(ii[hi]); r["lo"]=pb(ii[lo])
    rows.append(r)

def act(key,Wm,tfi):
    out=[]
    for r in rows:
        out.append(ulm(r[key][None,:],Wm)[0][tfi] if key in r else np.nan)
    return np.array(out,dtype=object)
H=tfs.index("HIF1A") if "HIF1A" in tfs else None
U=pd.DataFrame([{k:r[k] for k in ["patient","grp","n"]} for r in rows])
def onesided(x):
    x=np.asarray(x,float); x=x[~np.isnan(x)]
    if len(x)<5: return {"n":int(len(x))}
    return {"n":int(len(x)),"mean":float(x.mean()),"median":float(np.median(x)),"n_pos":int((x>0).sum()),
            "p_one_sided":float(stats.wilcoxon(x,alternative="greater").pvalue),"p_two_sided":float(stats.wilcoxon(x).pvalue)}
res={"plan_sha256":sha,"collectri_sha256":net_sha,"n_tf_eligible":len(elig),"HIF1A_coverage":{k:(float(v) if isinstance(v,(np.floating,float,int,np.integer)) and not isinstance(v,bool) else v) for k,v in hifcov.items()},"HIF1A_decirc_targets_detected":n_dc}

# all-TF activity matrices per unit
def allact(key):
    M=np.full((len(rows),len(tfs)),np.nan)
    for i,r in enumerate(rows):
        if key in r: M[i]=ulm(r[key][None,:],W)[0]
    return M
Apos,Aneg,Ahi,Alo,Aall=[allact(k) for k in ["pos","neg","hi","lo","all"]]
if H is None:
    res["verdict"]="HIF1A not eligible (coverage) -> R not evaluable"
else:
    U["R1"]=Apos[:,H]-Aneg[:,H]; U["R3"]=Ahi[:,H]-Alo[:,H]; U["all_HIF1A"]=Aall[:,H]
    if Wdc is not None:
        U["R1_decirc"]=[ (ulm(r["pos"][None,:],Wdc)[0][0]-ulm(r["neg"][None,:],Wdc)[0][0]) if "pos" in r else np.nan for r in rows]
    res["R1"]=onesided(U.groupby("patient").R1.mean().values)
    w=U.pivot(index="patient",columns="grp",values="all_HIF1A").dropna()
    res["R2"]=onesided((w["LUAD"]-w["Precursor"]).values) if {"LUAD","Precursor"}<=set(w.columns) else {"n":0}
    ps=[res[k].get("p_one_sided",np.nan) for k in ["R1","R2"]]
    order=np.argsort(ps); m=len(ps); holm=[None]*m; run=0
    for rank,i in enumerate(order): run=max(run,min(1,(m-rank)*ps[i])); holm[i]=run
    res["R1"]["p_holm"]=holm[0]; res["R2"]["p_holm"]=holm[1]
    res["R1_decirc"]=onesided(U.groupby("patient").R1_decirc.mean().values) if "R1_decirc" in U else {"note":f"not evaluable: {n_dc} targets"}
    res["R3"]=onesided(U.groupby("patient").R3.mean().values)
    for g in ["Precursor","LUAD"]: res[f"R1|{g}"]=onesided(U[U.grp==g].R1.values)
    # all-TF ranking (R1-type)
    D=pd.DataFrame(Apos-Aneg,columns=tfs); D["patient"]=U.patient.values; Dp=D.groupby("patient").mean()
    rk=[]
    for t in tfs:
        x=Dp[t].dropna().values
        if len(x)>=5: rk.append((t,float(np.median(x)),int((x>0).sum()),len(x),float(stats.wilcoxon(x,alternative="greater").pvalue)))
    R=pd.DataFrame(rk,columns=["TF","median_diff","n_pos","n","p_one_sided"]).sort_values("median_diff",ascending=False).reset_index(drop=True)
    from statsmodels.stats.multitest import multipletests
    R["q_BH"]=multipletests(R.p_one_sided,method="fdr_bh")[1]; R["rank"]=np.arange(1,len(R)+1)
    R.to_csv(OUT+"R_tf_rank.tsv",sep="\t",index=False)
    hr=int(R[R.TF=="HIF1A"]["rank"].iloc[0]); res["HIF1A_rank"]={"rank":hr,"of":len(R),"percentile_top":hr/len(R)}
    res["top20"]=R.head(20)[["TF","median_diff","n_pos","n","q_BH"]].to_dict("records")
U.to_csv(OUT+"R_units.tsv",sep="\t",index=False)
# optional cross-check with decoupler's own ULM on first pseudobulk
try:
    r0=next(r for r in rows if "all" in r); mat=pd.DataFrame(r0["all"][None,:],columns=genes,index=["s0"])
    try: est=dc.mt.ulm(data=mat,net=net_e,tmin=10); est=est[0] if isinstance(est,tuple) else est
    except AttributeError: est=dc.run_ulm(mat=mat,net=net_e,min_n=10,verbose=False)[0]
    mine=pd.Series(ulm(r0["all"][None,:],W)[0],index=tfs); common=[t for t in mine.index if t in est.columns]
    res["crosscheck_decoupler_r"]=float(np.corrcoef(mine[common],est.loc["s0",common].astype(float))[0,1])
except Exception as e: res["crosscheck_decoupler_r"]=f"skipped: {type(e).__name__}: {e}"
json.dump(res,open(OUT+"R_results.json","w"),indent=1,ensure_ascii=False,default=float)
for k,v in res.items():
    if k not in ("plan_sha256","collectri_sha256","top20"): print(k,v)
print("top10",[d["TF"] for d in res.get("top20",[])[:10]]); print("DONE")
