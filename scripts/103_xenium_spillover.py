"""docs/预设方案_预后稳健性与Xenium串扰校准_2026-10-12.md 第 S 节。细胞分型与 101 完全相同。环境 ggn_target。"""
import os, glob, json, hashlib, re
import numpy as np, pandas as pd, scanpy as sc, scipy.sparse as sp
from scipy import stats
from scipy.spatial import cKDTree
P=os.path.expanduser("~/全部数据库/GGN_target"); D=f"{P}/raw/Xenium_Takano2024/"; OUT=f"{P}/results/xenium/"; os.makedirs(OUT,exist_ok=True)
sha=hashlib.sha256(open(f"{P}/docs/预设方案_预后稳健性与Xenium串扰校准_2026-10-12.md","rb").read()).hexdigest(); assert sha.startswith("ec46b0c7"),sha
SETS={"Macrophage":["CD68","CD163","MARCO","MSR1","C1QC","MPEG1","MS4A4A"],"Epithelial":["KRT8","NAPSA","SFTPC","SFTPB","NKX2-1","CDH1","SCGB3A2"],
      "Fibroblast":["COL1A1","PDGFRA","FAP","POSTN"],"Endothelial":["PECAM1","KDR","PLVAP"],"T/NK":["CD3E","CD8A","GZMA","LCK"],"B/Plasma":["CD19","MZB1","BANK1"]}
REF=["CD68","C1QC","MARCO","MSR1"]
rows=[]
for d in sorted(glob.glob(D+"output-*")):
    sid=re.search(r"__(LUAD_No\d+|TSU-\d+)__",d).group(1); grp="invasive" if sid.startswith("LUAD") else "AIS/MIA"
    a=sc.read_10x_h5(d+"/cell_feature_matrix.h5"); a.var_names_make_unique()
    cells=pd.read_csv(d+"/cells.csv.gz").set_index("cell_id"); cells.index=cells.index.astype(str)
    a=a[a.obs_names.isin(cells.index)].copy(); xy=cells.loc[a.obs_names,["x_centroid","y_centroid"]].values
    X=sp.csr_matrix(a.X); gi=lambda g:a.var_names.get_loc(g)
    M=np.column_stack([np.asarray(X[:,[gi(g) for g in v if g in a.var_names]].sum(1)).ravel() for v in SETS.values()]); names=np.array(list(SETS))
    o=np.argsort(-M,axis=1); top=M[np.arange(len(M)),o[:,0]]; sec=M[np.arange(len(M)),o[:,1]]
    ctype=np.where((top>=2)&(top>sec),names[o[:,0]],"Unassigned")
    spp1=np.asarray(X[:,gi("SPP1")].todense()).ravel(); ref=np.asarray(X[:,[gi(g) for g in REF]].sum(1)).ravel()
    mac=ctype=="Macrophage"; epi=ctype=="Epithelial"
    L=ref[epi].sum()/ref[mac].sum(); exp_=L*spp1[mac].sum(); obs=spp1[epi].sum(); tot=spp1.sum()
    r={"section":sid,"group":grp,"leak_rate_L":float(L),"epi_SPP1_observed":float(obs),"epi_SPP1_expected_spill":float(exp_),
       "obs_over_exp":float(obs/exp_) if exp_>0 else np.nan,"epi_share_raw":float(obs/tot),"epi_share_corrected":float(max(0,obs-exp_)/tot),
       "mac_share_raw":float(spp1[mac].sum()/tot),"exceeds_2x":bool(obs>2*exp_)}
    pos=mac&(spp1>=2)
    if pos.sum()>=20:
        dist=cKDTree(xy[pos]).query(xy[epi])[0]; rho,p=stats.spearmanr(spp1[epi],dist); r["epi_SPP1_vs_dist_to_SPP1pos_mac_rho"]=float(rho); r["rho_p"]=float(p)
    rows.append(r); print(sid,grp,"L=%.3f obs/exp=%.1f epi_share raw=%.2f corrected=%.2f"%(L,r["obs_over_exp"],r["epi_share_raw"],r["epi_share_corrected"]),flush=True)
T=pd.DataFrame(rows); T.to_csv(OUT+"xenium_spillover_per_section.tsv",sep="\t",index=False)
res={"plan_sha256":sha}
for g,s in T.groupby("group"):
    res[g]={"median_epi_share_raw":float(s.epi_share_raw.median()),"median_epi_share_corrected":float(s.epi_share_corrected.median()),
            "median_obs_over_exp":float(s.obs_over_exp.median()),"n_exceeds_2x":int(s.exceeds_2x.sum()),"n":int(len(s))}
json.dump(res,open(OUT+"xenium_spillover_results.json","w"),indent=1,ensure_ascii=False); print(json.dumps(res,indent=1,ensure_ascii=False)); print("DONE")
