"""docs/预设方案_补充检查T_2026-10-12.md T1、T2。分型同 101/103。环境 ggn_target。"""
import os, glob, json, hashlib, re
import numpy as np, pandas as pd, scanpy as sc, scipy.sparse as sp
from scipy import stats
from scipy.spatial import cKDTree
P=os.path.expanduser("~/全部数据库/GGN_target"); D=f"{P}/raw/Xenium_Takano2024/"; OUT=f"{P}/results/xenium/"
sha=hashlib.sha256(open(f"{P}/docs/预设方案_补充检查T_2026-10-12.md","rb").read()).hexdigest(); assert sha.startswith("b29834f7"),sha
SETS={"Macrophage":["CD68","CD163","MARCO","MSR1","C1QC","MPEG1","MS4A4A"],"Epithelial":["KRT8","NAPSA","SFTPC","SFTPB","NKX2-1","CDH1","SCGB3A2"],
      "Fibroblast":["COL1A1","PDGFRA","FAP","POSTN"],"Endothelial":["PECAM1","KDR","PLVAP"],"T/NK":["CD3E","CD8A","GZMA","LCK"],"B/Plasma":["CD19","MZB1","BANK1"]}
REF=["CD68","C1QC","MARCO","MSR1"]; rows=[]
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
    mac=ctype=="Macrophage"; epi=ctype=="Epithelial"; tot=spp1.sum()
    L=ref[epi].sum()/ref[mac].sum(); exp_=L*spp1[mac].sum()
    dist=cKDTree(xy[mac]).query(xy[epi])[0]; rho,p=stats.spearmanr(ref[epi],dist)
    r={"section":sid,"group":grp,"T1_rho_epiRef_vs_dist_to_mac":float(rho),"T1_p":float(p),
       "T1_percell_leak":float(ref[epi].mean()/ref[mac].mean()),"L_total":float(L),"n_epi":int(epi.sum()),"n_mac":int(mac.sum()),
       "T2_mac_share_raw":float(spp1[mac].sum()/tot),"T2_mac_share_corrected":float(min(1,(spp1[mac].sum()+exp_)/tot)),
       "T2_epi_share_raw":float(spp1[epi].sum()/tot),"T2_epi_share_corrected":float(max(0,spp1[epi].sum()-exp_)/tot)}
    rows.append(r); print(sid,grp,"rho=%.3f percell_leak=%.4f mac %.2f->%.2f epi %.2f->%.2f"%(rho,r["T1_percell_leak"],r["T2_mac_share_raw"],r["T2_mac_share_corrected"],r["T2_epi_share_raw"],r["T2_epi_share_corrected"]),flush=True)
T=pd.DataFrame(rows); T.to_csv(OUT+"xenium_checkT_per_section.tsv",sep="\t",index=False)
res={"plan_sha256":sha,"T1_n_negative_rho":int((T.T1_rho_epiRef_vs_dist_to_mac<0).sum()),"T1_n_negative_and_p05":int(((T.T1_rho_epiRef_vs_dist_to_mac<0)&(T.T1_p<0.05)).sum()),"n":int(len(T))}
for g,s in T.groupby("group"):
    res[g]={k:[float(s[k].min()),float(s[k].median()),float(s[k].max())] for k in ["T2_mac_share_raw","T2_mac_share_corrected","T2_epi_share_raw","T2_epi_share_corrected","T1_rho_epiRef_vs_dist_to_mac"]}
json.dump(res,open(OUT+"xenium_checkT_results.json","w"),indent=1,ensure_ascii=False); print(json.dumps(res,indent=1,ensure_ascii=False)); print("DONE")
