"""docs/预设方案_Xenium独立验证_2026-10-12.md（SHA 818813d5）：Takano 2024 Xenium X1–X3。环境 ggn_target。"""
import os, glob, json, hashlib, re
import numpy as np, pandas as pd, scanpy as sc, scipy.sparse as sp
from scipy import stats
from scipy.spatial import cKDTree
P=os.path.expanduser("~/全部数据库/GGN_target"); D=f"{P}/raw/Xenium_Takano2024/"; OUT=f"{P}/results/xenium/"; os.makedirs(OUT,exist_ok=True)
sha=hashlib.sha256(open(f"{P}/docs/预设方案_Xenium独立验证_2026-10-12.md","rb").read()).hexdigest(); assert sha.startswith("818813d5"),sha
SETS={"Macrophage":["CD68","CD163","MARCO","MSR1","C1QC","MPEG1","MS4A4A"],"Epithelial":["KRT8","NAPSA","SFTPC","SFTPB","NKX2-1","CDH1","SCGB3A2"],
      "Fibroblast":["COL1A1","PDGFRA","FAP","POSTN"],"Endothelial":["PECAM1","KDR","PLVAP"],"T/NK":["CD3E","CD8A","GZMA","LCK"],"B/Plasma":["CD19","MZB1","BANK1"]}
FFPE={"LUAD_No2","LUAD_No3"}
rows=[]
for d in sorted(glob.glob(D+"output-*")):
    sid=re.search(r"__(LUAD_No\d+|TSU-\d+)__",d).group(1)
    grp="invasive" if sid.startswith("LUAD") else "AIS/MIA"
    a=sc.read_10x_h5(d+"/cell_feature_matrix.h5"); a.var_names_make_unique()
    cells=pd.read_csv(d+"/cells.csv.gz").set_index("cell_id"); cells.index=cells.index.astype(str)
    a=a[a.obs_names.isin(cells.index)].copy(); xy=cells.loc[a.obs_names,["x_centroid","y_centroid"]].values
    X=sp.csr_matrix(a.X)
    sc_={k:np.asarray(X[:,[a.var_names.get_loc(g) for g in v if g in a.var_names]].sum(1)).ravel() for k,v in SETS.items()}
    M=np.column_stack([sc_[k] for k in SETS]); names=np.array(list(SETS))
    o=np.argsort(-M,axis=1); top=M[np.arange(len(M)),o[:,0]]; sec=M[np.arange(len(M)),o[:,1]]
    ctype=np.where((top>=2)&(top>sec),names[o[:,0]],"Unassigned")
    spp1=np.asarray(X[:,a.var_names.get_loc("SPP1")].todense()).ravel()
    tot=spp1.sum()
    r={"section":sid,"group":grp,"fixation":"FFPE" if sid in FFPE else "FF","n_cells":int(len(ctype))}
    for k in list(SETS)+["Unassigned"]: r[f"n_{k}"]=int((ctype==k).sum())
    r["SPP1_total_counts"]=float(tot)
    r["X1_frac_SPP1_in_macrophage"]=float(spp1[ctype=="Macrophage"].sum()/tot) if tot>0 else np.nan
    r["X1_frac_SPP1_in_epithelial"]=float(spp1[ctype=="Epithelial"].sum()/tot) if tot>0 else np.nan
    r["X1_frac_SPP1_in_other"]=1-r["X1_frac_SPP1_in_macrophage"]-r["X1_frac_SPP1_in_epithelial"] if tot>0 else np.nan
    mac=ctype=="Macrophage"; pos=mac&(spp1>=2); neg=mac&(spp1==0)
    r["n_SPP1pos_mac"]=int(pos.sum()); r["n_SPP1neg_mac"]=int(neg.sum())
    r["X2_SPP1pos_frac_of_macrophages"]=float(pos.sum()/mac.sum()) if mac.sum()>0 else np.nan
    epi=ctype=="Epithelial"
    if pos.sum()>=20 and neg.sum()>=20 and epi.sum()>0:
        t=cKDTree(xy[epi]); dp=t.query(xy[pos])[0]; dn=t.query(xy[neg])[0]
        r["X3_median_dist_pos_um"]=float(np.median(dp)); r["X3_median_dist_neg_um"]=float(np.median(dn)); r["X3_diff_neg_minus_pos"]=float(np.median(dn)-np.median(dp)); r["X3_eligible"]=True
    else: r["X3_eligible"]=False
    rows.append(r); print(sid,grp,r["n_cells"],"mac",r["n_Macrophage"],"SPP1+mac",r["n_SPP1pos_mac"],"X1mac=%.2f"%r["X1_frac_SPP1_in_macrophage"],flush=True)
T=pd.DataFrame(rows); T.to_csv(OUT+"xenium_per_section.tsv",sep="\t",index=False)
res={"plan_sha256":sha,"n_sections":len(T)}
for g,s in T.groupby("group"):
    res[f"X1_{g}"]={k:float(s[k].median()) for k in ["X1_frac_SPP1_in_macrophage","X1_frac_SPP1_in_epithelial","X1_frac_SPP1_in_other"]}
    res[f"X2_{g}_median"]=float(s["X2_SPP1pos_frac_of_macrophages"].median())
inv=T[T.group=="invasive"]["X2_SPP1pos_frac_of_macrophages"]; pre=T[T.group=="AIS/MIA"]["X2_SPP1pos_frac_of_macrophages"]
res["X2_MW_p_one_sided"]=float(stats.mannwhitneyu(inv,pre,alternative="greater").pvalue)
invFF=T[(T.group=="invasive")&(T.fixation=="FF")]["X2_SPP1pos_frac_of_macrophages"]
res["X2_FFonly_MW_p_one_sided"]=float(stats.mannwhitneyu(invFF,pre,alternative="greater").pvalue); res["X2_FFonly_n"]=[int(len(invFF)),int(len(pre))]
E=T[T.X3_eligible]
res["X3_n_eligible"]=int(len(E))
if len(E)>=5:
    res["X3_wilcoxon_p_one_sided"]=float(stats.wilcoxon(E.X3_diff_neg_minus_pos,alternative="greater").pvalue)
    res["X3_n_positive"]=int((E.X3_diff_neg_minus_pos>0).sum()); res["X3_median_diff_um"]=float(E.X3_diff_neg_minus_pos.median())
    for g,s in E.groupby("group"): res[f"X3_{g}"]={"n":int(len(s)),"n_positive":int((s.X3_diff_neg_minus_pos>0).sum()),"median_diff_um":float(s.X3_diff_neg_minus_pos.median())}
json.dump(res,open(OUT+"xenium_results.json","w"),indent=1,ensure_ascii=False)
print(json.dumps(res,indent=1,ensure_ascii=False)); print("DONE")
