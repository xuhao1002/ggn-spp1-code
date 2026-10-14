"""按 docs/预设方案_SPP1表达来源诊断_2026-10-11.md 运行（SHA256 见 results/30_spatial/diagnostic_prereg_sha256.txt）。"""
import json, time, warnings
import numpy as np, pandas as pd, anndata as ad
from scipy import stats
warnings.filterwarnings("ignore")
t0=time.time()
def log(*a): print(f"[{time.time()-t0:5.0f}s]",*a,flush=True)
R="/Users/dudu/全部数据库/GGN_target/results/"
OUT=R+"30_spatial/exploratory/"
a=ad.read_h5ad(R+"01_snRNA/snRNA_final.h5ad")
obs=a.obs[["patient","stage","sample","final_label","total_counts"]].copy()
for c in ["patient","stage","sample","final_label"]: obs[c]=obs[c].astype(str)
j=a.var_names.get_loc("SPP1")
cnt=a.layers["counts"].tocsc()[:,j].toarray().ravel().astype(np.float32)
del a
log("loaded",len(obs),"SPP1 nonzero cells",int((cnt>0).sum()))
lab=pd.read_csv(R+"30_spatial/deconv/ref/reference_labels.tsv.gz",sep="\t",index_col=0)["label"]
fl=obs["final_label"].astype(str)
grp=pd.Series("Other",index=obs.index)
for k in ["KAC-high","CAS-high","AT2L-other"]:
    grp[lab.reindex(obs.index)==k]=k
grp[fl=="Macrophage"]="Macrophage"
coarse={"Monocyte":"Myeloid-other","cDC":"Myeloid-other","pDC":"Myeloid-other","Neutrophil":"Myeloid-other","Proliferating myeloid":"Myeloid-other",
 "AT1":"Epithelium-other","Ciliated":"Epithelium-other","Basal":"Epithelium-other","AT2-like (patient-specific)":"Epithelium-other",
 "Fibroblast":"Stroma","Endothelial":"Stroma","SMC/pericyte":"Stroma","Lymphatic EC":"Stroma",
 "T cell":"T cell","NK":"Lymphoid-other","B cell":"Lymphoid-other","Plasma":"Lymphoid-other","Proliferating T/NK":"Lymphoid-other","Proliferating B/Plasma":"Lymphoid-other","Mast":"Mast"}
m=(grp=="Other")
grp[m]=fl[m].map(coarse).fillna("Other")
# 保留 Ciliated 单独用于环境参照
grp[fl=="Ciliated"]="Ciliated"
obs["grp"]=grp.values
obs["spp1"]=cnt; obs["lg"]=np.log1p(cnt/obs["total_counts"].values*1e4); obs["det"]=(cnt>0).astype(float)
obs["group"]=obs.stage.map({"AAH":"precursor","AIS":"precursor","MIA":"precursor","LUAD":"LUAD","Normal":"Normal"})
log("group sizes",obs.grp.value_counts().to_dict())

g=obs.groupby(["sample","patient","group","grp"],observed=True).agg(n=("lg","size"),mean_lg=("lg","mean"),det=("det","mean"),spp1_sum=("spp1","sum")).reset_index()
tot=obs.groupby("sample",observed=True)["spp1"].sum().rename("sample_spp1_total")
g=g.merge(tot,on="sample"); g["share"]=g.spp1_sum/g.sample_spp1_total.replace(0,np.nan)
g.to_csv(OUT+"spp1_source_by_sample_group.tsv",sep="\t",index=False)

def patient_mat(metric,grp_name,gname):
    x=g[(g.grp==grp_name)&(g.group==gname)&(g.n>=50)]
    return x.groupby("patient")[metric].mean()
res={}
# D2 份额（样本中位数，按阶段）
sh={}
for gname in ["LUAD","precursor","Normal"]:
    s=g[g.group==gname]
    piv=s.pivot_table(index="sample",columns="grp",values="share",aggfunc="mean").fillna(0)
    sh[gname]={c:round(float(piv[c].median()),4) for c in piv.columns}
    sh[gname]["_n_samples"]=int(len(piv))
    # 上皮总份额
    ep=piv[[c for c in ["KAC-high","CAS-high","AT2L-other","Epithelium-other","Ciliated"] if c in piv.columns]].sum(axis=1)
    sh[gname]["_epithelium_total_median"]=round(float(ep.median()),4)
res["D2_share_median_by_stage"]=sh
# 检出率与表达（样本×组，≥50 细胞）
det={}
for gname in ["LUAD","precursor"]:
    det[gname]={c:{"mean_lg_median":round(float(g[(g.grp==c)&(g.group==gname)&(g.n>=50)].mean_lg.median()),4) if len(g[(g.grp==c)&(g.group==gname)&(g.n>=50)]) else None,
                   "det_median":round(float(g[(g.grp==c)&(g.group==gname)&(g.n>=50)].det.median()),4) if len(g[(g.grp==c)&(g.group==gname)&(g.n>=50)]) else None,
                   "n_samples":int(len(g[(g.grp==c)&(g.group==gname)&(g.n>=50)]))}
                for c in ["Macrophage","KAC-high","CAS-high","AT2L-other","T cell","Ciliated","Epithelium-other","Stroma","Myeloid-other"]}
res["descriptive_by_group"]=det
# D1, D3
def paired(g1,g2,gname):
    a=patient_mat("mean_lg",g1,gname); b=patient_mat("mean_lg",g2,gname)
    d=(a-b).dropna()
    if len(d)<5: return {"n":int(len(d)),"p":None}
    return {"n":int(len(d)),"median_diff":float(d.median()),"n_first_higher":int((d>0).sum()),"p_two_sided":float(stats.wilcoxon(d).pvalue)}
tests={"D1_LUAD_Mac_minus_KAC":paired("Macrophage","KAC-high","LUAD"),"D1_precursor_Mac_minus_KAC":paired("Macrophage","KAC-high","precursor"),
       "D3_LUAD_KAC_minus_AT2Lother":paired("KAC-high","AT2L-other","LUAD"),"D3_precursor_KAC_minus_AT2Lother":paired("KAC-high","AT2L-other","precursor")}
ps=[(k,v["p_two_sided"]) for k,v in tests.items() if v.get("p_two_sided") is not None]
order=sorted(ps,key=lambda x:x[1]);k=len(order);prev=0
for i,(n,p) in enumerate(order):
    prev=max(prev,min(1,(k-i)*p)); tests[n]["p_holm_4"]=prev
res["tests"]=tests
# 预设解读指标
luad=sh["LUAD"]
res["verdict_inputs"]={"LUAD_macrophage_share_median":luad.get("Macrophage"),"LUAD_epithelium_total_share_median":luad["_epithelium_total_median"],
 "LUAD_KAC_high_share_median":luad.get("KAC-high"),
 "ambient_det_T_cell_vs_Mac_LUAD":[det["LUAD"]["T cell"]["det_median"],det["LUAD"]["Ciliated"]["det_median"],det["LUAD"]["Macrophage"]["det_median"]]}
json.dump(res,open(OUT+"spp1_source_results.json","w"),ensure_ascii=False,indent=1)
log("DONE"); print(json.dumps(res,ensure_ascii=False,indent=1))
