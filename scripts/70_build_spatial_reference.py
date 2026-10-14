"""空间去卷积参考构建。严格按 docs/第三部分_空间分析预设方案.md 第2节与第10节确认记录。
环境: ggn_spatial。输入: results/01_snRNA/snRNA_final.h5ad (counts 层), results/50_screen/AT2L_cell_state_assignments.tsv.gz。
"""
import json, time, hashlib, warnings
import numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad
warnings.filterwarnings("ignore")
t0=time.time()
def log(*a): print(f"[{time.time()-t0:7.0f}s]",*a,flush=True)
R="/Users/dudu/全部数据库/GGN_target/results/"
OUT=R+"30_spatial/deconv/ref/"
rng=np.random.default_rng(0)

a=ad.read_h5ad(R+"01_snRNA/snRNA_final.h5ad")
log("loaded",a.shape)
counts=a.layers["counts"].tocsr().astype(np.float32)
obs=a.obs.copy()
lab=obs["final_label"].astype(str).copy()

# 1) 上皮亚型：KAC-high / CAS-high / AT2L-other（阈值冻结，来自 state_thresholds_frozen.json，与方案第2节一致）
thr=json.load(open(R+"50_screen/state_thresholds_frozen.json"))["thresholds"]
log("thresholds",thr)
asg=pd.read_csv(R+"50_screen/AT2L_cell_state_assignments.tsv.gz",sep="\t").set_index("cell_id")
asg=asg[~asg.index.duplicated()]
inter=asg.index.intersection(obs.index)
log("AT2L assignment cells",len(asg),"in h5ad",len(inter))
asg=asg.loc[inter]
kr=asg["KAC"].rank(pct=True); cr=asg["CAS"].rank(pct=True)
khi=asg["KAC"]>thr["KAC"]; chi=asg["CAS"]>thr["CAS"]
new=pd.Series("AT2L-other",index=asg.index)
new[khi&~chi]="KAC-high"; new[chi&~khi]="CAS-high"
both=khi&chi
new[both&(kr>=cr)]="KAC-high"; new[both&(kr<cr)]="CAS-high"
log("AT2L subtype counts",new.value_counts().to_dict(),"both-high",int(both.sum()))
lab.loc[new.index]=new.values

# 2) 巨噬细胞亚型：SPP1+ / IL1B+ / Mac-other（Normal 髓系 SPP1 第90百分位；IL1B 用冻结 0.5904 并复核）
genes=a.var_names
tot=obs["total_counts"].values.astype(np.float32)
def lognorm(g):
    j=genes.get_loc(g); x=counts[:,j].toarray().ravel()
    return np.log1p(x/tot*1e4)
spp1=lognorm("SPP1"); il1b=lognorm("IL1B")
myeloid=["Macrophage","Monocyte","cDC","pDC","Neutrophil","Proliferating myeloid"]
is_norm=(obs["stage"]=="Normal").values
mye_norm=is_norm&obs["final_label"].isin(myeloid).values
p90_spp1=float(np.quantile(spp1[mye_norm],0.9)); log("Normal myeloid SPP1 p90",p90_spp1,"n",int(mye_norm.sum()))
mac_norm=is_norm&(obs["final_label"]=="Macrophage").values
p90_il1b=float(np.quantile(il1b[mac_norm],0.9)); log("Normal macrophage IL1B p90 (recomputed)",p90_il1b,"(frozen 0.5904)")
IL1B_THR=0.5904
mac=(obs["final_label"]=="Macrophage").values
sp_pos=(spp1>p90_spp1) if p90_spp1>0 else (spp1>0)
il_pos=(il1b>IL1B_THR)
drop=mac&sp_pos&il_pos
newlab=np.where(mac&sp_pos&~il_pos,"Mac-SPP1",np.where(mac&il_pos&~sp_pos,"Mac-IL1B",np.where(mac,"Mac-other",lab.values)))
lab=pd.Series(newlab,index=obs.index)
keep=~drop
log("Macrophage subtype counts",pd.Series(newlab[mac]).value_counts().to_dict(),"dropped double-positive",int(drop.sum()))

# 3) 最小细胞数 >=200
vc=lab[keep].value_counts(); log("label counts\n"+vc.to_string())
small=vc[vc<200].index.tolist(); log("labels <200 cells (excluded from reference):",small)
keep&=~lab.isin(small).values
lab_df=pd.DataFrame({"label":lab,"sample":obs["sample"],"keep":keep}); lab_df.to_csv(OUT+"reference_labels.tsv.gz",sep="\t")
json.dump({"p90_spp1_normal_myeloid":p90_spp1,"p90_il1b_recomputed":p90_il1b,"il1b_frozen_used":IL1B_THR,"thresholds":thr,
 "labels_kept_counts":vc[vc>=200].to_dict(),"labels_excluded_lt200":small,"dropped_double_positive_mac":int(drop.sum())},open(OUT+"reference_summary.json","w"),ensure_ascii=False,indent=1)

# 4) 每标签抽 <=2000 个细胞（种子0），保存参考 h5ad（counts）
idx=[]
for L,g in lab_df[lab_df.keep].groupby("label"):
    ii=np.where(obs.index.isin(g.index))[0] if False else None
pos=pd.Series(np.arange(len(obs)),index=obs.index)
for L,g in lab_df[lab_df.keep].groupby("label"):
    p=pos.loc[g.index].values
    idx.extend(rng.choice(p,min(2000,len(p)),replace=False).tolist())
idx=np.sort(np.array(idx))
ref=ad.AnnData(X=counts[idx],obs=pd.DataFrame({"label":lab.values[idx],"sample":obs["sample"].values[idx]},index=obs.index[idx]),var=pd.DataFrame(index=genes))
ref.write_h5ad(OUT+"reference_subsample_counts.h5ad")
log("reference subsample",ref.shape,ref.obs.label.value_counts().to_dict())
log("DONE step A (labels + subsample). Next: fit RegressionModel.")
