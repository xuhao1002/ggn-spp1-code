"""方案 §12.1：合并 KAC-high+CAS-high→KACCAS-high、Mac-SPP1+Mac-IL1B→Mac-act，重抽 ≤2000/标签（种子0），重拟合参考回归模型（参数同 §2）。环境 ggn_spatial。"""
import matplotlib; matplotlib.use("Agg")
import time, json, warnings
import numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad
warnings.filterwarnings("ignore")
import torch
from cell2location.utils.filtering import filter_genes
from cell2location.models import RegressionModel
t0=time.time()
def log(*a): print(f"[{time.time()-t0:7.0f}s]",*a,flush=True)
R="/Users/dudu/全部数据库/GGN_target/results/"
REF=R+"30_spatial/deconv/ref/"; OUT=R+"30_spatial/deconv/ref_merged/"
import os; os.makedirs(OUT,exist_ok=True)
torch.manual_seed(0); rng=np.random.default_rng(0)
L=pd.read_csv(REF+"reference_labels.tsv.gz",sep="\t",index_col=0)
L["label"]=L["label"].replace({"KAC-high":"KACCAS-high","CAS-high":"KACCAS-high","Mac-SPP1":"Mac-act","Mac-IL1B":"Mac-act"})
log("merged label counts (keep)",L[L.keep]["label"].value_counts().to_dict())
a=ad.read_h5ad(R+"01_snRNA/snRNA_final.h5ad"); assert (a.obs_names==L.index).all()
counts=a.layers["counts"].tocsr().astype(np.float32)
pos=pd.Series(np.arange(len(L)),index=L.index); idx=[]
for lab,g in L[L.keep].groupby("label"):
    p=pos.loc[g.index].values; idx.extend(rng.choice(p,min(2000,len(p)),replace=False).tolist())
idx=np.sort(np.array(idx))
ref=ad.AnnData(X=sp.csr_matrix(counts[idx]),obs=pd.DataFrame({"label":L["label"].values[idx],"sample":L["sample"].values[idx]},index=L.index[idx]),var=pd.DataFrame(index=a.var_names))
del a,counts
pd.Series(ref.obs_names).to_csv(OUT+"reference_subsample_cells.tsv.gz",index=False)
log("ref",ref.shape,ref.obs.label.value_counts().to_dict())
sel=filter_genes(ref,cell_count_cutoff=5,cell_percentage_cutoff2=0.03,nonz_mean_cutoff=1.12)
sel=np.array(sel)
ref=ref[:,sel].copy() if sel.dtype==bool else ref[:,ref.var_names.isin(sel)].copy()
log("genes after filter",ref.shape[1])
RegressionModel.setup_anndata(adata=ref,batch_key="sample",labels_key="label")
mod=RegressionModel(ref)
mod.train(max_epochs=250,batch_size=2500,train_size=1,lr=0.002,accelerator="cpu")
log("trained; final ELBO",float(mod.history["elbo_train"].iloc[-1,0]))
ref=mod.export_posterior(ref,sample_kwargs={"num_samples":1000,"batch_size":2500,"accelerator":"cpu"})
fn=list(ref.uns["mod"]["factor_names"])
sig=ref.varm["means_per_cluster_mu_fg"][[f"means_per_cluster_mu_fg_{i}" for i in fn]].copy(); sig.columns=fn
sig.to_csv(OUT+"signatures_per_label.csv"); mod.history["elbo_train"].to_csv(OUT+"elbo_history.csv")
c=pd.DataFrame(np.corrcoef(np.log1p(sig.values.T)),index=fn,columns=fn); c.to_csv(OUT+"signature_correlations.csv")
log("DONE merged reference;",len(fn),"labels")
