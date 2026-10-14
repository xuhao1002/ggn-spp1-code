"""拟合 cell2location 参考回归模型（方案第2节：基因选择 cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12；250 epoch，lr=0.002）。"""
import time, json, warnings
import numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp
warnings.filterwarnings("ignore")
import torch
from cell2location.utils.filtering import filter_genes
from cell2location.models import RegressionModel
t0=time.time()
def log(*a): print(f"[{time.time()-t0:7.0f}s]",*a,flush=True)
OUT="/Users/dudu/全部数据库/GGN_target/results/30_spatial/deconv/ref/"
torch.manual_seed(0); np.random.seed(0)
ref=ad.read_h5ad(OUT+"reference_subsample_counts.h5ad")
ref.X=sp.csr_matrix(ref.X)
log("ref",ref.shape,"torch threads",torch.get_num_threads())
sel=filter_genes(ref,cell_count_cutoff=5,cell_percentage_cutoff2=0.03,nonz_mean_cutoff=1.12)
sel=np.array(sel) if not isinstance(sel,np.ndarray) else sel
if sel.dtype==bool: ref=ref[:,sel].copy()
else: ref=ref[:,ref.var_names.isin(sel)].copy()
log("genes after filter",ref.shape[1])
RegressionModel.setup_anndata(adata=ref,batch_key="sample",labels_key="label")
mod=RegressionModel(ref)
mod.train(max_epochs=250,batch_size=2500,train_size=1,lr=0.002,accelerator="cpu")
log("trained; final ELBO",float(mod.history["elbo_train"].iloc[-1,0]))
ref=mod.export_posterior(ref,sample_kwargs={"num_samples":1000,"batch_size":2500,"accelerator":"cpu"})
fn=list(ref.uns["mod"]["factor_names"])
sig=ref.varm["means_per_cluster_mu_fg"][[f"means_per_cluster_mu_fg_{i}" for i in fn]].copy(); sig.columns=fn
sig.to_csv(OUT+"signatures_per_label.csv")
mod.history["elbo_train"].to_csv(OUT+"elbo_history.csv")
# 方案：签名相关 >0.9 则合并并重拟合——这里只计算并报告
c=np.corrcoef(np.log1p(sig.values.T))
cd=pd.DataFrame(c,index=fn,columns=fn)
cd.to_csv(OUT+"signature_correlations.csv")
for a,b in [("KAC-high","CAS-high"),("Mac-SPP1","Mac-IL1B"),("Mac-SPP1","Mac-other"),("KAC-high","AT2L-other"),("CAS-high","AT2L-other")]:
    log("corr(log1p sig)",a,b,round(cd.loc[a,b],3))
log("DONE step B")
