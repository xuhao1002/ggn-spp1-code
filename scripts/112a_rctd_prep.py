"""D：导出 RCTD 参考（不合并参考，同 cell2location 参考细胞）。环境 ggn_spatial。"""
import os, hashlib, numpy as np, pandas as pd, anndata as ad, scipy.io as sio, scipy.sparse as sp
P="/Users/dudu/全部数据库/GGN_target/"; O=P+"results/mech/rctd/ref/"; os.makedirs(O,exist_ok=True)
sha=hashlib.sha256(open(P+"docs/预设方案_调控互作RCTD_2026-10-13.md","rb").read()).hexdigest()
assert sha==open(P+"results/mech/prereg_sha256_D.txt").read().split()[0],sha
a=ad.read_h5ad(P+"results/30_spatial/deconv/ref/reference_subsample_counts.h5ad")
X=a.layers["counts"] if "counts" in a.layers else a.X
X=sp.csr_matrix(X); assert np.allclose(X.data,np.round(X.data)),"reference not integer counts"
if "label" in a.obs: lab=a.obs["label"].astype(str)
else:
    L=pd.read_csv(P+"results/30_spatial/deconv/ref/reference_labels.tsv.gz",sep="\t",index_col=0); lab=L.loc[a.obs_names,"label"].astype(str)
safe=lab.str.replace("/","_").str.replace(" ","_").str.replace("(","").str.replace(")","")
pd.DataFrame({"label":lab.unique(),"safe":safe.unique()}).to_csv(O+"label_map.tsv",sep="\t",index=False)
sio.mmwrite(O+"counts.mtx",X.T.tocoo().astype(np.int32))   # genes x cells
pd.Series(a.var_names).to_csv(O+"genes.tsv",index=False,header=False)
pd.DataFrame({"cell":a.obs_names,"type":safe.values}).to_csv(O+"cells.tsv",sep="\t",index=False)
print("ref",a.shape,safe.value_counts().to_dict()); print("DONE")
