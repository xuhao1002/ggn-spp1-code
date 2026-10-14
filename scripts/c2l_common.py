"""cell2location 空间映射的共用函数（方案 §2、§12.1 参数）。"""
import time, glob, os, re
import numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad
R="/Users/dudu/全部数据库/GGN_target/"

def read_st(d):
    import scanpy as sc
    sub=glob.glob(d+"/*/filtered_feature_bc_matrix")[0]
    a=sc.read_10x_mtx(sub,var_names="gene_symbols",make_unique=True)
    pos=pd.read_csv(os.path.dirname(sub)+"/spatial/tissue_positions.csv").set_index("barcode")
    pos=pos.loc[pos.index.intersection(a.obs_names)]; pos=pos[pos.in_tissue==1]
    a=a[pos.index].copy(); a.obsm["rc"]=pos[["array_row","array_col"]].values.astype(int)
    a.X=sp.csr_matrix(a.X).astype(np.float32)
    return a

def st_samples():
    out=[]
    for d in sorted(glob.glob(R+"raw/ST_GSE307534/extracted/GSM*")):
        b=os.path.basename(d); m=re.match(r"GSM\d+_(P\d+)_([A-Za-z]+)(-\d+)?$",b)
        out.append((b,d,m.group(1),m.group(2)))
    return out

def load_sig(path):
    return pd.read_csv(path,index_col=0)

def run_c2l(st,sig,ncells=6,alpha=20,epochs=30000,seed=0):
    import torch
    from cell2location.models import Cell2location
    torch.manual_seed(seed); np.random.seed(seed)
    st=st[:,~st.var_names.str.upper().str.startswith("MT-")].copy()
    shared=st.var_names.intersection(sig.index)
    st=st[:,shared].copy(); s=sig.loc[shared].copy()
    keep=np.asarray(st.X.sum(1)).ravel()>0
    st=st[keep].copy()
    Cell2location.setup_anndata(adata=st)
    s=s.astype(np.float32)
    mod=Cell2location(st,cell_state_df=s,N_cells_per_location=ncells,detection_alpha=alpha)
    if os.environ.get("C2L_ACCEL","cpu")=="mps": mod.module=mod.module.float()  # MPS 不支持 float64
    t=time.time()
    ACC=os.environ.get("C2L_ACCEL","cpu")
    mod.train(max_epochs=epochs,batch_size=None,train_size=1,accelerator=ACC)
    st=mod.export_posterior(st,sample_kwargs={"num_samples":1000,"batch_size":mod.adata.n_obs,"accelerator":ACC})
    q=st.obsm["q05_cell_abundance_w_sf"].copy()
    q.columns=[c.replace("q05cell_abundance_w_sf_","") for c in q.columns]
    elbo=mod.history["elbo_train"].iloc[:,0].values
    return q, {"n_spots":int(st.n_obs),"n_genes":int(len(shared)),"train_sec":round(time.time()-t,1),
               "accelerator":ACC,"elbo_final":float(elbo[-1]),"elbo_at_90pct":float(elbo[int(len(elbo)*0.9)]),
               "umi_vs_abund_r":float(np.corrcoef(np.asarray(st.X.sum(1)).ravel(),q.sum(1).values)[0,1])}
