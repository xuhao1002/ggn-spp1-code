"""方案 §12.2：可区分性模拟。用未合并 26 标签签名，对 3000 个由留出细胞合成的 spot 做 cell2location 映射。环境 ggn_spatial。"""
import matplotlib; matplotlib.use("Agg")
import sys, time, json, warnings
sys.path.insert(0,"/Users/dudu/全部数据库/GGN_target/scripts")
import numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad
warnings.filterwarnings("ignore")
from c2l_common import R, read_st, st_samples, load_sig, run_c2l
t0=time.time()
def log(*a): print(f"[{time.time()-t0:7.0f}s]",*a,flush=True)
REF=R+"results/30_spatial/deconv/ref/"; OUT=R+"results/30_spatial/deconv/identifiability/"
import os; os.makedirs(OUT,exist_ok=True)
rng=np.random.default_rng(0)
sig=load_sig(REF+"signatures_per_label.csv"); labels=list(sig.columns); assert len(labels)==26
L=pd.read_csv(REF+"reference_labels.tsv.gz",sep="\t",index_col=0)
sub=ad.read_h5ad(REF+"reference_subsample_counts.h5ad",backed="r"); in_ref=set(sub.obs_names); sub.file.close()
pool=L[L.keep & ~L.index.isin(in_ref)]
log("held-out pool",len(pool),pool.label.value_counts().to_dict())
FOCAL=["KAC-high","CAS-high","Mac-SPP1","Mac-IL1B"]
w=np.array([4.0 if l in FOCAL else 1.0 for l in labels])
npool=pool.label.value_counts()
w=np.array([0.0 if npool.get(l,0)==0 else x for l,x in zip(labels,w)]); w/=w.sum()  # §12.5
log("labels with no held-out cells (weight 0):",[l for l in labels if npool.get(l,0)==0])
# depth distribution from first 5 ST samples (shared genes)
samp=st_samples()[:5]; tots=[]
genes=None
for b,d,p,s in samp:
    st=read_st(d); g=st.var_names.intersection(sig.index); g=g[~g.str.upper().str.startswith("MT-")]
    tots.append(np.asarray(st[:,g].X.sum(1)).ravel()); log("depth from",b,len(tots[-1]))
tots=np.concatenate(tots); tots=tots[tots>0]
# composition
N=3000; comp=np.zeros((N,len(labels)),dtype=int); cell_lists=[]
pools={l:pool.index[pool.label==l].values for l in labels}
for i in range(N):
    k=rng.choice([2,3,4]); chosen=rng.choice(len(labels),size=k,replace=False,p=w)
    n=1+rng.poisson(5); props=rng.dirichlet(np.ones(k)); cnt=rng.multinomial(n,props)
    cells=[]
    for j,c in zip(chosen,cnt):
        comp[i,j]+=c
        if c>0: cells.extend(rng.choice(pools[labels[j]],size=c,replace=True).tolist())
    cell_lists.append(cells)
need=sorted(set(c for cl in cell_lists for c in cl)); log("unique cells needed",len(need))
a=ad.read_h5ad(R+"results/01_snRNA/snRNA_final.h5ad")
gi=a.var_names.get_indexer(sig.index); assert (gi>=0).all()
ci=pd.Series(np.arange(a.n_obs),index=a.obs_names)
X=a.layers["counts"][ci.loc[need].values][:,gi].tocsr().astype(np.float32); del a
row=pd.Series(np.arange(len(need)),index=need)
rows=[]
for i,cl in enumerate(cell_lists):
    v=np.asarray(X[row.loc[cl].values].sum(0)).ravel()
    tgt=rng.choice(tots); tot=v.sum()
    if tot>tgt: v=rng.binomial(v.astype(np.int64),tgt/tot).astype(np.float32)
    rows.append(v)
syn=ad.AnnData(X=sp.csr_matrix(np.vstack(rows)),var=pd.DataFrame(index=sig.index),obs=pd.DataFrame(index=[f"syn{i}" for i in range(N)]))
truth=pd.DataFrame(comp,index=syn.obs_names,columns=labels); truth.to_csv(OUT+"truth_counts.tsv.gz",sep="\t")
log("synthetic built; median UMI",float(np.median(np.asarray(syn.X.sum(1)))))
q,qc=run_c2l(syn,sig,ncells=6,alpha=20,epochs=30000); log("mapped",qc)
q.to_csv(OUT+"q05_abundance.tsv.gz",sep="\t")
truth=truth.loc[q.index]
est=q.div(q.sum(1),axis=0); tru=truth.div(truth.sum(1),axis=0)
def r(x,y): return float(np.corrcoef(x,y)[0,1])
res={"qc":qc,"r_own":{l:r(est[l],tru[l]) for l in labels}}
def pair(A,B):
    ro_a=r(est[A],tru[A]); ro_b=r(est[B],tru[B]); xab=r(est[A],tru[B]); xba=r(est[B],tru[A])
    ok=(ro_a>=0.6) and (ro_b>=0.6) and (ro_a-xab>=0.2) and (ro_b-xba>=0.2)
    return {"r_own_A":ro_a,"r_own_B":ro_b,"r_cross_A_to_B":xab,"r_cross_B_to_A":xba,"pass":bool(ok)}
res["pairs"]={"KAC-high|CAS-high":pair("KAC-high","CAS-high"),"Mac-SPP1|Mac-IL1B":pair("Mac-SPP1","Mac-IL1B"),
  "control: Mac-SPP1|Mac-other":pair("Mac-SPP1","Mac-other"),"control: KAC-high|AT2L-other":pair("KAC-high","AT2L-other")}
json.dump(res,open(OUT+"results.json","w"),ensure_ascii=False,indent=1)
log("DONE"); print(json.dumps(res["pairs"],ensure_ascii=False,indent=1))
