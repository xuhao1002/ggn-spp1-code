"""按 docs/预设方案_空间基因层面探索性分析_2026-10-11.md 运行（SHA256 见 results/30_spatial/exploratory_prereg_sha256.txt）。
全部为探索性。环境 ggn_spatial。"""
import glob, os, re, json, time, warnings
import numpy as np, pandas as pd, scipy.sparse as sp
from scipy import stats
import scanpy as sc
warnings.filterwarnings("ignore")
t0=time.time()
def log(*a): print(f"[{time.time()-t0:6.0f}s]",*a,flush=True)
R="/Users/dudu/全部数据库/GGN_target/"
OUT=R+"results/30_spatial/exploratory/"
CTRL=["CD68","C1QA","C1QB","C1QC","MS4A6A","CD14"]
EPI=["EPCAM","NAPSA","SFTPB","SFTPC","SFTPA1","KRT19"]
kac_all=pd.read_csv(R+"results/20_stage2/peng_TableS3_KAC_human.txt",sep="\t",header=None)[0].astype(str).tolist()
OFFS=[(0,2),(0,-2),(1,1),(1,-1),(-1,1),(-1,-1)]

def adjacency(rc):
    idx={(r,c):i for i,(r,c) in enumerate(rc)}
    rows=[];cols=[]
    for i,(r,c) in enumerate(rc):
        rows.append(i);cols.append(i)
        for dr,dc in OFFS:
            j=idx.get((r+dr,c+dc))
            if j is not None: rows.append(i);cols.append(j)
    A=sp.csr_matrix((np.ones(len(rows)),(rows,cols)),shape=(len(rc),len(rc)))
    d=np.asarray(A.sum(1)).ravel()
    return sp.diags(1/d)@A

def zs(v):
    s=v.std(); return (v-v.mean())/s if s>0 else v*0

def spearman(a,b):
    if np.std(a)==0 or np.std(b)==0: return np.nan
    return float(stats.spearmanr(a,b)[0])

def analyse(adata,name):
    sc.pp.filter_cells(adata,min_counts=1)
    n=adata.n_obs
    out={"sample":name,"n_spots":int(n)}
    if n<500:
        out["eligible"]=False; out["reason"]="spots<500"; return out
    sc.pp.normalize_total(adata,target_sum=1e4); sc.pp.log1p(adata)
    genes=set(adata.var_names)
    kac=[g for g in kac_all if g in genes]
    sc.tl.score_genes(adata,kac,ctrl_size=50,n_bins=25,random_state=0,use_raw=False,score_name="KAC")
    X=lambda g: np.asarray(adata[:,g].X.todense()).ravel() if sp.issparse(adata.X) else np.asarray(adata[:,g].X).ravel()
    epi=[g for g in EPI if g in genes]
    ei=np.mean([zs(X(g)) for g in epi],axis=0)
    kacv=adata.obs["KAC"].values.astype(float)
    A=np.column_stack([np.ones(n),ei]); beta=np.linalg.lstsq(A,kacv,rcond=None)[0]
    kres=kacv-A@beta
    W=adjacency(adata.obsm["rc"])
    sm=lambda v: W@v
    kres_s=sm(kres)
    spp1=X("SPP1"); il1b=X("IL1B"); il1r1=X("IL1R1")
    ctrl=[g for g in CTRL if g in genes]
    spp1_s=sm(spp1); il1b_s=sm(il1b); il1r1_s=sm(il1r1)
    ctrl_v={g:X(g) for g in ctrl}; ctrl_s={g:sm(v) for g,v in ctrl_v.items()}
    out.update(n_kac_genes=len(kac),n_epi_genes=len(epi),n_ctrl_genes=len(ctrl),
        spp1_nonzero_smoothed=int((spp1_s>0).sum()),il1b_nonzero_smoothed=int((il1b_s>0).sum()))
    R_spp1=spearman(spp1_s,kres_s); R_ctrl=[spearman(ctrl_s[g],kres_s) for g in ctrl]
    out["R_SPP1_KACres"]=R_spp1; out["R_ctrl_KACres_mean"]=float(np.nanmean(R_ctrl)); out["dR1"]=R_spp1-float(np.nanmean(R_ctrl))
    R_il=spearman(il1b_s,il1r1_s); R_ctrl3=[spearman(ctrl_s[g],il1r1_s) for g in ctrl]
    out["R_IL1B_IL1R1"]=R_il; out["R_ctrl_IL1R1_mean"]=float(np.nanmean(R_ctrl3)); out["dR3"]=R_il-float(np.nanmean(R_ctrl3))
    out["SPP1_mean"]=float(spp1.mean()); out["ctrl_mean"]=float(np.mean([v.mean() for v in ctrl_v.values()]))
    out["M"]=out["SPP1_mean"]/out["ctrl_mean"] if out["ctrl_mean"]>0 else np.nan
    out["eligible_E1E2"]=bool(out["spp1_nonzero_smoothed"]>=100)
    out["eligible_E3"]=bool(out["il1b_nonzero_smoothed"]>=100)
    out["eligible"]=True
    return out

# ---------- GSE307534 ----------
rows=[]
for d in sorted(glob.glob(R+"raw/ST_GSE307534/extracted/GSM*")):
    base=os.path.basename(d); m=re.match(r"GSM\d+_(P\d+)_([A-Za-z]+)(-\d+)?$",base)
    pat,stage=m.group(1),m.group(2)
    sub=glob.glob(d+"/*/filtered_feature_bc_matrix")[0]
    ad_=sc.read_10x_mtx(sub,var_names="gene_symbols",make_unique=True)
    pos=pd.read_csv(os.path.dirname(sub)+"/spatial/tissue_positions.csv").set_index("barcode")
    pos=pos.loc[pos.index.intersection(ad_.obs_names)]; pos=pos[pos.in_tissue==1]
    ad_=ad_[pos.index].copy(); ad_.obsm["rc"]=pos[["array_row","array_col"]].values.astype(int)
    o=analyse(ad_,base); o.update(patient=pat,stage=stage); rows.append(o); log(base,o.get("n_spots"),o.get("dR1"))
S=pd.DataFrame(rows)
S["group"]=S.stage.map({"AAH":"precursor","AIS":"precursor","MIA":"precursor","LUAD":"LUAD","Normal":"Normal"})
S.to_csv(OUT+"sample_level_GSE307534.tsv",sep="\t",index=False)

res={"n_samples":int(len(S)),"n_eligible":int(S.eligible.fillna(False).sum()),
     "excluded":S.loc[~S.eligible.fillna(False).astype(bool),["sample","reason"]].values.tolist() if "reason" in S else []}
E12=S[S.eligible.fillna(False).astype(bool)&S.get("eligible_E1E2",pd.Series(False,index=S.index)).fillna(False).astype(bool)]
E3s=S[S.eligible.fillna(False).astype(bool)&S.get("eligible_E3",pd.Series(False,index=S.index)).fillna(False).astype(bool)&(S.group=="precursor")]
res["n_E1E2_eligible_samples"]=int(len(E12)); res["n_E3_eligible_samples"]=int(len(E3s))
def wilc(x,alt="greater"):
    x=np.asarray(x,float);x=x[~np.isnan(x)]
    if len(x)<5 or np.all(x==0): return {"n":int(len(x)),"p":np.nan}
    return {"n":int(len(x)),"median":float(np.median(x)),"mean":float(np.mean(x)),"n_pos":int((x>0).sum()),"p_one_sided":float(stats.wilcoxon(x,alternative=alt).pvalue)}
# E1
pm=E12.groupby("patient")["dR1"].mean()
res["E1_all_main"]=wilc(pm.values)
res["E1_precursor"]=wilc(E12[E12.group=="precursor"].groupby("patient")["dR1"].mean().values)
res["E1_LUAD"]=wilc(E12[E12.group=="LUAD"].groupby("patient")["dR1"].mean().values)
res["E1_all_excl_Normal_sens"]=wilc(E12[E12.group!="Normal"].groupby("patient")["dR1"].mean().values)
res["E1_mean_R_SPP1_vs_ctrl_by_patient"]={"R_SPP1_mean":float(E12.groupby("patient")["R_SPP1_KACres"].mean().mean()),"R_ctrl_mean":float(E12.groupby("patient")["R_ctrl_KACres_mean"].mean().mean())}
# E2
pre=E12[E12.group=="precursor"].groupby("patient")[["M","SPP1_mean"]].mean(); lu=E12[E12.group=="LUAD"].groupby("patient")[["M","SPP1_mean"]].mean()
pairs=pre.join(lu,lsuffix="_pre",rsuffix="_lu",how="inner").dropna()
for v in ["M","SPP1_mean"]:
    d=pairs[v+"_lu"]-pairs[v+"_pre"]
    res[f"E2_{v}"]=wilc(d.values)
res["E2_n_pairs"]=int(len(pairs))
pairs.to_csv(OUT+"E2_patient_pairs.tsv",sep="\t")
# E3
res["E3_precursor"]=wilc(E3s.groupby("patient")["dR3"].mean().values)
# Holm
ps=[res["E1_all_main"].get("p_one_sided"),res["E2_M"].get("p_one_sided"),res["E3_precursor"].get("p_one_sided")]
names=["E1","E2(M)","E3"]
ok=[(n,p) for n,p in zip(names,ps) if p==p]
order=sorted(ok,key=lambda x:x[1]); k=len(order); holm={}; prev=0
for i,(n,p) in enumerate(order):
    adj=min(1,(k-i)*p); prev=max(prev,adj); holm[n]=prev
res["Holm_three_main"]=holm
# ---------- GSE189487 ----------
stg={"TD1":3,"TD2":3,"TD3":2,"TD6":2,"TD5":1,"TD8":1}
rows=[]
for mt in sorted(glob.glob(R+"raw/ST_GSE189487/*_matrix.mtx.gz")):
    pref=mt.replace("_matrix.mtx.gz",""); td=re.search(r"(TD\d+)",pref).group(1)
    import scipy.io, gzip
    M=scipy.io.mmread(mt).T.tocsr()
    bc=pd.read_csv(pref+"_barcodes.tsv.gz",header=None)[0].values
    ft=pd.read_csv(pref+"_features.tsv.gz",sep="\t",header=None)
    ad_=sc.AnnData(M); ad_.obs_names=bc; ad_.var_names=pd.Index(ft[1].astype(str).values); ad_.var_names_make_unique()
    pos=pd.read_csv(pref+"_tissue_positions_list.csv.gz",header=None,names=["barcode","in_tissue","array_row","array_col","pr","pc"]).set_index("barcode")
    pos=pos.loc[pos.index.intersection(ad_.obs_names)]; pos=pos[pos.in_tissue==1]
    ad_=ad_[pos.index].copy(); ad_.obsm["rc"]=pos[["array_row","array_col"]].values.astype(int)
    o=analyse(ad_,td); o["stage_order"]=stg.get(td); rows.append(o); log("GSE189487",td,o.get("n_spots"),o.get("dR1"))
G=pd.DataFrame(rows); G.to_csv(OUT+"sample_level_GSE189487.tsv",sep="\t",index=False)
res["GSE189487_direction_only"]=G[["sample","stage_order","n_spots","dR1","dR3","M","SPP1_mean"]].round(4).to_dict("records")
if G["M"].notna().sum()>=4:
    r,_=stats.spearmanr(G["stage_order"],G["M"]); res["GSE189487_E2_spearman_M_vs_stage"]=float(r)
json.dump(res,open(OUT+"results.json","w"),ensure_ascii=False,indent=1,default=float)
log("DONE"); print(json.dumps(res,ensure_ascii=False,indent=1,default=float))
