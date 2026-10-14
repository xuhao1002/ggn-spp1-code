"""docs/预设方案_缺氧H1H2_2026-10-12.md H2（Visium，处理同 scripts/80）。环境 ggn_spatial。"""
import glob, os, re, json, hashlib, warnings
import numpy as np, pandas as pd, scipy.sparse as sp, scanpy as sc
from scipy import stats
warnings.filterwarnings("ignore")
P=os.path.expanduser("~/全部数据库/GGN_target"); OUT=f"{P}/results/hypoxia/"
sha=hashlib.sha256(open(f"{P}/docs/预设方案_缺氧H1H2_2026-10-12.md","rb").read()).hexdigest(); assert sha.startswith("d1839e14"),sha
G=json.load(open(OUT+"gene_sets_used.json"))
CTRL=["CD68","C1QA","C1QB","C1QC","MS4A6A","CD14"]; OFFS=[(0,2),(0,-2),(1,1),(1,-1),(-1,1),(-1,-1)]
def adjacency(rc):
    idx={(r,c):i for i,(r,c) in enumerate(rc)}; rows=[];cols=[]
    for i,(r,c) in enumerate(rc):
        rows.append(i);cols.append(i)
        for dr,dc in OFFS:
            j=idx.get((r+dr,c+dc))
            if j is not None: rows.append(i);cols.append(j)
    A=sp.csr_matrix((np.ones(len(rows)),(rows,cols)),shape=(len(rc),len(rc))); d=np.asarray(A.sum(1)).ravel(); return sp.diags(1/d)@A
def spear(a,b): return np.nan if (np.std(a)==0 or np.std(b)==0) else float(stats.spearmanr(a,b)[0])
rows=[]
for d in sorted(glob.glob(P+"/raw/ST_GSE307534/extracted/GSM*")):
    base=os.path.basename(d); m=re.match(r"GSM\d+_(P\d+)_([A-Za-z]+)(-\d+)?$",base); pat,stage=m.group(1),m.group(2)
    sub=glob.glob(d+"/*/filtered_feature_bc_matrix")[0]
    ad=sc.read_10x_mtx(sub,var_names="gene_symbols",make_unique=True)
    pos=pd.read_csv(os.path.dirname(sub)+"/spatial/tissue_positions.csv").set_index("barcode")
    pos=pos.loc[pos.index.intersection(ad.obs_names)]; pos=pos[pos.in_tissue==1]
    ad=ad[pos.index].copy(); rc=pos[["array_row","array_col"]].values.astype(int)
    sc.pp.filter_cells(ad,min_counts=1); rc=pos.loc[ad.obs_names,["array_row","array_col"]].values.astype(int)
    r={"sample":base,"patient":pat,"stage":stage,"n_spots":int(ad.n_obs)}
    if ad.n_obs<500: r["eligible"]=False; rows.append(r); continue
    sc.pp.normalize_total(ad,target_sum=1e4); sc.pp.log1p(ad)
    W=adjacency(rc); X=lambda g: np.asarray(ad[:,g].X.todense()).ravel()
    spp1_s=W@X("SPP1")
    for nm,gs in [("HIF",G["HIF14_used"]),("HALL",G["HALLMARK_used"])]:
        gs=[g for g in gs if g in ad.var_names]
        sc.tl.score_genes(ad,gs,ctrl_size=50,n_bins=25,random_state=0,use_raw=False,score_name=nm)
        h=ad.obs[nm].values.astype(float); hs=W@h
        ctrl=[spear(W@X(g),hs) for g in CTRL if g in ad.var_names]
        r[f"dH_{nm}"]=spear(spp1_s,hs)-float(np.nanmean(ctrl)); r[f"mean_{nm}"]=float(h.mean())
    r["eligible"]=bool((spp1_s>0).sum()>=100); rows.append(r); print(base,r.get("dH_HIF"),flush=True)
S=pd.DataFrame(rows); S.to_csv(OUT+"H2_sections.tsv",sep="\t",index=False)
S=S[S.eligible.fillna(False).astype(bool)]
S["grp"]=S.stage.map({"AAH":"Precursor","AIS":"Precursor","MIA":"Precursor","LUAD":"LUAD","Normal":"Normal"})
def onesided(x):
    x=np.asarray(x,float); x=x[~np.isnan(x)]
    return {"n":int(len(x)),"mean":float(x.mean()),"n_pos":int((x>0).sum()),"p_one_sided":float(stats.wilcoxon(x,alternative="greater").pvalue)} if len(x)>=5 else {"n":int(len(x))}
res={"plan_sha256":sha}
for nm in ["HIF","HALL"]:
    res[f"H2a_{nm}"]=onesided(S.groupby("patient")[f"dH_{nm}"].mean().values)
    for g in ["Precursor","LUAD"]: res[f"H2a_{nm}|{g}"]=onesided(S[S.grp==g].groupby("patient")[f"dH_{nm}"].mean().values)
    pre=S[S.grp=="Precursor"].groupby("patient")[f"mean_{nm}"].mean(); lu=S[S.grp=="LUAD"].groupby("patient")[f"mean_{nm}"].mean()
    pr=pd.concat([pre,lu],axis=1,keys=["pre","lu"]).dropna(); res[f"H2b_{nm}"]=onesided((pr.lu-pr.pre).values)
json.dump(res,open(OUT+"H2_results.json","w"),indent=1,ensure_ascii=False)
for k,v in res.items():
    if k!="plan_sha256": print(k,v)
print("DONE")
