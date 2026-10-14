"""第三部分 §7（图 7）：SPP1+ 生态位基因集在 7 个 I 期队列中的预后。严格按 docs/第三部分_空间分析预设方案.md §7、§10。环境 ggn_target。
实现细节（方案未明写，运行前写定）：
- Wilcoxon 用 scanpy rank_genes_groups（Mac-SPP1 vs Mac-other，log1p(CP10k)，BH）；检出比例 = 计数>0 的细胞比例。
- 反向来源过滤：第一阶段 ratio_peng.tsv 中 med_Myeloid / epi ≥0.2（epi=0 视为通过）；不在该表中的基因剔除。
- 平台可测 = 在 7 个队列的基因×样本矩阵中存在（≥5/7）。
- 排序后取前 30；SPP1 强制纳入：若不在前 30，则取前 29 + SPP1（方案"包含 SPP1 本身"）。
- 髓系标志基因与第一阶段微环境第 4 层相同（scripts/15 的 MARK['Myeloid']）；基因集与标志基因的重叠只报告，不剔除。
"""
import matplotlib; matplotlib.use("Agg")
import os, sys, json, importlib.util, time
import numpy as np, pandas as pd, anndata as ad, scanpy as sc, scipy.sparse as sp
from scipy.stats import norm
from statsmodels.duration.hazard_regression import PHReg
t0=time.time()
def log(*a): print(f"[{time.time()-t0:7.0f}s]",*a,flush=True)
P=os.path.expanduser("~/全部数据库/GGN_target"); OUT=f"{P}/results/30_spatial/clinical/"; os.makedirs(OUT,exist_ok=True)
sys.path.insert(0,f"{P}/scripts")
spec=importlib.util.spec_from_file_location("l4",f"{P}/scripts/09_layer4_clinical.py"); l4=importlib.util.module_from_spec(spec); spec.loader.exec_module(l4)
from prolif_genes import TIROSH
MYE="CD68 C1QA C1QB MARCO APOE CD14 VCAN S100A8 CLEC9A CD1C FCER1A CLEC10A LILRA4 IL3RA CLEC4C CSF3R FCGR3B CXCR2".split()

# 1. 队列（先读，用于平台可测）
cohorts={}
for n in ["GSE31210","GSE50081","GSE30219","GSE37745","GSE68465","GSE72094","TCGA"]:
    cohorts[n]=l4.tcga() if n=="TCGA" else l4.cohort_def(n); log("cohort",n,len(cohorts[n][0]))
cov=pd.Series(0,dtype=int)
for n,c in cohorts.items(): cov=cov.add(pd.Series(1,index=c[1].index),fill_value=0)

# 2. 基因集
L=pd.read_csv(f"{P}/results/30_spatial/deconv/ref/reference_labels.tsv.gz",sep="\t",index_col=0)
mac=L[L.keep & L.label.isin(["Mac-SPP1","Mac-other"])]
a=ad.read_h5ad(f"{P}/results/01_snRNA/snRNA_final.h5ad")
a=a[mac.index].copy(); X=a.layers["counts"].tocsr().astype(np.float32)
tot=np.asarray(X.sum(1)).ravel(); X=sp.diags(1e4/np.maximum(tot,1))@X; X.data=np.log1p(X.data)
b=ad.AnnData(X=X.tocsr(),obs=pd.DataFrame({"grp":pd.Categorical(mac.label.values)},index=mac.index),var=pd.DataFrame(index=a.var_names))
det=pd.Series(np.asarray((a.layers["counts"][(mac.label=="Mac-SPP1").values]>0).mean(0)).ravel(),index=a.var_names); del a
sc.tl.rank_genes_groups(b,"grp",groups=["Mac-SPP1"],reference="Mac-other",method="wilcoxon",corr_method="benjamini-hochberg",pts=False)
de=sc.get.rank_genes_groups_df(b,"Mac-SPP1").set_index("names")
rat=pd.read_csv(f"{P}/results/04_source_filter/ratio_peng.tsv",sep="\t",index_col=0)
src=(rat["med_Myeloid"]/rat["epi"].replace(0,np.nan)).fillna(np.inf)
de["det_SPP1grp"]=det.reindex(de.index); de["src_ratio"]=src.reindex(de.index); de["n_cohorts_measured"]=cov.reindex(de.index).fillna(0).astype(int)
ok=(de.pvals_adj<0.05)&(de.det_SPP1grp>=0.3)&(de.src_ratio>=0.2)&(de.n_cohorts_measured>=5)&(de.logfoldchanges>0)
cand=de[ok].sort_values("logfoldchanges",ascending=False)
top=list(cand.index[:30])
if "SPP1" not in top: top=list(cand.drop(index="SPP1",errors="ignore").index[:29])+["SPP1"]
de.to_csv(OUT+"MacSPP1_vs_Macother_DE.tsv.gz",sep="\t")
pd.Series(top,name="gene").to_csv(OUT+"SPP1_niche_geneset.tsv",sep="\t",index=False)
log("gene set",len(top),top); log("overlap with myeloid markers",sorted(set(top)&set(MYE)))

# 3. Cox
def zscore_rows(m): return m.sub(m.mean(1),axis=0).div(m.std(1)+1e-9,axis=0)
def score(g,genes):
    pres=[x for x in genes if x in g.index]; s=zscore_rows(g.loc[pres].astype(float)).mean(0); return (s-s.mean())/(s.std()+1e-9), len(pres)
def cox(d,x,extra):
    X0=d[["age","male","stage"]].copy(); X0["age"]=(X0.age-X0.age.mean())/(X0.age.std()+1e-9); X0=X0.loc[:,X0.std()>0]
    for k,v in extra.items(): X0[k]=v.reindex(d.index).values
    X=pd.concat([x.reindex(d.index).rename("s"),X0],axis=1).values
    r=PHReg(d.time.values,X,status=d.event.values,ties="breslow").fit(disp=False); return float(r.params[0]),float(r.bse[0])
def meta(lh,se):
    lh,se=np.array(lh),np.array(se); w=1/se**2; fe=(w*lh).sum()/w.sum(); Q=(w*(lh-fe)**2).sum(); k=len(lh)
    tau2=max(0,(Q-(k-1))/(w.sum()-(w**2).sum()/w.sum())); w2=1/(se**2+tau2); mu=(w2*lh).sum()/w2.sum(); s=np.sqrt(1/w2.sum())
    return {"HR":float(np.exp(mu)),"lo":float(np.exp(mu-1.96*s)),"hi":float(np.exp(mu+1.96*s)),"p_one_sided":float(norm.sf(mu/s)),"I2":float(max(0,(Q-(k-1))/Q)) if Q>0 else 0.0,"tau2":float(tau2)}
models={"main":(top,["mye"]),"sens_a_plus_prolif":(top,["mye","prolif"]),"sens_b_without_SPP1":([x for x in top if x!="SPP1"],["mye"]),"sens_c_SPP1_only":(["SPP1"],["mye"]),"unadjusted_for_myeloid":(top,[])}
res={}; per=[]
for mname,(gs,adj) in models.items():
    lh=[];se=[]
    for n,(d,g,plat,ep,_) in cohorts.items():
        s,np_=score(g,gs); extra={}
        if "mye" in adj: extra["mye"]=score(g,MYE)[0]
        if "prolif" in adj: extra["prolif"]=score(g,TIROSH)[0]
        b_,s_=cox(d,s,extra); lh.append(b_); se.append(s_)
        per.append({"model":mname,"cohort":n,"endpoint":ep,"n":len(d),"events":int(d.event.sum()),"genes_present":np_,"HR":float(np.exp(b_)),"lo":float(np.exp(b_-1.96*s_)),"hi":float(np.exp(b_+1.96*s_)),"p_one_sided":float(norm.sf(b_/s_))})
    res[mname]=meta(lh,se); log(mname,res[mname])
pd.DataFrame(per).to_csv(OUT+"per_cohort.tsv",sep="\t",index=False)
res["main_supported"]=bool(res["main"]["p_one_sided"]<0.05 and res["main"]["HR"]>1)
res["gene_set"]=top; res["overlap_myeloid_markers"]=sorted(set(top)&set(MYE))
json.dump(res,open(OUT+"results.json","w"),ensure_ascii=False,indent=1)
log("DONE; main supported:",res["main_supported"])
