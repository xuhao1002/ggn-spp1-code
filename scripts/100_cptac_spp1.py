"""预设方案_敏感性与蛋白验证_2026-10-12.md 第 C 节：CPTAC-LUAD SPP1 蛋白层验证。"""
import re, sys, json, hashlib, numpy as np, pandas as pd
from scipy import stats
D=sys.argv[1]; PLAN=sys.argv[2]; OUT=sys.argv[3]
sha=hashlib.sha256(open(PLAN,"rb").read()).hexdigest(); assert sha.startswith("b15055f8"), sha
m=pd.read_csv(D+"aliquot_map_labeled151.tsv",sep="\t")
MAC=["CD68","CD163","C1QA","CD14","MARCO"]
def load(f):
    d=pd.read_csv(D+f,sep="\t",index_col=0)
    d=d.drop([i for i in ["Mean","Median","StdDev"] if i in d.index])
    keep={}
    for c in d.columns:
        if "Unshared" in c: continue
        mm=re.search(r"CPT\d{10}",c)
        if mm: keep[mm.group(0)]=c
    d=d[list(keep.values())]; d.columns=list(keep.keys()); return d
def zrow(x): return (x-x.mean())/x.std()
res={"plan_sha256":sha}
for name,f in [("Discovery","PDC000153_CPTAC3_Lung_Adeno_Carcinoma_Proteome.tmt10.tsv"),("Confirmatory","PDC000489_CPTAC3_LUAD_Confirmatory_Proteome.tmt11.tsv")]:
    d=load(f); mm=m[m.aliquot.isin(d.columns)]
    T=mm[mm.sample_type=="Primary Tumor"].drop_duplicates("case").set_index("case").aliquot
    N=mm[mm.sample_type=="Solid Tissue Normal"].drop_duplicates("case").set_index("case").aliquot
    cases=sorted(set(T.index)&set(N.index))
    r={"n_pairs_mapped":len(cases),"macrophage_proteins_present":[g for g in MAC if g in d.index]}
    sp_t=d.loc["SPP1",T[cases]].values.astype(float); sp_n=d.loc["SPP1",N[cases]].values.astype(float)
    ok=~np.isnan(sp_t)&~np.isnan(sp_n); diff=sp_t[ok]-sp_n[ok]
    r["C1"]={"n":int(ok.sum()),"median_diff_log2":float(np.median(diff)),"mean_diff_log2":float(diff.mean()),"n_tumor_higher":int((diff>0).sum()),
             "p_one_sided":float(stats.wilcoxon(diff,alternative="greater").pvalue)}
    # macrophage score: z across all samples (tumor+normal) per protein, mean
    allcols=list(T[cases])+list(N[cases])
    mac=[g for g in MAC if g in d.index]
    Z=pd.DataFrame({g:zrow(d.loc[g,allcols].astype(float)) for g in mac})
    macscore=Z.mean(axis=1,skipna=True)
    st=d.loc["SPP1",T[cases]].astype(float); ms=macscore[T[cases]]; ep=d.loc["EPCAM",T[cases]].astype(float) if "EPCAM" in d.index else None
    k=st.notna()&ms.notna(); rho,p=stats.spearmanr(st[k],ms[k]); p1=p/2 if rho>0 else 1-p/2
    r["C2"]={"n":int(k.sum()),"rho_SPP1_macscore_tumor":float(rho),"p_one_sided":float(p1)}
    if ep is not None:
        k2=st.notna()&ep.notna(); rho2,p2=stats.spearmanr(st[k2],ep[k2]); r["C2"]["rho_SPP1_EPCAM_tumor_control"]=float(rho2); r["C2"]["p_two_sided_EPCAM"]=float(p2)
    r["C2"]["per_protein_rho"]={g:float(stats.spearmanr(st,d.loc[g,T[cases]].astype(float),nan_policy="omit")[0]) for g in mac}
    # C3
    ds=pd.Series(sp_t-sp_n,index=cases); dm=pd.Series(macscore[T[cases]].values-macscore[N[cases]].values,index=cases)
    k3=ds.notna()&dm.notna(); x=dm[k3].values; y=ds[k3].values
    X=np.column_stack([np.ones(len(x)),x]); b,res_,rk,sv=np.linalg.lstsq(X,y,rcond=None)
    resid=y-X@b; s2=resid@resid/(len(y)-2); cov=s2*np.linalg.inv(X.T@X); t=b[0]/np.sqrt(cov[0,0])
    r["C3"]={"n":int(k3.sum()),"intercept":float(b[0]),"slope":float(b[1]),"t_intercept":float(t),"p_one_sided_intercept_gt0":float(stats.t.sf(t,len(y)-2)),
             "mean_delta_macscore":float(x.mean())}
    res[name]=r
res["C1_supported"]=all(res[n]["C1"]["p_one_sided"]<0.05 for n in ["Discovery","Confirmatory"])
res["C2_supported"]=all(res[n]["C2"]["p_one_sided"]<0.05 and res[n]["C2"]["rho_SPP1_macscore_tumor"]>0 for n in ["Discovery","Confirmatory"])
json.dump(res,open(OUT,"w"),indent=1,ensure_ascii=False); print(json.dumps(res,indent=1,ensure_ascii=False))
