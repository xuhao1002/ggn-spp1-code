"""第五部分 §8/§8.1：SPP1 用巨噬细胞 eQTL（Nedelec_2016 macrophage_naive, QTD000384）做工具变量。由 95 派生，GWAS 处理、对齐、Wald 比、共定位全部相同。
原 95 说明：Wald 比 MR + 近似 Bayes 因子共定位。严格按 docs/第五部分_遗传学证据预设方案.md 第 3、7 节（SHA 见 results/50_genetics_prereg_sha256.txt）。
环境 ggn_target（需 pyliftover）。输出 results/50_genetics/。"""
import matplotlib; matplotlib.use("Agg")
import os, sys, json, time, gzip, hashlib
import numpy as np, pandas as pd
from scipy import stats
t0=time.time()
def log(*a): print(f"[{time.time()-t0:7.0f}s]",*a,flush=True)
R="/Users/dudu/全部数据库/GGN_target/"; RAW=R+"raw/"; OUT=R+"results/50_genetics/macrophage_SPP1/"; os.makedirs(OUT,exist_ok=True)

# ---------- 0. 预设核对 ----------
sha=hashlib.sha256(open(R+"docs/第五部分_遗传学证据预设方案.md","rb").read()).hexdigest()
log("plan sha256",sha); assert sha.startswith("31da180a"), "方案文件已改动，停止"

FAM={"primary":["SPP1"]}
genes=["SPP1"]
fam_of={g:f for f,v in FAM.items() for g in v}

# ---------- 1. coloc 实现与自检 ----------
def labf(beta,se,W):
    z=beta/se; V=se**2; r=W/(W+V); return 0.5*(np.log(1-r)+r*z**2)
def lsum(x): m=np.max(x); return m+np.log(np.sum(np.exp(x-m)))
def ldiff(a,b): return a+np.log1p(-np.exp(b-a)) if a>b else -np.inf
def coloc(b1,s1,b2,s2,W1=0.15**2,W2=0.2**2,p1=1e-4,p2=1e-4,p12=1e-5):
    l1=labf(b1,s1,W1); l2=labf(b2,s2,W2); s1_=lsum(l1); s2_=lsum(l2); s12=lsum(l1+l2)
    lh=np.array([0.0,np.log(p1)+s1_,np.log(p2)+s2_,np.log(p1)+np.log(p2)+ldiff(s1_+s2_,s12),np.log(p12)+s12])
    pp=np.exp(lh-lsum(lh)); return dict(zip(["PP.H0","PP.H1","PP.H2","PP.H3","PP.H4"],pp.round(4).tolist()))
rng=np.random.default_rng(0); n=300; se=np.full(n,0.05)
b1=rng.normal(0,0.05,n); b2=rng.normal(0,0.05,n); b1[100]=0.5; b2[100]=0.5
same=coloc(b1,se,b2,se)
b2b=rng.normal(0,0.05,n); b2b[200]=0.5
diff=coloc(b1,se,b2b,se)
log("coloc self-check same",same,"diff",diff)
assert max(same,key=same.get)=="PP.H4" and max(diff,key=diff.get)=="PP.H3", "coloc 自检未通过，停止"
json.dump({"same_causal":same,"distinct_causal":diff},open(OUT+"coloc_selfcheck.json","w"),indent=1)

# ---------- 2. 基因 ID（GTEx v8 表，§7.9） ----------
g=pd.read_csv(RAW+"safety_druggability/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz",sep="\t",skiprows=2,usecols=[0,1])
g.columns=["ens","sym"]; g["ens"]=g.ens.str.split(".").str[0]
gm=g[g.sym.isin(genes)].groupby("sym").ens.unique()
ens={s:v[0] for s,v in gm.items() if len(v)==1}; log("mapped",len(ens),"of",len(genes))

# ---------- 3. eQTL ----------
E=RAW+"eQTL_macrophage/"
perm=pd.read_csv(E+"QTD000384.permuted.tsv.gz",sep="\t")
log("permuted cols",list(perm.columns),perm.shape)
perm=perm.dropna(subset=["p_beta"]).copy()
from statsmodels.stats.multitest import multipletests
perm["q_beta"]=multipletests(perm.p_beta,method="fdr_bh")[1]
cache=OUT+"eqtl_selected_genes.tsv.gz"
want=set(ens.values())
if not os.path.exists(cache):
    keep=[]
    for ch in pd.read_csv(E+"QTD000384.all.tsv.gz",sep="\t",chunksize=2_000_000,dtype={"chromosome":str}):
        k=ch[ch.gene_id.isin(want)]
        if len(k): keep.append(k)
    eq=pd.concat(keep); eq.to_csv(cache,sep="\t",index=False)
eq=pd.read_csv(cache,sep="\t",dtype={"chromosome":str})
log("eQTL rows for genes",len(eq),"cols",list(eq.columns))

# ---------- 4. GWAS 读入、质控、λ、阳性对照、荟萃 ----------
def read_gwas(fn):
    d=pd.read_csv(RAW+"GWAS_NCC/"+fn,sep="\t",usecols=["MARKER","ALLELES","FREQ1","RSQR","EFFECT1","STDERR","PVALUE"],dtype=str)
    n0=len(d)
    a=d.ALLELES.str.split(",",expand=True); d["a1"]=a[0]; d["a2"]=a[1]
    d["f1"]=pd.to_numeric(d.FREQ1.str.split("(").str[0],errors="coerce")
    for c in ["RSQR","EFFECT1","STDERR","PVALUE"]: d[c]=pd.to_numeric(d[c],errors="coerce")
    m=d.MARKER.str.split("_",expand=True); d["chr"]=m[0]; d["pos"]=pd.to_numeric(m[1],errors="coerce")
    ok=(d.RSQR>=0.8)&d.a1.str.fullmatch("[ACGT]")&d.a2.str.fullmatch("[ACGT]")&(d.f1>0.01)&(d.f1<0.99)&d.STDERR.gt(0)&d.EFFECT1.notna()&d.pos.notna()
    d=d[ok][["chr","pos","a1","a2","f1","EFFECT1","STDERR","PVALUE"]].copy()
    lam=float(np.median(stats.chi2.isf(d.PVALUE.clip(1e-300,1),1))/stats.chi2.ppf(0.5,1))
    log(fn,"rows",n0,"after QC",len(d),"lambdaGC",round(lam,3)); return d,lam,n0
g1,lam1,n1=read_gwas("NCCandBBJ_mach2dat_5PCs_20150824.txt")
g2,lam2,n2=read_gwas("Kyodai_mach2dat_5PCs_20150508.txt")
m=g1.merge(g2,on=["chr","pos"],how="outer",suffixes=("_1","_2"))
same=(m.a1_1==m.a1_2)&(m.a2_1==m.a2_2); flip=(m.a1_1==m.a2_2)&(m.a2_1==m.a1_2)
b2=np.where(same,m.EFFECT1_2,np.where(flip,-m.EFFECT1_2,np.nan))
both=m.a1_1.notna()&~np.isnan(b2)
only1=m.a1_1.notna()&m.a1_2.isna(); only2=m.a1_1.isna()&m.a1_2.notna()
w1=1/m.STDERR_1**2; w2=1/m.STDERR_2**2
meta=pd.DataFrame({"chr":m.chr,"pos":m.pos})
meta["a1"]=np.where(m.a1_1.notna(),m.a1_1,m.a1_2); meta["a2"]=np.where(m.a1_1.notna(),m.a2_1,m.a2_2)
meta["f1"]=np.where(m.a1_1.notna(),m.f1_1,m.f1_2)
bm=np.where(both,(m.EFFECT1_1*w1+b2*w2)/(w1+w2),np.where(only1,m.EFFECT1_1,np.where(only2,m.EFFECT1_2,np.nan)))
sm=np.where(both,1/np.sqrt(w1+w2),np.where(only1,m.STDERR_1,np.where(only2,m.STDERR_2,np.nan)))
meta["beta"]=bm; meta["se"]=sm
meta["b1"]=m.EFFECT1_1.values; meta["s1"]=m.STDERR_1.values; meta["b2"]=b2; meta["s2"]=m.STDERR_2.values
meta=meta.dropna(subset=["beta","se"]); meta["p"]=2*stats.norm.sf(np.abs(meta.beta/meta.se))
log("meta variants",len(meta),"both",int(both.sum()),"allele-mismatch dropped",int((m.a1_1.notna()&m.a1_2.notna()&~same&~flip).sum()))
PC={"TERT":("5",1253282-500000,1295183+500000),"TP63":("3",189349205-500000,189615068+500000)}
pc={}
for k,(c,a,b) in PC.items():
    s=meta[(meta.chr==c)&(meta.pos>=a)&(meta.pos<=b)]
    pc[k]={"n_variants":int(len(s)),"min_p_meta":float(s.p.min()) if len(s) else None,"pass":bool(len(s) and s.p.min()<5e-8)}
log("positive controls",pc,"lambda",lam1,lam2)
qc={"lambdaGC_GWAS1":lam1,"lambdaGC_GWAS2":lam2,"rows_raw":[n1,n2],"rows_qc":[int(len(g1)),int(len(g2))],"meta_variants":int(len(meta)),"positive_controls":pc}
json.dump(qc,open(OUT+"gwas_qc.json","w"),indent=1)
if not any(v["pass"] for v in pc.values()):
    log("两个阳性对照都未达标：按方案不解读 MR 结果，停止"); sys.exit(0)

# ---------- 5. 坐标转换 hg19 -> hg38 ----------
from pyliftover import LiftOver
lo=LiftOver(RAW+"GWAS_NCC/hg19ToHg38.over.chain.gz")
chrs=set(eq.chromosome.astype(str))
meta=meta[meta.chr.isin(chrs)].copy()
p38=[]; 
for c,p in zip(meta.chr.values,meta.pos.values.astype(int)):
    r=lo.convert_coordinate("chr"+c,p-1)
    p38.append(r[0][1]+1 if (r and len(r)==1 and r[0][0]=="chr"+c) else -1)
meta["pos38"]=p38; nfail=int((meta.pos38<0).sum()); meta=meta[meta.pos38>0]
log("lifted",len(meta),"failed/ambiguous",nfail)

# ---------- 6. 对齐 ----------
comp={"A":"T","T":"A","C":"G","G":"C"}
eq=eq.rename(columns={"chromosome":"chr","position":"pos38"})
eq["chr"]=eq.chr.astype(str)
J=eq.merge(meta[["chr","pos38","a1","a2","f1","beta","se","b1","s1","b2","s2"]],on=["chr","pos38"],how="inner",suffixes=("_e","_g"))
def align(r):
    ref,alt,a1,a2=r["ref"],r["alt"],r["a1"],r["a2"]
    if len(ref)!=1 or len(alt)!=1: return np.nan
    pal={ref,alt} in ({"A","T"},{"C","G"})
    if pal:
        if 0.42<=r["maf"]<=0.5: return np.nan
        # 用频率定方向：eQTL 的 alt 频率需与 GWAS 同向（均 <0.5 或均 >0.5）
        alt_f=r["ac"]/r["an"] if ("ac" in r and r["an"]>0) else np.nan
        if np.isnan(alt_f): return np.nan
        if a1==alt and (alt_f<0.5)==(r["f1"]<0.5): return 1.0
        if a1==ref and (alt_f<0.5)==((1-r["f1"])<0.5): return -1.0
        return np.nan
    if a1==alt and a2==ref: return 1.0
    if a1==ref and a2==alt: return -1.0
    if comp.get(a1)==alt and comp.get(a2)==ref: return 1.0
    if comp.get(a1)==ref and comp.get(a2)==alt: return -1.0
    return np.nan
J["sign"]=J.apply(align,axis=1); J=J.dropna(subset=["sign"])
for b in ["beta_g","b1","b2"]: J[b]=J[b]*J["sign"]
J=J.rename(columns={"beta_e":"beta_eqtl","se_e":"se_eqtl","beta_g":"beta_gwas","se_g":"se_gwas"}) if "beta_e" in J else J.rename(columns={"beta_g":"beta_gwas","se_g":"se_gwas"})
log("aligned eQTL×GWAS rows",len(J))

# ---------- 7. 每个基因 ----------
def wald(bg,sg,be):
    b=bg/be; s=sg/abs(be); z=b/s
    return {"logOR":b,"se":s,"OR":float(np.exp(b)),"OR_lo":float(np.exp(b-1.96*s)),"OR_hi":float(np.exp(b+1.96*s)),"p_one_sided_risk":float(stats.norm.sf(z)),"p_two_sided":float(2*stats.norm.sf(abs(z)))}
rows=[]
for sym in genes:
    r={"gene":sym,"family":fam_of[sym],"ensembl":ens.get(sym)}
    if sym not in ens: r["status"]="无法检验：ID 映射失败"; rows.append(r); continue
    e=ens[sym]; pr=perm[perm.molecular_trait_id==e] if "molecular_trait_id" in perm else perm[perm.gene_id==e]
    r["eGene_q"]=float(pr.q_beta.iloc[0]) if len(pr) else None
    r["is_eGene"]=bool(len(pr) and pr.q_beta.iloc[0]<0.05)
    lead=pr.variant.iloc[0] if len(pr) else None; r["lead_variant"]=lead
    Jg=J[J.gene_id==e]
    if len(Jg)==0: r["status"]="无法检验：无可对齐变异"; rows.append(r); continue
    cand=Jg[Jg.variant==lead] if lead is not None else Jg.iloc[0:0]
    if len(cand) and r["is_eGene"]:
        iv=cand.iloc[0]; r["instrument_is_lead"]=True
    else:
        Jt=Jg[Jg.pvalue<1e-5].sort_values("pvalue")
        if not r["is_eGene"] and len(Jt)==0: r["status"]="无法检验：非 eGene 且无 P<1e-5 的 cis 变异"; rows.append(r); continue
        if len(Jt)==0: r["status"]="无法检验：主变异缺失且无 P<1e-5 的替代变异"; rows.append(r); continue
        iv=Jt.iloc[0]; r["instrument_is_lead"]=False
    r.update(instrument=iv.variant,rsid=iv.get("rsid"),eqtl_beta=float(iv.beta_eqtl),eqtl_p=float(iv.pvalue),gwas_beta=float(iv.beta_gwas),gwas_se=float(iv.se_gwas))
    r.update({"wald_"+k:v for k,v in wald(iv.beta_gwas,iv.se_gwas,iv.beta_eqtl).items()})
    for tag,bb,ss in [("GWAS1",iv.b1,iv.s1),("GWAS2",iv.b2,iv.s2)]:
        if pd.notna(bb) and pd.notna(ss): w=wald(bb,ss,iv.beta_eqtl); r[f"{tag}_OR"]=w["OR"]; r[f"{tag}_p_one_sided"]=w["p_one_sided_risk"]
    reg=Jg[(Jg.pos38>=iv.pos38-500000)&(Jg.pos38<=iv.pos38+500000)].drop_duplicates("variant")
    r["coloc_n_variants"]=int(len(reg))
    if len(reg)>=50: r.update(coloc(reg.beta_eqtl.values,reg.se_eqtl.values,reg.beta_gwas.values,reg.se_gwas.values))
    else: r["coloc_note"]="区域变异 <50，未做"
    r["status"]="已检验"; rows.append(r)
T=pd.DataFrame(rows)
# ---------- 8. 多重检验与判定 ----------
T["holm"]=np.nan
for f in ["A_stage2","B_literature","C_microenv"]:
    idx=T[(T.family==f)&T.wald_p_one_sided_risk.notna()].index
    if len(idx): T.loc[idx,"holm"]=multipletests(T.loc[idx,"wald_p_one_sided_risk"],method="holm")[1]
pri=T.family=="primary"; T.loc[pri,"holm"]=T.loc[pri,"wald_p_one_sided_risk"]
# 保护方向（OR<1）同样标准，单侧反向 P
T["p_one_sided_protective"]=1-T.wald_p_one_sided_risk
T["holm_protective"]=np.nan
for f in ["A_stage2","B_literature","C_microenv"]:
    idx=T[(T.family==f)&T.p_one_sided_protective.notna()].index
    if len(idx): T.loc[idx,"holm_protective"]=multipletests(T.loc[idx,"p_one_sided_protective"],method="holm")[1]
T.loc[pri,"holm_protective"]=T.loc[pri,"p_one_sided_protective"]
def call(r):
    if r.get("status")!="已检验": return r.get("status")
    if r.family=="descriptive": return "描述性（不判定）"
    h4=r.get("PP.H4"); h4=h4 if pd.notna(h4) else 0
    if r.wald_p_one_sided_risk<0.05 and r.holm<0.05 and h4>=0.8: return "遗传支持（风险方向）"
    if r.p_one_sided_protective<0.05 and r.holm_protective<0.05 and h4>=0.8: return "遗传支持（保护方向）"
    return "无遗传证据"
T["call"]=T.apply(call,axis=1)
T.to_csv(OUT+"mr_coloc_results.tsv",sep="\t",index=False)
cols=["gene","family","is_eGene","instrument_is_lead","wald_OR","wald_OR_lo","wald_OR_hi","wald_p_one_sided_risk","holm","PP.H3","PP.H4","call"]
log("\n"+T[[c for c in cols if c in T]].to_string())
log("DONE")
