"""按 docs/预设方案_KAC比例跨队列合并_2026-10-11.md 运行（SHA256 见 results/60_kac_pooled/prereg_sha256.txt）。"""
import json, numpy as np, pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
R="/Users/dudu/全部数据库/GGN_target/results/"
OUT=R+"60_kac_pooled/"
rng=np.random.default_rng(20261011)

p=pd.read_csv(R+"40_drivers/state_matched/step0_KAC_counts.tsv",sep="\t")
p=p[(p.normal_AT2L>=50)&(p.precursor_AT2L>=50)&(p.patient!="P7")]
peng=pd.concat([
 pd.DataFrame({"cohort":"Peng","unit":p.patient,"side":0,"kac":p.normal_KAC,"n":p.normal_AT2L}),
 pd.DataFrame({"cohort":"Peng","unit":p.patient,"side":1,"kac":p.precursor_KAC,"n":p.precursor_AT2L})])
g=pd.read_csv(R+"20_stage2/validation/GSE164789_cell_fractions.tsv",sep="\t")
g=g.dropna(subset=["KAC"]).copy()
g["pat"]=g.unit; g["kac"]=(g.KAC*g.n_at2).round().astype(int)
a=g[g.tissue=="A"].set_index("unit"); t=g[g.tissue=="T"].set_index("unit")
ids=[i for i in a.index if i in t.index and a.loc[i,"n_at2"]>=50 and t.loc[i,"n_at2"]>=50]
gse=pd.concat([
 pd.DataFrame({"cohort":"GSE164789","unit":ids,"side":0,"kac":a.loc[ids,"kac"].values,"n":a.loc[ids,"n_at2"].values}),
 pd.DataFrame({"cohort":"GSE164789","unit":ids,"side":1,"kac":t.loc[ids,"kac"].values,"n":t.loc[ids,"n_at2"].values})])
d=pd.concat([peng,gse],ignore_index=True)
d["unit_id"]=d.cohort+"_"+d.unit.astype(str)
d["other"]=d.n-d.kac
d["prop"]=d.kac/d.n
res={"n_units":{c:int(d[d.cohort==c].unit.nunique()) for c in d.cohort.unique()}}

def gee(df):
    m=smf.gee("kac + other ~ side"+(" + C(cohort)" if df.cohort.nunique()>1 else ""),groups="unit_id",data=df,family=sm.families.Binomial(),cov_struct=sm.cov_struct.Exchangeable())
    # 用 Binomial 的二列响应：kac, other
    f=m.fit()
    return f
def gee2(df):
    y=np.column_stack([df.kac,df.other])
    X=pd.DataFrame({"const":1.0,"side":df.side.values.astype(float)})
    if df.cohort.nunique()>1: X["GSE"]=(df.cohort=="GSE164789").astype(float).values
    f=sm.GEE(y,X,groups=df.unit_id.values,family=sm.families.Binomial(),cov_struct=sm.cov_struct.Exchangeable()).fit()
    b=f.params["side"];se=f.bse["side"]
    return {"OR":float(np.exp(b)),"CI95":[float(np.exp(b-1.96*se)),float(np.exp(b+1.96*se))],"p":float(f.pvalues["side"]),"b":float(b),"se":float(se)}
res["primary_pooled_GEE"]=gee2(d)
per={c:gee2(d[d.cohort==c]) for c in d.cohort.unique()}
res["a_per_cohort"]=per
# (a) DL random effects
bs=np.array([per[c]["b"] for c in per]);vs=np.array([per[c]["se"]**2 for c in per])
w=1/vs;fe=(w*bs).sum()/w.sum();Q=(w*(bs-fe)**2).sum();k=len(bs)
tau2=max(0,(Q-(k-1))/(w.sum()-(w**2).sum()/w.sum()))
wr=1/(vs+tau2);re=(wr*bs).sum()/wr.sum();sre=np.sqrt(1/wr.sum())
res["a_DL_random"]={"OR":float(np.exp(re)),"CI95":[float(np.exp(re-1.96*sre)),float(np.exp(re+1.96*sre))],"tau2":float(tau2),"I2":float(max(0,(Q-(k-1))/Q)) if Q>0 else 0.0}
# (b) 成对比例差 + bootstrap
w_=d.pivot_table(index=["cohort","unit_id"],columns="side",values="prop")
diff=(w_[1]-w_[0]).dropna()
def boot(x,B=10000):
    m=[rng.choice(x,len(x),replace=True).mean() for _ in range(B)]
    return [float(np.percentile(m,2.5)),float(np.percentile(m,97.5))]
res["b_pooled_paired_diff"]={"n":int(len(diff)),"mean_pp":float(diff.mean()*100),"CI95_pp":[x*100 for x in boot(diff.values)],
  "n_positive":int((diff>0).sum()),"n_negative":int((diff<0).sum())}
res["b_per_cohort"]={c:{"n":int(len(x)),"mean_pp":float(x.mean()*100),"CI95_pp":[v*100 for v in boot(x.values)],"n_positive":int((x>0).sum())} for c,x in diff.groupby(level=0)}
# (c) 符号翻转置换
obs=diff.mean();cnt=0;B=10000
for _ in range(B):
    s=rng.choice([-1,1],len(diff));cnt+= abs((diff.values*s).mean())>=abs(obs)
res["c_signflip_p_two_sided"]=float((cnt+1)/(B+1))
# (d) 分母>=100
d100=d.copy();ok=d100.groupby("unit_id").n.transform("min")>=100;d100=d100[ok]
res["d_min_denominator_100"]={"n_units":int(d100.unit_id.nunique()),**gee2(d100)}
# (e) 留一单位
loo=[]
for u in d.unit_id.unique():
    r=gee2(d[d.unit_id!=u]);loo.append((u,r["OR"],r["CI95"][0]))
loo=pd.DataFrame(loo,columns=["left_out","OR","CI_low"]).sort_values("CI_low")
loo.to_csv(OUT+"e_leave_one_out.tsv",sep="\t",index=False)
res["e_loo_min_CI_low"]=loo.head(3).values.tolist()
res["e_loo_all_CI_low_gt1"]=bool((loo.CI_low>1).all())
# (f) Peng 内分层
pf=pd.read_csv(R+"40_drivers/state_matched/step0_KAC_counts.tsv",sep="\t")
res["f_peng_AAH_pairs"]=int(((pf.AAH_KAC+pf.AAH_other)>=50).sum())
d.to_csv(OUT+"analysis_input.tsv",sep="\t",index=False);diff.to_csv(OUT+"paired_diffs.tsv",sep="\t")
# 判定（预设）
pooled=res["primary_pooled_GEE"]
c1=pooled["CI95"][0]>1
c2=all(per[c]["OR"]>1 for c in per)
c3=res["b_pooled_paired_diff"]["CI95_pp"][0]>0
res["verdict"]={"criterion1_pooled_CI_low_gt1":bool(c1),"criterion2_both_cohorts_OR_gt1":bool(c2),"criterion3_paired_diff_CI_low_gt0":bool(c3),
 "label":"复现成立" if (c1 and c2 and c3) else ("合并支持但队列间不一致" if c1 else "合并后不支持")}
json.dump(res,open(OUT+"kac_pooled_results.json","w"),ensure_ascii=False,indent=1)
print(json.dumps(res,ensure_ascii=False,indent=1))
