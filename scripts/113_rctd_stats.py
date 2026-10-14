"""D1/D2 统计。环境 ggn_target 或 ggn_spatial。"""
import os, re, glob, json, hashlib, numpy as np, pandas as pd
from scipy import stats
P="/Users/dudu/全部数据库/GGN_target/"; O=P+"results/mech/rctd/"
sha=hashlib.sha256(open(P+"docs/预设方案_调控互作RCTD_2026-10-13.md","rb").read()).hexdigest()
assert sha==open(P+"results/mech/prereg_sha256_D.txt").read().split()[0],sha
mp=pd.read_csv(O+"ref/label_map.tsv",sep="\t"); back=dict(zip(mp.safe,mp.label))
MAC=["Mac-SPP1","Mac-IL1B","Mac-other"]; rows=[]; d2=[]
for f in sorted(glob.glob(O+"GSM*_rctd.tsv.gz")):
    b=os.path.basename(f).replace("_rctd.tsv.gz",""); m=re.match(r"GSM\d+_(P\d+)_([A-Za-z]+)(-\d+)?$",b)
    w=pd.read_csv(f,sep="\t",index_col=0).rename(columns=back)
    mac=w[MAC].sum(1)
    rows.append({"sample":b,"patient":m.group(1),"stage":m.group(2),"n_spots":len(w),"MacSPP1":w["Mac-SPP1"].mean(),
                 "MacSPP1_of_mac":float(w["Mac-SPP1"].sum()/mac.sum()) if mac.sum()>0 else np.nan})
    c=P+f"results/30_spatial/deconv/map_unmerged_n6/{b}_q05.tsv.gz"
    if os.path.exists(c):
        q=pd.read_csv(c,sep="\t",index_col=0); q=q.div(q.sum(1),axis=0); sh=w.index.intersection(q.index)
        d2.append({"sample":b,"spot_rho":float(stats.spearmanr(w.loc[sh,"Mac-SPP1"],q.loc[sh,"Mac-SPP1"])[0]),
                   "rctd_mean":float(w.loc[sh,"Mac-SPP1"].mean()),"c2l_mean":float(q.loc[sh,"Mac-SPP1"].mean())})
S=pd.DataFrame(rows); S.to_csv(O+"D_sections.tsv",sep="\t",index=False)
S["grp"]=np.where(S.stage=="LUAD","LUAD",np.where(S.stage.isin(["AAH","AIS","MIA"]),"Precursor","Normal"))
def paired(col):
    w=S[S.grp!="Normal"].groupby(["patient","grp"])[col].mean().unstack().dropna(); d=(w.LUAD-w.Precursor).values
    return {"n_pairs":int(len(d)),"median_diff":float(np.median(d)),"n_pos":int((d>0).sum()),
            "p_one_sided":float(stats.wilcoxon(d,alternative="greater").pvalue),"p_two_sided":float(stats.wilcoxon(d).pvalue)} if len(d)>=5 else {"n_pairs":int(len(d))}
res={"plan_sha256":sha,"n_sections":len(S),"D1_MacSPP1":paired("MacSPP1"),"desc_MacSPP1_share_of_macrophages":paired("MacSPP1_of_mac")}
if d2:
    D=pd.DataFrame(d2); D.to_csv(O+"D2_agreement.tsv",sep="\t",index=False)
    res["D2"]={"n_sections":len(D),"spot_rho_median":float(D.spot_rho.median()),"spot_rho_range":[float(D.spot_rho.min()),float(D.spot_rho.max())],
               "section_mean_rho":float(stats.spearmanr(D.rctd_mean,D.c2l_mean)[0]) if len(D)>=4 else None}
json.dump(res,open(O+"D_results.json","w"),indent=1); [print(k,v) for k,v in res.items()]; print("DONE")
