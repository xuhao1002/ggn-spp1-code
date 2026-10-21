"""偏离记录 D7（post hoc，2026-10-21 14:24 写于运行前）：外部队列 HIF1A 程序（A5c）的检出深度匹配。环境 ggn_target。
数据、细胞、单位、检出基因、覆盖、伪批量与 ULM 与 A5c 完全相同（复用 A5 脚本逻辑并先复现 A5c 的患者层面值）；每侧伪批量只由匹配后细胞构成。"""
import os, sys, json, hashlib
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import P, SEED, check_plan, log, outdir, dump, holm_fixed
import numpy as np, pandas as pd, scanpy as sc
from scipy import stats

sha = check_plan(); log("方案哈希校验通过", sha[:16])
dev = open(f"{P}/docs/偏离与修订记录_外审回应A1-A6.md", encoding="utf-8").read()
assert "## D7（2026-10-21" in dev, "偏离记录中未找到 D7，拒绝运行"
dev_sha = hashlib.sha256(dev.encode("utf-8")).hexdigest(); log("偏离记录 SHA", dev_sha[:16])
A5RUN = f"{P}/results/review_A5/run_20261021_140428/"
OUT = outdir("review_D7")
net = pd.read_csv(f"{P}/results/mech/collectri_human.tsv", sep="\t")
assert hashlib.sha256(open(f"{P}/results/mech/collectri_human.tsv", "rb").read()).hexdigest() == json.load(open(f"{P}/results/mech/R_results.json"))["collectri_sha256"]
hif = net[net.source == "HIF1A"]
CFG = {"GSE164789": "T", "GSE189357": None}
def ulm(Y, W):
    n = W.shape[0]; Wc = W - W.mean(0); ssw = (Wc ** 2).sum(0); Yc = Y - Y.mean(1, keepdims=True)
    beta = (Yc @ Wc) / ssw; resid_ss = (Yc ** 2).sum(1, keepdims=True) - beta ** 2 * ssw
    return beta / np.sqrt(np.maximum(resid_ss, 1e-12) / (n - 2) / ssw)
rng = np.random.default_rng(SEED)
units = []
for ds, tf in CFG.items():
    a = sc.read_h5ad(f"{P}/results/05_validation/{ds}/{ds}_annotated.h5ad")
    ann = pd.read_csv(f"{P}/results/05_validation/{ds}/cluster_annotation_auto.tsv", sep="\t", index_col=0)
    lab = a.obs["leiden_1.0"].astype(int).map(ann["final_label"]).astype(str); assert lab.ne("nan").all()
    C = a.layers["counts"].tocsr(); vn = list(a.var_names)
    ep = np.asarray(C[:, vn.index("EPCAM")].todense()).ravel()
    keep = (lab == "Macrophage").values & (ep < 2)
    if tf: keep &= (a.obs["tissue"].astype(str) == tf).values
    o = a.obs.loc[keep].copy(); Cm = C[np.where(keep)[0]]; del a, C
    o["patient"] = o["patient"].astype(str); h = o["histology"].astype(str)
    o["grp"] = np.where(h == "IAC", "Invasive", np.where(h.isin(["AIS", "MIA", "AAH"]), "Precursor", "Other"))
    o["r"] = np.arange(len(o)); o = o[o.grp != "Other"]
    o["spp1"] = np.asarray(Cm[o.r.values][:, vn.index("SPP1")].todense()).ravel()
    o["ngenes"] = np.asarray((Cm[o.r.values] > 0).sum(1)).ravel()
    tot = np.asarray(Cm[o.r.values].sum(0)).ravel(); det = tot >= 10       # 与 A5c 相同
    genes = np.array(vn)[det]; gset = set(genes); Cd = Cm[:, det]
    nT = int(hif.target.nunique()); nD = int(hif.target[hif.target.isin(gset)].nunique())
    assert nD >= 10 and nD / nT >= 0.5
    gi = {g: i for i, g in enumerate(genes)}; W = np.zeros((len(genes), 1))
    for t_, w_ in hif[hif.target.isin(gset)][["target", "weight"]].itertuples(index=False): W[gi[t_], 0] = w_
    def pb(rows): v = np.asarray(Cd[rows].sum(0)).ravel(); return np.log1p(v / v.sum() * 1e6)[None, :]
    for (p, g), u in sorted(o.groupby(["patient", "grp"]), key=lambda kv: kv[0]):
        pos, neg = u[u.spp1 > 0], u[u.spp1 == 0]
        r = {"cohort": ds, "patient": p, "grp": g, "n_pos": len(pos), "n_neg": len(neg)}
        if len(pos) >= 50 and len(neg) >= 50:
            r["unmatched"] = float(ulm(pb(pos.r.values), W)[0, 0] - ulm(pb(neg.r.values), W)[0, 0])
            r["ng_pos"], r["ng_neg"] = float(pos.ngenes.median()), float(neg.ngenes.median())
            st = pd.qcut(u.ngenes, q=5, duplicates="drop", labels=False); mp, mn = [], []
            for s_ in sorted(st.unique()):
                uu = u[st == s_]; P_, N_ = uu[uu.spp1 > 0], uu[uu.spp1 == 0]; k = min(len(P_), len(N_))
                if k == 0: continue
                if len(P_) > k: P_ = P_.iloc[np.sort(rng.choice(len(P_), k, replace=False))]
                if len(N_) > k: N_ = N_.iloc[np.sort(rng.choice(len(N_), k, replace=False))]
                mp.append(P_); mn.append(N_)
            MP, MN = pd.concat(mp), pd.concat(mn); r["pairs"] = len(MP); r["strata"] = int(st.nunique())
            if len(MP) >= 50:
                r["matched"] = float(ulm(pb(MP.r.values), W)[0, 0] - ulm(pb(MN.r.values), W)[0, 0])
                r["ng_pos_m"], r["ng_neg_m"] = float(MP.ngenes.median()), float(MN.ngenes.median())
        units.append(r)
    log(ds, "完成")
U = pd.DataFrame(units); U.to_csv(OUT + "D7_units.tsv", sep="\t", index=False)
out, pv = [], {}
for ds in CFG:
    d = U[U.cohort == ds]; pu = d.groupby("patient").unmatched.mean().dropna(); pm = d.groupby("patient").matched.mean().dropna()
    ref = pd.read_csv(A5RUN + f"A5_patients_{ds}.tsv", sep="\t", index_col=0)["A5c_HIF1A"].dropna()
    assert set(ref.index) == set(pu.index) and np.allclose(ref.loc[pu.index].values, pu.values, atol=1e-8), f"{ds}: 未复现 A5c"
    same = pu.loc[pm.index]
    r = {"cohort": ds, "unmatched_n": int(len(pu)), "unmatched_npos": int((pu > 0).sum()), "unmatched_median": float(pu.median()),
         "matched_n": int(len(pm)), "matched_npos": int((pm > 0).sum()), "matched_median": float(pm.median()) if len(pm) else None,
         "unmatched_median_same_patients": float(same.median()) if len(pm) else None,
         "pairs_median": float(d.pairs.median()), "units_lost_after_matching": int(d.unmatched.notna().sum() - d.matched.notna().sum()),
         "ng_pos_before": float(d.ng_pos.median()), "ng_neg_before": float(d.ng_neg.median()),
         "ng_pos_after": float(d.ng_pos_m.median()), "ng_neg_after": float(d.ng_neg_m.median())}
    r["retention"] = r["matched_median"] / r["unmatched_median_same_patients"] if len(pm) and r["unmatched_median_same_patients"] else None
    r["p_one_sided"] = float(stats.wilcoxon(pm.values, alternative="greater").pvalue) if len(pm) >= 5 else None
    pv[ds] = r["p_one_sided"]; out.append(r)
adj = holm_fixed(pv, 2)
for r in out:
    q = adj[r["cohort"]]; r["holm_p"] = q; med, ret = r["matched_median"], r["retention"]
    if med is None or q is None: r["category_A5"] = r["category_D6"] = "not evaluable"; continue
    r["category_A5"] = "opposite direction" if med < 0 else ("meets replication criterion" if (med > 0 and q < 0.05) else "does not meet replication criterion")
    r["category_D6"] = "robust" if (med > 0 and q < 0.05) else ("attenuated but retained" if (med > 0 and ret is not None and ret >= 0.5) else "largely explained by detection depth")
R = pd.DataFrame(out); R.to_csv(OUT + "D7_summary.tsv", sep="\t", index=False)
dump({"plan_sha256": sha, "deviation_record_sha256_at_run": dev_sha, "deviation": "D7 (post hoc)", "results": out}, OUT + "D7_results.json")
log("D7 完成"); print(R.to_string())
