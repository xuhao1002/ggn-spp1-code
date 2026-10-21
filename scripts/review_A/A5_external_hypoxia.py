"""A5（方案 v1.3 §2-A5、§6.3）外部单细胞队列中缺氧评分与 HIF1A 程序的复算。环境 ggn_target。
实现细节见 docs/执行补充_A5_2026-10-21.md（执行前写定，启动时校验 SHA-256）。
A5a 缺氧评分（族 1，Holm m=4）；A5b 检出深度匹配后（族 2，m=4）；A5c HIF1A ULM（族 3，m=2）。四分类判定见执行补充 §6。"""
import os, sys, json, hashlib, inspect
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import P, SEED, check_plan, log, outdir, dump, holm_fixed
import numpy as np, pandas as pd, scanpy as sc, scipy.sparse as sp
from scipy import stats

sha = check_plan(); log("方案哈希校验通过", sha[:16])
NOTE = f"{P}/docs/执行补充_A5_2026-10-21.md"; NOTE_SHA = "74ef4e5cc0e4f9bc316573a75d7125ad4c525051f8def87b5182d76c7f0a9038"
ns = hashlib.sha256(open(NOTE, "rb").read()).hexdigest()
if ns != NOTE_SHA: sys.exit(f"执行补充记录哈希不符，拒绝运行：{ns}")
log("执行补充记录哈希校验通过", ns[:16])
# 评分调用与默认参数须与发现队列一致（§6.3）
src106 = open(f"{P}/scripts/106_hypoxia_H1.py", encoding="utf-8").read()
assert "sc.tl.score_genes(b,[g for g in gs if g in b.var_names],score_name=nm,random_state=0)" in src106, "106 的评分调用形式与预期不符"
sig = inspect.signature(sc.tl.score_genes).parameters
assert sig["ctrl_size"].default == 50 and sig["n_bins"].default == 25 and sig["gene_pool"].default is None, "scanpy score_genes 默认参数不符"
OUT = outdir("review_A5")

G = json.load(open(f"{P}/results/hypoxia/gene_sets_used.json"))
GS = {"HIF13": G["HIF14_used"], "HALLMARK": G["HALLMARK_used"]}
assert len(GS["HIF13"]) == 13 and len(GS["HALLMARK"]) == 197
net = pd.read_csv(f"{P}/results/mech/collectri_human.tsv", sep="\t")
net_sha = hashlib.sha256(open(f"{P}/results/mech/collectri_human.tsv", "rb").read()).hexdigest()
assert net_sha == json.load(open(f"{P}/results/mech/R_results.json"))["collectri_sha256"], "CollecTRI 快照与发现队列不一致"
hif = net[net.source == "HIF1A"]
CFG = {"GSE164789": "T", "GSE189357": None}

def ulm(Y, W):   # 与 110_tf_activity.py 相同
    n = W.shape[0]; Wc = W - W.mean(0); ssw = (Wc ** 2).sum(0); Yc = Y - Y.mean(1, keepdims=True)
    beta = (Yc @ Wc) / ssw; resid_ss = (Yc ** 2).sum(1, keepdims=True) - beta ** 2 * ssw
    return beta / np.sqrt(np.maximum(resid_ss, 1e-12) / (n - 2) / ssw)
def onesided(x):
    x = np.asarray(x, float); x = x[~np.isnan(x)]
    if len(x) < 5: return {"n_patients": int(len(x)), "evaluable": False, "p_one_sided": None}
    return {"n_patients": int(len(x)), "evaluable": True, "mean": float(x.mean()), "median": float(np.median(x)), "n_pos": int((x > 0).sum()),
            "p_one_sided": float(stats.wilcoxon(x, alternative="greater").pvalue), "p_two_sided": float(stats.wilcoxon(x).pvalue)}

res = {"plan_sha256": sha, "execution_note_sha256": ns, "collectri_sha256": net_sha, "scanpy": sc.__version__,
       "hallmark_source": "Enrichr MSigDB_Hallmark_2020 'Hypoxia' via results/hypoxia/gene_sets_used.json (D1 closed)", "cohorts": {}}
TESTS = {}
rng = np.random.default_rng(SEED)
for ds, tf in CFG.items():
    a = sc.read_h5ad(f"{P}/results/05_validation/{ds}/{ds}_annotated.h5ad")
    ann = pd.read_csv(f"{P}/results/05_validation/{ds}/cluster_annotation_auto.tsv", sep="\t", index_col=0)
    lab = a.obs["leiden_1.0"].astype(int).map(ann["final_label"]).astype(str); assert lab.ne("nan").all()
    C = a.layers["counts"].tocsr(); vn = list(a.var_names)
    ep = np.asarray(C[:, vn.index("EPCAM")].todense()).ravel()
    keep = (lab == "Macrophage").values & (ep < 2)
    if tf: keep &= (a.obs["tissue"].astype(str) == tf).values
    b = sc.AnnData(X=C[np.where(keep)[0]].astype(np.float32), obs=a.obs.loc[keep].copy(), var=pd.DataFrame(index=a.var_names))
    del a, C
    b.layers["counts"] = b.X.copy()
    sc.pp.normalize_total(b, target_sum=1e4); sc.pp.log1p(b)
    cr = {"n_macrophages": int(b.n_obs), "missing_genes": {}}
    for nm, gs in GS.items():
        pres = [g for g in gs if g in b.var_names]; cr["missing_genes"][nm] = sorted(set(gs) - set(pres))
        sc.tl.score_genes(b, pres, score_name=nm, random_state=0)
    Cm = b.layers["counts"].tocsr()
    o = b.obs.copy(); o["patient"] = o["patient"].astype(str); h = o["histology"].astype(str)
    o["grp"] = np.where(h == "IAC", "Invasive", np.where(h.isin(["AIS", "MIA", "AAH"]), "Precursor", "Other"))
    cr["histology_counts"] = h.value_counts().to_dict(); o = o[o.grp != "Other"]
    o["spp1"] = np.asarray(Cm[:, vn.index("SPP1")].todense()).ravel()[b.obs.index.get_indexer(o.index)]
    o["ngenes"] = np.asarray((Cm > 0).sum(1)).ravel()[b.obs.index.get_indexer(o.index)]
    # A5c 检出基因与 HIF1A 覆盖
    rowi = b.obs.index.get_indexer(o.index); tot = np.asarray(Cm[rowi].sum(0)).ravel(); det = tot >= 10
    genes = np.array(vn)[det]; gset = set(genes); Cd = Cm[rowi][:, det]
    nT = int(hif.target.nunique()); nD = int(hif.target[hif.target.isin(gset)].nunique())
    cov_ok = nD >= 10 and nD / nT >= 0.5
    cr["HIF1A_coverage"] = {"n_targets": nT, "n_detected": nD, "frac": nD / nT, "eligible": bool(cov_ok)}
    gi = {g: i for i, g in enumerate(genes)}; W = np.zeros((len(genes), 1))
    for s_, t_, w_ in hif[hif.target.isin(gset)][["source", "target", "weight"]].itertuples(index=False): W[gi[t_], 0] = w_
    pos_in_o = {ix: k for k, ix in enumerate(o.index)}
    def pb(ix): v = np.asarray(Cd[[pos_in_o[i] for i in ix]].sum(0)).ravel(); return np.log1p(v / v.sum() * 1e6)
    units, match_log = [], []
    for (p, g), u in sorted(o.groupby(["patient", "grp"]), key=lambda kv: kv[0]):
        pos, neg = u[u.spp1 > 0], u[u.spp1 == 0]
        r = {"cohort": ds, "patient": p, "grp": g, "n": len(u), "n_pos": len(pos), "n_neg": len(neg),
             "ngenes_med_pos": float(pos.ngenes.median()) if len(pos) else np.nan, "ngenes_med_neg": float(neg.ngenes.median()) if len(neg) else np.nan}
        okA = len(pos) >= 20 and len(neg) >= 20 and len(u) >= 50
        r["A5a_reason"] = "" if okA else ("unit <50" if len(u) < 50 else "a side <20")
        for nm in GS: r[f"A5a_{nm}"] = float(pos[nm].mean() - neg[nm].mean()) if okA else np.nan
        # A5b 匹配
        npairs, nstr = 0, 0; mp, mn = [], []
        if len(u) >= 50:
            st = pd.qcut(u.ngenes, q=5, duplicates="drop", labels=False); nstr = int(st.nunique())
            for s_ in sorted(st.unique()):
                uu = u[st == s_]; P_, N_ = uu[uu.spp1 > 0], uu[uu.spp1 == 0]; k = min(len(P_), len(N_))
                if k == 0: continue
                P_ = P_.iloc[np.sort(rng.choice(len(P_), k, replace=False))] if len(P_) > k else P_
                N_ = N_.iloc[np.sort(rng.choice(len(N_), k, replace=False))] if len(N_) > k else N_
                mp.append(P_); mn.append(N_); npairs += k
        okB = len(u) >= 50 and npairs >= 20
        r["A5b_pairs"], r["A5b_strata"] = npairs, nstr
        r["A5b_reason"] = "" if okB else ("unit <50" if len(u) < 50 else "matched pairs <20")
        for nm in GS: r[f"A5b_{nm}"] = float(pd.concat(mp)[nm].mean() - pd.concat(mn)[nm].mean()) if okB else np.nan
        if okB:
            MP, MN = pd.concat(mp), pd.concat(mn); r["ngenes_med_pos_matched"] = float(MP.ngenes.median()); r["ngenes_med_neg_matched"] = float(MN.ngenes.median())
        okC = cov_ok and len(pos) >= 50 and len(neg) >= 50
        r["A5c_reason"] = "" if okC else ("HIF1A coverage" if not cov_ok else "a side <50")
        r["A5c_HIF1A"] = float(ulm(pb(pos.index)[None, :], W)[0, 0] - ulm(pb(neg.index)[None, :], W)[0, 0]) if okC else np.nan
        units.append(r)
    U = pd.DataFrame(units); U.to_csv(OUT + f"A5_units_{ds}.tsv", sep="\t", index=False)
    PT = U.groupby("patient")[[c for c in U.columns if c.startswith(("A5a_", "A5b_", "A5c_")) and not c.endswith(("reason", "pairs", "strata"))]].mean()
    PT.to_csv(OUT + f"A5_patients_{ds}.tsv", sep="\t")
    d = (U.ngenes_med_pos - U.ngenes_med_neg)
    cr["detected_genes_unit_median_diff_pos_minus_neg"] = {"median": float(d.median()), "q1": float(d.quantile(.25)), "q3": float(d.quantile(.75)), "n_units": int(d.notna().sum())}
    for col in PT.columns:
        TESTS[(ds, col)] = onesided(PT[col].values)
        TESTS[(ds, col)]["n_units_evaluable"] = int(U[col].notna().sum())
        rc = col.split("_")[0] + "_reason"
        TESTS[(ds, col)]["units_not_evaluable"] = U.loc[U[col].isna(), rc].value_counts().to_dict()
    res["cohorts"][ds] = cr
    log(ds, "巨噬细胞", cr["n_macrophages"], "HIF1A 覆盖", cr["HIF1A_coverage"])

FAM = {"family1_A5a": (["A5a_HIF13", "A5a_HALLMARK"], 4), "family2_A5b": (["A5b_HIF13", "A5b_HALLMARK"], 4), "family3_A5c": (["A5c_HIF1A"], 2)}
out = []
for fam, (cols, m) in FAM.items():
    keys = [(ds, c) for ds in CFG for c in cols]
    adj_s = holm_fixed({f"{k[0]}|{k[1]}": TESTS[k]["p_one_sided"] for k in keys}, m)   # holm_fixed 需字符串键
    adj = {k: adj_s[f"{k[0]}|{k[1]}"] for k in keys}
    for k in keys:
        t = TESTS[k]; q = adj[k]
        if not t["evaluable"]: cat = "not evaluable"
        elif t["median"] < 0: cat = "opposite direction"
        elif t["median"] > 0 and q is not None and q < 0.05: cat = "meets replication criterion"
        else: cat = "does not meet replication criterion"
        out.append({"family": fam, "m": m, "cohort": k[0], "metric": k[1], **{kk: vv for kk, vv in t.items()}, "holm_p": q, "category": cat})
R = pd.DataFrame(out); R.to_csv(OUT + "A5_tests.tsv", sep="\t", index=False)
res["tests"] = out
dump(res, OUT + "A5_results.json")
log("A5 完成"); print(R[["family", "cohort", "metric", "n_patients", "median", "p_one_sided", "holm_p", "category"]].to_string())
