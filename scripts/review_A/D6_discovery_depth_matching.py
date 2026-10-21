"""偏离记录 D6（post hoc，2026-10-21 14:10 写于运行前）：发现队列缺氧/HIF1A 比较的检出深度匹配。环境 ggn_target。
细胞、单位、评分、门槛与 106/108/110 相同；只增加按检出基因数五分位层内 1:1 匹配。先复现原报告（不符即中止），再匹配。"""
import os, sys, json, hashlib
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import P, SEED, check_plan, log, outdir, dump, holm_fixed
import numpy as np, pandas as pd, scanpy as sc, scipy.sparse as sp
from scipy import stats

sha = check_plan(); log("方案哈希校验通过", sha[:16])
dev = open(f"{P}/docs/偏离与修订记录_外审回应A1-A6.md", encoding="utf-8").read()
assert "## D6（2026-10-21" in dev, "偏离记录中未找到 D6，拒绝运行"
dev_sha = hashlib.sha256(dev.encode("utf-8")).hexdigest(); log("偏离记录 SHA", dev_sha[:16])
OUT = outdir("review_D6")
G = json.load(open(f"{P}/results/hypoxia/gene_sets_used.json"))
H1 = json.load(open(f"{P}/results/hypoxia/H1_results.json")); H1c = json.load(open(f"{P}/results/hypoxia/H1c_results.json"))
RR = json.load(open(f"{P}/results/mech/R_results.json"))
net = pd.read_csv(f"{P}/results/mech/collectri_human.tsv", sep="\t")
assert hashlib.sha256(open(f"{P}/results/mech/collectri_human.tsv", "rb").read()).hexdigest() == RR["collectri_sha256"]
cov = pd.read_csv(f"{P}/results/mech/R_tf_coverage.tsv", sep="\t"); elig = cov[cov.eligible].TF.tolist()
assert "HIF1A" in elig
HIF14 = ["CA9", "VEGFA", "SLC2A1", "PGK1", "LDHA", "BNIP3", "NDRG1", "ADM", "ANKRD37", "P4HA1", "ENO1", "ALDOA", "PDK1", "EGLN3"]
HALL200 = set(G["HALLMARK_used"]) | {"ANXA2", "ENO1", "MIF"}   # D1：与 Enrichr/MSigDB 200 基因相同

a = sc.read_h5ad(f"{P}/results/01_snRNA/snRNA_final.h5ad")
epcam = np.asarray(a[:, "EPCAM"].layers["counts"].todense()).ravel()
b = a[((a.obs["final_label"].astype(str) == "Macrophage").values) & (epcam < 2)].copy(); del a
C = sp.csr_matrix(b.layers["counts"]); vn = list(b.var_names)
for nm, gs in [("HIF", G["HIF14_used"]), ("HALL", G["HALLMARK_used"]), ("STATE", G["state_genes"])]:
    sc.tl.score_genes(b, [g for g in gs if g in b.var_names], score_name=nm, random_state=0)   # 与 106/108 相同
o = b.obs.copy(); o["patient"] = o.patient.astype(str); o["stage"] = o.stage.astype(str)
o["spp1"] = np.asarray(C[:, vn.index("SPP1")].todense()).ravel(); o["ngenes"] = np.asarray((C > 0).sum(1)).ravel()
o["grp"] = np.where(o.stage == "LUAD", "LUAD", np.where(o.stage.isin(["AAH", "AIS", "MIA"]), "Precursor", "Other"))
o["row"] = np.arange(len(o)); o = o[o.grp != "Other"]
# 110 的检出基因与 ULM
# 与 110 相同：检出基因按全部巨噬细胞（含非 LUAD/前体组）计
tot = np.asarray(C.sum(0)).ravel(); det = tot >= 10
genes = np.array(vn)[det]; gset = set(genes); Cd = C[:, det]
gi = {g: i for i, g in enumerate(genes)}
def W_of(tg):
    W = np.zeros((len(genes), 1))
    for t_, w_ in tg[["target", "weight"]].itertuples(index=False):
        if t_ in gi: W[gi[t_], 0] = w_
    return W
hifT = net[net.source == "HIF1A"]; W_h = W_of(hifT); W_dc = W_of(hifT[~hifT.target.isin(set(HIF14) | HALL200)])
assert int(hifT[~hifT.target.isin(set(HIF14) | HALL200)].target.isin(gset).sum()) == RR["HIF1A_decirc_targets_detected"], "去循环靶基因数与 110 不符"
def ulm(Y, W):
    n = W.shape[0]; Wc = W - W.mean(0); ssw = (Wc ** 2).sum(0); Yc = Y - Y.mean(1, keepdims=True)
    beta = (Yc @ Wc) / ssw; rs = (Yc ** 2).sum(1, keepdims=True) - beta ** 2 * ssw
    return float((beta / np.sqrt(np.maximum(rs, 1e-12) / (n - 2) / ssw))[0, 0])
def pb(rows): v = np.asarray(Cd[rows].sum(0)).ravel(); return np.log1p(v / v.sum() * 1e6)[None, :]
rng = np.random.default_rng(SEED)
def match(A, B):
    """A、B 为两侧细胞（DataFrame）；按合并后检出基因数五分位分层，层内 1:1 无放回随机匹配。"""
    allc = pd.concat([A.assign(_s=0), B.assign(_s=1)]); st = pd.qcut(allc.ngenes, q=5, duplicates="drop", labels=False)
    ma, mb = [], []
    for s_ in sorted(st.unique()):
        x = allc[st == s_]; xa, xb = x[x._s == 0], x[x._s == 1]; k = min(len(xa), len(xb))
        if k == 0: continue
        if len(xa) > k: xa = xa.iloc[np.sort(rng.choice(len(xa), k, replace=False))]
        if len(xb) > k: xb = xb.iloc[np.sort(rng.choice(len(xb), k, replace=False))]
        ma.append(xa); mb.append(xb)
    return (pd.concat(ma), pd.concat(mb), int(st.nunique())) if ma else (A.iloc[:0], B.iloc[:0], int(st.nunique()))
def val(A, B, kind):
    if kind in ("HIF", "HALL"): return float(A[kind].mean() - B[kind].mean())
    W = W_h if kind == "TF" else W_dc
    return ulm(pb(A.row.values), W) - ulm(pb(B.row.values), W)

rec = []   # 每单位/患者一条：比较、指标、不匹配值、匹配值、配对数、检出基因中位数
def add(comp, kind, key, A, B, thr):
    r = {"comparison": comp, "metric": kind, "key": key, "nA": len(A), "nB": len(B)}
    if len(A) >= thr and len(B) >= thr:
        r["unmatched"] = val(A, B, kind); r["ngA"] = float(A.ngenes.median()); r["ngB"] = float(B.ngenes.median())
        mA, mB, ns = match(A, B); r["pairs"] = len(mA); r["strata"] = ns
        if len(mA) >= thr:
            r["matched"] = val(mA, mB, kind); r["ngA_m"] = float(mA.ngenes.median()); r["ngB_m"] = float(mB.ngenes.median())
    rec.append(r)
units = sorted(o.groupby(["patient", "grp"]), key=lambda kv: kv[0])
for (p, g), u in units:
    pos, neg = u[u.spp1 > 0], u[u.spp1 == 0]
    q1, q3 = u.STATE.quantile([.25, .75]); hi, lo = u[u.STATE >= q3], u[u.STATE <= q1]
    for k in ["HIF", "HALL"]:
        add("H1a-i", k, (p, g), pos, neg, 20); add("H1a-ii", k, (p, g), hi, lo, 20)
    add("R1", "TF", (p, g), pos, neg, 20); add("R1-decirc", "TFdc", (p, g), pos, neg, 20)
byp = {}
for (p, g), u in units: byp.setdefault(p, {})[g] = u
for p in sorted(byp):
    if not {"LUAD", "Precursor"} <= set(byp[p]): continue
    L, Pr = byp[p]["LUAD"], byp[p]["Precursor"]
    topL, topP = L[L.STATE >= L.STATE.quantile(.75)], Pr[Pr.STATE >= Pr.STATE.quantile(.75)]
    for k in ["HIF", "HALL"]:
        add("H1b", k, (p,), L, Pr, 50); add("H1c-ii", k, (p,), topL, topP, 20)
    add("R2", "TF", (p,), L, Pr, 50)
D = pd.DataFrame(rec); D["patient"] = D.key.map(lambda t: t[0]); D["key"] = D.key.astype(str)
D.to_csv(OUT + "D6_units.tsv", sep="\t", index=False)

REPRO = {("H1a-i", "HIF"): H1["H1a_i_HIF"], ("H1a-i", "HALL"): H1["H1a_i_HALL"], ("H1a-ii", "HIF"): H1["H1a_ii_HIF"], ("H1a-ii", "HALL"): H1["H1a_ii_HALL"],
         ("H1b", "HIF"): H1["H1b_HIF"], ("H1b", "HALL"): H1["H1b_HALL"], ("H1c-ii", "HIF"): H1c["H1c_ii_HIF"], ("H1c-ii", "HALL"): H1c["H1c_ii_HALL"],
         ("R1", "TF"): RR["R1"], ("R1-decirc", "TFdc"): RR["R1_decirc"], ("R2", "TF"): RR["R2"]}
out, pv = [], {}
for (comp, kind), ref in REPRO.items():
    d = D[(D.comparison == comp) & (D.metric == kind)]
    pu = d.groupby("patient").unmatched.mean().dropna()
    got = (int(len(pu)), int((pu > 0).sum())); exp = (int(ref["n"]), int(ref["n_pos"]))
    if got != exp: sys.exit(f"复现失败 {comp} {kind}: 得 {got}，原报告 {exp}")
    pm = d.groupby("patient").matched.mean().dropna(); pu_same = pu.loc[pm.index]
    r = {"comparison": comp, "metric": {"HIF": "HIF13", "HALL": "HALLMARK", "TF": "HIF1A", "TFdc": "HIF1A (hypoxia genes removed)"}[kind],
         "reproduced_n_npos": got, "unmatched_n": got[0], "unmatched_npos": got[1], "unmatched_median": float(pu.median()),
         "unmatched_p_one_sided_recomputed": float(stats.wilcoxon(pu.values, alternative="greater").pvalue),
         "matched_n": int(len(pm)), "matched_npos": int((pm > 0).sum()), "matched_median": float(pm.median()) if len(pm) else None,
         "unmatched_median_same_patients": float(pu_same.median()) if len(pm) else None,
         "pairs_median": float(d.pairs.median()), "units_not_evaluable_after_matching": int(d.unmatched.notna().sum() - d.matched.notna().sum()),
         "ngenes_median_A_before": float(d.ngA.median()), "ngenes_median_B_before": float(d.ngB.median()),
         "ngenes_median_A_after": float(d.ngA_m.median()), "ngenes_median_B_after": float(d.ngB_m.median())}
    r["retention"] = (r["matched_median"] / r["unmatched_median_same_patients"]) if len(pm) and r["unmatched_median_same_patients"] else None
    r["matched_p_one_sided"] = float(stats.wilcoxon(pm.values, alternative="greater").pvalue) if len(pm) >= 5 else None
    pv[f"{comp}|{kind}"] = r["matched_p_one_sided"]; out.append(r)
assert len(pv) == 11
adj = holm_fixed(pv, 11)   # D6 更正：检验数为 11（原文误写 13）
for r in out:
    kind = {"HIF13": "HIF", "HALLMARK": "HALL", "HIF1A": "TF", "HIF1A (hypoxia genes removed)": "TFdc"}[r["metric"]]
    q = adj[f"{r['comparison']}|{kind}"]; r["matched_holm_p_m11"] = q
    med, ret = r["matched_median"], r["retention"]
    if med is None or q is None: r["category"] = "not evaluable"
    elif med > 0 and q < 0.05: r["category"] = "robust"
    elif med > 0 and ret is not None and ret >= 0.5: r["category"] = "attenuated but retained"
    else: r["category"] = "largely explained by detection depth"
R = pd.DataFrame(out); R.to_csv(OUT + "D6_summary.tsv", sep="\t", index=False)
dump({"plan_sha256": sha, "deviation_record_sha256_at_run": dev_sha, "deviation": "D6 (post hoc)", "seed": SEED, "results": out}, OUT + "D6_results.json")
log("D6 完成"); print(R[["comparison", "metric", "unmatched_npos", "unmatched_n", "unmatched_median", "matched_npos", "matched_n", "matched_median",
                       "retention", "matched_p_one_sided", "matched_holm_p_m11", "category"]].to_string())
