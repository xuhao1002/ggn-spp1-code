"""A1–A6 公用头：方案哈希校验与日志。预设方案 v1.3，锁定于 2026-10-19。"""
import os, sys, json, time, hashlib
import numpy as np, pandas as pd

P = os.path.expanduser("~/全部数据库/GGN_target")
PLAN = f"{P}/docs/预设方案_外审回应A1-A6_v1.3_2026-10-19.md"
PLAN_SHA = "1cde317f3bccac30050ddb8239ef6d623090f8aab8203765396548f9a136f94d"
SEED = 20261019

def check_plan():
    """方案哈希不符即拒绝运行。方案锁定后不得修改；偏离另建记录。"""
    sha = hashlib.sha256(open(PLAN, "rb").read()).hexdigest()
    if sha != PLAN_SHA:
        sys.exit(f"方案哈希不符，拒绝运行。\n  期望 {PLAN_SHA}\n  实际 {sha}\n"
                 f"方案锁定后不应被修改；若确需修订，请另建版本并同步更新脚本中的 PLAN_SHA。")
    return sha

_t0 = time.time()
def log(*a):
    print(f"[{time.time()-_t0:7.1f}s]", *a, flush=True)

def outdir(name):
    """每次运行写入独立目录 results/<name>/run_<本地时间>/，目录已存在即拒绝运行（v1.3 0.3：原始结果不得覆盖）。
    2026-10-20 起生效；此前的输出位于 results/<name>/ 根目录，保持原样。"""
    import datetime
    stamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    d = f"{P}/results/{name}/run_{stamp}/"
    if os.path.exists(d):
        sys.exit(f"输出目录已存在，拒绝覆盖：{d}")
    os.makedirs(d)
    with open(d + "RUN_INFO.txt", "w", encoding="utf-8") as fh:
        fh.write(f"script\t{os.path.abspath(sys.argv[0])}\nstarted\t{stamp}\npython\t{sys.executable}\nplan_sha256\t{PLAN_SHA}\n")
    log("输出目录", d)
    return d

def dump(obj, path):
    json.dump(obj, open(path, "w"), indent=1, ensure_ascii=False, default=float)
    log("写出", path)

def holm_fixed(pvals, m):
    """固定族大小 m 的 Holm 校正（方案 6.3）。
    pvals 为 dict：键 -> P 值或 None。None 表示无法评估，以 P=1 占位占满一个族位置，
    但返回值中该键仍为 None，绝不报告为 P=1。"""
    keys = list(pvals)
    filled = [(k, (1.0 if pvals[k] is None else float(pvals[k]))) for k in keys]
    assert len(filled) <= m, f"族成员数 {len(filled)} 超过预设 m={m}"
    while len(filled) < m:                      # 族成员不足时同样以 P=1 占位
        filled.append((f"__placeholder{len(filled)}", 1.0))
    order = sorted(range(len(filled)), key=lambda i: filled[i][1])
    adj = [None]*len(filled); run = 0.0
    for rank, i in enumerate(order):
        v = min(1.0, (m - rank) * filled[i][1]); run = max(run, v); adj[i] = run
    out = {}
    for i, (k, _) in enumerate(filled):
        if k.startswith("__placeholder"): continue
        out[k] = None if pvals[k] is None else adj[i]
    return out

def boot_ci(x, n=10000, seed=SEED):
    """患者级差值的自助法 95% CI（方案 6.2）。"""
    x = np.asarray([v for v in x if np.isfinite(v)], float)
    if len(x) < 3: return (np.nan, np.nan)
    rng = np.random.default_rng(seed)
    bs = rng.choice(x, size=(n, len(x)), replace=True).mean(1)
    return float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))
