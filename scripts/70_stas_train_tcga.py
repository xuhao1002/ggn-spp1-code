"""STAS 分子评分 P1 第1步：仅用 TCGA-LUAD 训练并冻结系数（不接触 CPTAC）。
预设方案：docs/STAS分子评分_预设方案_冻结_2026-10-11.md
实现细节（预设方案未规定的）：外层交叉验证 5 折×5 次重复；样本取每名患者首个 -01 肿瘤样本；
log2(TPM+1)；蛋白编码基因中方差前 5000；TCGA 内 z 标准化；弹性网络 l1_ratio=0.5，C 用 10 折内层 CV（AUC）。
"""
import gzip, hashlib, json, re, sys
import numpy as np, pandas as pd
from sklearn.linear_model import LogisticRegressionCV
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.metrics import roc_auc_score

R = '/Users/dudu/全部数据库/GGN_target/'
OUT = R + 'results/70_stas/'
import os; os.makedirs(OUT, exist_ok=True)
SEED = 20261011

lab = pd.read_excel(R + 'raw/STAS_labels/STAS_labels_TCGA_from_PanLR_20261011.xlsx')
y_pt = lab.groupby('患者编号').label.max()

hdr = pd.read_csv(R + 'raw/TCGA_LUAD/TCGA-LUAD.star_tpm.tsv.gz', sep='\t', nrows=0).columns
pick = {}
for c in hdr[1:]:
    p = c[:12]
    if p in y_pt.index and c.split('-')[3].startswith('01') and p not in pick:
        pick[p] = c
print('患者', len(y_pt), '有 01 肿瘤 RNA', len(pick))
pts = sorted(pick)
tpm = pd.read_csv(R + 'raw/TCGA_LUAD/TCGA-LUAD.star_tpm.tsv.gz', sep='\t', usecols=['Ensembl_ID'] + [pick[p] for p in pts], index_col=0)
tpm.columns = pts
tpm.index = tpm.index.str.split('.').str[0]
tpm = tpm[~tpm.index.duplicated()]

pc = set()
with gzip.open(R + 'raw/annotation/gencode.v44.basic.annotation.gtf.gz', 'rt') as f:
    for line in f:
        if '\tgene\t' in line and 'gene_type "protein_coding"' in line:
            pc.add(re.search(r'gene_id "([^.]+)', line).group(1))
names = {}
with gzip.open(R + 'raw/annotation/gencode.v44.basic.annotation.gtf.gz', 'rt') as f:
    for line in f:
        if '\tgene\t' in line:
            m = re.search(r'gene_id "([^.]+)\.[^"]*"; gene_type "[^"]*"; gene_name "([^"]+)"', line)
            if m: names[m.group(1)] = m.group(2)
tpm = tpm.loc[tpm.index.isin(pc)]
X = np.log2(tpm.T + 1)
top = X.var().sort_values(ascending=False).index[:5000]
X = X[top]
mu, sd = X.mean(), X.std(ddof=0).replace(0, 1)
Z = ((X - mu) / sd).values
y = y_pt.loc[pts].values.astype(int)
print('X', Z.shape, '阳性', y.sum(), '阴性', (1 - y).sum())


def fit(Zm, ym, seed):
    return LogisticRegressionCV(Cs=10, cv=10, penalty='elasticnet', solver='saga', l1_ratios=[0.5],
                                scoring='roc_auc', max_iter=3000, tol=1e-3, class_weight='balanced',
                                random_state=seed, n_jobs=-1, refit=True).fit(Zm, ym)

# 内部嵌套 CV（外层 5 折 × 5 次）
rs = RepeatedStratifiedKFold(n_splits=5, n_repeats=5, random_state=SEED)
aucs = []
for i, (tr, te) in enumerate(rs.split(Z, y)):
    m = fit(Z[tr], y[tr], SEED + i)
    aucs.append(roc_auc_score(y[te], m.decision_function(Z[te])))
    print('outer', i, round(aucs[-1], 3), flush=True)
aucs = np.array(aucs)
rep = aucs.reshape(5, 5).mean(1)
print('TCGA 内部嵌套 CV AUC 均值', aucs.mean().round(3), '重复均值范围', rep.min().round(3), rep.max().round(3))

# 最终模型（全部 TCGA）
m = fit(Z, y, SEED)
coef = pd.Series(m.coef_[0], index=top)
nz = coef[coef != 0]
tab = pd.DataFrame({'ensembl': nz.index, 'symbol': [names.get(g, '') for g in nz.index], 'coef': nz.values})
tab = tab.reindex(tab.coef.abs().sort_values(ascending=False).index)
tab.to_csv(OUT + 'score_full_coefficients.tsv', sep='\t', index=False)
print('完整模型非零基因', len(nz), 'C', m.C_[0])

# 精简版（≤20 基因）
t20 = tab.head(20).ensembl.tolist()
Z20 = ((X[t20] - mu[t20]) / sd[t20]).values
m20 = fit(Z20, y, SEED)
c20 = pd.DataFrame({'ensembl': t20, 'symbol': [names.get(g, '') for g in t20], 'coef': m20.coef_[0]})
c20.to_csv(OUT + 'score_compact20_coefficients.tsv', sep='\t', index=False)
json.dump({'intercept_full': float(m.intercept_[0]), 'intercept_compact20': float(m20.intercept_[0]),
           'C_full': float(m.C_[0]), 'n_nonzero_full': int(len(nz)),
           'tcga_nested_cv_auc_mean': float(aucs.mean()), 'tcga_nested_cv_auc_repeat_means': rep.tolist(),
           'top5000_genes': list(top), 'seed': SEED}, open(OUT + 'train_meta.json', 'w'))

with open(OUT + 'frozen_sha256.txt', 'w') as f:
    for fn in ['score_full_coefficients.tsv', 'score_compact20_coefficients.tsv', 'train_meta.json']:
        f.write(hashlib.sha256(open(OUT + fn, 'rb').read()).hexdigest() + '  ' + fn + '\n')
print(open(OUT + 'frozen_sha256.txt').read())
