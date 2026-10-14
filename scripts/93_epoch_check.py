"""§12.7 轮数验证：比较 P1_LUAD 合并参考 10,000 轮 vs 30,000 轮。输出 deconv/epoch_check.json，并打印 USE_EPOCHS=..."""
import json, numpy as np, pandas as pd
D="/Users/dudu/全部数据库/GGN_target/results/30_spatial/deconv/"
a=pd.read_csv(D+"map_merged_n6/GSM9226169_P1_LUAD_q05.tsv.gz",sep="\t",index_col=0)
b=pd.read_csv(D+"map_merged_n6_e10000/GSM9226169_P1_LUAD_q05.tsv.gz",sep="\t",index_col=0).loc[a.index,a.columns]
pa=a.div(a.sum(1),axis=0); pb=b.div(b.sum(1),axis=0)
share=a.sum(0)/a.values.sum()
big=share[share>=0.01].index
r={l:float(np.corrcoef(pa[l],pb[l])[0,1]) for l in big}
key=["KACCAS-high","Mac-act","Mac-other","AT2L-other"]
sb=b.sum(0)/b.values.sum()
rel={l:float(abs(sb[l]-share[l])/share[l]) for l in key}
c1=all(v>=0.95 for v in r.values()); c2=all(v<0.10 for v in rel.values())
out={"corr_by_label":r,"rel_diff_key":rel,"cond1":c1,"cond2":c2,"use_epochs":10000 if (c1 and c2) else 30000}
json.dump(out,open(D+"epoch_check.json","w"),indent=1); print(json.dumps(out,indent=1)); print(f"USE_EPOCHS={out['use_epochs']}")
