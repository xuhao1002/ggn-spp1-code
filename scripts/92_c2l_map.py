"""方案 §2/§12.1：GSE307534 逐样本 cell2location 映射。可续跑（已有结果的样本跳过）。
用法: python 92_c2l_map.py merged|unmerged NCELLS [样本名 ...|all]   环境 ggn_spatial。"""
import matplotlib; matplotlib.use("Agg")
import sys, time, json, os, warnings
sys.path.insert(0,"/Users/dudu/全部数据库/GGN_target/scripts")
warnings.filterwarnings("ignore")
import pandas as pd
from c2l_common import R, read_st, st_samples, load_sig, run_c2l
tag=sys.argv[1]; ncells=int(sys.argv[2]); which=sys.argv[3:] or ["all"]
EPOCHS=int(os.environ.get("C2L_EPOCHS","30000")); suf="" if EPOCHS==30000 else f"_e{EPOCHS}"
sigp={"merged":"ref_merged","unmerged":"ref"}[tag]
sig=load_sig(R+f"results/30_spatial/deconv/{sigp}/signatures_per_label.csv")
OUT=R+f"results/30_spatial/deconv/map_{tag}_n{ncells}{suf}/"; os.makedirs(OUT,exist_ok=True)
t0=time.time()
for b,d,p,s in st_samples():
    if which!=["all"] and b not in which: continue
    f=OUT+f"{b}_q05.tsv.gz"
    if os.path.exists(f): print("skip",b,flush=True); continue
    st=read_st(d)
    q,qc=run_c2l(st,sig,ncells=ncells,alpha=20,epochs=EPOCHS)
    q.to_csv(f,sep="\t")
    qc.update(epochs=EPOCHS,sample=b,patient=p,stage=s,signature=tag,ncells=ncells)
    with open(OUT+"qc.jsonl","a") as fh: fh.write(json.dumps(qc,ensure_ascii=False)+"\n")
    print(f"[{time.time()-t0:7.0f}s]",b,qc,flush=True)
print("DONE",flush=True)
