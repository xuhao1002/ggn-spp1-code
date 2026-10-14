import json, re, numpy as np, pandas as pd, matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from scipy import stats
# ---------------------------------------------------------------------------
# Paths. Nothing below is hard-coded to the authors' machine.
#   SPP1_RESULTS   directory holding the original results tree (optional)
#   SPP1_RAW       directory holding the raw inputs (optional)
#   SPP1_SOURCE_DATA  directory holding the released source-data files
#                     (default: ../source_data next to this script)
#   SPP1_OUT       output directory (default: ../figures_out)
# Any input is looked up first under SPP1_RESULTS / SPP1_RAW using its original
# relative path, then in SPP1_SOURCE_DATA under the flattened name in which it
# is released (path separators replaced by "__").
# ---------------------------------------------------------------------------
import os as _os
_HERE = _os.path.dirname(_os.path.abspath(__file__))
R   = _os.environ.get("SPP1_RESULTS", "")
RAW = _os.environ.get("SPP1_RAW", "")
SRC = _os.environ.get("SPP1_SOURCE_DATA", _os.path.join(_HERE, _os.pardir, "source_data"))
OUT = _os.environ.get("SPP1_OUT", _os.path.join(_HERE, _os.pardir, "figures_out"))
if not OUT.endswith(_os.sep): OUT += _os.sep
_os.makedirs(OUT, exist_ok=True)

def P(rel, base=""):
    """Resolve an input file: original tree first, released source_data second."""
    if base:
        cand = _os.path.join(base, rel)
        if _os.path.exists(cand): return cand
    flat = _os.path.join(SRC, rel.replace("/", "__"))
    if _os.path.exists(flat): return flat
    raise FileNotFoundError(
        "Could not find %r. Set SPP1_RESULTS/SPP1_RAW to the original directories, "
        "or SPP1_SOURCE_DATA to the released source_data folder." % rel)
mpl.rcParams.update({"font.family":"DejaVu Sans","font.size":7,"axes.titlesize":7.5,"axes.labelsize":7,"xtick.labelsize":6.5,"ytick.labelsize":6.5,
 "legend.fontsize":6.5,"axes.spines.top":False,"axes.spines.right":False,"axes.linewidth":0.6,"xtick.major.width":0.6,"ytick.major.width":0.6,
 "xtick.major.size":2.5,"ytick.major.size":2.5,"pdf.fonttype":42,"svg.fonttype":"none","savefig.dpi":300})
PRE="#2a78d6"; INV="#eb6834"; NOR="#8a8984"; OK="#1baf7a"; NO="#b9b8b3"; INK="#0b0b0b"; MUT="#52514e"
W=7.09  # 180 mm double column
def lab(ax,s,dx=-0.18,dy=1.06): ax.text(dx,dy,s.lower()+")",transform=ax.transAxes,fontsize=9,fontweight="bold",va="bottom",ha="left")
def save(fig,name):
    fig.savefig(OUT+name+".png",dpi=300,bbox_inches="tight"); fig.savefig(OUT+name+".pdf",bbox_inches="tight"); fig.savefig(OUT+name+".tiff",dpi=300,bbox_inches="tight"); plt.close(fig)
def paired(ax,a,b,c1,c2,labels,ylab,title):
    for x,y in zip(a,b): ax.plot([0,1],[x,y],color="#c9c8c3",lw=0.6,zorder=1)
    ax.scatter(np.zeros(len(a)),a,s=12,color=c1,zorder=2,edgecolor="white",lw=0.4); ax.scatter(np.ones(len(b)),b,s=12,color=c2,zorder=2,edgecolor="white",lw=0.4)
    for x,v,c in [(0,a,c1),(1,b,c2)]: ax.plot([x-0.18,x+0.18],[np.median(v)]*2,color=INK,lw=1.2,zorder=3)
    ax.set_xticks([0,1]); ax.set_xticklabels(labels); ax.set_xlim(-0.45,1.45); ax.set_ylabel(ylab); ax.set_title(title,loc="left")
def ptxt(ax,s,y=0.98): ax.text(0.5,y,s,transform=ax.transAxes,ha="center",va="top",fontsize=6.3,color=MUT)

# ---------------- Figure 1: design ----------------
fig=plt.figure(figsize=(W,2.9)); ax=fig.add_axes([0,0,1,1]); ax.set_xlim(0,100); ax.set_ylim(0,42); ax.axis("off")
def box(x,y,w,h,title,body,fc="#f4f3ef",ec="#c9c8c3",tc=INK):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.4,rounding_size=1.2",fc=fc,ec=ec,lw=0.7))
    ax.text(x+1.2,y+h-1.6,title,fontsize=7.2,fontweight="bold",va="top",color=tc); ax.text(x+1.2,y+h-4.6,body,fontsize=6.2,va="top",color=MUT,linespacing=1.35)
# progression strip
stages=[("Normal",NOR),("AAH",PRE),("AIS",PRE),("MIA",PRE),("Invasive LUAD",INV)]
for i,(s,c) in enumerate(stages):
    x=4+i*18.5; ax.add_patch(FancyBboxPatch((x,36),15.5,4.4,boxstyle="round,pad=0.3,rounding_size=1",fc=c,ec="none",alpha=0.9))
    ax.text(x+7.75,38.2,s,ha="center",va="center",fontsize=7,color="white",fontweight="bold")
    if i<4: ax.annotate("",xy=(x+18.2,38.2),xytext=(x+15.9,38.2),arrowprops=dict(arrowstyle="-|>",color=MUT,lw=0.7))
ax.text(50,32.6,"Same-patient pairs: precursor (AAH/AIS/MIA) vs invasive LUAD, compared cross-sectionally  |  staged analysis plans, hash-locked before each analysis",ha="center",fontsize=6.5,color=MUT)
box(2,16.5,30,13,"Discovery: snRNA-seq (GSE308103)","23 patients, 75 samples, 700,199 nuclei\n22 precursor–invasive pairs\nPaired pseudobulk (~ patient + condition)")
box(35,16.5,30,13,"Validation: RNA and space","scRNA-seq: GSE164789 (62), GSE189357 (9)\nVisium (GSE307534): 25 patients, 56 sections\nXenium (Takano 2024): 15 lesions")
box(68,16.5,30,13,"Validation: protein and outcome","CPTAC-LUAD tumour vs adjacent:\n  Discovery 65 + Confirmatory 68 pairs\n7 stage I cohorts (n = 1,163)")
box(2,1.5,30,11,"1  Reproducibility audit","KAC/CAS expansion, IL-1β niche,\nTIM-3, CD47 (Fig. 2)",fc="#eef4fc",ec="#9cc0ea")
box(35,1.5,30,11,"2  SPP1+ macrophage state","State vs marker score, source,\nhypoxia program (Figs. 3–5)",fc="#fdf0ea",ec="#f2b393")
box(68,1.5,30,11,"3  Clinical association & limits","Stage I outcome, Mendelian randomization,\nrenal expression (Fig. 6)",fc="#ecf7f2",ec="#97d7bd")
save(fig,"Figure1")

# ---------------- Figure 2: audit ----------------
fr=pd.read_csv(P("40_drivers/stepwise/stepwise_fractions_per_patient.tsv", R),sep="\t")
B=json.load(open(P("sens_protein/B_il1b_power.json", R)))
fig=plt.figure(figsize=(W,5.0)); gs=fig.add_gridspec(2,6,height_ratios=[1,1.05],hspace=0.95,wspace=1.6)
axA=[fig.add_subplot(gs[0,2*k:2*k+2]) for k in range(3)]
for ax,ind,t,key in [(axA[0],"KAC_frac","KAC-high","AISMIA|KAC_frac"),(axA[1],"CAS_frac","CAS-high","AISMIA|CAS_frac"),(axA[2],"IL1Bpos_macrophage_frac","IL1B+ macrophages","AISMIA|IL1Bpos_macrophage_frac")]:
    s_=fr[(fr.step=="AISMIA")&(fr.indicator==ind)]
    paired(ax,s_.normal*100,s_.precursor*100,NOR,PRE,["Normal","AIS/MIA"],"% of lineage",t)
    v=B[key]; ptxt(ax,f"Δ {v['mean_diff']*100:+.1f} pp ({v['ci95'][0]*100:.1f} to {v['ci95'][1]*100:.1f})\n{int((s_['diff']>0).sum())}/{len(s_)} patients higher",y=-0.2)
for ax,l in zip(axA,"ABC"): lab(ax,l,dx=-0.3)
ax=fig.add_subplot(gs[1,1:3]); rows=[("KAC-high, AIS/MIA","AISMIA|KAC_frac"),("CAS-high, AIS/MIA","AISMIA|CAS_frac"),("IL1B+ mac., AIS/MIA","AISMIA|IL1Bpos_macrophage_frac"),("KAC-high, AAH","AAH|KAC_frac"),("CAS-high, AAH","AAH|CAS_frac"),("IL1B+ mac., AAH","AAH|IL1Bpos_macrophage_frac")]
for i,(n,k) in enumerate(rows):
    v=B[k]; y=len(rows)-1-i; c=OK if v["ci95"][0]>0 else (NO)
    ax.plot([v["ci95"][0]*100,v["ci95"][1]*100],[y,y],color=c,lw=1.6,solid_capstyle="round"); ax.scatter(v["mean_diff"]*100,y,s=18,color=c,zorder=3)
    if "IL1B" in k: ax.scatter(v["upper95_one_sided"]*100,y,marker="|",s=60,color=INK,zorder=4)
ax.axvline(0,color=MUT,lw=0.6,ls=(0,(2,2))); ax.set_yticks(range(len(rows))); ax.set_yticklabels([r[0] for r in rows][::-1]); ax.set_xlabel("Change vs paired normal (pp, 95% CI)\n| = one-sided 95% upper bound")
ax.set_title("Effect sizes by stage",loc="left"); lab(ax,"D",dx=-0.55)
claims=["KAC expansion","CAS expansion","CD47 ↑ (AT2 lineage)","IL-1β+ macrophage ↑","IL1R1 ↑ (epithelium)","HAVCR2/TIM-3 ↑","CEACAM1 ↑ (ligand)"]
cols=["Peng\nsnRNA","GSE\n164789","Bulk\ncohorts","Peng\nVisium"]
G=[[1,0,1,9],[1,0,1,9],[1,9,2,9],[0,0,0,0],[0,0,0,9],[0,0,0,9],[1,1,1,9]]
cmap={1:OK,2:"#a6dcc7",0:"#e9e8e4",9:"white"}; txt={1:"✓",2:"~",0:"✗",9:"–"}
ax=fig.add_subplot(gs[1,4:6])
for i,r in enumerate(G):
    for j,v in enumerate(r):
        ax.add_patch(plt.Rectangle((j,len(G)-1-i),0.94,0.9,fc=cmap[v],ec="#d8d7d2" if v==9 else "none",lw=0.5)); ax.text(j+0.47,len(G)-1-i+0.45,txt[v],ha="center",va="center",fontsize=7.5,color=INK if v!=1 else "white")
ax.set_xlim(0,len(cols)); ax.set_ylim(0,len(G)); ax.set_xticks(np.arange(len(cols))+0.47); ax.set_xticklabels(cols,fontsize=5.8); ax.xaxis.tick_top()
ax.set_yticks(np.arange(len(G))+0.45); ax.set_yticklabels(claims[::-1]); [sp.set_visible(False) for sp in ax.spines.values()]; ax.tick_params(length=0)
ax.text(-0.2,-0.9,"✓ supported   ~ direction only   ✗ not supported   – not tested",fontsize=5.8,color=MUT,transform=ax.transData,ha="center" if False else "left"); lab(ax,"E",dx=-1.05,dy=1.2)
save(fig,"Figure2")

# ---------------- Figure 3: SPP1 state (RNA) ----------------
mc=pd.read_csv(P("14_microenv/microenv_layer1_candidates.tsv", R),sep="\t"); m=mc[(mc.gene=="SPP1")&(mc.lineage=="Myeloid")].iloc[0]
def se_from(l,p1): z=stats.norm.isf(p1); return l/z
est=[("Peng snRNA\n22 paired patients",m.log2FoldChange,m.lfcSE),("GSE164789\nscRNA",m.GSE164789_lfc,se_from(m.GSE164789_lfc,m.GSE164789_p1)),("GSE189357\nscRNA",m.GSE189357_lfc,se_from(m.GSE189357_lfc,m.GSE189357_p1))]
X=pd.read_csv(P("xenium/xenium_per_section.tsv", R),sep="\t")
fig,axs=plt.subplots(1,2,figsize=(W*0.72,2.3),gridspec_kw={"width_ratios":[1.4,1],"wspace":0.6})
ax=axs[0]
for i,(n,l,s) in enumerate(est):
    y=2-i; ax.plot([l-1.96*s,l+1.96*s],[y,y],color=INV,lw=1.6,solid_capstyle="round"); ax.scatter(l,y,s=24,color=INV,zorder=3); ax.text(l,y+0.2,f"{l:.2f}",ha="center",fontsize=6.3,color=MUT)
ax.axvline(0,color=MUT,lw=0.6,ls=(0,(2,2))); ax.set_yticks([2,1,0]); ax.set_yticklabels([e[0] for e in est]); ax.set_xlabel("Myeloid SPP1 log2 fold change\ninvasive vs precursor (95% CI)")
ax.set_ylim(-0.5,2.7); ax.set_title("SPP1 induction in myeloid cells",loc="left"); lab(ax,"A",dx=-0.55)
ax=axs[1]; rng=np.random.default_rng(1)
for j,(g,c) in enumerate([("AIS/MIA",PRE),("invasive",INV)]):
    s=X[X.group==g]; xs=j+rng.uniform(-0.12,0.12,len(s))
    for x,(_,r) in zip(xs,s.iterrows()): ax.scatter(x,r.X2_SPP1pos_frac_of_macrophages*100,s=22,color=c if r.fixation=="FF" else "white",edgecolor=c,lw=0.9,zorder=3)
    ax.plot([j-0.22,j+0.22],[s.X2_SPP1pos_frac_of_macrophages.median()*100]*2,color=INK,lw=1.2)
ax.set_xticks([0,1]); ax.set_xticklabels(["AIS/MIA\n(n=10)","Invasive\n(n=5)"]); ax.set_xlim(-0.5,1.5); ax.set_ylabel("SPP1+ among macrophages (%)")
ax.set_title("Independent Xenium cohort",loc="left"); ptxt(ax,"median 8.5% vs 30%\nP = 0.08 one-sided (0.16 two-sided)",y=-0.3); ax.text(1.02,0.02,"○ FFPE\n● FF",transform=ax.transAxes,fontsize=6,color=MUT); lab(ax,"B",dx=-0.4)
save(fig,"Figure3")

# ---------------- Figure 4: state not abundance; source ----------------
E2=pd.read_csv(P("30_spatial/exploratory/E2_patient_pairs.tsv", R),sep="\t")
# CPTAC per pair
CPTAC=lambda f: P("CPTAC_LUAD_PDC/"+f, RAW); mp=pd.read_csv(CPTAC("aliquot_map_labeled151.tsv"),sep="\t"); MAC=["CD68","CD163","C1QA","CD14","MARCO"]
def load(f):
    d=pd.read_csv(CPTAC(f),sep="\t",index_col=0); d=d.drop([i for i in ["Mean","Median","StdDev"] if i in d.index]); keep={}
    for c in d.columns:
        if "Unshared" in c: continue
        mm=re.search(r"CPT\d{10}",c)
        if mm: keep[mm.group(0)]=c
    d=d[list(keep.values())]; d.columns=list(keep.keys()); return d
cpt={}
for nm,f in [("Discovery","PDC000153_CPTAC3_Lung_Adeno_Carcinoma_Proteome.tmt10.tsv"),("Confirmatory","PDC000489_CPTAC3_LUAD_Confirmatory_Proteome.tmt11.tsv")]:
    d=load(f); mm=mp[mp.aliquot.isin(d.columns)]
    T=mm[mm.sample_type=="Primary Tumor"].drop_duplicates("case").set_index("case").aliquot; N=mm[mm.sample_type=="Solid Tissue Normal"].drop_duplicates("case").set_index("case").aliquot
    cs=sorted(set(T.index)&set(N.index)); cols=list(T[cs])+list(N[cs])
    Z=pd.DataFrame({g:(d.loc[g,cols].astype(float)-d.loc[g,cols].astype(float).mean())/d.loc[g,cols].astype(float).std() for g in MAC}).mean(axis=1)
    cpt[nm]=pd.DataFrame({"t":d.loc["SPP1",T[cs]].astype(float).values,"n":d.loc["SPP1",N[cs]].astype(float).values,"dm":Z[T[cs]].values-Z[N[cs]].values})
XS=pd.read_csv(P("xenium/xenium_checkT_per_section.tsv", R),sep="\t")
fig=plt.figure(figsize=(W,5.0)); gs=fig.add_gridspec(2,6,height_ratios=[1,1.1],hspace=0.75,wspace=1.5)
axs=[fig.add_subplot(gs[0,0:2]),fig.add_subplot(gs[0,2:4]),fig.add_subplot(gs[0,4:6]),fig.add_subplot(gs[1,0:3]),fig.add_subplot(gs[1,4:6])]
paired(axs[0],E2.M_pre,E2.M_lu,PRE,INV,["Precursor","Invasive"],"SPP1 / control mac. genes","Visium, 24 patients"); ptxt(axs[0],"21/24 higher; Holm P = 4.5×10⁻⁵",y=-0.2); lab(axs[0],"A",dx=-0.35)
for k,(nm,ax) in enumerate(zip(["Discovery","Confirmatory"],axs[1:3])):
    c=cpt[nm]; paired(ax,c.n,c.t,NOR,INV,["Adjacent","Tumour"],"SPP1 protein (log2 ratio)" if k==0 else "",f"CPTAC {nm}, {len(c)} pairs")
    ptxt(ax,"P = 8×10⁻¹⁰" if k==0 else "P = 1×10⁻⁷",y=-0.2); lab(ax,"BC"[k],dx=-0.35)
ax=axs[3]
for nm,c,mk in [("Discovery",INV,"o"),("Confirmatory","#6250d6","s")]:
    q=cpt[nm]; y=q.t-q.n; ax.scatter(q.dm,y,s=9,color=c,marker=mk,alpha=0.8,lw=0,label=nm)
    b=np.polyfit(q.dm,y,1); xx=np.linspace(q.dm.min(),q.dm.max(),10); ax.plot(xx,np.polyval(b,xx),color=c,lw=1)
ax.axhline(0,color=MUT,lw=0.5,ls=(0,(2,2))); ax.axvline(0,color=MUT,lw=0.5,ls=(0,(2,2)))
ax.set_xlabel("Δ macrophage marker score (SD)"); ax.set_ylabel("Δ SPP1 (tumour − adjacent)"); ax.set_title("Rise beyond the macrophage marker score",loc="left")
ax.legend(frameon=False,loc="upper left",handletextpad=0.2); ptxt(ax,"Intercept > 0: P = 1×10⁻⁹ (Discovery), 2×10⁻⁷ (Confirmatory)",y=-0.24); lab(ax,"D",dx=-0.22)
ax=axs[4]; XF=pd.read_csv(P("xenium/xenium_source_fixed.tsv"),sep="\t")
XF["g"]=pd.Categorical(XF.group,["AIS/MIA","invasive"]); XF=XF.sort_values(["g","mac_c"],ascending=[True,False]).reset_index(drop=True)
xs=[];x=0;prev=None
for r in XF.itertuples():
    if prev is not None and r.group!=prev: x+=1.5
    xs.append(x); x+=1; prev=r.group
xs=np.array(xs); oth_c=1-XF.mac_c-XF.epi_c
ax.bar(xs,XF.mac_c*100,color=INV,width=0.82,lw=0); ax.bar(xs,XF.epi_c*100,bottom=XF.mac_c*100,color=PRE,width=0.82,lw=0); ax.bar(xs,oth_c*100,bottom=(XF.mac_c+XF.epi_c)*100,color="#d8d7d2",width=0.82,lw=0)
ax.scatter(xs,XF.mac_r*100,s=10,facecolor="white",edgecolor=INK,lw=0.6,zorder=4)
ax.set_xticks(xs); ax.set_xticklabels([f"{sec}{' (F)' if fx=='FFPE' else ''}" for sec,fx in zip(XF.section,XF.fixation)],fontsize=4.6,rotation=90,ha="center"); ax.tick_params(axis="x",length=0,pad=1)
for g,gl in [("AIS/MIA","AIS/MIA (10)"),("invasive","Invasive (5)")]:
    xx=xs[XF.group.values==g]; ax.text(xx.mean(),-40,gl,ha="center",va="top",fontsize=6.5)
ax.set_ylim(0,100); ax.set_ylabel("Share of SPP1 transcripts (%)"); ax.set_title("Xenium: corrected SPP1 source per section",loc="left")
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
ax.legend(handles=[Patch(color=INV,label="Macrophage"),Patch(color=PRE,label="Epithelium"),Patch(color="#d8d7d2",label="Other"),Line2D([],[],marker="o",ls="",markerfacecolor="white",markeredgecolor=INK,markersize=4,label="Macrophage, uncorrected")],frameon=False,fontsize=6,loc="upper left",bbox_to_anchor=(1.0,1.0),ncol=1,handlelength=1)
ax.text(1.02,0.02,"F = FFPE",transform=ax.transAxes,fontsize=5.8,color=MUT)
lab(ax,"E",dx=-0.4)
save(fig,"Figure4")

# ---------------- Figure 5: hypoxia ----------------
H=pd.read_csv(P("hypoxia/H1_units.tsv", R),sep="\t"); H2=pd.read_csv(P("hypoxia/H2_sections.tsv", R),sep="\t")
fig,axs=plt.subplots(2,3,figsize=(W,4.6),gridspec_kw={"wspace":0.7,"hspace":0.85}); axs=axs.ravel()
def strip(ax,vals_by,title,ylab,note,letter):
    for j,(name,v,c) in enumerate(vals_by):
        v=np.asarray(v,float); v=v[~np.isnan(v)]; xs=j+np.random.default_rng(j).uniform(-0.13,0.13,len(v))
        ax.scatter(xs,v,s=12,color=c,edgecolor="white",lw=0.4,zorder=3); ax.plot([j-0.22,j+0.22],[np.median(v)]*2,color=INK,lw=1.2,zorder=4)
    ax.axhline(0,color=MUT,lw=0.6,ls=(0,(2,2))); ax.set_xticks(range(len(vals_by))); ax.set_xticklabels([v[0] for v in vals_by]); ax.set_xlim(-0.5,len(vals_by)-0.5)
    ax.set_title(title,loc="left"); ax.set_ylabel(ylab); ptxt(ax,note,y=-0.18); lab(ax,letter,dx=-0.25)
pm=H.groupby("patient")[["H1a_i_HIF","H1a_i_HALL","H1a_ii_HIF","H1a_ii_HALL"]].mean()
strip(axs[0],[("HIF13",pm.H1a_i_HIF,INV),("HALLMARK",pm.H1a_i_HALL,"#f2b393")],"SPP1+ vs SPP1− macrophages","Δ hypoxia score","18/19 patients; Holm P = 1.5×10⁻⁵","A")
strip(axs[1],[("HIF13",pm.H1a_ii_HIF,INV),("HALLMARK",pm.H1a_ii_HALL,"#f2b393")],"SPP1-free state score\n(top vs bottom quartile)","Δ hypoxia score","19/23 patients; P = 3×10⁻⁴","B")
w=H.pivot(index="patient",columns="grp",values="mean_HIF").dropna()
paired(axs[2],w.Precursor,w.LUAD,PRE,INV,["Precursor","Invasive"],"Hypoxia score (HIF13)","Macrophages, 22 pairs"); ptxt(axs[2],"17/22 higher; Holm P = 0.003",y=-0.18); lab(axs[2],"C",dx=-0.25)
H2e=H2[H2.eligible==True]; pp=H2e.groupby("patient")[["dH_HIF","dH_HALL"]].mean()
strip(axs[3],[("HIF13",pp.dH_HIF,"#8a8984"),("HALLMARK",pp.dH_HALL,"#c9c8c3")],"Spatial co-localization\n(Visium)","Δ Spearman ρ","Holm P = 0.055;\nsignature-dependent","D")
RU=pd.read_csv(P("mech/R_units.tsv", R),sep="\t"); rp=RU.groupby("patient")[["R1","R1_decirc"]].mean()
strip(axs[4],[("All targets",rp.R1,INV),("Hypoxia genes\nremoved",rp.R1_decirc,"#f2b393")],"HIF1A activity:\nSPP1+ vs SPP1− macrophages","Δ HIF1A activity (ULM t)","","E")
ptxt(axs[4],"15/19 and 16/19 patients;\nP = 4×10⁻⁴ and 6×10⁻⁴",y=-0.3)
w2=RU.pivot(index="patient",columns="grp",values="all_HIF1A").dropna()
paired(axs[5],w2.Precursor,w2.LUAD,PRE,INV,["Precursor","Invasive"],"Macrophage HIF1A activity","HIF1A activity,\n22 pairs"); ptxt(axs[5],"20/22 higher; Holm P = 3×10⁻⁵",y=-0.18); lab(axs[5],"F",dx=-0.25)
save(fig,"Figure5")

# ---------------- Figure 6: outcome and limits ----------------
pc=pd.read_csv(P("30_spatial/clinical/per_cohort.tsv", R),sep="\t"); pc=pc[pc.model=="main"]
A=json.load(open(P("sens_protein/A_niche_sensitivity.json", R))); Rr=json.load(open(P("sens_protein/R_prognosis_robustness.json", R))); T=json.load(open(P("sens_protein/T_hypoxia.json", R)))
mr=pd.read_csv(P("50_genetics/macrophage_SPP1/mr_coloc_results.tsv", R),sep="\t").iloc[0]
g=pd.read_csv(P("safety_druggability/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz", RAW),sep="\t",skiprows=2); sp1=g[g.Description=="SPP1"].iloc[0,2:].astype(float)
fig=plt.figure(figsize=(W,5.4)); gs=fig.add_gridspec(2,2,width_ratios=[1,1],height_ratios=[1.25,1],wspace=0.75,hspace=0.6)
ax=fig.add_subplot(gs[0,0]); main=Rr["main"]
rows=[(f"{r.cohort} ({r.endpoint}, n={r.n})",r.HR,r.lo,r.hi,INV) for r in pc.itertuples()]
for i,(n,h,l,u,c) in enumerate(rows):
    y=len(rows)+1-i; ax.plot([l,u],[y,y],color=c,lw=1.3); ax.scatter(h,y,s=16,color=c,marker="s",zorder=3)
y=0.7; ax.add_patch(plt.Polygon([[main["lo"],y],[main["HR"],y+0.38],[main["hi"],y],[main["HR"],y-0.38]],color=INK))
ax.plot(main["PI95"],[y-0.9,y-0.9],color=MUT,lw=1); ax.set_ylim(-0.6,len(rows)+1.6)
ax.set_xscale("log"); ax.axvline(1,color=MUT,lw=0.6,ls=(0,(2,2)))
ax.set_yticks([len(rows)+1-i for i in range(len(rows))]+[y,y-0.9]); ax.set_yticklabels([r[0] for r in rows]+[f"Random effects  {main['HR']:.2f} ({main['lo']:.2f}–{main['hi']:.2f})","95% prediction interval"])
ax.get_yticklabels()[-2].set_fontweight("bold"); ax.get_yticklabels()[-1].set_color(MUT)
ax.set_xticks([0.25,0.5,1,2,4]); ax.set_xticklabels(["0.25","0.5","1","2","4"]); ax.xaxis.set_minor_formatter(mpl.ticker.NullFormatter()); ax.set_xlabel("HR per SD of niche score (95% CI)")
ax.set_title("Stage I outcome (RFS, DFS or OS)",loc="left"); ptxt(ax,f"Two-sided P = {main['p_two_sided']:.3f}; I² = {main["I2"]:.2f}",y=-0.25); lab(ax,"A",dx=-1.05,dy=1.02)
ax=fig.add_subplot(gs[0,1]); loo=[v["HR"] for v in A["A1_leave_one_out"].values()]
mods=[("Main model",main),("+ purity proxy",Rr["R1_purity"]),("+ CAF score",A["A2_plus_CAF"]),("Myeloid-source genes only",A["A3_myeloid_source_genes_only"]),("+ hypoxia (HALLMARK)",T["T3_main_plus_HALLMARK"]),("+ hypoxia (HIF13)",Rr["R2_hypoxia"]),("Fully adjusted",Rr["R3_all"])]
for i,(n,v) in enumerate(mods):
    y=len(mods)-i; c=INK if v["lo"]>1 else NO; ax.plot([v["lo"],v["hi"]],[y,y],color=c,lw=1.3); ax.scatter(v["HR"],y,s=16,color=c,zorder=3)
ax.plot([min(loo),max(loo)],[0,0],color=INV,lw=4,alpha=0.5,solid_capstyle="round")
ax.set_yticks(list(range(len(mods),0,-1))+[0]); ax.set_yticklabels([m[0] for m in mods]+["Leave-one-cohort-out range"]); ax.set_xscale("log"); ax.axvline(1,color=MUT,lw=0.6,ls=(0,(2,2)))
ax.set_xlim(0.8,2.8); ax.set_xticks([0.8,1,1.5,2,2.5]); ax.set_xticklabels(["0.8","1","1.5","2","2.5"]); ax.xaxis.set_minor_formatter(mpl.ticker.NullFormatter()); ax.set_ylim(-0.7,len(mods)+0.7)
ax.set_xlabel("Pooled HR per SD (95% CI)"); ax.set_title("Robustness of the association",loc="left"); ptxt(ax,"grey: 95% CI includes 1",y=-0.25); lab(ax,"B",dx=-0.72,dy=1.02)
ax=fig.add_subplot(gs[1,0]); pps=[mr[f"PP.H{i}"] for i in range(5)]
ax.bar(range(5),pps,color=[NO,INV,NO,NO,"#6250d6"],width=0.65)
for i,p in enumerate(pps): ax.text(i,p+0.02,f"{p:.2f}",ha="center",fontsize=6,color=MUT)
ax.set_xticks(range(5)); ax.set_xticklabels(["H0","H1\neQTL only","H2","H3","H4\nshared"]); ax.set_ylim(0,1.05); ax.set_ylabel("Posterior probability")
ax.set_title("SPP1 macrophage eQTL × LUAD GWAS",loc="left"); ax.text(0.97,0.95,f"Wald-ratio MR\nOR {mr.wald_OR:.2f} ({mr.wald_OR_lo:.2f}–{mr.wald_OR_hi:.2f})",transform=ax.transAxes,ha="right",va="top",fontsize=6.3,color=MUT); lab(ax,"C",dx=-0.25)
ax=fig.add_subplot(gs[1,1]); top=sp1.sort_values(ascending=False); sel=list(top.index[:6])+["Lung","Whole Blood"]; sel=list(dict.fromkeys(sel))
vals=sp1[sel]; cols=[INV if "Kidney" in t else (PRE if t=="Lung" else NO) for t in sel]
short=lambda t: t.replace("Brain - Spinal cord (cervical c-1)","Spinal cord").replace(" - "," – ")
ax.barh(range(len(sel))[::-1],vals,color=cols,height=0.7); ax.set_yticks(range(len(sel))[::-1]); ax.set_yticklabels([short(t) for t in sel],fontsize=6.3); ax.set_xscale("log"); ax.xaxis.set_minor_formatter(mpl.ticker.NullFormatter()); ax.set_xlabel("SPP1 median TPM (GTEx v8)")
ax.set_title("SPP1 in normal tissues",loc="left"); lab(ax,"D",dx=-0.55)
save(fig,"Figure6")
print("ok")
