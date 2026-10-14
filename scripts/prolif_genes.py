"""增殖基因集（2026-10-08 写定）：Tirosh 2016 的 S 期 + G2/M 期 97 个基因（基因名按新注释更新）；GO 细胞周期相关基因集（Enrichr GO_Biological_Process_2023，名称关键词写定）。"""
import os
TIROSH_S = "MCM5 PCNA TYMS FEN1 MCM2 MCM4 RRM1 UNG GINS2 MCM6 CDCA7 DTL PRIM1 UHRF1 CENPU HELLS RFC2 RPA2 NASP RAD51AP1 GMNN WDR76 SLBP CCNE2 UBR7 POLD3 MSH2 ATAD2 RAD51 RRM2 CDC45 CDC6 EXO1 TIPIN DSCC1 BLM CASP8AP2 USP1 CLSPN POLA1 CHAF1B BRIP1 E2F8".split()
TIROSH_G2M = "HMGB2 CDK1 NUSAP1 UBE2C BIRC5 TPX2 TOP2A NDC80 CKS2 NUF2 CKS1B MKI67 TMPO CENPK TACC3 PIMREG SMC4 CCNB2 CKAP2L CKAP2 AURKB BUB1 KIF11 ANP32E TUBB4B GTSE1 KIF20B HJURP CDCA3 JPT1 CDC20 TTK CDC25C KIF2C RANGAP1 NCAPD2 DLGAP5 CDCA2 CDCA8 ECT2 KIF23 HMMR AURKA PSRC1 ANLN LBR CKAP5 CENPE CTCF NEK2 GAS2L3 CBX5 CENPA".split()
TIROSH = TIROSH_S + TIROSH_G2M
GO_KEYWORDS = ["cell cycle", "mitotic", "cell division", "dna replication", "chromosome segregation", "spindle", "cytokinesis", "kinetochore"]

def go_cell_cycle_genes(path):
    genes = set()
    for line in open(path):
        parts = line.rstrip("\n").split("\t")
        name = parts[0].lower()
        if any(k in name for k in GO_KEYWORDS):
            genes.update(g for g in parts[2:] if g)
    return genes

def prolif_set():
    root = os.path.expanduser("~/全部数据库/GGN_target/raw/lincs")
    return set(TIROSH) | go_cell_cycle_genes(f"{root}/Enrichr_GO_Biological_Process_2023.txt")

def prolif_set_stage2():
    """第二阶段专用增殖基因集（2026-10-08 用户确认）：Tirosh 96 ∪ GO「Mitotic Nuclear Division (GO:0140014)」∪ GO「DNA Replication (GO:0006260)」，
    剔除属于 GO「Negative Regulation Of Cell Cycle (GO:0045786)」的基因（Enrichr GO_Biological_Process_2023）。与第一阶段的宽定义不同。"""
    root = os.path.expanduser("~/全部数据库/GGN_target/raw/lincs")
    terms = {}
    for line in open(f"{root}/Enrichr_GO_Biological_Process_2023.txt"):
        p = line.rstrip("\n").split("\t"); terms[p[0]] = {g for g in p[2:] if g}
    get = lambda name: set().union(*[v for k, v in terms.items() if k.lower().startswith(name.lower() + " (go:")])
    return (set(TIROSH) | get("Mitotic Nuclear Division") | get("DNA Replication")) - get("Negative Regulation Of Cell Cycle")
