#!/bin/bash
# 第一步（靶点筛选）剩余全部数据：验证队列 + 正常组织/依赖性/药物数据库
# 用法：bash ~/全部数据库/GGN_target/scripts/download_rest.sh   （可重复运行，已完成的会跳过）
BASE="$HOME/全部数据库/GGN_target/raw"
source "$(dirname "$0")/geo_lib.sh"
get () {  # $1=URL  $2=保存目录  $3=文件名(可选)
  mkdir -p "$2" && cd "$2" || return
  local name=${3:-$(basename "${1%%\?*}")}
  if [ -s "$name" ] && [ ! -f "$name.part" ]; then echo "  已存在 $name"; return; fi
  touch "$name.part"
  if curl -sS -L --fail --connect-timeout 30 --retry 10 --retry-delay 10 --retry-all-errors -C - -o "$name" "$1" 2>/dev/null; then
    rm -f "$name.part"; echo "  完成 $name"
  else echo "  !! 失败 $name —— 重新运行脚本；仍失败请手动下载：$1"; fi
}

echo "== 1. 浏览器下载的两个单细胞/空间数据集（重新完整下载到项目文件夹）"
get_series GSE189357 scRNA_GSE189357
get_series GSE189487 ST_GSE189487

echo "== 2. I 期复发/生存验证队列（GEO 芯片，series matrix）"
for g in GSE31210 GSE50081 GSE72094 GSE37745 GSE30219 GSE68465; do
  get "$FTP/series/${g:0:${#g}-3}nnn/$g/matrix/${g}_series_matrix.txt.gz" "$BASE/stageI_cohorts"
done

echo "== 3. TCGA-LUAD（UCSC Xena GDC hub）"
for f in TCGA-LUAD.star_tpm.tsv.gz TCGA-LUAD.clinical.tsv.gz TCGA-LUAD.survival.tsv.gz TCGA-LUAD.somaticmutation_wxs.tsv.gz; do
  get "https://gdc-hub.s3.us-east-1.amazonaws.com/download/$f" "$BASE/TCGA_LUAD"
done

echo "== 4. DepMap（CRISPR 基因依赖性，最新版本）"
mkdir -p "$BASE/DepMap"
curl -sS -L --retry 5 "https://depmap.org/portal/api/download/files" -o "$BASE/DepMap/_index.csv" 2>/dev/null
for f in CRISPRGeneEffect.csv Model.csv OmicsExpressionProteinCodingGenesTPMLogp1.csv; do
  url=$(grep ",$f," "$BASE/DepMap/_index.csv" 2>/dev/null | head -1 | grep -o 'https://[^",]*' | head -1)
  if [ -n "$url" ]; then get "$url" "$BASE/DepMap" "$f"; else echo "  !! 未找到 $f 下载地址，请到 https://depmap.org/portal/data_page/ 手动下载"; fi
done

echo "== 5. 正常组织表达与药物-靶点数据库（第 5 层：长期用药安全性与可成药性）"
get "https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz" "$BASE/safety_druggability"
get "https://www.proteinatlas.org/download/tsv/normal_tissue.tsv.zip" "$BASE/safety_druggability"
get "https://dgidb.org/data/latest/interactions.tsv" "$BASE/safety_druggability" dgidb_interactions.tsv
echo "== 全部处理完毕；有 !! 标记的项目请重新运行或手动下载"
