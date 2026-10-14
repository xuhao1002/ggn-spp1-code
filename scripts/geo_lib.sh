#!/bin/bash
# GEO 下载公共函数：限并发、失败自动重试、按服务器大小核对，只有大小一致才算完成
FTP="https://ftp.ncbi.nlm.nih.gov/geo"
NJOBS=${NJOBS:-3}        # NCBI 建议并发不超过 3
ROUNDS=${ROUNDS:-8}      # 每轮结束后对未完成文件再试，最多 8 轮

fetch () {   # 带重试地下载一个 URL 到当前目录（支持续传）
  curl -sS --fail --connect-timeout 30 --retry 10 --retry-delay 10 --retry-all-errors -C - -O "$1" 2>/dev/null
}

sample_url () {  # 由样本文件名推出 GSM 文件地址
  local f=$1 g=${1%%_*}
  echo "$FTP/samples/${g:0:${#g}-3}nnn/$g/suppl/$f"
}
export -f fetch sample_url; export FTP

pending_list () {  # 列出大小不符的文件
  while IFS=$'\t' read f s; do
    c=$(stat -f %z "$f" 2>/dev/null || echo 0)
    # 超过应有大小说明文件损坏（例如混入了报错页），删掉重下
    if [ "$c" -gt "$s" ]; then rm -f "$f"; c=0; fi
    [ "$c" = "$s" ] || echo "$f"
  done < files_expected.tsv
}

get_series () {   # $1=GSE  $2=保存目录（相对 BASE）
  local gse=$1 dir="$BASE/$2" sub="${1:0:${#1}-3}nnn"
  mkdir -p "$dir" && cd "$dir" || return
  echo "== $gse -> $dir"
  # 1) 系列级补充文件（RAW.tar 除外）
  local listing; listing=$(curl -sS --fail --retry 10 --retry-delay 10 --retry-all-errors "$FTP/series/$sub/$gse/suppl/" 2>/dev/null)
  if [ -z "$listing" ]; then echo "  !! 连不上 NCBI，请稍后重新运行脚本"; return; fi
  echo "$listing" | grep -o 'href="[^"]*"' | cut -d'"' -f2 \
    | grep -v -E '^\.\.|/$|_RAW\.tar$|^filelist\.txt$|^https?:|^/' | while read f; do
      fetch "$FTP/series/$sub/$gse/suppl/$f" && echo "  完成 $f" || echo "  !! 失败 $f（重新运行脚本即可）"; done
  # 2) 样本级文件
  if ! echo "$listing" | grep -q 'filelist.txt'; then echo "== $gse 无样本级文件"; return; fi
  fetch "$FTP/series/$sub/$gse/suppl/filelist.txt" || { echo "  !! 取不到文件清单，请稍后重新运行"; return; }
  awk -F'\t' '$1=="File"{print $2"\t"$4}' filelist.txt > files_expected.tsv
  local total; total=$(wc -l < files_expected.tsv | tr -d ' ')
  for r in $(seq 1 $ROUNDS); do
    local todo; todo=$(pending_list)
    [ -z "$todo" ] && break
    echo "  第 $r 轮：还有 $(echo "$todo" | wc -l | tr -d ' ') / $total 个文件待下载"
    echo "$todo" | xargs -P $NJOBS -I{} bash -c 'fetch "$(sample_url {})" && echo "  完成 {}"'
    sleep 5
  done
  local left; left=$(pending_list)
  if [ -z "$left" ]; then echo "== $gse 全部 $total 个文件下载完整"
  else echo "== $gse 仍有 $(echo "$left" | wc -l | tr -d ' ') 个文件未完成，请重新运行脚本续传"; fi
}
