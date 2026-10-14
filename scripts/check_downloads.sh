#!/bin/bash
# 对照 GEO 服务器上的文件大小，检查「下载」文件夹中的下载进度
D="$HOME/mnt/Downloads"
while read f s; do
  if [ -f "$D/$f" ]; then c=$(stat -c %s "$D/$f" 2>/dev/null || stat -f %z "$D/$f"); else c=0; fi
  pct=$(awk -v c=$c -v s=$s 'BEGIN{printf "%.1f", 100*c/s}')
  st=$([ "$c" = "$s" ] && echo "完成" || echo "下载中")
  printf "%-48s %8.1f / %8.1f MB  %5s%%  %s\n" "$f" $(awk -v c=$c 'BEGIN{print c/1e6}') $(awk -v s=$s 'BEGIN{print s/1e6}') "$pct" "$st"
done <<LIST
GSE102511_Smruthy_etal_allsamples_TPM.txt.gz 3857920
GSE166720_NIHAD_FPKM_normalized.txt.gz 7783875
GSE193725_GGO_pNS_Solid_RNAseq_FPKM.txt.gz 2175549
GSE277206_RAW.tar 53985280
GSE189487_RAW.tar 219719680
GSE189357_RAW.tar 654202880
GSE308103_RAW.tar 1561128960
GSE307534_RAW.tar 10090321920
LIST
