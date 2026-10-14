#!/bin/bash
# Takano 2024 Xenium（SSBD ssbd-000458）：每张切片只下载表达矩阵和细胞坐标，共约 223 MB
B=https://ssbd.riken.jp/data/ssbd-000458/repos
D=~/全部数据库/GGN_target/raw/Xenium_Takano2024; mkdir -p $D
for s in output-XETG00042__0003418__LUAD_No2__20230222__093449 output-XETG00042__0003418__LUAD_No3__20230222__093449 output-XETG00042__0003419__TSU-24__20230222__093449 output-XETG00042__0003516__LUAD_No14__20230324__075904 output-XETG00042__0003516__LUAD_No16__20230324__075904 output-XETG00042__0003516__LUAD_No17__20230324__075904 output-XETG00130__0017208__TSU-23__20240305__025626 output-XETG00130__0017208__TSU-36__20240305__025626 output-XETG00130__0017208__TSU-41__20240305__025626 output-XETG00130__0018133__TSU-28__20240229__043945 output-XETG00130__0018133__TSU-30__20240229__043945 output-XETG00130__0018133__TSU-33__20240229__043945 output-XETG00130__0018153__TSU-25__20240229__043945 output-XETG00130__0018153__TSU-27__20240229__043945 output-XETG00130__0018153__TSU-35__20240229__043945; do
  mkdir -p $D/$s
  for f in cell_feature_matrix.h5 cells.csv.gz gene_panel.json metrics_summary.csv; do
    curl -s -L -C - --retry 5 -o $D/$s/$f $B/$s/$f
  done
  echo "done $s"
done
du -sh $D
