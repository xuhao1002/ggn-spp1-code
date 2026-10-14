# D：RCTD 逐切片（docs/预设方案_调控互作RCTD_2026-10-13.md）。可续跑。
suppressMessages({library(spacexr); library(Matrix)})
P <- "/Users/dudu/全部数据库/GGN_target/"; O <- paste0(P,"results/mech/rctd/")
g <- readLines(paste0(O,"ref/genes.tsv")); cl <- read.delim(paste0(O,"ref/cells.tsv"))
M <- as(readMM(paste0(O,"ref/counts.mtx")),"CsparseMatrix"); rownames(M) <- make.unique(g); colnames(M) <- cl$cell
ct <- factor(setNames(cl$type,cl$cell)); nU <- setNames(colSums(M),cl$cell)
ref <- Reference(M, ct, nU)
cat("reference built:", nrow(M), "genes,", ncol(M), "cells\n")
dirs <- sort(list.dirs(paste0(P,"raw/ST_GSE307534/extracted"),recursive=FALSE))
for (d in dirs) {
  b <- basename(d); f <- paste0(O,b,"_rctd.tsv.gz")
  if (file.exists(f)) { cat("skip",b,"\n"); next }
  t0 <- Sys.time()
  sub <- list.dirs(d,recursive=FALSE)[1]; mx <- paste0(sub,"/filtered_feature_bc_matrix/")
  C <- as(readMM(gzfile(paste0(mx,"matrix.mtx.gz"))),"CsparseMatrix")
  fe <- read.delim(gzfile(paste0(mx,"features.tsv.gz")),header=FALSE); bc <- readLines(gzfile(paste0(mx,"barcodes.tsv.gz")))
  rownames(C) <- make.unique(fe$V2); colnames(C) <- bc
  pos <- read.csv(paste0(sub,"/spatial/tissue_positions.csv")); pos <- pos[pos$in_tissue==1 & pos$barcode %in% bc,]
  C <- C[!grepl("^MT-",rownames(C),ignore.case=TRUE), pos$barcode]
  coords <- data.frame(x=pos$pxl_col_in_fullres,y=pos$pxl_row_in_fullres,row.names=pos$barcode)
  sp_ <- SpatialRNA(coords, C, colSums(C))
  r <- create.RCTD(sp_, ref, max_cores=6, CELL_MIN_INSTANCE=25)
  r <- run.RCTD(r, doublet_mode="full")
  w <- as.matrix(normalize_weights(r@results$weights))
  out <- gzfile(f,"w"); write.table(data.frame(barcode=rownames(w),w,check.names=FALSE),out,sep="\t",quote=FALSE,row.names=FALSE); close(out)
  cat(format(Sys.time()), b, "spots", nrow(w), "min", round(as.numeric(difftime(Sys.time(),t0,units="mins")),1), "\n")
}
cat("DONE\n")
