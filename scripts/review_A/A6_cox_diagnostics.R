# A6（方案 v1.3 §2-A6、§6.5）第 2 步：队列内 Cox 的 PH 诊断、收敛与自由度；小事件队列的简化模型与 Firth。
# 用法：Rscript scripts/review_A/A6_cox_diagnostics.R <A6_prepare 的输出目录>
# - 每个队列内模型：survival::coxph(ties = "breslow")，设计矩阵与 Python 完全相同；记录有效人数、事件数、参数个数、
#   是否收敛（警告即记录）、Schoenfeld 残差 PH 诊断 cox.zph（全局与逐变量）。
# - 事件数 < 20 的队列（主模型）：(i) 仅含评分的简化 Cox；(ii) Firth 惩罚 Cox（coxphf），公式执行前定死：
#   Surv(time, event) ~ niche_score + age + sex + stage + myeloid_marker_score
#   不收敛或报错即记为拟合失败，不更换实现、不调参重试、不删协变量重试。
# - 方案哈希在启动时校验。
suppressPackageStartupMessages({ library(survival); library(coxphf) })
args <- commandArgs(trailingOnly = TRUE); stopifnot(length(args) == 1)
RUN <- normalizePath(args[1]); P <- path.expand("~/全部数据库/GGN_target")
PLAN <- file.path(P, "docs/预设方案_外审回应A1-A6_v1.3_2026-10-19.md")
PLAN_SHA <- "1cde317f3bccac30050ddb8239ef6d623090f8aab8203765396548f9a136f94d"
sha <- tolower(trimws(strsplit(system2("shasum", c("-a", "256", shQuote(PLAN)), stdout = TRUE), " ")[[1]][1]))
if (!identical(sha, PLAN_SHA)) stop("方案哈希不符，拒绝运行：", sha)
cat("方案哈希校验通过", substr(sha, 1, 16), "\n")
OUT <- file.path(RUN, "R"); if (dir.exists(OUT)) stop("R 输出目录已存在，拒绝覆盖：", OUT); dir.create(OUT)
man <- read.delim(file.path(RUN, "designs", "manifest.tsv"), check.names = FALSE)

fit_one <- function(df, cols) {
  f <- as.formula(paste("Surv(time, event) ~", paste(sprintf("`%s`", cols), collapse = " + ")))
  warn <- character(0)
  fit <- withCallingHandlers(tryCatch(coxph(f, data = df, ties = "breslow"), error = function(e) e),
                             warning = function(w) { warn <<- c(warn, conditionMessage(w)); invokeRestart("muffleWarning") })
  list(fit = fit, warn = warn)
}
rows <- list(); zrows <- list()
for (i in seq_len(nrow(man))) {
  m <- man[i, ]; df <- read.delim(file.path(RUN, "designs", m$file), check.names = FALSE)
  cols <- setdiff(names(df), c("sample", "time", "event"))
  r <- fit_one(df, cols)
  rec <- list(model = m$model, cohort = m$cohort, n = nrow(df), events = sum(df$event), n_params = length(cols))
  if (inherits(r$fit, "error")) {
    rec$status <- paste("failed:", conditionMessage(r$fit))
  } else {
    fit <- r$fit; cf <- summary(fit)$coefficients
    rec$beta_R <- unname(coef(fit)[1]); rec$se_R <- cf[1, "se(coef)"]
    rec$n_used <- fit$n; rec$events_used <- fit$nevent; rec$iter <- fit$iter
    rec$converged <- length(r$warn) == 0 && all(is.finite(coef(fit)))
    rec$warnings <- paste(r$warn, collapse = " | ")
    z <- tryCatch(cox.zph(fit, transform = "km"), error = function(e) e)
    if (inherits(z, "error")) { rec$zph_status <- paste("failed:", conditionMessage(z)) } else {
      tb <- z$table; rec$zph_status <- "ok"; rec$zph_global_p <- tb["GLOBAL", "p"]
      rec$zph_score_p <- tb[1, "p"]
      for (k in rownames(tb)) zrows[[length(zrows) + 1]] <- data.frame(model = m$model, cohort = m$cohort, term = gsub("`", "", k),
                                                                       chisq = tb[k, "chisq"], df = tb[k, "df"], p = tb[k, "p"])
    }
    rec$status <- if (rec$converged) "ok" else "warning"
  }
  rows[[i]] <- as.data.frame(rec[!sapply(rec, is.null)], stringsAsFactors = FALSE)
}
all_cols <- unique(unlist(lapply(rows, names)))
rows <- lapply(rows, function(x) { for (k in setdiff(all_cols, names(x))) x[[k]] <- NA; x[all_cols] })
write.table(do.call(rbind, rows), file.path(OUT, "A6_R_cox.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(do.call(rbind, zrows), file.path(OUT, "A6_R_zph_terms.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

# ---- 小事件队列（主模型）：简化 Cox 与 Firth ----
small <- unique(man$cohort[man$model == "main" & man$events < 20])
frows <- list()
for (cn in small) {
  df <- read.delim(file.path(RUN, "designs", paste0("main__", cn, ".tsv")), check.names = FALSE)
  r <- fit_one(df, "score")
  rec <- list(cohort = cn, n = nrow(df), events = sum(df$event))
  if (inherits(r$fit, "error")) rec$simplified_status <- paste("failed:", conditionMessage(r$fit)) else {
    cf <- summary(r$fit)$coefficients; rec$simplified_beta <- cf[1, "coef"]; rec$simplified_se <- cf[1, "se(coef)"]
    rec$simplified_HR <- exp(cf[1, "coef"]); rec$simplified_p_two_sided <- cf[1, "Pr(>|z|)"]
    rec$simplified_status <- if (length(r$warn) == 0) "ok" else paste("warning:", paste(r$warn, collapse = " | "))
  }
  need <- c("score", "age", "male", "stage", "mye")
  if (!all(need %in% names(df))) {
    rec$firth_status <- paste("failed: covariate(s) absent (constant and dropped in the main design):", paste(setdiff(need, names(df)), collapse = ","))
  } else {
    d2 <- data.frame(time = df$time, event = df$event, niche_score = df$score, age = df$age, sex = df$male, stage = df$stage,
                     myeloid_marker_score = df$mye)
    warn <- character(0)
    ff <- withCallingHandlers(tryCatch(coxphf(Surv(time, event) ~ niche_score + age + sex + stage + myeloid_marker_score, data = d2),
                                       error = function(e) e),
                              warning = function(w) { warn <<- c(warn, conditionMessage(w)); invokeRestart("muffleWarning") })
    if (inherits(ff, "error")) rec$firth_status <- paste("failed:", conditionMessage(ff)) else {
      j <- which(names(ff$coefficients) == "niche_score")
      rec$firth_beta <- unname(ff$coefficients[j]); rec$firth_HR <- exp(rec$firth_beta)
      rec$firth_lo <- unname(ff$ci.lower[j]); rec$firth_hi <- unname(ff$ci.upper[j]); rec$firth_p_two_sided <- unname(ff$prob[j])
      rec$firth_iter <- if (!is.null(ff$iter)) ff$iter else NA
      conv_ok <- length(warn) == 0 && is.finite(rec$firth_beta)
      rec$firth_status <- if (conv_ok) "ok" else paste("failed (not converged or warning):", paste(warn, collapse = " | "))
    }
  }
  frows[[length(frows) + 1]] <- as.data.frame(rec, stringsAsFactors = FALSE)
}
if (length(frows)) {
  ac <- unique(unlist(lapply(frows, names)))
  frows <- lapply(frows, function(x) { for (k in setdiff(ac, names(x))) x[[k]] <- NA; x[ac] })
  write.table(do.call(rbind, frows), file.path(OUT, "A6_R_small_event.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
}
writeLines(c(R.version.string, paste("survival", as.character(packageVersion("survival"))),
             paste("coxphf", as.character(packageVersion("coxphf"))), paste("plan_sha256", sha)), file.path(OUT, "R_versions.txt"))
cat("A6 第 2 步完成：", OUT, "\n")
