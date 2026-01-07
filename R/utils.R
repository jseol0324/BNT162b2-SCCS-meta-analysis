## R/utils.R
suppressPackageStartupMessages({
  library(metafor)
})

ensure_dirs <- function() {
  if (!dir.exists("output")) dir.create("output", recursive = TRUE)
}

read_inputs <- function(dat_path = "data/dat.csv", meta_path = "data/meta_info.csv") {
  dat <- read.csv(dat_path, stringsAsFactors = FALSE)
  meta <- read.csv(meta_path, stringsAsFactors = FALSE)
  
  # 최소 필드 체크
  req_dat <- c("study","IRR","SE")
  if (!all(req_dat %in% names(dat))) stop("dat.csv must contain: ", paste(req_dat, collapse=", "))
  
  # 효과크기
  dat$yi <- log(dat$IRR)
  dat$vi <- dat$SE^2
  
  # 병합
  dat2 <- merge(dat, meta, by="study", all.x=TRUE, sort=FALSE)
  
  # dose 정규화('mixed->booster')
  dose_chr <- tolower(as.character(dat2$dose))
  dose_chr[dose_chr %in% c("mixed","booster_any","3","3rd","third")] <- "booster"
  dat2$dose <- factor(dose_chr, levels=c("primary","booster"))
  
  dat2$continent <- factor(dat2$continent, levels=c("Europe","Asia","North America"))
  dat2$risk_window <- factor(dat2$risk_window, levels=c("<=21d",">21d"))
  dat2$continent_bin <- factor(ifelse(dat2$continent == "Europe","Europe","NonEurope"),
                               levels=c("Europe","NonEurope"))
  dat2
}

fit_reml_kh <- function(yi, vi, slab) {
  rma.uni(yi=yi, vi=vi, method="REML", test="knha", slab=slab)
}

pseudo_R2 <- function(m0, m1) {
  if (is.null(m0) || is.null(m1)) return(NA_real_)
  if (is.na(m0$tau2) || is.na(m1$tau2) || m0$tau2 <= 0) return(NA_real_)
  max(0, (m0$tau2 - m1$tau2) / m0$tau2)
}

fixed_Q_I2 <- function(yi, vi) {
  w <- 1/vi
  theta <- sum(w*yi)/sum(w)
  Q <- sum(w*(yi-theta)^2)
  I2 <- if (Q > 0) max(0, (Q - (length(yi)-1))/Q)*100 else 0
  list(Q=Q, df=length(yi)-1, I2=I2, theta=theta, w=w)
}

safe_write_csv <- function(x, path) {
  ensure_dirs()
  write.csv(x, path, row.names = FALSE)
}

save_baujat_png <- function(m, file = "output/Fig_S_Baujat.png", res = 300) {
  ensure_dirs()
  png(file, width=1600, height=1200, res=res)
  baujat(m, main="Baujat plot — contribution to heterogeneity vs influence")
  dev.off()
}

print_re_summary <- function(m, label="REML+KH") {
  pooled <- exp(as.numeric(m$b))
  lcl <- exp(m$ci.lb)
  ucl <- exp(m$ci.ub)
  cat(sprintf("\n[%s] Pooled IRR = %.3f [%.3f, %.3f] | I^2 = %.1f%% | tau^2 = %.6f\n",
              label, pooled, lcl, ucl, as.numeric(m$I2), as.numeric(m$tau2)))
}
