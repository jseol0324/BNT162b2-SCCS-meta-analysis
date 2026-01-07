## R/02_meta_regression.R
source("R/utils.R")

dat2 <- read_inputs()
ensure_dirs()

m0 <- fit_reml_kh(dat2$yi, dat2$vi, dat2$study)
cat("\n[Base model]\n")
print(summary(m0), digits=6)
print_re_summary(m0, "Base (REML+KH)")

fit_mod <- function(formula, name) {
  fit <- try(rma.uni(yi=dat2$yi, vi=dat2$vi, mods=formula,
                     method="REML", test="knha", data=dat2), silent=TRUE)
  if (inherits(fit,"try-error")) return(NULL)
  s <- summary(fit)
  cat("\n====================\n")
  cat(sprintf("[%s]\n", name))
  print(s, digits=6)
  cat(sprintf("tau^2=%.6f | pseudo-R^2=%.1f%% | QM=%.6f (p=%.4f)\n",
              as.numeric(fit$tau2), 100*pseudo_R2(m0, fit),
              as.numeric(fit$QM), as.numeric(fit$QMp)))
  fit
}

m_cont <- fit_mod(~ continent, "continent")
m_rw   <- fit_mod(~ risk_window, "risk_window")
m_win  <- fit_mod(~ window_days, "window_days")
m_dose <- fit_mod(~ dose, "dose (primary vs booster)")

# 다중(k>p 확인)
idx <- complete.cases(dat2$yi, dat2$vi, dat2$dose, dat2$risk_window, dat2$continent_bin)
X <- model.matrix(~ dose + risk_window + continent_bin, data=dat2[idx, ])
cat(sprintf("\n[check MULTI] k=%d, p=%d (need k>p)\n", nrow(X), ncol(X)))

m_multi <- try(
  rma.uni(yi=dat2$yi, vi=dat2$vi,
          mods=~ dose + risk_window + continent_bin,
          method="REML", test="knha",
          data=dat2, subset=idx),
  silent=TRUE
)
if (!inherits(m_multi,"try-error")) {
  cat("\n====================\n[MULTI: dose + risk_window + continent_bin]\n")
  print(summary(m_multi), digits=6)
  cat(sprintf("tau^2=%.6f | pseudo-R^2=%.1f%% | QM=%.6f (p=%.4f)\n",
              as.numeric(m_multi$tau2), 100*pseudo_R2(m0, m_multi),
              as.numeric(m_multi$QM), as.numeric(m_multi$QMp)))
} else {
  m_multi <- NULL
}

# 민감도: AbRahman 제외 dose 모형
m_dose_noAR <- try(
  rma.uni(yi=dat2$yi, vi=dat2$vi, mods=~ dose,
          method="REML", test="knha", data=dat2,
          subset = (study != "AbRahman 2024")),
  silent=TRUE
)
if (!inherits(m_dose_noAR,"try-error")) {
  cat("\n====================\n[dose exclude AbRahman 2024]\n")
  print(summary(m_dose_noAR), digits=6)
} else {
  m_dose_noAR <- NULL
}

# 모델 수준 요약표(Table S1 스타일)
model_row <- function(fit, name) {
  if (is.null(fit)) return(NULL)
  data.frame(
    model = name,
    k = fit$k,
    tau2 = as.numeric(fit$tau2),
    I2 = as.numeric(fit$I2),
    QM = as.numeric(fit$QM),
    QMp = as.numeric(fit$QMp),
    pseudoR2 = 100*pseudo_R2(m0, fit),
    stringsAsFactors = FALSE
  )
}

tbl_model <- do.call(rbind, Filter(Negate(is.null), list(
  model_row(m0, "null"),
  model_row(m_cont, "continent"),
  model_row(m_rw, "risk_window"),
  model_row(m_win, "window_days"),
  model_row(m_dose, "dose"),
  model_row(m_multi, "MULTI"),
  model_row(m_dose_noAR, "dose_noAR")
)))

safe_write_csv(tbl_model, "output/Table_S_model_summary.csv")
print(tbl_model, digits=6, row.names=FALSE)
