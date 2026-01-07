## R/01_meta_analysis.R
source("R/utils.R")

dat2 <- read_inputs()
ensure_dirs()

# 1) 주모형: REML + KH
m <- fit_reml_kh(dat2$yi, dat2$vi, dat2$study)
print(summary(m), digits=6)
print_re_summary(m, "Random-effects (REML + KH)")

# 2) 이질성: Q-based vs tau^2-based(I^2는 metafor가 tau^2 기반)
fe <- fixed_Q_I2(dat2$yi, dat2$vi)
cat(sprintf("\n[Fixed-effects] Q = %.6f, df = %d, I^2(Q-based)=%.1f%%\n",
            fe$Q, fe$df, fe$I2))

# 3) LOO 이질성(Leave-One-out)
loo <- metafor::leave1out(m)
loo_tab <- data.frame(
  study = m$slab,
  tau2_loo = as.numeric(loo$tau2),
  I2_loo   = as.numeric(loo$I2)
)
loo_tab$delta_tau2 <- loo_tab$tau2_loo - as.numeric(m$tau2)
loo_tab$delta_I2   <- loo_tab$I2_loo   - as.numeric(m$I2)
safe_write_csv(loo_tab, "output/Table_S_LOO_heterogeneity.csv")
print(loo_tab, digits=6, row.names=FALSE)

# 4) Baujat plot 저장
save_baujat_png(m, "output/Fig_S_Baujat.png", res=300)

# 5) 영향진단(표로 저장)
inf <- influence(m)
rs <- rstudent(m)

diag_tab <- data.frame(
  study = dat2$study,
  std_resid = if (!is.null(rs$z)) as.numeric(rs$z) else NA_real_,
  cook_d    = if (!is.null(inf$cook.d)) as.numeric(inf$cook.d) else NA_real_,
  hat       = if (!is.null(inf$hat)) as.numeric(inf$hat) else NA_real_
)
safe_write_csv(diag_tab, "output/Table_S_influence_diagnostics.csv")
print(diag_tab, digits=6, row.names=FALSE)

# 6) “고정효과 기준 Q 기여도” (Baujat x-axis)
Q_i <- fe$w * (dat2$yi - fe$theta)^2
Q_pct <- 100 * Q_i / sum(Q_i)
Q_tab <- data.frame(study=dat2$study, Q_i_FE=Q_i, Q_i_FE_pct=Q_pct)
Q_tab <- Q_tab[order(Q_tab$Q_i_FE, decreasing=TRUE), ]
safe_write_csv(Q_tab, "output/Table_S_Q_contribution_FE.csv")
print(Q_tab, digits=6, row.names=FALSE)
