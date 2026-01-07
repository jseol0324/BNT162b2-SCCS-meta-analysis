## R/03_subgroup_analysis.R
source("R/utils.R")

dat2 <- read_inputs()
ensure_dirs()

m_null <- fit_reml_kh(dat2$yi, dat2$vi, dat2$study)

# 2수준 대비용 factor 생성
dat2$continent_AsiaBin <- factor(ifelse(dat2$continent=="Asia","Asia","Non-Asia"),
                                 levels=c("Asia","Non-Asia"))
dat2$risk_bin <- factor(dat2$risk_window, levels=c("<=21d",">21d"))
dat2$dose_bin <- factor(ifelse(dat2$dose=="primary","Primary","Booster"),
                        levels=c("Primary","Booster"))

subgroup_table_2lvl <- function(dat, var, label) {
  lv <- levels(dat[[var]])
  stopifnot(length(lv)==2)
  
  # 그룹별 적합
  rows <- lapply(lv, function(g) {
    idx <- dat[[var]] == g
    fit <- fit_reml_kh(dat$yi[idx], dat$vi[idx], dat$study[idx])
    data.frame(
      Subgroup = label,
      Group = g,
      k = fit$k,
      IRR = exp(as.numeric(fit$b)),
      LCL = exp(fit$ci.lb),
      UCL = exp(fit$ci.ub),
      I2 = as.numeric(fit$I2),
      tau2 = as.numeric(fit$tau2),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  
  # 대비(조절변수) 검정: QM, QMp (Q_between, p)
  mod_fit <- rma.uni(yi=dat$yi, vi=dat$vi, mods=~ dat[[var]],
                     method="REML", test="knha", data=dat)
  out$Q_between <- as.numeric(mod_fit$QM)
  out$p <- as.numeric(mod_fit$QMp)
  out
}

tab_cont <- subgroup_table_2lvl(dat2, "continent_AsiaBin", "Continent (Asia vs Non-Asia)")
tab_risk <- subgroup_table_2lvl(dat2, "risk_bin", "Risk window (≤21 d vs >21 d)")
tab_dose <- subgroup_table_2lvl(dat2, "dose_bin", "Dose (Primary vs Booster)")

tab_all <- rbind(tab_cont, tab_risk, tab_dose)

# 보기용 반올림
tab_all$IRR <- round(tab_all$IRR, 2)
tab_all$LCL <- round(tab_all$LCL, 2)
tab_all$UCL <- round(tab_all$UCL, 2)
tab_all$I2  <- round(tab_all$I2, 2)
tab_all$tau2 <- round(tab_all$tau2, 4)
tab_all$Q_between <- round(tab_all$Q_between, 3)
tab_all$p <- signif(tab_all$p, 3)

safe_write_csv(tab_all, "output/Table_subgroups_2level.csv")
print(tab_all, row.names=FALSE)
