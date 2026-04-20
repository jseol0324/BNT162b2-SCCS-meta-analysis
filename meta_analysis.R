## ============================================================
## Meta-analysis: BNT162b2 & Stroke Risk
## 8 SCCS studies | REML + Knapp-Hartung adjustment
##
## Citation:
##   Jang S, Pyo S, Kim EJ, Lee W, Kwon S-S.
##   Stroke risk following BNT162b2 vaccination: A systematic review
##   and meta-analysis of self-controlled case series studies.
##   Scientific Reports (submitted 2026).
##
## Repository: https://github.com/jseol0324/BNT162b2-SCCS-meta-analysis
## R version:  4.4.1
## metafor:    4.8.0
## Last updated: 2026-03-16
## ============================================================
## Included studies (chronological):
##   1. Hippisley-Cox et al. 2021 (UK)
##   2. Chui et al. 2022          (Hong Kong)
##   3. Jabagi et al. 2022        (France, ≥75y)
##   4. Botton et al. 2022        (France, 18-74y)
##   5. Ab Rahman et al. 2024     (Malaysia)
##   6. Xu et al. 2024            (USA)
##   7. Makkar et al. 2025        (USA)
##   8. Salmaggi et al. 2025      (Italy)
## ============================================================


## ── 0) Packages ────────────────────────────────────────────
for (pkg in c("metafor", "ggplot2", "ggrepel")) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, repos = "https://cloud.r-project.org")
}
library(metafor)
library(ggplot2)
library(ggrepel)


## ============================================================
## SECTION A: Within-study inverse-variance pooling
##
## Applied to studies that did not report a single combined IRR
## but provided stratum-specific estimates (by dose / outcome
## subtype). Pooling is performed on the log scale.
## ============================================================

iv_pool <- function(IRR, LCL, UCL) {
  yi_vec <- log(IRR)
  se_vec <- (log(UCL) - log(LCL)) / (2 * 1.96)
  vi_vec <- se_vec^2
  w      <- 1 / vi_vec
  yi_p   <- sum(w * yi_vec) / sum(w)
  se_p   <- sqrt(1 / sum(w))
  list(IRR = exp(yi_p),
       LCL = exp(yi_p - 1.96 * se_p),
       UCL = exp(yi_p + 1.96 * se_p),
       yi  = yi_p,
       sei = se_p)
}

## A-1) Jabagi 2022 ------------------------------------------
## Source: JAMA 2022, Table 2 — BNT162b2, ≥75y
## IS Dose1 1-14d: 0.90 (0.84-0.98) | IS Dose2 1-14d: 0.92 (0.84-1.02)
## HS Dose1 1-14d: 0.90 (0.78-1.04) | HS Dose2 1-14d: 0.97 (0.81-1.15)
jabagi_pool <- iv_pool(
  IRR = c(0.90, 0.92, 0.90, 0.97),
  LCL = c(0.84, 0.84, 0.78, 0.81),
  UCL = c(0.98, 1.02, 1.04, 1.15)
)
cat(sprintf("[Jabagi 2022]   IRR=%.4f [%.4f, %.4f]  log(IRR)=%.5f  SE=%.5f\n",
            jabagi_pool$IRR, jabagi_pool$LCL, jabagi_pool$UCL,
            jabagi_pool$yi, jabagi_pool$sei))

## A-2) Botton 2022 ------------------------------------------
## Source: Ann Intern Med 2022, Table 3 — Pfizer-BioNTech, 18-74y
## IS: Dose1(Wk1=0.84,Wk2=0.95), Dose2(Wk1=0.93,Wk2=1.09)
## HS: Dose1(Wk1=0.97,Wk2=1.07), Dose2(Wk1=0.98,Wk2=0.86)
botton_pool <- iv_pool(
  IRR = c(0.84, 0.95, 0.93, 1.09, 0.97, 1.07, 0.98, 0.86),
  LCL = c(0.74, 0.85, 0.81, 0.96, 0.80, 0.88, 0.77, 0.67),
  UCL = c(0.94, 1.06, 1.06, 1.23, 1.19, 1.30, 1.25, 1.11)
)
cat(sprintf("[Botton 2022]   IRR=%.4f [%.4f, %.4f]  log(IRR)=%.5f  SE=%.5f\n",
            botton_pool$IRR, botton_pool$LCL, botton_pool$UCL,
            botton_pool$yi, botton_pool$sei))

## A-3) Makkar 2025 ------------------------------------------
## Source: J Neurol Sci 2025, Table 3 — Robust SCCS, BNT162b2
## Dose1-4 x 1-14d/15-28d; age ≥16
makkar_pool <- iv_pool(
  IRR = c(0.723, 0.870, 0.631, 0.814, 0.964, 0.833, 0.989, 1.136),
  LCL = c(0.360, 0.466, 0.237, 0.402, 0.339, 0.328, 0.001, 0.002),
  UCL = c(1.141, 1.358, 1.050, 1.297, 1.675, 1.582, 2.512, 2.898)
)
cat(sprintf("[Makkar 2025]   IRR=%.4f [%.4f, %.4f]  log(IRR)=%.5f  SE=%.5f\n",
            makkar_pool$IRR, makkar_pool$LCL, makkar_pool$UCL,
            makkar_pool$yi, makkar_pool$sei))

## A-4) Salmaggi 2025 ----------------------------------------
## Source: Neurol Sci 2025, Table 2 — BNT162b2, Italy (Lombardia), ≥12y
## IS: 0.98 (0.91-1.06) | HS: 0.95 (0.83-1.08)
## Note: reference period = pre-vaccination baseline (differs from other studies)
salmaggi_pool <- iv_pool(
  IRR = c(0.98, 0.95),
  LCL = c(0.91, 0.83),
  UCL = c(1.06, 1.08)
)
cat(sprintf("[Salmaggi 2025] IRR=%.4f [%.4f, %.4f]  log(IRR)=%.5f  SE=%.5f\n\n",
            salmaggi_pool$IRR, salmaggi_pool$LCL, salmaggi_pool$UCL,
            salmaggi_pool$yi, salmaggi_pool$sei))


## ============================================================
## SECTION B: Dataset
## ============================================================

dat2 <- data.frame(
  study = c("Hippisley-Cox 2021",
            "Chui 2022",      "Jabagi 2022",   "Botton 2022",
            "AbRahman 2024",  "Xu 2024",
            "Makkar 2025",    "Salmaggi 2025"),

  IRR   = c(1.04162890,
            1.67447740, jabagi_pool$IRR,  botton_pool$IRR,
            0.93130640, 0.97000000,
            makkar_pool$IRR, salmaggi_pool$IRR),

  SE    = c(0.01963612,
            0.18003846, jabagi_pool$sei, botton_pool$sei,
            0.02841848, 0.08940900,
            makkar_pool$sei, salmaggi_pool$sei),

  events = c(6439, 117, 17014, 28717, 74700, 1057, 119275, NA),

  continent = c("Europe", "Asia", "Europe", "Europe",
                "Asia", "North America", "North America", "Europe"),

  window_days = c(28, 27, 14, 21, 21, 42, 28, 28),

  risk_window = c(">21d", ">21d", "<=21d", "<=21d",
                  "<=21d", ">21d", ">21d", ">21d"),

  ## primary = doses 1-2 only; booster = includes dose >=3
  dose = c("primary", "primary", "primary", "primary",
           "booster", "booster", "booster", "booster"),

  sccs_method = c("Standard", "Modified", "Modified", "Modified",
                  "Standard", "Modified", "Robust",   "Modified"),

  sccs_bin = c("Standard",     "Non-standard", "Non-standard", "Non-standard",
               "Standard",     "Non-standard", "Non-standard", "Non-standard"),

  ## Both = IS+HS combined, IS = ischemic only, HS = hemorrhagic only
  outcome_type = c("Both", "HS", "Both", "Both",
                   "Both", "IS", "IS",   "Both"),

  year = c(2021, 2022, 2022, 2022, 2024, 2024, 2025, 2025),

  ## Post-vax = post-vaccination control; Mixed = pre+post; Pre-vax = pre-vaccination only
  baseline_type = c("Post-vax", "Post-vax", "Mixed",    "Mixed",
                    "Post-vax", "Post-vax", "Post-vax", "Pre-vax"),

  min_age = c(16, 16, 75, 18, 18, 12, 16, 12),

  stringsAsFactors = FALSE
)

## Derived variables
dat2$yi <- log(dat2$IRR)
dat2$vi <- dat2$SE^2

## Age group variables
## Adult_only (≥18):           Jabagi(75), Botton(18), AbRahman(18)  k=3
## Includes_adolescents (<18): HC(16), Chui(16), Xu(12), Makkar(16), Salmaggi(12) k=5
dat2$age_group <- factor(
  ifelse(dat2$min_age < 18, "Includes_adolescents", "Adult_only"),
  levels = c("Adult_only", "Includes_adolescents"))

dat2$elderly_only <- factor(
  ifelse(dat2$min_age >= 75, "Elderly_only", "Non-elderly"),
  levels = c("Non-elderly", "Elderly_only"))

## Factor conversions
dat2$continent_bin <- factor(
  ifelse(dat2$continent == "Asia", "Asia", "Non-Asia"),
  levels = c("Asia", "Non-Asia"))
dat2$risk_window   <- factor(dat2$risk_window,  levels = c("<=21d", ">21d"))
dat2$dose          <- factor(dat2$dose,         levels = c("primary", "booster"))
dat2$sccs_method   <- factor(dat2$sccs_method,
                             levels = c("Standard", "Partial", "Modified", "Robust"))
dat2$sccs_bin      <- factor(dat2$sccs_bin,     levels = c("Standard", "Non-standard"))
dat2$outcome_type  <- factor(dat2$outcome_type, levels = c("Both", "IS", "HS"))
dat2$baseline_type <- factor(dat2$baseline_type,levels = c("Post-vax", "Mixed", "Pre-vax"))

cat("=== Cross-tab checks ===\n")
print(table(dat2$continent_bin,  useNA = "ifany"))
print(table(dat2$risk_window,    useNA = "ifany"))
print(table(dat2$dose,           useNA = "ifany"))
print(table(dat2$sccs_method,    useNA = "ifany"))
print(table(dat2$sccs_bin,       useNA = "ifany"))
print(table(dat2$outcome_type,   useNA = "ifany"))
print(table(dat2$baseline_type,  useNA = "ifany"))
print(table(dat2$age_group,      useNA = "ifany"))
print(table(dat2$elderly_only,   useNA = "ifany"))


## ============================================================
## SECTION C: Primary meta-analysis
## ============================================================

m <- rma.uni(yi = dat2$yi, vi = dat2$vi,
             method = "REML", test = "knha", slab = dat2$study)
cat("\n=== Primary meta-analysis (k=8, REML+KH) ===\n")
print(summary(m))

tau2_show <- as.numeric(m$tau2)
I2_show   <- as.numeric(m$I2)
cat(sprintf("\n[Pooled] IRR=%.3f [%.3f, %.3f] | I2=%.1f%% | tau2=%.4f\n",
            exp(as.numeric(m$b)), exp(m$ci.lb), exp(m$ci.ub),
            I2_show, tau2_show))

## Study-level + pooled summary table
LCL_i <- exp(dat2$yi - 1.96 * dat2$SE)
UCL_i <- exp(dat2$yi + 1.96 * dat2$SE)
wi_   <- 1 / (dat2$vi + tau2_show)
wt_   <- round(100 * wi_ / sum(wi_), 2)
wt_[which.max(wt_)] <- wt_[which.max(wt_)] + (100 - sum(wt_))

tab_main <- data.frame(
  Study  = c(as.character(dat2$study), "Pooled (REML+KH)"),
  Events = c(dat2$events, sum(dat2$events, na.rm = TRUE)),
  IRR    = round(c(dat2$IRR,        exp(as.numeric(m$b))), 3),
  LCL    = round(c(LCL_i,           exp(m$ci.lb)),         3),
  UCL    = round(c(UCL_i,           exp(m$ci.ub)),         3),
  Weight = c(paste0(wt_, "%"), "—"),
  stringsAsFactors = FALSE
)
cat("\n=== Study-level + pooled results ===\n")
print(tab_main, row.names = FALSE)


## ============================================================
## SECTION D: Forest plot
## ============================================================

metafor::forest(m,
                atransf  = exp,
                at       = log(c(0.5, 0.75, 1.0, 1.5, 2.0, 2.5)),
                xlim     = c(-5.0, 3.0),
                xlab     = "Incidence Rate Ratio (IRR)",
                header   = c("Study", "IRR [95% CI]"),
                ilab     = cbind(formatC(wt_, format = "f", digits = 2)),
                ilab.xpos = c(-3.3),
                ilab.lab  = c("Weight (%)"),
                refline  = 0,
                mlab     = "Random-Effects Model (REML + KH)",
                cex      = 0.9)
title("Forest plot — BNT162b2 & Stroke")


## ============================================================
## SECTION E: Leave-one-out analysis
## ============================================================

loo <- leave1out(m)
loo_tab <- data.frame(
  study      = m$slab,
  IRR_loo    = round(exp(as.numeric(loo$estimate)), 3),
  LCL_loo    = round(exp(as.numeric(loo$ci.lb)),    3),
  UCL_loo    = round(exp(as.numeric(loo$ci.ub)),    3),
  I2_loo     = round(as.numeric(loo$I2),   1),
  tau2_loo   = round(as.numeric(loo$tau2), 4),
  delta_I2   = round(as.numeric(loo$I2)   - I2_show,   1),
  delta_tau2 = round(as.numeric(loo$tau2) - tau2_show, 4)
)
cat("\n=== Leave-one-out analysis ===\n")
print(loo_tab, row.names = FALSE)


## ============================================================
## SECTION F: Subgroup analyses
## ============================================================

## F-0a) Safe rma wrapper (fallback estimation methods)
safe_rma <- function(yi, vi, slab) {
  for (meth in c("REML", "ML", "DL")) {
    kh  <- (meth != "DL")
    fit <- tryCatch(
      rma.uni(yi = yi, vi = vi, method = meth,
              test = if (kh) "knha" else "z",
              slab = slab),
      error   = function(e) NULL,
      warning = function(w) NULL
    )
    if (!is.null(fit)) { attr(fit, "method_used") <- meth; return(fit) }
  }
  stop("All estimation methods failed (k=", length(yi), ")")
}

## F-0b) Two-level subgroup helper
subgroup_2lvl <- function(data, var, label) {
  lv  <- levels(data[[var]])
  mod <- tryCatch(
    rma.uni(yi = data$yi, vi = data$vi,
            mods = ~data[[var]], method = "REML", test = "knha"),
    error = function(e) NULL, warning = function(w) NULL)
  QM  <- if (!is.null(mod)) round(as.numeric(mod$QM),  3) else NA_real_
  QMp <- if (!is.null(mod)) signif(as.numeric(mod$QMp), 3) else NA_real_

  rows <- lapply(lv, function(g) {
    idx <- which(data[[var]] == g)
    k_g <- length(idx)
    if (k_g < 2) {
      return(data.frame(
        Subgroup = label, Group = g, k = k_g,
        IRR  = round(data$IRR[idx], 3),
        LCL  = round(exp(data$yi[idx] - 1.96 * data$SE[idx]), 3),
        UCL  = round(exp(data$yi[idx] + 1.96 * data$SE[idx]), 3),
        I2 = NA_real_, tau2 = NA_real_,
        Qbetween = QM, p = QMp,
        note = "k=1; single study", stringsAsFactors = FALSE))
    }
    fit  <- safe_rma(data$yi[idx], data$vi[idx], data$study[idx])
    used <- attr(fit, "method_used")
    data.frame(
      Subgroup = label, Group = g, k = fit$k,
      IRR  = round(exp(as.numeric(fit$b)), 3),
      LCL  = round(exp(fit$ci.lb),         3),
      UCL  = round(exp(fit$ci.ub),         3),
      I2   = round(as.numeric(fit$I2),     1),
      tau2 = round(as.numeric(fit$tau2),   4),
      Qbetween = QM, p = QMp,
      note = if (used != "REML") paste0("fallback: ", used) else "",
      stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

## F-1) Pre-specified subgroups
tab_continent <- subgroup_2lvl(dat2, "continent_bin", "Continent (Asia vs Non-Asia)")
tab_risk      <- subgroup_2lvl(dat2, "risk_window",   "Risk window (<=21d vs >21d)")
tab_dose      <- subgroup_2lvl(dat2, "dose",          "Dose (Primary vs Booster)")
tab_sccs      <- subgroup_2lvl(dat2, "sccs_bin",      "SCCS method (Standard vs Non-standard)")
tab_age       <- subgroup_2lvl(dat2, "age_group",     "Age group (Adult-only [>=18] vs Includes adolescents [<18])")
tab_elderly   <- subgroup_2lvl(dat2, "elderly_only",  "Elderly-only (>=75 [Jabagi] vs Non-elderly)")

## F-2) Outcome subtype (descriptive — groups not mutually exclusive)
idx_IS <- which(dat2$outcome_type %in% c("Both", "IS"))
idx_HS <- which(dat2$outcome_type %in% c("Both", "HS"))
fit_IS <- safe_rma(dat2$yi[idx_IS], dat2$vi[idx_IS], dat2$study[idx_IS])
fit_HS <- safe_rma(dat2$yi[idx_HS], dat2$vi[idx_HS], dat2$study[idx_HS])
tab_outcome <- data.frame(
  Subgroup = rep("Outcome subtype (descriptive only)", 2),
  Group    = c(paste0("Includes IS (k=", fit_IS$k, ")"),
               paste0("Includes HS (k=", fit_HS$k, ")")),
  k        = c(fit_IS$k,  fit_HS$k),
  IRR      = round(c(exp(as.numeric(fit_IS$b)), exp(as.numeric(fit_HS$b))), 3),
  LCL      = round(c(exp(fit_IS$ci.lb), exp(fit_HS$ci.lb)), 3),
  UCL      = round(c(exp(fit_IS$ci.ub), exp(fit_HS$ci.ub)), 3),
  I2       = round(c(as.numeric(fit_IS$I2),  as.numeric(fit_HS$I2)),  1),
  tau2     = round(c(as.numeric(fit_IS$tau2), as.numeric(fit_HS$tau2)), 4),
  Qbetween = c(NA, NA), p = c(NA, NA),
  note     = rep("Non-mutually exclusive; no between-group comparison", 2),
  stringsAsFactors = FALSE)

## F-3) Baseline period type (descriptive)
baseline_rows <- lapply(levels(dat2$baseline_type), function(g) {
  idx <- which(dat2$baseline_type == g)
  k_g <- length(idx)
  if (k_g < 2) {
    return(data.frame(
      Subgroup = "Baseline type (descriptive only)", Group = g, k = k_g,
      IRR  = round(dat2$IRR[idx], 3),
      LCL  = round(LCL_i[idx],   3),
      UCL  = round(UCL_i[idx],   3),
      I2 = NA_real_, tau2 = NA_real_,
      Qbetween = NA_real_, p = NA_real_,
      note = "k=1; single study", stringsAsFactors = FALSE))
  }
  fit  <- safe_rma(dat2$yi[idx], dat2$vi[idx], dat2$study[idx])
  used <- attr(fit, "method_used")
  data.frame(
    Subgroup = "Baseline type (descriptive only)", Group = g, k = fit$k,
    IRR  = round(exp(as.numeric(fit$b)), 3),
    LCL  = round(exp(fit$ci.lb),         3),
    UCL  = round(exp(fit$ci.ub),         3),
    I2   = round(as.numeric(fit$I2),     1),
    tau2 = round(as.numeric(fit$tau2),   4),
    Qbetween = NA_real_, p = NA_real_,
    note = if (used != "REML")
             paste0("fallback: ", used, "; descriptive only")
           else "Descriptive only; no between-group comparison",
    stringsAsFactors = FALSE)
})
tab_baseline <- do.call(rbind, baseline_rows)

tab_subgroups <- rbind(tab_continent, tab_risk, tab_dose, tab_sccs,
                       tab_age, tab_elderly, tab_outcome, tab_baseline)

cat("\n=== Subgroup analyses ===\n")
print(tab_subgroups, row.names = FALSE)


## ============================================================
## SECTION G: Meta-regression
## ============================================================

m0 <- m

print_mreg <- function(fit, name) {
  if (inherits(fit, "try-error")) {
    cat(sprintf("\n[Meta-reg: %s] Model fitting failed\n", name))
    return(invisible(NULL))
  }
  pR2 <- if (!is.na(m0$tau2) && m0$tau2 > 0)
    round(max(0, (m0$tau2 - as.numeric(fit$tau2)) / m0$tau2) * 100, 1) else NA
  cat(sprintf(
    "\n[Meta-reg: %s]\n  QM(df=%d)=%.3f, p=%.4f | tau2=%.4f | pseudo-R2=%.1f%%\n",
    name, fit$QMdf[2], fit$QM, fit$QMp, as.numeric(fit$tau2), pR2))
}

m_cont <- try(rma.uni(yi=dat2$yi, vi=dat2$vi, mods=~continent_bin,
                      method="REML", test="knha", data=dat2), silent=TRUE)
m_rw   <- try(rma.uni(yi=dat2$yi, vi=dat2$vi, mods=~risk_window,
                      method="REML", test="knha", data=dat2), silent=TRUE)
m_dose <- try(rma.uni(yi=dat2$yi, vi=dat2$vi, mods=~dose,
                      method="REML", test="knha", data=dat2), silent=TRUE)
m_sccs <- try(rma.uni(yi=dat2$yi, vi=dat2$vi, mods=~sccs_bin,
                      method="REML", test="knha", data=dat2), silent=TRUE)
m_win  <- try(rma.uni(yi=dat2$yi, vi=dat2$vi, mods=~window_days,
                      method="REML", test="knha", data=dat2), silent=TRUE)
m_base <- try(rma.uni(yi=dat2$yi, vi=dat2$vi, mods=~baseline_type,
                      method="REML", test="knha", data=dat2), silent=TRUE)
m_age_cont <- try(rma.uni(yi=dat2$yi, vi=dat2$vi, mods=~min_age,
                          method="REML", test="knha", data=dat2), silent=TRUE)
m_age_grp  <- try(rma.uni(yi=dat2$yi, vi=dat2$vi, mods=~age_group,
                          method="REML", test="knha", data=dat2), silent=TRUE)

print_mreg(m_cont,     "Continent (Asia vs Non-Asia)")
print_mreg(m_rw,       "Risk window (<=21d vs >21d)")
print_mreg(m_dose,     "Dose (Primary vs Booster)")
print_mreg(m_sccs,     "SCCS method (Standard vs Non-standard)")
print_mreg(m_win,      "Window days (continuous)")
print_mreg(m_base,     "Baseline type (Post-vax vs Mixed vs Pre-vax)")
print_mreg(m_age_cont, "Min inclusion age (continuous)")
print_mreg(m_age_grp,  "Age group (Adult-only vs Includes adolescents)")


## ============================================================
## SECTION H: Sensitivity analyses
## ============================================================

sens_run <- function(idx, label) {
  k_s  <- length(idx)
  fit  <- safe_rma(dat2$yi[idx], dat2$vi[idx], dat2$study[idx])
  used <- attr(fit, "method_used")
  note <- if (used != "REML") sprintf(" [fallback: %s]", used) else ""
  cat(sprintf("%-65s k=%d  IRR=%.3f [%.3f, %.3f]  I2=%.1f%%%s\n",
              label, k_s,
              exp(as.numeric(fit$b)),
              exp(fit$ci.lb), exp(fit$ci.ub),
              as.numeric(fit$I2), note))
  invisible(fit)
}

cat("\n=== Sensitivity analyses ===\n")

## (i)   Exclude Hippisley-Cox (influential study)
m_noHC   <- sens_run(which(dat2$study != "Hippisley-Cox 2021"),
                     "[Sens i]    Exclude Hippisley-Cox")
## (ii)  Standard SCCS only
m_std    <- sens_run(which(dat2$sccs_bin == "Standard"),
                     "[Sens ii]   Standard SCCS only")
## (iii) Non-standard SCCS only
m_nonstd <- sens_run(which(dat2$sccs_bin == "Non-standard"),
                     "[Sens iii]  Non-standard SCCS only")
## (iv)  Studies including booster doses
m_boost  <- sens_run(which(dat2$dose == "booster"),
                     "[Sens iv]   Booster only")
## (v)   Risk window <=21d only
m_short  <- sens_run(which(dat2$risk_window == "<=21d"),
                     "[Sens v]    Risk window <=21d only")
## (vi)  Exclude Salmaggi (pre-vaccination baseline)
m_noSal  <- sens_run(which(dat2$study != "Salmaggi 2025"),
                     "[Sens vi]   Exclude Salmaggi (pre-vax baseline)")


## ============================================================
## SECTION I: Publication bias
## ============================================================

cat("\n=== Egger's test ===\n")
cat("Note: k=8; statistical power is limited (recommended k>=10)\n")
egg <- regtest(m, model = "lm")
print(egg)

## Funnel plot (ggplot2)
se_seq <- seq(0, max(dat2$SE) * 1.1, length.out = 100)
mu_hat <- as.numeric(m$b)
df_tri <- data.frame(
  se = c(se_seq, rev(se_seq)),
  x  = c(mu_hat - 1.96 * se_seq, rev(mu_hat + 1.96 * se_seq)))
df_pts <- data.frame(
  x     = dat2$yi,
  se    = dat2$SE,
  label = dat2$study)

ggplot() +
  geom_polygon(data = df_tri, aes(x = x, y = se),
               fill = "grey88", colour = NA) +
  geom_vline(xintercept = mu_hat, linetype = "dashed", colour = "grey40") +
  geom_point(data = df_pts, aes(x = x, y = se), size = 2.5) +
  geom_text_repel(data = df_pts, aes(x = x, y = se, label = label),
                  size = 3.2, max.overlaps = 20,
                  box.padding = 0.4, point.padding = 0.3,
                  segment.colour = "grey60", segment.size = 0.35) +
  scale_y_reverse(name = "Standard Error") +
  scale_x_continuous(name = "Observed Outcome (log IRR)") +
  ggtitle("Funnel plot — BNT162b2 & Stroke") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))


## ============================================================
## SECTION J: Influence diagnostics
## ============================================================

inf <- influence(m)
cat("\n=== Influence diagnostics ===\n")
print(inf)
baujat(m, main = "Baujat plot")


## ============================================================
## SECTION K: Save outputs
## ============================================================

## Set output directory (edit as needed)
out_dir <- "output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

write.csv(tab_main,      file.path(out_dir, "Table_Main.csv"),      row.names = FALSE)
write.csv(tab_subgroups, file.path(out_dir, "Table_Subgroups.csv"), row.names = FALSE)
write.csv(loo_tab,       file.path(out_dir, "Table_LOO.csv"),       row.names = FALSE)

cat("\n=== Analysis complete (k=8) ===\n")
cat(sprintf("Output saved to: %s/\n", out_dir))
