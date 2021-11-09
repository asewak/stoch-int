require(data.table)
source("est_stoch_int.R")
source("plot_cuminc.R")
source("plot_rrr.R")
options("digits" = 4)

# Read in data
dat <- readRDS("dat.RDS")
dat[,Acum:=cumsum(A),keyby=id]

# Fit propensity score model
form <- as.formula("A ~ cavgL1 + L2 + t0")
(mod_ps <- est_lr(dat[Acum<=1], form))

# Estimate risk via a plug-in IPW estimator in the g-formula
dat_est <- copy(dat)
(mu <- est_risk(dat_est, mod_ps=mod_ps, dta=0.4, cuminc=F, comp=F))
(cuminc <- est_risk(dat_est, mod_ps=mod_ps, dta=0.4, cuminc=T, comp=F))

# Deterministic
# Note: when det=T, dta can only take either 0 or 1
# Here dta=0 indicates "never treat" and dta=1 indicates "always treat" for the full population
(cuminc_trt <- est_risk(dat_est, mod_ps=mod_ps, dta=1, det=T, cuminc=F, comp=F))
(cuminc_ctrl <- est_risk(dat_est, mod_ps=mod_ps, dta=0, det=T, cuminc=F, comp=F))

# Plots
# Plot cumulative incidence
plot_cuminc(cuminc)

# Plot relative risk reduction
K <- 10
n_reps <- 1000
dtas <- seq(0,20)/20
boot_T0 <- lapply(dtas, function(d) est_risk(dat=dat_est, mod_ps=mod_ps, dta=d, cuminc=T, comp=F))
boot_dta <- readRDS("boot_dta.RDS")
plot_rrr(boot_dta, boot_T0, level=0.95)


