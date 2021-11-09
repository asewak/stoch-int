require(data.table)
require(parallel)
source("est_stoch_int.R")
set.seed(100)
n_cores <- 6

# Read in data
dat <- readRDS("dat.RDS")
setkey(dat, "id")

# Specify propensity score model
form <- as.formula("A ~ cavgL1 + L2 + t0")

# Bootstrap IDs
fn_boot_resamp <- function(dat, ids, dtas, K){
  ids_star <- sample(ids, replace = T)
  dat_id <- data.table(id=ids_star, samp=seq(length(ids)), key = "id")
  dat_tmp <- copy(dat[dat_id, allow.cartesian=TRUE])
  mod_ps <- est_lr(dat_tmp[Acum<=1], form)
  cuminc <- lapply(dtas, function(d) est_risk(dat=dat_tmp, mod_ps=mod_ps, dta=d, cuminc=T, comp=F))
  cuminc_dta <- rbindlist(cuminc)
  cuminc_dta[,dta:=rep(dtas, each=K)]
  return(cuminc_dta)
}

ids <- dat[,unique(id)]
n_reps <- 1000
reps <- seq(1,n_reps)
dtas <- seq(0,20)/20
K <- 10

boot_resamp_dta <- mclapply(reps, function(x) fn_boot_resamp(dat, ids, dtas, K), mc.cores = n_cores)
saveRDS(boot_resamp_dta, "boot_dta.RDS")

