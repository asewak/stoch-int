require(data.table)
require(ggplot2)
require(Epi)

plot_rrr <- function(boot_dta, boot_T0, level=0.95){
  # Clean up
  boot_dta <- rbindlist(boot_dta)
  n_dtas <- length(boot_dta[,unique(dta)])
  boot_dta[,rep:=rep(1:n_reps, each=K*n_dtas)]
  
  # Add risk ratio
  boot_T0 <- rbindlist(boot_T0)
  boot_T0[,dta:=rep(dtas, each=K)]
  boot_T0[,Y_rr:=1-(Y_int/Y_nc)]
  boot_T0[,Y_rd:=Y_nc-Y_int]
  boot_dta[,Y_rr:=fifelse(Y_nc!=0,1-(Y_int/Y_nc),0)]
  boot_dta[,Y_rd:=Y_nc-Y_int]
  
  # Add quantiles
  
  boot_ci <- boot_dta[,.(ql_rr=quantile(Y_rr,(1-level)/2),
                         qu_rr=quantile(Y_rr,(1+level)/2),
                         mean_rr=mean(Y_rr),
                         ql_rd=quantile(Y_rd,(1-level)/2),
                         qu_rd=quantile(Y_rd,(1+level)/2)),by=.(t0,dta)]
  boot_ci[boot_T0, `:=`(Y_rr=i.Y_rr, Y_rd=i.Y_rd,
                        A_int=i.A_int), on=.(t0,dta)]
  
  # Make final plot
  res <- boot_ci[t0==K-1]
  dtas <- res$dta
  rr <- res$Y_rr
  ainc <- res$A_int
  xo <- seq(min(ainc), max(ainc), by=0.01)
  s <- spline(dtas, ainc, method = "hyman")
  a <- approx(x=s$y, y=s$x, xout = xo)
  d <- c(a$y,1)
  l <- c(round(xo,2), round(max(ainc),3))
  
  par(mai = c(1.02,0.82,0.82,0.42) * c(1, 1.5, 1, 2))
  plot(dtas, rr, type = "l", axes = FALSE, 
       xlab = expression("Shift in propensity score,"~delta),
       ylab="Relative risk reduction at end of follow-up",
       ylim = c(-0.1,1),
       panel.first = c(abline(v=seq(0,1,by=0.2), col="gray90", lwd=0.5),
                       abline(h=seq(0,1,by=0.2), col="gray90", lwd=0.5)))
  matshade(res$dta, res[, .(Y_rr,ql_rr,qu_rr)], )
  axis(1)
  axis(2)
  axis(3, at = d, label = l)
  mtext("Population treatment uptake at end of follow-up", side = 3, line = 2.5)
}

