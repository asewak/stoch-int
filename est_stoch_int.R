require(data.table)

# Function to estimate the propensity score
est_lr <- function(dat, form){
  # Propensity score model and predictions
  fit_lr <- glm(form, family=binomial, data=dat)
  return(fit_lr)
}

# Function to estimate the risk
est_risk <- function(dat, mod_ps, dta=0.1, det=F, cuminc=F, comp=F){
  # Add sample ID from bootstrap
  if(!"samp" %in% colnames(dat)){dat[,samp:=id]}
  
  # Add propensity scores (initiation)
  dat[Acum<=1,ps:=mod_ps$fitted.values]
  dat[Acum>1,qs:=1]
  dat[Acum>1,ps:=1]
  
  if(det){
    if( !dta %in% c(0,1) ) stop('dta is not 0 or 1')
    dat[Acum<=1,qs:=dta]
  } else {
    dat[Acum<=1,qs:=pmin(1,ps+dta)*rf + ps*(1-rf)]
  }

  # Create numerator and denominator for the weights
  dat[,numer:=1]
  dat[,denom:=1]
  dat[Acum==1,`:=`(numer=qs, denom=ps)]
  dat[Acum==0, `:=`(numer=1-qs, denom=1-ps)]
  dat[,ratio:=numer/denom]
  
  # Censoring
  dat[,censratio:=1]
  
  # Weights for cumulative incidence of outcome
  dat[,`:=`(wt=cumprod(ratio),
            censwt=cumprod(censratio),
            ccwt=cumprod(ratio*censratio)), by=.(id,samp)]
  
  # Weights for cumulative incidence of treatment
  # "s": shifted by one interval, "h": hazard formulation
  dat[,scwth:=shift(ccwt,1,1),keyby=.(id,samp)]
  
  # Not accounting for competing risks
  dat_l <- dat[,.(lambda=sum(fifelse(Y==1,ccwt,0)/sum(ccwt)),
                  lambda0=sum(fifelse(Y==1,censwt,0)/sum(censwt))), keyby=t0]
  dat_l[,`:=`(lamsurr=shift(lambda,1,fill=0),
              lamsurr0=shift(lambda0,1,fill=0))]
  dat_l[,`:=`(cumsurv=cumprod(1-lamsurr),
              cumsurv0=cumprod(1-lamsurr0))]
  mu <- dat_l[,.(Y_int=sum(lambda*cumsurv),
                 Y_nc=sum(lambda0*cumsurv0))]
  
  # Cumulative incidence of outcome (no competing risks)
  cuminc_base <- dat_l[,.(t0,
                          Y_int=cumsum(lambda*cumsurv), 
                          Y_nc=cumsum(lambda0*cumsurv0))]
  
  # Accounting for competing risks
  dat_comp_l <- dat[,.(lambda1=sum(fifelse(Y==1 & D==0,ccwt,0)/sum(fifelse(D==0,ccwt,0))),
                       lambda2=sum(fifelse(D==1,ccwt,0)/sum(ccwt)),
                       lambda01=sum(fifelse(Y==1 & D==0,censwt,0)/sum(fifelse(D==0,censwt,0))),
                       lambda02=sum(fifelse(D==1,censwt,0)/sum(censwt))),keyby=t0]
  dat_comp_l[,`:=`(lamsurr1=shift(lambda1,1,fill=0),
                   lamsurr2=shift(lambda2,1,fill=0),
                   lamsurr01=shift(lambda01,1,fill=0),
                   lamsurr02=shift(lambda02,1,fill=0))]
  dat_comp_l[,`:=`(cumsurv=cumprod((1-lamsurr1)*(1-lamsurr2)),
                   cumsurv0=cumprod((1-lamsurr01)*(1-lamsurr02)))]
  
  # Risk at the end of follow-up
  mu_comp <- dat_comp_l[,.(Y_int=sum(lambda1*(1-lambda2)*cumsurv),
                           Y_nc=sum(lambda01*(1-lambda02)*cumsurv0))]
  
  # Cumulative incidence of outcome (with competing risks)
  cuminc_comp <- dat_comp_l[,.(t0,
                               Y_int=cumsum(lambda1*(1-lambda2)*cumsurv), 
                               Y_nc=cumsum(lambda01*(1-lambda02)*cumsurv0),
                               D_int=cumsum(lambda2*cumsurv),
                               D_nc=cumsum(lambda02*cumsurv0))]
  
  # Cumulative incidence of treatment
  n <- dat[,length(unique(samp))]
  cuminc_A <- dat[Acum<=1, .(h_int=mean(qs*scwth)*.N/n, 
                             h_nc=mean(ps)*.N/n), by=t0]
  cuminc_A[,`:=`(A_int=cumsum(h_int),
                 A_nc=cumsum(h_nc))]
  mu_A <- cuminc_A[,.(A_int=sum(h_int),
                      A_nc=sum(h_nc))]
  
  # Final output table
  if(comp) {
    cuminc_all <- cuminc_comp[cuminc_A,]
    mu_all <- unlist(c(mu_comp,mu_A))
  } else {
    cuminc_all <- cuminc_base[cuminc_A,]
    mu_all <- unlist(c(mu,mu_A))
  }
  
  if(cuminc) cuminc_all else mu_all
}

# Example
if(FALSE){
  mod_ps <- est_lr(dat[Acum<=1], form)
  est_risk(dat, mod_ps=mod_ps, dta=0.1)
}


