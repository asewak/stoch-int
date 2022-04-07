source("sim.R")
set.seed(100)

dat <- mydata(n=10000, K=10, sigma=2,
              alpha0=-3, alpha1=0.5, alpha2=-0.3, alpha3=0.001,
              beta1_0=-0.1, beta1_1=0.3, beta1_2=-0.7, beta1_3=-0.1, beta1_4=0.5, beta1_5=0.1, beta1_6=0.001,
              beta2_0=-4, beta2_1=1, beta2_2=0.2, beta2_3=0.5, beta2_4=-0.7, beta2_5=-0.1, beta2_6=0.001, 
              theta0=-5, theta1=-0.3, theta2=-2, theta3=-2, theta4=0.5, theta5=3, theta6=0.01,
              cens0=-100, cens2=0, cens3=0)
dat[,Acum := cumsum(A), by = id]
dat[,`:=`(rf = L2, 
          D = 0)]
saveRDS(dat, "dat.RDS")
