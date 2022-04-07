require(data.table)
###################################################################################
##data generation for one individual

#alpha's: coefficents in treatment model
#beta1's: coefficients in model for generating the continuous covariate
#beta2's: coefficients in model for generating the binary covariates
#theta's: coefficents in outcome model
#cens's: coefficents in censoring model
datagen <- function(i, K, sigma,
                    alpha0, alpha1, alpha2, alpha3,
                    beta1_0, beta1_1, beta1_2, beta1_3, beta1_4,beta1_5, beta1_6, 
                    beta2_0, beta2_1, beta2_2, beta2_3, beta2_4, beta2_5, beta2_6,
                    theta0, theta1, theta2, theta3, theta4, theta5,theta6,
                    cens0,cens2,cens3){
  id <- as.numeric(i) # Define individual id number
  
  # Define baseline time
  t0 <- 0
  
  # Generate baseline common cause of time-varying covariates
  U <- rbinom(1, 1, 0.7) 
  
  # Generate baseline data
  
  # Continuous covariate L1
  L1 <- rnorm(1, mean=beta1_0+beta1_1*U, sd=sigma)
  cavgL1 <- cumsum(L1)[1]/1 # Calculate cumavg(L1)
  # Binary covariate L2
  L2 <- rbinom(1, 1, plogis(beta2_0+beta2_1*U+beta2_2*L1))
  
  # Binary treatment exposure indicator A
  A <- rbinom(1, 1, plogis(alpha0+alpha1*cavgL1+
                             alpha2*L2+alpha3*t0))
  
  
  
  Cen <- rbinom(1, 1, plogis(cens0+cens2*L1+cens3*L2))
  
  # Binary outcome indicator Y; write Y=2 if censored but this Y value won't be used in estimation
  
  Y <- ifelse(Cen==0,rbinom(1, 1, plogis(theta0+theta1*U+theta3*A+theta4*L1+theta5*L2)),2)
  
  
  # Coerce NA to num for future data
  enay <- 0
  enay <- NA
  
  # Generate vectors to build data.frame for data at initial time point
  id_ <- c(id)
  t0_ <- c(t0)
  U_  <- c(U)
  L1_ <- c(L1)
  cavgL1_ <- c(cavgL1)  
  
  
  L2_ = c(L2)
  A_ = c(A)
  Y_ = c(Y)
  Cen_=c(Cen)
  
  
  
  # If data for only one time point are to be generated, no more computations 
  # need to be performed
  
  # If data for multiple time points are to be generated and individual is alive
  # at end of inital measured interval, continue generating data
  if ((K > 1) && (Y==0)){
    # Generate data beyond baseline interval
    for (j in 2:K){
      # Define time interval for which data are generated
      t0 <- j-1
      # Simplified function when no data at A_[j-2]
      if (j==2){
        Ustar=rbinom(1, 1, 1/(1+exp(-0.6+A-0.2*L1+0.4*L2))) 
        
        L1star <- rnorm(1, mean=beta1_0+beta1_1*Ustar+beta1_2*A_[j-1]+
                          beta1_4*cumsum(L1_)[j-1]/(j-1)+beta1_5*L2_[j-1]+
                          beta1_6*t0, sd=sigma)
        temp_L1 <- c(L1_, L1star) # Store new and prior L1
        cavgL1 <- cumsum(temp_L1)[j]/j # Calculate cumavg(L1)
        
        L2star <- rbinom(1, 1, plogis(beta2_0+beta2_1*Ustar+beta2_2*cumsum(temp_L1)[j]/j+
                                        beta2_3*L2_[j-1]+
                                        beta2_4*A_[j-1]+beta2_6*t0))
        
        temp_L2 <- c(L2_, L2star) # Store new nad prior L2
      }
      else{
        # Function when data for A_[j-2]
        Ustar=rbinom(1, 1, 1/(1+exp(-0.6+A_[j-1]-0.2*L1_[j-1]+0.4*L2_[j-1]))) 
        
        
        L1star <- rnorm(1, mean=beta1_0+beta1_1*Ustar+beta1_2*A_[j-1]+
                          beta1_3*cumsum(A_)[j-2]/(j-2)+
                          beta1_4*cumsum(L1_)[j-1]/(j-1)+beta1_5*L2_[j-1]+
                          beta1_6*t0, sd=sigma)
        temp_L1 <- c(L1_, L1star) # Store new and prior L1
        cavgL1 <- cumsum(temp_L1)[j]/j # Calculate cumavg(L1)
        
        L2star <- rbinom(1, 1, plogis(beta2_0+beta2_1*Ustar+beta2_2*cumsum(temp_L1)[j]/j+
                                        beta2_3*L2_[j-1]+
                                        beta2_4*A_[j-1]+beta2_5*cumsum(A_)[j-2]/(j-2)+
                                        beta2_6*t0))
        temp_L2 <- c(L2_, L2star) # Store new and prior L2
      }
      
      Astar <- ifelse(A_[j-1]==0, rbinom(1, 1, plogis(alpha0+alpha1*cavgL1+
                                                        alpha2*L2star+alpha3*t0)),1)
      
      
      temp_A <- c(A_, Astar) # Store new and prior A
      
      
      Censtar <- rbinom(1, 1, plogis(cens0+cens2*L1star+cens3*L2star))
      
      Ystar <- ifelse(Censtar==0,rbinom(1, 1, plogis(theta0+theta1*Ustar+theta2*cumsum(temp_A)[j-1]/(j-1)+theta3*Astar+
                                                       theta4*cumsum(temp_L1)[j]/j+theta5*L2star+
                                                       theta6*t0)),2)
      
      # Finalize new data to add to longitudinal data frame
      
      id_[j]        <- id
      t0_[j]        <- t0
      U_[j]         <- Ustar
      L1_[j]        <- L1star
      L2_[j]        <- L2star
      
      cavgL1_[j]    <- cavgL1
      A_[j]         <- Astar
      
      Y_[j]         <- Ystar
      
      Cen_[j]       <- Censtar
      
      if(Ystar>0) break
      
    }
    
  }
  
  # Consolidate data in a single data frame
  temp_data <- data.frame(id = id_, t0 = t0_, U = U_, L1 = L1_, 
                          L2 = L2_, cavgL1 = cavgL1_, A = A_, 
                          Y = Y_, Cen=Cen_)
  return(temp_data)
}


##generate a sample of size n
mydata<-function(n, K, sigma,
                 alpha0, alpha1, alpha2, alpha3,
                 beta1_0, beta1_1, beta1_2, beta1_3, beta1_4,beta1_5, beta1_6, 
                 beta2_0, beta2_1, beta2_2, beta2_3, beta2_4, beta2_5, beta2_6,
                 theta0, theta1, theta2, theta3, theta4, theta5,theta6,
                 cens0,cens2,cens3){
  
  # Generate data with datagen function
  df <- lapply(as.list(1:n), FUN=function(ind){
    datagen(ind, K=K, sigma=sigma,
            alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, alpha3=alpha3, 
            beta1_0=beta1_0, beta1_1=beta1_1,beta1_2=beta1_2, beta1_3=beta1_3, beta1_4=beta1_4,beta1_5=beta1_5,beta1_6=beta1_6, 
            beta2_0=beta2_0, beta2_1=beta2_1, beta2_2=beta2_2,beta2_3=beta2_3, beta2_4=beta2_4, beta2_5=beta2_5,beta2_6=beta2_6, 
            theta0=theta0,theta1=theta1, theta2=theta2, theta3=theta3, theta4=theta4, theta5=theta5,theta6=theta6,
            cens0=cens0,cens2=cens2,cens3=cens3)
  })
  
  # Store final dataframe
  simdat <- data.table(rbindlist(df))
  simdat
}

