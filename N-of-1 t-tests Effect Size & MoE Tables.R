

library(pwr)

#=========================================== Std Err Functions ============================================
pst.se <- function(SD, RHO, N){
  c <- (N + 2*RHO^(N+1) - N*RHO^2 - 2*RHO) / (N^2*(RHO-1)^2)
  b <- ( N - (N+2*RHO^(N+1)-N*RHO^2-2*RHO)/(N*(RHO-1)^2) ) / (N-1)
  nprime <- 1 / (1 - b*(N-1)/N)
  se <- SD * sqrt( c/b )
  cbnse <- cbind(nprime=nprime, se=se, c=c, b=b)
  return(cbnse)
}

pstr.se <- function(SD, RHO, N){
  c <- 12/(N^2-1)^2 * ( - 6*RHO*(RHO+1)^2*(RHO^N-1)/(N^2*(RHO-1)^4) +
                               2*RHO*(6*RHO^(N+1)+6*RHO^N+RHO^2-2*RHO+1)/(N*(RHO-1)^3) -
                               6*RHO*(RHO^N+1)/(RHO-1)^2 - 2*N*RHO/(RHO-1) + (N^2-1)/N )
  b <- (N - (2 + 2*RHO*(RHO^N-N*RHO+N-1)/(N*(RHO-1)^2) +
                    (-2*RHO/(N*(N^2-1)*(RHO-1)^4) *
                       (N^3*(RHO-1)^3 + 3*N^2*(RHO-1)^2*(RHO^N+1) +
                          3*(RHO+1)^2*(RHO^N-1) -
                          N*(RHO-1)*(6*RHO^(N+1)+6*RHO^N+RHO^2-2*RHO+1))))) / (N-2)
  nprime <- 2 / (1-b*(N-2)/N)
  se <- SD * sqrt( c/b )
  cbnse <- cbind(nprime=nprime, se=se, c=c, b=b)
  return(cbnse)
}

st.se <- function(SD, RHO, N1, N2=N1){
  N <- N1
  c1 <- (N + 2*RHO^(N+1) - N*RHO^2 - 2*RHO) / (N^2*(RHO-1)^2)
  b1 <- ( N - (N+2*RHO^(N+1)-N*RHO^2-2*RHO)/(N*(RHO-1)^2) ) / (N-1)
  nprime1 <- 1 / (1 - b1*(N-1)/N)
  
  N <- N2
  c2 <- (N + 2*RHO^(N+1) - N*RHO^2 - 2*RHO) / (N^2*(RHO-1)^2)
  b2 <- ( N - (N+2*RHO^(N+1)-N*RHO^2-2*RHO)/(N*(RHO-1)^2) ) / (N-1)
  nprime2 <- 1 / (1 - b2*(N-1)/N)
  
  nprime <- nprime1 + nprime2
  se <- SD * sqrt( c1/b1 + c2/b2 )
  cbnse <- cbind(nprime=nprime, se=se, c1=c1, b1=b1, c2=c2, b2=b2)
  return(cbnse)
}

str.se <- function(SD, RHO, N1, N2=N1){
  N <- N1
  c1 <- 12/(N^2-1)^2 * ( - 6*RHO*(RHO+1)^2*(RHO^N-1)/(N^2*(RHO-1)^4) +
                          2*RHO*(6*RHO^(N+1)+6*RHO^N+RHO^2-2*RHO+1)/(N*(RHO-1)^3) -
                          6*RHO*(RHO^N+1)/(RHO-1)^2 - 2*N*RHO/(RHO-1) + (N^2-1)/N )
  b1 <- (N - (2 + 2*RHO*(RHO^N-N*RHO+N-1)/(N*(RHO-1)^2) +
               (-2*RHO/(N*(N^2-1)*(RHO-1)^4) *
                  (N^3*(RHO-1)^3 + 3*N^2*(RHO-1)^2*(RHO^N+1) +
                     3*(RHO+1)^2*(RHO^N-1) -
                     N*(RHO-1)*(6*RHO^(N+1)+6*RHO^N+RHO^2-2*RHO+1))))) / (N-2)
  nprime1 <- 2 / (1-b1*(N-2)/N)
  
  N <- N2
  c2 <- 12/(N^2-1)^2 * ( - 6*RHO*(RHO+1)^2*(RHO^N-1)/(N^2*(RHO-1)^4) +
                           2*RHO*(6*RHO^(N+1)+6*RHO^N+RHO^2-2*RHO+1)/(N*(RHO-1)^3) -
                           6*RHO*(RHO^N+1)/(RHO-1)^2 - 2*N*RHO/(RHO-1) + (N^2-1)/N )
  b2 <- (N - (2 + 2*RHO*(RHO^N-N*RHO+N-1)/(N*(RHO-1)^2) +
                (-2*RHO/(N*(N^2-1)*(RHO-1)^4) *
                   (N^3*(RHO-1)^3 + 3*N^2*(RHO-1)^2*(RHO^N+1) +
                      3*(RHO+1)^2*(RHO^N-1) -
                      N*(RHO-1)*(6*RHO^(N+1)+6*RHO^N+RHO^2-2*RHO+1))))) / (N-2)
  nprime2 <- 2 / (1-b2*(N-2)/N)
  
  
  nprime <- nprime1 + nprime2
  se <- SD * sqrt( c1/b1 + c2/b2 )
  cbnse <- cbind(nprime=nprime, se=se, c1=c1, b1=b1, c2=c2, b2=b2)
  return(cbnse)
}


#=================================== Computing Margin of Error Tables =====================================

#----- Setting the parameters -----
ALPHA <-  0.05
ONE.SIDED <- FALSE
RHO <- seq(0,0.8, 0.2)
N <- c(4,8,12,100)
SIGFIG <- 2  # Rounds output

#----- Margin of Error Table for  pst -----
MOE.TABLE.pst <- data.frame()
for(rho in RHO){
  
  n.MOE <- NULL
  column.name <- NULL
  for(n in N){
    
    # Getting nprime for pst        
    se.parts <- pst.se(SD=1, RHO=rho, N=n)
    nprime <- se.parts[1]
    se <- se.parts[2]
    df <- nprime - 1
    .alpha <- ifelse(ONE.SIDED, ALPHA, ALPHA/2)
    tcrit <- qt(1-.alpha, df, lower.tail=TRUE)
    MOE <- tcrit*se
    n.MOE <- c(n.MOE, round(MOE, SIGFIG))
    column.name <- c(column.name, paste0("m=",n))
  }
  MOE.TABLE.pst <- rbind(MOE.TABLE.pst, c(ALPHA, ONE.SIDED, rho, n.MOE) )
  colnames(MOE.TABLE.pst) <- c('alpha', 'one-sided', 'rho', column.name)
} 

#----- Margin of Error Table for st -----
MOE.TABLE.st <- data.frame()
for(rho in RHO){
  
  n.MOE <- NULL
  column.name <- NULL
  for(n in N){
    
    # Getting nprime for pst
    se.parts <- st.se(SD=1, RHO=rho, N1=n, N2=n)
    nprime <- se.parts[1]
    se <- se.parts[2]
    df <- nprime - 2
    .alpha <- ifelse(ONE.SIDED, ALPHA, ALPHA/2)
    tcrit <- qt(1-.alpha, df, lower.tail=TRUE)
    MOE <- tcrit*se
    n.MOE <- c(n.MOE, round(MOE,SIGFIG))
    column.name <- c(column.name, paste0("m=",n))
    
  }
  MOE.TABLE.st <- rbind(MOE.TABLE.st, c(ALPHA, ONE.SIDED, rho, n.MOE) )
  colnames(MOE.TABLE.st) <- c('alpha', 'one-sided', 'rho', column.name)
}

#----- Margin of Error Table for str -----
MOE.TABLE.str <- data.frame()
for(rho in RHO){
  
  n.MOE <- NULL
  column.name <- NULL
  for(n in N){
    
    # Getting nprime for pst
    se.parts <- str.se(SD=1, RHO=rho, N1=n, N2=n)
    nprime <- se.parts[1]
    se <- se.parts[2]
    df <- nprime - 4
    .alpha <- ifelse(ONE.SIDED, ALPHA, ALPHA/2)
    tcrit <- qt(1-.alpha, df, lower.tail=TRUE)
    MOE <- tcrit*se
    n.MOE <- c(n.MOE, round(MOE,SIGFIG))
    column.name <- c(column.name, paste0("m=",n))
    
  }
  MOE.TABLE.str <- rbind(MOE.TABLE.str, c(ALPHA, ONE.SIDED, rho, n.MOE) )
  colnames(MOE.TABLE.str) <- c('alpha', 'one-sided', 'rho', column.name)
}

#----- Margin of Error Table for pstr -----
MOE.TABLE.pstr <- data.frame()
for(rho in RHO){
  
  n.MOE <- NULL
  column.name <- NULL
  for(n in N){
    
    # Getting nprime for pst            
    se.parts <- pstr.se(SD=1, RHO=rho, N=n)
    nprime <- se.parts[1]
    se <- se.parts[2]
    df <- nprime - 2
    .alpha <- ifelse(ONE.SIDED, ALPHA, ALPHA/2)
    tcrit <- qt(1-.alpha, df, lower.tail=TRUE)
    MOE <- tcrit*se
    n.MOE <- c(n.MOE, round(MOE,SIGFIG))
    column.name <- c(column.name, paste0("m=",n))
    
  }
  MOE.TABLE.pstr <- rbind(MOE.TABLE.pstr, c(ALPHA, ONE.SIDED, rho, n.MOE) )
  colnames(MOE.TABLE.pstr) <- c('alpha', 'one-sided', 'rho', column.name)
}
#----- Printing Margin of Error Tables -----
MOE.TABLE.pst
MOE.TABLE.st
MOE.TABLE.str
MOE.TABLE.pstr



#================================== Computing Eff Size Tables =====================================
POWER <- 0.80
ALPHA <- 0.05
ONE.SIDED <- TRUE
ALTERNATIVE <- ifelse(ONE.SIDED, "greater","two.sided")
RHO <- seq(0, 0.8, 0.2)
N.level <- c(4:12)   # Minimum allowed is 4
N.rate <- c(5:12)    # Minimum allowed is 5
SIGFIG <- 2
#------- Effect size for pst test ---------
ES.TABLE.pst <- data.frame()
for(rho in RHO)   {
  row.tmp <- NULL
  column.name <- NULL
  for(n in N.level){
    se.parts <- pst.se(SD=1, RHO=rho, N=n)
    c <- se.parts[3] #---------------------------- Variance of delta estimator
    nprime <- se.parts[1]    #-------------------- Effective sample size
    n.a <- nprime    #---------------------------- Assumed sample size for assumed one-sample t-test
    d.a <- ifelse(n.a >= 1.8,  #------------------ Effect size from assumed one-sample t-test. It seems  
                  pwr.t.test(n=n.a,  #             n.a needs to be greater than 1.8 for pwr.t.test to work.
                             d=NULL, #             We set the effect size to "inf" if n.a <= 1.8.
                             sig.level=ALPHA, 
                             power = POWER,
                             type="one.sample", 
                             alternative=ALTERNATIVE)$d,
                  Inf)
    d <- d.a*sqrt(n.a*c) #------------------ Effect size for pst test
    row.tmp <- c(row.tmp, round(d,SIGFIG))
    column.name <- c(column.name, paste0("m=",n))
    rm(se.parts, nprime,n.a,d.a,d, c)
  }
  ES.TABLE.pst <- rbind(ES.TABLE.pst, c(POWER, ALPHA, ONE.SIDED, rho,row.tmp))
  colnames(ES.TABLE.pst) <- c('power', 'alpha', 'one-sided', 'rho', column.name)
}


#------- Effect size for st test ---------
ES.TABLE.st <- data.frame()
for(rho in RHO)   {
  row.tmp <- NULL
  column.name <- NULL
  for(n in N.level){
    se.parts <- st.se(SD=1, RHO=rho, N1=n, N2=n)
    c <- se.parts[3] + se.parts[5]
    nprime <- se.parts[1]/2   #-------------------- Effective sample size in 1 group; 
    #                                               both groups assumed to have same sample size
    n.a <- 2*nprime -1  #-------------------------- Assumed sample size for assumed one-sample t-test
    d.a <- ifelse(n.a >= 1.8,  #------------------- Effect size from assumed one-sample t-test. It seems
                  pwr.t.test(n=n.a, #               n.a needs to be greater than 1.8 for pwr.t.test to work.
                             d=NULL, #              We set the effect size to "inf" if n.a <= 1.8.
                             sig.level=ALPHA, 
                             power = POWER,
                             type="one.sample", 
                             alternative=ALTERNATIVE)$d,
                  Inf)
    d <- d.a*sqrt(n.a*c ) #---------------- Effect size for st test
    row.tmp <- c(row.tmp, round(d,SIGFIG))
    column.name <- c(column.name, paste0("m=",n))
    rm(se.parts, nprime,n.a,d.a,d, c)
  }
  ES.TABLE.st <- rbind(ES.TABLE.st, c(POWER, ALPHA, ONE.SIDED, rho,row.tmp))
  colnames(ES.TABLE.st) <- c('power', 'alpha', 'one-sided', 'rho', column.name)
}





#------- Effect size for pstr test ---------
ES.TABLE.pstr <- data.frame()
for(rho in RHO)   {
  row.tmp <- NULL
  column.name <- NULL
  for(n in N.rate){
    se.parts <- pstr.se(SD=1, RHO=rho, N=n)
    c <- se.parts[3]
    nprime <- se.parts[1]    #-------------------- Effective sample size
    sig.x <- sqrt( (n^2 - 1) / 12 ) #------------- SD of X vector 
    n.a <- nprime -1  #--------------------------- Assumed sample size for assumed one-sample t-test
    d.a <- ifelse(n.a >= 1.8,  #------------------ Effect size from assumed one-sample t-test. It seems
                  pwr.t.test(n=n.a, #              n.a needs to be greater than 1.8 for pwr.t.test to work.
                             d=NULL, #             We set the effect size to "inf" if n.a <= 1.8.
                             sig.level=ALPHA, 
                             power = POWER,
                             type="one.sample", 
                             alternative=ALTERNATIVE)$d,
                  Inf)
    d <- d.a*sqrt(n.a*c) #-------------- Effect size for pstr test
    row.tmp <- c(row.tmp, round(d,SIGFIG))
    column.name <- c(column.name, paste0("m=",n))
    rm(se.parts, nprime,n.a,d.a,d, c)
  }
  ES.TABLE.pstr <- rbind(ES.TABLE.pstr, c(POWER, ALPHA, ONE.SIDED, rho,row.tmp))
  colnames(ES.TABLE.pstr) <- c('power', 'alpha', 'one-sided', 'rho', column.name)
}


#------- Effect size for str test ---------
ES.TABLE.str <- data.frame()
for(rho in RHO)   {
    row.tmp <- NULL
    column.name <- NULL
  for(n in N.rate){
    se.parts <- str.se(SD=1, RHO=rho, N1=n, N2=n)
    c <- se.parts[3] + se.parts[5]
    nprime <- se.parts[1]/2  #-------------------- Effective sample size in 1 group; 
    #                                              both groups assumed to have same sample size
    sig.x <- sqrt( (n^2 - 1) / 12 ) #------------- SD of X vector 
    n.a <- 2*nprime -3  #------------------------- Assumed sample size for assumed one-sample t-test
    d.a <- ifelse(n.a >= 1.8,  #------------------ Effect size from assumed one-sample t-test. It seems
                  pwr.t.test(n=n.a, #              n.a needs to be greater than 1.8 for pwr.t.test to work.
                             d=NULL, #             We set the effect size to "inf" if n.a <= 1.8.
                             sig.level=ALPHA, 
                             power = POWER,
                             type="one.sample", 
                             alternative=ALTERNATIVE)$d,
                  Inf)
    d <- d.a*sqrt(n.a*c) #----------- Effect size for str test-
    row.tmp <- c(row.tmp, round(d, SIGFIG))
    column.name <- c(column.name, paste0("m=",n))
    rm(se.parts, nprime,n.a,d.a,d, c)
    }
    ES.TABLE.str <- rbind(ES.TABLE.str, c(POWER, ALPHA, ONE.SIDED, rho,row.tmp))
    colnames(ES.TABLE.str) <- c('power', 'alpha', 'one-sided', 'rho', column.name)
}

#----- Printing Effect Size Tables -----
ES.TABLE.pst
ES.TABLE.st
ES.TABLE.pstr
ES.TABLE.str
    
    
#========================= Computing Individual Eff Size via Grid-Search ==========================

#----- Function computing P(T >= t.crit | delta) -----
est.power <- function(DELTA, DF, DIFF.VEC, SE.VEC, SIGNIF){
  t.stat <- (DIFF.VEC + DELTA) / SE.VEC
  prob.t <- pt(t.stat, df=DF, lower.tail=FALSE)
  pow <- mean(prob.t<=SIGNIF)
  return(power = pow)
}
#----- Grid-search function for single effect size  -----
ES.single.case <- function(test, n, rho, power, alpha, one.sided, sigfig){
  # Sets correct degrees of freedom for serial t-tests
  # and gets the variance ("c") of estimator.
  if(test=="st"){
    se.parts.truerho <- st.se(SD=0, RHO=rho, N1=n, N2=n)  
    truerho.df <- unname(se.parts.truerho[,"nprime"]) - 2
    truerho.c <- 2 * unname(se.parts.truerho[,"c1"])
  } else if(test=="str"){
    se.parts.truerho <- str.se(SD=0, RHO=rho, N1=n, N2=n)  
    truerho.df <- unname(se.parts.truerho[,"nprime"]) - 4
    truerho.c <- 2 * unname(se.parts.truerho[,"c1"])
  } else if(test=="pst"){
    se.parts.truerho <- pst.se(SD=0, RHO=rho, N=n)  
    truerho.df <- unname(se.parts.truerho[,"nprime"]) - 1
    truerho.c <- unname(se.parts.truerho[,"c"])
  } else if(test=="pstr"){
    se.parts.truerho <- pstr.se(SD=0, RHO=rho, N=n)  
    truerho.df <- unname(se.parts.truerho[,"nprime"]) - 2
    truerho.c <- unname(se.parts.truerho[,"c"])
  } else{
    stop("test not found!")
  }
  
  t.crit <- qt(1-alpha, df=truerho.df)  # critical value for testing H0
  
  # Determine effect size to achieve given power by searching from 0:.dmax
  notFound <- TRUE  # changes to false when correct effect size for given power is found
  upper <- .dmax
  lower <- 0
  increment <- 1
  while(notFound){
    upperNotFound <- TRUE  # changes to false when upper bound for effect size is found
    while(upperNotFound){
      estpow.up <- round(pt(t.crit, df=truerho.df, ncp=upper / sqrt(truerho.c),
                            lower.tail=FALSE), digits=sigfig)
      if(estpow.up == power){
        upperNotFound <- FALSE
      } else if(estpow.up < power){
        upperNotFound <- FALSE
        upper <- upper + increment
      } else{
        upper <- upper - increment
      }
    } 
    
    lowerNotFound <- TRUE  # changes to false when lower bound for effect size is found
    while(lowerNotFound){
      estpow.low <- round(pt(t.crit, df=truerho.df, ncp=lower / sqrt(truerho.c),
                             lower.tail=FALSE), digits=sigfig)
      if(estpow.low == power){
        lowerNotFound <- FALSE
      } else if(estpow.low > power){
        lowerNotFound <- FALSE
        lower <- lower - increment
      } else{
        lower <- lower + increment
      }
    }
    
    if(upper == lower){
      notFound <- FALSE
      effectsize <- upper
    } else if(estpow.low==power && estpow.up==power){
      notFound <- FALSE
      effectsize <- (upper + lower) / 2
    } else if(estpow.low==power && estpow.up!=power){
      lower <- lower - increment
      increment <- increment / 10
    } else if(estpow.up==power && estpow.low!=power){
      upper <- upper + increment
      increment <- increment / 10
    } else{
      increment <- increment / 10   # If effect size not found, check
      # smaller range with greater precision
    }
  } 
  output <- c(power, ALPHA, ONE.SIDED, rho, n, effectsize)
  names(output) <- c("power", "alpha", "one-sided", "rho", "m", "effect.size")
  return(output)
}

#----- Setting the parameters -----
        N <- 4
      RHO <- 0.8
     TEST <- "pst"  
    D.MAX <- 10000    # Upper bound of grid search 
    POWER <- 0.80
    ALPHA <- 0.05
ONE.SIDED <- TRUE     # alpha <- ifelse(ONE.SIDED, ALPHA, ALPHA/2)
   SIGFIG <- 2

#----- Single case effect size calculation -----
ES.single.case(test = TEST, n = N, rho = RHO, power = POWER, alpha = ALPHA, 
               one.sided = ONE.SIDED, sigfig = SIGFIG)





