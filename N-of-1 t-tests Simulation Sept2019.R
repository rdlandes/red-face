#==========================================================================================================
#**********************************************************************************************************
# N-of-1 t-tests: Simulation datasets generator (19 August 2019) ----
# Run time: approx 75 minutes
#           on a Dell Precision M6800 with an 
#           Intel Core i7-8650U CPU @ 1.90 GHz
#           with 32.0 GB of RAM
#           R version 3.5.3
#**********************************************************************************************************
#==========================================================================================================


#==========================================================================================================
#------------------------------------------------- Set-up -------------------------------------------------
#==========================================================================================================
# Change to correct folder
# setwd("<insert location of datasets outputted from simulation>")
setwd("C:/Temp for simulation data")

# Load libraries
library(MASS)
library(data.table)
library(plotrix)

# Parameters
#  .n is the number of observations from a treatment; referred to as "m" in article text. 
#     Can take values {4,5,6,...}. Default is {4:12,30,50,100}.
#  .rho.w is the serial correlation. Can take values between (-1,1).
#     Default is {-0.33, 0, 0.33, 0.67}.
#  .paired indicates whether the data from 2 treatments are paired (see ey in sim loop).
#     Can take values {FALSE, TRUE}. Default has both.
#     **NOTE** If .paired includes FALSE ..rho.a must include 0. Likewise if .paired
#              includes TRUE ..rho.a will include only non-zero values.
#  ..rho.a is the correlation in observations from paired treatments. 
#     Can take values between (-1,1) but see note under .paired. Default is {0, 0.33, 0.67}.
#  .sigsqr is the variance. Default variance is 1.
#  .M is the number of iterations for each combination of test, .n, rho.a, and rho.w.
#     Default is 10000 to maintain Monte Carlo error within 0.01.
#  .delta is the true difference between the 2 treatments. Default is 0, the null value. 
.n <- c(4:12,30,50,100)
.M <- 10000  # Results based on .M=10000.
.paired <- c(FALSE, TRUE)
..rho.a <- c(0, 0.33, 0.67)
.rho.w <- c(-0.33, 0, 0.33, 0.67)
.sigsqr <- 1
.slope <- 1
.delta <- 0

# List to hold generated data sets, counter for list
SUMSTAT <- list()
COMBO <- 1


#==========================================================================================================
#----------------------------------------------- Functions ------------------------------------------------
#==========================================================================================================
#--------------------- Function computing Wayne Fuller's serial correlation estimator ---------------------
WF.r <- function(resid){
  .n <- length(resid)
  r.bar <- mean(resid) 
  .c0 <-  (resid[1] - r.bar)^2 
  .c1 <- 0
  for(k in 2:.n){
    .c0 <- .c0 +  (resid[k] - r.bar)^2  
    .c1 <- .c1 +  (resid[k-1] - r.bar) * (resid[k] - r.bar) 
  } 
  mle.r <- .c1 / .c0
  WF.r <- mle.r + (1-mle.r^2)/(.n-1)
  return(WF.r)
}

#--------------------- Summary statistics for 2-sample serial t-test for rate-change ----------------------
# Returns difference in intercepts, difference in slopes, standard deviation, serial correlation estimate.
#
# Parameters
#  vA is the vector of data under condition A.
#  vB is the vector of data under condition B.
#  rhotype is the type of serial correlation estimate to use. If "WF" (Wayne Fuller's) is not
#    specified, rhow must be inputted and will be used as the serial correlation estimate.
#  rhow is an inputted serial correlation. Only used if rhotype is not specified.
strfxn <- function(vA, vB, rhotype, rhow){
  # Counting length of sequences from 2 "independent" conditions
  nA <- length(vA)
  nB <- length(vB)
  
  # Constructing the variables for the linear model
  oneA <- c(rep(1, nA), rep(0, nB))
  oneB <- c(rep(0, nA), rep(1, nB))
  xA <- c( 1:nA - mean(1:nA), 
           rep(0, nB) )
  xB <- c( rep(0, nA),
           1:nB - mean(1:nB) )
  y <- c(vA, vB)
  
  # Fitting the linear model & obtaining summary stats
  m <- summary( lm(y ~ oneA + oneB + xA + xB -1))
  intA <- m$coefficients[1,1]
  intB <- m$coefficients[2,1]
  bA   <- m$coefficients[3,1]
  bB   <- m$coefficients[4,1]
  sd.yr <- m$sigma 
  names(intA) <- NULL
  names(intB) <- NULL
  names(bA)   <- NULL
  names(bB)   <- NULL
  names(sd.yr) <- NULL
  
  # Obtaining residuals & computing correlation
  resA <- m$residuals[1:nA]
  resB <- m$residuals[-(1:nA)]
  if(rhotype=="WF") {
    rhoB <- WF.r(resB)
    rhoA <- WF.r(resA)
    rho <- (nA*rhoA + nB*rhoB) / (nA + nB)
  } else{
    rho <- rhow
  }
  
  # In case serial correlation fails to compute, set all to NA
  if(is.na(rho)){
    unp.int <- NA
    unp.b <- NA
    sd.yr <- NA
  } else {
    unp.int <- intB - intA
    unp.b <- bB - bA
  }
  
  # Gathering summary statistics
  strSum <- cbind(unp.int, unp.b, sd.yr, rho)
  colnames(strSum) <- c("unp.int","unp.b", "sd.yr", paste0("rho.unp.b.", rhotype))
  return(strSum)
} 

#--------------------- Summary statistics for 2-sample serial t-test for level-change ---------------------
# Returns difference in means, standard deviation, serial correlation estimate.
#
# Parameters
#  vA is the vector of data under condition A.
#  vB is the vector of data under condition B.
#  rhotype is the type of serial correlation estimate to use. If "WF" (Wayne Fuller's) is not
#    specified, rhow must be inputted and will be used as the serial correlation estimate.
#  rhow is an inputted serial correlation. Only used if rhotype is not specified.
stfxn <- function(vA, vB, rhotype, rhow) {
  # Counting length of sequences from 2 "independent" conditions
  nA <- length(vA)
  nB <- length(vB)
  
  # Constructing the variables for the linear model
  oneA <- c(rep(1, nA), rep(0, nB))
  oneB <- c(rep(0, nA), rep(1, nB))
  y <- c(vA, vB)
  
  # Fitting the linear model & obtaining summary stats
  m <- summary( lm(y ~ oneA + oneB -1 ))
  intA <- m$coefficients[1,1]
  intB <- m$coefficients[2,1]
  sd.y <- m$sigma 
  names(intA) <- NULL
  names(intB) <- NULL
  names(sd.y) <- NULL
  
  # Obtaining residuals & computing correlation
  resA <- m$residuals[1:nA]
  resB <- m$residuals[-(1:nA)]
  if(rhotype=="WF") {
    rhoB <- WF.r(resB)
    rhoA <- WF.r(resA)
    rho <- (nA*rhoA + nB*rhoB) / (nA + nB)
  } else{
    rho <- rhow
  }
  
  # In case serial correlation fails to compute, set all to NA
  if(is.na(rho)){
    unp.d <- NA
    sd.y <- NA
  } else {
    unp.d <- intB - intA
  }
  
  # Gathering summary statistics
  stSum <- cbind(unp.d, sd.y, rho)
  colnames(stSum) <- c("unp.d", "sd.y", paste0("rho.unp.d.", rhotype))
  return(stSum)
}

#---------------------- Summary statistics for paired serial t-test for level-change ----------------------
# Returns mean difference, standard deviation, serial correlation estimate.
#
# Parameters
#  vA is the vector of data under condition A.
#  vB is the vector of data under condition B.
#  rhotype is the type of serial correlation estimate to use. If "WF" (Wayne Fuller's) is not
#    specified, rhow must be inputted and will be used as the serial correlation estimate.
#  rhow is an inputted serial correlation. Only used if rhotype is not specified.
pstfxn <- function(vA, vB, rhotype, rhow) {
  # Calculating vector of differences
  diff <- vB - vA
  
  # Counting length of the sequence of paired differences between conditions
  n <- length(diff)
  
  # Constructing the variables for the linear model
  one <- rep(1, n)
  y <- diff
  
  # Fitting the linear model & obtaining summary stats
  m <- summary( lm(y ~ one -1 ))
  d <- m$coefficients[1,1]
  sd.d <- m$sigma 
  names(d) <- NULL
  names(sd.d) <- NULL
  
  # Obtaining residuals & computing correlation
  res <- m$residuals
  if(rhotype=="WF"){
    rho <- WF.r(res)
  } else{
    rho <- rhow
  }
  
  # In case serial correlation fails to compute, set all to NA
  if(is.na(rho)){
    d <- NA
    sd.d <- NA
  }
  
  # Gathering summary statistics
  pstSum <- cbind(d, sd.d, rho)
  colnames(pstSum) <- c("d", "sd.d", paste0("rho.d.", rhotype))
  return(pstSum)
}

#---------------------- Summary statistics for paired serial t-test for rate-change -----------------------
# Returns intercept of differences, slope of differences, standard deviation, serial correlation estimate.
#
# Parameters
#  vA is the vector of data under condition A.
#  vB is the vector of data under condition B.
#  rhotype is the type of serial correlation estimate to use. If "WF" (Wayne Fuller's) is not
#    specified, rhow must be inputted and will be used as the serial correlation estimate.
#  rhow is an inputted serial correlation. Only used if rhotype is not specified.
pstrfxn <- function(vA, vB, rhotype, rhow){
  # Calculating vector of differences
  diff <- vB - vA
  
  # Counting length of sequences from 2 "independent" conditions
  n <- length(diff)
  
  # Constructing the variables for the linear model
  one <- rep(1, n)
  x <- 1:n - mean(1:n)
  y <- diff
  
  # Fitting the linear model & obtaining summary stats
  m <- summary( lm(y ~ one + x - 1))
  int  <- m$coefficients[1,1]
  b    <- m$coefficients[2,1]
  sd.db <- m$sigma 
  names(int)  <- NULL
  names(b)    <- NULL
  names(sd.db) <- NULL
  
  # Obtaining residuals & computing correlation
  res <- m$residuals
  if(rhotype=="WF"){
    rho <- WF.r(res)
  } else{
    rho <- rhow
  }
  
  # In case serial correlation fails to compute, set all to NA
  if(is.na(rho)){
    b <- NA
    sd.db <- NA
  }
  
  # Gathering summary statistics
  pstrSum <- cbind(int, b, sd.db, rho)
  colnames(pstrSum) <- c("int","b", "sd.db", paste0("rho.b.", rhotype))
  return(pstrSum)
}


#==========================================================================================================
#----------------------------------------------- Simulation -----------------------------------------------
#==========================================================================================================
# Loop through delta
for(delta in .delta) {
  
  # Loop through paired/unpaired
  for(paired in .paired) {
    if(paired) .rho.a <- ..rho.a[..rho.a!=0] else .rho.a <- 0  # Sets paired correlations accordingly
    
    # Loop through paired correlations
    for(rho.a in .rho.a) {
      
      # Loop through serial correlations
      for(rho.w in .rho.w) {
        
        nmax <- max(.n) # Finds largest inputted n. For the article, this is 100. 
        
        # Exponents for 1st-order serial var-covar matrix
        #   Used in creating "R.matrix" below.
        e1 <- matrix(rep(1:nmax,nmax),nmax,nmax,byrow=TRUE) #-- a temporary step 
        e2 <- matrix(sort(rep(1:nmax,nmax)),nmax,nmax,byrow=TRUE) #-- a temporary step 
        e <- abs(e1-e2)  #-- rho exponents
        
        # x and ey 
        x <- 1:nmax
        ey <- .slope*x   
        x.matrix <- as.matrix(cbind(c(rep(0,nmax),rep(1,nmax)),c(x,x)))
        colnames(x.matrix) <- c('treatment','x')
        
        # Empty data frame to store data
        SUMSTAT[[COMBO]] <- data.frame()
        
        # Loop through M individuals
        for(ITER in 1:.M) {
          
          #----- Generating simulation data -----
          # Set seed
          seed.val <- (201706+ITER)*COMBO
          set.seed(seed.val)
          # Generates matrix based on error structure
          R.matrix <- ((diag(rep(1,nmax))+rho.w)-diag(rep(rho.w,nmax)))^e
          # Var-covar matrix
          tmp.errsigma <- diag(rep(sqrt(.sigsqr),nmax))%*%R.matrix%*%diag(rep(sqrt(.sigsqr),nmax))
          errsigma <- cbind(rbind(tmp.errsigma, rho.a*tmp.errsigma),rbind(rho.a*tmp.errsigma, tmp.errsigma))
          # Subject ID
          subj <- rep(ITER,2*nmax)
          # Generating serially correlated errors
          epsilon <- mvrnorm(1, mu=rep(0,2*nmax), Sigma=errsigma)
          # Data values
          y <- rep(1, 2*length(ey)) + epsilon   # data for level-change tests
          yr <- c(ey, ey) + epsilon             # data for rate-change tests
          
          #----- Estimating summary stats -----
          # Set-up
          ss.unp.levl <- list()   # summary stats from unpaired (2-sample) tests for level change
          ss.unp.rate <- list()   # summary stats from unpaired (2-sample) tests for rate change   
          ss.levl <- list()       # summary stats from paired tests for level change
          ss.rate <- list()       # summary stats from paired tests for rate change

          # Loop through n
          for(n in .n){
            y0 <- y[1:n]
            y1 <- y[(nmax+1):(nmax+n)]
            yr0 <- yr[1:n]
            yr1 <- yr[(nmax+1):(nmax+n)]
            
            if(paired){
              
              # Summary stats for paired tests for level change
              ss.levl[[n]] <- pstfxn(y0, y1, "WF", rho.w)
              
              # Summary stats for paired tests for rate change 
              ss.rate[[n]] <- pstrfxn(yr0, yr1, "WF", rho.w)
              
              # Adding results to SUMSTAT
              if(n==.n[1]){
                dat <- cbind(delta, paired, rho.a, rho.w, ITER,
                             ss.levl[[n]], ss.rate[[n]])
                colnames(dat) <- c('delta', 'paired', 'rho.a', 'rho.w', 'subj',
                                   paste0('diff.pst', n),       # mean of differences in paired tests for level-change
                                   paste0('sd.diff.pst.', n),   # SD of differences in paired tests for level-change
                                   paste0('rho.pst.WF.', n),    # serial correlation in paired tests for level-change   
                                   paste0('intd.pstr.', n),     # intercept of diffs in paired tests for rate-change   
                                   paste0('diff.pstr.', n),     # slope of diffs in paired tests for rate-change
                                   paste0('sd.diff.pstr.', n),  # SD of diffs in paired tests for rate-change
                                   paste0('rho.pstr.WF.', n))   # serial correlation in paired tests for rate-change  
              } else{
                tmp.dat <- cbind(ss.levl[[n]], ss.rate[[n]])
                colnames(tmp.dat) <- c(paste0('diff.pst', n), 
                                       paste0('sd.diff.pst.', n), 
                                       paste0('rho.pst.WF.', n),        
                                       paste0('intd.pstr.', n),        
                                       paste0('diff.pstr.', n), 
                                       paste0('sd.diff.pstr.', n),
                                       paste0('rho.pstr.WF.', n))     
                dat <- cbind(dat, tmp.dat)
              }
            } else{
              
              # Summary stats for unpaired tests for level change
              ss.unp.levl[[n]] <- stfxn(y0, y1, "WF", rho.w)
              
              # Summary stats for unpaired tests for rate change       
              ss.unp.rate[[n]] <- strfxn(yr0, yr1, "WF", rho.w)     
              
              # Adding results to SUMSTAT
              if(n==.n[1]){
                dat <- cbind(delta, paired, rho.a, rho.w, ITER,
                             ss.unp.levl[[n]], ss.unp.rate[[n]])
                colnames(dat) <- c('delta', 'paired', 'rho.a', 'rho.w', 'subj',
                                   paste0('diff.st.', n),     # diff in means for 2-sample tests for level-change
                                   paste0('sd.y.st', n),      # SD of y for 2-sample tests for level-change
                                   paste0('rho.st.WF.', n),   # serial correlation for 2-sample tests for level-change
                                   paste0('intd.str.', n),    # diff in intercepts for 2-sample tests for rate-change  
                                   paste0('diff.str.', n),    # diff in slopes for 2-sample tests for rate-change 
                                   paste0('sd.y.str', n),     # SD of y for 2-sample tests for rate-change
                                   paste0('rho.str.WF.', n))  # serial correlation for 2-sample tests for rate-change     
              } else{
                tmp.dat <- cbind(ss.unp.levl[[n]], ss.unp.rate[[n]])
                colnames(tmp.dat) <- c(paste0('diff.st.', n), paste0('sd.y.st', n),
                                       paste0('rho.st.WF.', n),
                                       paste0('intd.str.', n),  
                                       paste0('diff.str.', n), paste0('sd.y.str', n),
                                       paste0('rho.str.WF.', n))      
                dat <- cbind(dat, tmp.dat)
              }
            }
          }
          
          # Append new data
          SUMSTAT[[COMBO]] <- rbind(SUMSTAT[[COMBO]], dat)
          rownames(SUMSTAT[[COMBO]]) <- NULL
        }
        
        # Write datasets
        fwrite(SUMSTAT[[COMBO]], file=paste("SUMSTAT", COMBO, ".csv", sep=""))  # Writes one dataset for each combo
                                                                                # of paired and rho.w.
        if(paired){  # Writes one dataset for all paired data
          fwrite(SUMSTAT[[COMBO]],
                 file="pairedSUMSTAT.csv", append=TRUE, col.names=!file.exists("pairedSUMSTAT.csv"))
        } else{        # Writes one dataset for all unpaired data
          fwrite(SUMSTAT[[COMBO]],
                 file="unpairedSUMSTAT.csv", append=TRUE, col.names=!file.exists("unpairedSUMSTAT.csv"))
        }
        
        # Increment counter
        COMBO <- COMBO + 1
      }
    }
  }
}

#==========================================================================================================
#**********************************************************************************************************
# N-of-1 t-tests: Datasets combiner for simulation ----
  # Run time: approx 85 minutes
  #           on a Dell Precision M6800 with an 
  #           Intel Core i7-8650U CPU @ 1.90 GHz
  #           with 32.0 GB of RAM
  #           R version 3.5.3
#**********************************************************************************************************
#==========================================================================================================


#==========================================================================================================
#------------------------------------------------- Set-up -------------------------------------------------
#==========================================================================================================
# Read in datasets
dat1 <- fread("pairedSUMSTAT.csv")
dat2 <- fread("unpairedSUMSTAT.csv")

# Parameters
#  .n is the number of observations from a treatment; referred to as "m" in article text. 
#     Can take values {4,5,6,...}, as long as they exist in the generated dataset, though
#     setting it to the same as the simulation is recommended. Default is {4:12,30,50,100}
#  .tests1 indicates the tests in dat1 (the paired serial tests). Takes values {'pst', 'pstr'}.
#     -pst: paired serial t-test for level-change
#     -pstr: paired serial t-test for rate-change
#  .tests2 indicates the tests in dat2 (the 2-sample serial tests). Takes values {'st', 'str'}.
#     -st: 2-sample serial t-test for level-change
#     -str: 2-sample serial t-test for rate-change
.n <- c(4:12,30,50,100)
.tests1 <- c('pst', 'pstr')
.tests2 <- c('st', 'str')

# List to hold rearranged data sets, counter for list
tmp.dat <- list()
counter <- 1

#==========================================================================================================
#--------------------------------------------- Combining sets ---------------------------------------------
#==========================================================================================================
#----------------------------------------- Paired serial t-tests ------------------------------------------
headcolnames <- c('rho.a', 'rho.w', 'subj')
head <- dat1[, ..headcolnames]
# Loop through tests1
for(test in .tests1){
  if(test=='pstr') tstnm <- 'pstr.' else tstnm <- test
  
  # Loop through n
  for(n in .n){
    maincolnames <- c(paste0('diff.',tstnm,n), paste0('sd.diff.',test,'.',n))
    main <- dat1[, ..maincolnames]
    if(test=='pstr'){
      intd <- dat1[[paste0("intd.pstr.",n)]]
    } else{
      intd <- NA # sets intd to NA for level-change
    }
    rhovals <- dat1[[paste0('rho.',test,'.WF.',n)]]
    tmp.dat[[counter]] <- cbind(head, n, test, main, intd, rhovals)
    colnames(tmp.dat[[counter]]) <- c('rho_a','rho_w','subj','N','TEST','DIFF','SD','INTD','wf_RHO')
    counter <- counter+1
  }
}

#---------------------------------------- 2-sample serial t-tests -----------------------------------------
headcolnames <- c('rho.a', 'rho.w', 'subj')
head <- dat2[, ..headcolnames]
# Loop through tests2
for(test in .tests2){
  
  # Loop through n
  for(n in .n){
    maincolnames <- c(paste0('diff.',test,'.',n), paste0('sd.y.',test,n))
    main <- dat2[, ..maincolnames]
    if(test=='str'){
      intd <- dat2[[paste0("intd.str.",n)]]
    } else{
      intd <- NA # sets intd to NA for level-change
    }
    rhovals <- dat2[[paste0('rho.',test,'.WF.',n)]]
    tmp.dat[[counter]] <- cbind(head, n, test, main, intd, rhovals)
    colnames(tmp.dat[[counter]]) <- c('rho_a','rho_w','subj','N','TEST','DIFF','SD','INTD','wf_RHO')
    counter <- counter+1
  }
}

# Combined rearranged datasets
allSUMSTAT <- rbindlist(tmp.dat)

# Writes rearranged datasets (allSUMSTAT)
fwrite(allSUMSTAT, "allSUMSTAT.csv")

#==========================================================================================================
#**********************************************************************************************************
# N-of-1 t-tests: Type I errors for simulation ----
# Run time: less than 1 minute
#           on a Dell Precision M6800 with an 
#           Intel Core i7-8650U CPU @ 1.90 GHz
#           with 32.0 GB of RAM
#           R version 3.5.3
#**********************************************************************************************************
#==========================================================================================================


#==========================================================================================================
#------------------------------------------------- Set-up -------------------------------------------------
#==========================================================================================================
# Read in dataset
dat <- fread("allSUMSTAT.csv")
colnames(dat) <- c("RHO_A", "RHO_W", "SUBJ", "N", "TEST", "DIFF", "SD", "INTD", "WF_RHO")
dat.tab <- as.data.table(dat)
setkey(dat.tab, TEST, N, RHO_W)

# Parameters
#  alpha is the nominal Type I error level.
#  .rho is the serial correlation. Can take values between (-1,1), as long as they exist in
#     the generated dataset. Setting it to the same as the simulation is recommended.
#     Default is {-.33, 0, .33, .67}.
#  .tests is the four serial t-tests. Can take values {'st', 'str', 'pst', 'pstr'}, as long
#     as they exist in the generated dataset. Default has all 4.
#       -st: 2-sample serial t-test for level-change
#       -str: 2-sample serial t-test for rate-change
#       -pst: paired serial t-test for level-change
#       -pstr: paired serial t-test for rate-change
#  .nrate is the number of observations from a treatment for rate-change tests;
#     referred to as "m" in article text. Can take values {5,6,7,...}, as long as they
#     exist in the generated dataset. Setting it to the same as the simulation
#     (with the exception that the minimum is 5) is recommended. Default is {5:12,30,50,100}.
#  .nlevel is the number of observations from a treatment for level-change tests;
#     referred to as "m" in article text. Can take values {4,5,6,...}, as long
#     as they exist in the generated dataset. Setting it to the same as the simulation is
#     recommended. Default is {4:12,30,50,100}.
alpha <- .05
.rho <- c(-.33, 0, 0.33, 0.67)
.tests <- c("st", "str", "pst", "pstr")
.nrate <- c(5:12, 30, 50, 100)
.nlevel <- c(4:12, 30, 50, 100)

# Dataset to hold Type I errors
ERRORS <- data.table()

#==========================================================================================================
#----------------------------------------------- Functions ------------------------------------------------
#==========================================================================================================
pst.se <- function(N, SD, RHO){
  level.c <- (N + 2*RHO^(N+1) - N*RHO^2 - 2*RHO) / (N^2*(RHO-1)^2)
  level.b <- ( N - (N+2*RHO^(N+1)-N*RHO^2-2*RHO)/(N*(RHO-1)^2) ) / (N-1)
  level.nprime <- 1 / (1 - (N-1)*level.b/N)
  pst.se <- SD * sqrt( level.c/level.b )
  cbnse <- cbind(c=level.c, b=level.b, nprime=level.nprime, se=pst.se) 
  return(cbnse)
}

pstr.se <- function(N, SD, RHO){
  rate.c <- 12/(N^2-1)^2 * ( - 6*RHO*(RHO+1)^2*(RHO^N-1)/(N^2*(RHO-1)^4) +
                               2*RHO*(6*RHO^(N+1)+6*RHO^N+RHO^2-2*RHO+1)/(N*(RHO-1)^3) -
                               6*RHO*(RHO^N+1)/(RHO-1)^2 - 2*N*RHO/(RHO-1) + (N^2-1)/N )
  rate.b <- (N - (2 + 2*RHO*(RHO^N-N*RHO+N-1)/(N*(RHO-1)^2) +
                    (-2*RHO/(N*(N^2-1)*(RHO-1)^4) *
                       (N^3*(RHO-1)^3 + 3*N^2*(RHO-1)^2*(RHO^N+1) +
                          3*(RHO+1)^2*(RHO^N-1) -
                          N*(RHO-1)*(6*RHO^(N+1)+6*RHO^N+RHO^2-2*RHO+1))))) / (N-2)
  rate.nprime <- 2 / (1-rate.b*(N-2)/N)
  pstr.se <- SD * sqrt( rate.c/rate.b )
  cbnse <- cbind(c=rate.c, b=rate.b, nprime=rate.nprime, se=pstr.se)
  return(cbnse)
}

st.se <- function(N, SD, RHO){
  level.c <- (N + 2*RHO^(N+1) - N*RHO^2 - 2*RHO) / (N^2*(RHO-1)^2)
  level.b <- ( N - (N+2*RHO^(N+1)-N*RHO^2-2*RHO)/(N*(RHO-1)^2) ) / (N-1)
  level.nprime <- 2 / (1 - (N-1)*level.b/N)
  st.se <- SD * sqrt( 2*level.c/level.b )
  cbnse <- cbind(c=level.c, b=level.b, nprime=level.nprime, se=st.se)
  return(cbnse)
}

str.se <- function(N, SD, RHO){
  rate.c <- 12/(N^2-1)^2 * ( - 6*RHO*(RHO+1)^2*(RHO^N-1)/(N^2*(RHO-1)^4) +
                               2*RHO*(6*RHO^(N+1)+6*RHO^N+RHO^2-2*RHO+1)/(N*(RHO-1)^3) -
                               6*RHO*(RHO^N+1)/(RHO-1)^2 - 2*N*RHO/(RHO-1) + (N^2-1)/N )
  rate.b <- (N - (2 + 2*RHO*(RHO^N-N*RHO+N-1)/(N*(RHO-1)^2) +
                    (-2*RHO/(N*(N^2-1)*(RHO-1)^4) *
                       (N^3*(RHO-1)^3 + 3*N^2*(RHO-1)^2*(RHO^N+1) +
                          3*(RHO+1)^2*(RHO^N-1) -
                          N*(RHO-1)*(6*RHO^(N+1)+6*RHO^N+RHO^2-2*RHO+1))))) / (N-2)
  rate.nprime <- 4 / (1-rate.b*(N-2)/N)
  str.se <- SD * sqrt( 2*rate.c/rate.b )
  
  cbnse <- cbind(c=rate.c, b=rate.b, nprime=rate.nprime, se=str.se)
  return(cbnse)
}

#==========================================================================================================
#------------------------------- Type I error & Margin of Error calculation -------------------------------
#==========================================================================================================
# Loop through tests
for(test in .tests){
  
  # Sets correct .n vector based on whether test is rate- or level-change
  if(test=="pstr" || test=="str"){
    .n <- .nrate
  } else{
    .n <- .nlevel
  }
  
  # Loop through n
  for(n in .n){
    
    # Loop through serial correlation values
    for(rho in .rho){
      
      # Subsets data
      tmp.dat1 <- dat.tab[.(test, n, rho)]
      
      sd <- tmp.dat1$SD  # vector of SD's
      est.rho <- tmp.dat1$WF_RHO  # vector of estimated serial correlations
      
      # Calculates standard errors (se) and degrees of freedom (df) for serial and
      # analogous usual (usu) t-tests
      if(test=="st") {
        se.parts <- st.se(n, sd, est.rho)
        se <- unname(se.parts[,"se"])
        b <- unname(se.parts[,"b"])
        c <- unname(se.parts[,"c"])  
        n.eff <- unname(se.parts[,"nprime"])        # st.se multiplies n.eff by 2, assuming both
                                                    # samples are the same size
        df <- n.eff - 2
        
        usu.test <- "t"
        usu.c <- sqrt( 2/n )
        usu.se <- sd * usu.c
        usu.df <- 2*n-2
        
        true.se.parts <- st.se(n, 1, rho); 
        true.n.eff <- unname(true.se.parts[,"nprime"])
        true.se <- unname(true.se.parts[,"se"])
        true.MoE <- qt(alpha, df=true.n.eff-2, lower.tail=FALSE)*true.se  #MoE for (1-2*alpha)100% CI
      } else if(test=="str") {
        se.parts <- str.se(n, sd, est.rho)
        se <- unname(se.parts[,"se"])
        b <- unname(se.parts[,"b"])
        c <- unname(se.parts[,"c"])  
        n.eff <- unname(se.parts[,"nprime"])         # str.se multiplies n.eff by 2, assuming both
                                                     # samples are the same size
        df <- n.eff - 4
        
        usu.test <- "tr"
        usu.c <- unname(str.se(n, sd, 0)[,"c"])
        usu.se <- unname(str.se(n, sd, 0)[,"se"])
        usu.df <- 2*n-4
        
        true.se.parts <- str.se(n, 1, rho)
        true.n.eff <- unname(true.se.parts[,"nprime"])
        true.se <- unname(true.se.parts[,"se"])
        true.MoE <- qt(alpha, df=true.n.eff-4, lower.tail=FALSE)*true.se  #MoE for (1-2*alpha)100% CI
      } else if(test=="pst") {
        se.parts <- pst.se(n, sd, est.rho)
        se <- unname(se.parts[,"se"])
        b <- unname(se.parts[,"b"])
        c <- unname(se.parts[,"c"])  
        n.eff <- unname(se.parts[,"nprime"])
        df <- n.eff - 1
          
        usu.test <- "pt"
        usu.c <- sqrt( 1/n )
        usu.se <- sd*usu.c
        usu.df <- n-1
        
        true.se.parts <- pst.se(n, 1, rho)    
        true.n.eff <- unname(true.se.parts[,"nprime"])
        true.se <- unname(true.se.parts[,"se"])
        true.MoE <- qt(alpha, df=true.n.eff-1, lower.tail=FALSE)*true.se  #MoE for (1-2*alpha)100% CI
        
      } else if(test=="pstr") {
        se.parts <- pstr.se(n, sd, est.rho)
        se <- unname(se.parts[,"se"])
        b <- unname(se.parts[,"b"])
        c <- unname(se.parts[,"c"])  
        n.eff <- unname(se.parts[,"nprime"])
        df <- n.eff - 2
        
        usu.test <- "ptr"
        usu.c <- unname(pstr.se(n, sd, 0)[,"c"])
        usu.se <- unname(pstr.se(n, sd, 0)[,"se"])
        usu.df <- n-2
        
        true.se.parts <- pstr.se(n, 1, rho)
        true.n.eff <- unname(true.se.parts[,"nprime"])
        true.se <- unname(true.se.parts[,"se"])
        true.MoE <- qt(alpha, df=true.n.eff-2, lower.tail=FALSE)*true.se  #MoE for (1-2*alpha)100% CI
      } else {
        stop("test not found!")
      }
      
      t.stat <- tmp.dat1$DIFF / se   # t-statistic for serial
      prob.t <- pt(t.stat, df=df, lower.tail=FALSE)   # one-sided p-value of t-statistic for serial
      typeIerror <- mean(prob.t<=alpha)   # actual Type I error for serial
      ser.MoE <- mean( qt(alpha, df=df, lower.tail=FALSE)*se )
      
      usu.t.stat <- tmp.dat1$DIFF / usu.se   # t-statistic for usual
      usu.prob.t <- pt(usu.t.stat, df=usu.df, lower.tail=FALSE)   # one-sided p-value of t-statistic for usual
      usu.typeIerror <- mean(usu.prob.t<=alpha)   # actual Type I error for usual
      usu.MoE <- mean( qt(alpha, df=usu.df, lower.tail=FALSE)*usu.se )
      
      trueMoE <- mean(true.MoE*sd)     # expected Margin of Error
      ser.tru.MoE <- ser.MoE/trueMoE   # ratio serial observed MoE to expected MoE
      usu.tru.MoE <- usu.MoE/trueMoE   # ratio serial observed MoE to expected MoE
      
      # Append calculated Type I errors to dataset
      tmp.dat <- cbind(alpha, n, test, rho, typeIerror, usu.typeIerror, ser.MoE, usu.MoE, 
                       trueMoE, ser.tru.MoE, usu.tru.MoE)
      colnames(tmp.dat) <- c("nominaltypeI", "n", "test", "rho", "serialtypeI", "usualtypeI", 
                             "serialMoE","usualMoE", "trueMoE", "ser.tru.MoE", "usu.tru.MoE")
      ERRORS <- rbind(ERRORS, tmp.dat)
      
    }
  }
}

# Writes Type I errors dataset
fwrite(ERRORS, "TYPEIERRORS for allSUMSTAT.csv")

#----- Table 3 Margin of Error results ----
tmp1 <- subset(ERRORS, n==8 | n==12 | n==100)
tmp2 <- subset(ERRORS, (test=='st' | test=='pst') & n==4)
tmp3 <- subset(ERRORS, (test=='str' | test=='pstr') & n==5)
tmp4 <- rbind(tmp1,tmp2, tmp3)
tmp5 <- data.frame()

for(iter1 in 1:4){
  for(iter2 in 1:4){
    
  }}
fwrite(tmp4, "Margin of Errors from allSUMSTAT.csv")
#==========================================================================================================
#**********************************************************************************************************
# N-of-1 t-tests: Effect sizes for simulation ----
# Run time: approx 4 minutes
#           on a Dell Precision M6800 with an 
#           Intel Core i7-8650U CPU @ 1.90 GHz
#           with 32.0 GB of RAM
#           R version 3.5.3
#**********************************************************************************************************
#==========================================================================================================

#==========================================================================================================
#------------------------------------------------- Set-up -------------------------------------------------
#==========================================================================================================
# Read in dataset
dat <- fread("allSUMSTAT.csv")
colnames(dat) <- c("RHO_A", "RHO_W", "SUBJ", "N", "TEST", "DIFF", "SD", "INTD", "WF_RHO")
dat.tab <- as.data.table(dat)
setkey(dat.tab, TEST, N, RHO_W)

# Parameters
#  alpha is the nominal Type I error level.
#  power is the value for power.
#  .dmax is the maximum delta to be checked, where delta is the true difference
#     between the two treatments.
#  .rho is the serial correlation. Can take values between (-1,1), as long as they exist in
#     the generated dataset. Setting it to the same as the simulation is recommended.
#     Default is {-.33, 0, .33, .67}.
#  .tests is the four serial t-tests. Can take values {'st', 'str', 'pst', 'pstr'}, as long
#     as they exist in the generated dataset. Default has all 4.
#       -st: 2-sample serial t-test for level-change
#       -str: 2-sample serial t-test for rate-change
#       -pst: paired serial t-test for level-change
#       -pstr: paired serial t-test for rate-change
#  .types is the type of calculated effect size. Can take values {'U', 'S', 'TH'}, which
#     stands for usual, serial, and theoretical, respectively. Default has all 3.
#  .sigfig is the number of decimals against which power is checked.
#     i.e. if .sigfig=3 and power=.8, the power of returned effect sizes will
#     be .800+/-5e-4. Default is 3.
#  .nrate is the number of observations from a treatment for rate-change tests;
#     referred to as "m" in article text. Can take values {5,6,7,...}, as long as they
#     exist in the generated dataset. Setting it to the same as the simulation
#     (with the exception that the minimum is 5) is recommended. Default is {5:12,30,50,100}.
#  .nlevel is the number of observations from a treatment for level-change tests;
#     referred to as "m" in article text. Can take values {4,5,6,...}, as long as they 
#     exist in the generated dataset. Setting it to the same as the simulation is
#     recommended. Default is {4:12,30,50,100}.
alpha <- .05
power <- .8
.dmax <- 130
.types <- c("U", "S", "TH")
.sigfig <- 3
.rho <- c(-.33, 0, 0.33, 0.67)
.tests <- c("st", "str", "pst", "pstr")
.nrate <- c(5:12, 30, 50, 100)
.nlevel <- c(4:12, 30, 50, 100)

# Dataset to hold effect sizes
EFFECTSIZES <- data.table()

#==========================================================================================================
#----------------------------------------------- Functions ------------------------------------------------
#==========================================================================================================
#------------------- Functions computing standard errors of 4 tests: st, str, pst, pstr -------------------
pst.se <- function(N, SD, RHO){
  level.c <- (N + 2*RHO^(N+1) - N*RHO^2 - 2*RHO) / (N^2*(RHO-1)^2)
  level.b <- ( N - (N+2*RHO^(N+1)-N*RHO^2-2*RHO)/(N*(RHO-1)^2) ) / (N-1)
  level.nprime <- 1 / (1 - (N-1)*level.b/N)
  pst.se <- SD * sqrt( level.c/level.b )
  cbnse <- cbind(level.c=level.c, level.b=level.b, level.nprime=level.nprime, pst.se=pst.se)
  return(cbnse)
}

pstr.se <- function(N, SD, RHO){
  rate.c <- 12/(N^2-1)^2 * ( - 6*RHO*(RHO+1)^2*(RHO^N-1)/(N^2*(RHO-1)^4) +
                               2*RHO*(6*RHO^(N+1)+6*RHO^N+RHO^2-2*RHO+1)/(N*(RHO-1)^3) -
                               6*RHO*(RHO^N+1)/(RHO-1)^2 - 2*N*RHO/(RHO-1) + (N^2-1)/N )
  rate.b <- (N - (2 + 2*RHO*(RHO^N-N*RHO+N-1)/(N*(RHO-1)^2) +
                    (-2*RHO/(N*(N^2-1)*(RHO-1)^4) *
                       (N^3*(RHO-1)^3 + 3*N^2*(RHO-1)^2*(RHO^N+1) +
                          3*(RHO+1)^2*(RHO^N-1) -
                          N*(RHO-1)*(6*RHO^(N+1)+6*RHO^N+RHO^2-2*RHO+1))))) / (N-2)
  rate.nprime <- 2 / (1-rate.b*(N-2)/N)
  pstr.se <- SD * sqrt( rate.c/rate.b )
  cbnse <- cbind(rate.c=rate.c, rate.b=rate.b, rate.nprime=rate.nprime, pstr.se=pstr.se)
  return(cbnse)
}

st.se <- function(N, SD, RHO){
  level.c <- (N + 2*RHO^(N+1) - N*RHO^2 - 2*RHO) / (N^2*(RHO-1)^2)
  level.b <- ( N - (N+2*RHO^(N+1)-N*RHO^2-2*RHO)/(N*(RHO-1)^2) ) / (N-1)
  level.nprime <- 2 / (1 - (N-1)*level.b/N)
  st.se <- SD * sqrt( 2*level.c/level.b )
  cbnse <- cbind(level.c=level.c, level.b=level.b, level.nprime=level.nprime, st.se=st.se)
  return(cbnse)
}

str.se <- function(N, SD, RHO){
  rate.c <- 12/(N^2-1)^2 * ( - 6*RHO*(RHO+1)^2*(RHO^N-1)/(N^2*(RHO-1)^4) +
                               2*RHO*(6*RHO^(N+1)+6*RHO^N+RHO^2-2*RHO+1)/(N*(RHO-1)^3) -
                               6*RHO*(RHO^N+1)/(RHO-1)^2 - 2*N*RHO/(RHO-1) + (N^2-1)/N )
  rate.b <- (N - (2 + 2*RHO*(RHO^N-N*RHO+N-1)/(N*(RHO-1)^2) +
                    (-2*RHO/(N*(N^2-1)*(RHO-1)^4) *
                       (N^3*(RHO-1)^3 + 3*N^2*(RHO-1)^2*(RHO^N+1) +
                          3*(RHO+1)^2*(RHO^N-1) -
                          N*(RHO-1)*(6*RHO^(N+1)+6*RHO^N+RHO^2-2*RHO+1))))) / (N-2)
  rate.nprime <- 4 / (1-rate.b*(N-2)/N)
  str.se <- SD * sqrt( 2*rate.c/rate.b )
  
  cbnse <- cbind(rate.c=rate.c, rate.b=rate.b, rate.nprime=rate.nprime, str.se=str.se)
  return(cbnse)
}
#--------------------------- Function computing P(T >= t.crit | delta) ------------------------------------
est.power <- function(DELTA, DF, DIFF.VEC, SE.VEC, SIGNIF){
  t.stat <- (DIFF.VEC + DELTA) / SE.VEC
  prob.t <- pt(abs(t.stat), df=DF, lower.tail=FALSE)
  pow <- mean(prob.t<=SIGNIF)
  return(power = pow)
}


#==========================================================================================================
#---------------------------------------- Effect size calculation -----------------------------------------   
#==========================================================================================================
# Loop through types
for(type in .types){
  
  # Loop through tests
  for(test in .tests){
    
    # Sets correct .n vector based on whether test is rate- or level-change
    if(test=="pstr" || test=="str"){
      .n <- .nrate
    } else{
      .n <- .nlevel
    }
    
    # Loop through n
    for (n in .n){
      
      # Loop through serial correlation values
      for (iter in 1:length(.rho)){
        rho <- .rho[iter]
        
        # Subsets data
        tmp.dat1 <- dat.tab[.(test, n, rho)]
        
        sd <- tmp.dat1$SD  # vector of SD's
        est.rho <- tmp.dat1$WF_RHO  # vector of estimated serial correlations
        
        # Sets correct degrees of freedom for serial t-tests
        if(test=="st"){
          se.parts.truerho <- st.se(n, 0, rho)  # sd isn't used
          truerho.df <- unname(se.parts.truerho[,"level.nprime"]) - 2
          truerho.c <- 2 * unname(se.parts.truerho[,"level.c"])
          
          se.parts <- st.se(n, sd, est.rho)
          df <- unname(se.parts[,"level.nprime"]) - 2
          se <- unname(se.parts[,"st.se"])
        } else if(test=="str"){
          se.parts.truerho <- str.se(n, 0, rho)  # sd isn't used
          truerho.df <- unname(se.parts.truerho[,"rate.nprime"]) - 4
          truerho.c <- 2 * unname(se.parts.truerho[,"rate.c"])
          
          se.parts <- str.se(n, sd, est.rho)
          df <- unname(se.parts[,"rate.nprime"]) - 4
          se <- unname(se.parts[,"str.se"])
        } else if(test=="pst"){
          se.parts.truerho <- pst.se(n, 0, rho)  # sd isn't used
          truerho.df <- unname(se.parts.truerho[,"level.nprime"]) - 1
          truerho.c <- unname(se.parts.truerho[,"level.c"])
          
          se.parts <- pst.se(n, sd, est.rho)
          df <- unname(se.parts[,"level.nprime"]) - 1
          se <- unname(se.parts[,"pst.se"])
        } else if(test=="pstr"){
          se.parts.truerho <- pstr.se(n, 0, rho)  # sd isn't used
          truerho.df <- unname(se.parts.truerho[,"rate.nprime"]) - 2
          truerho.c <- unname(se.parts.truerho[,"rate.c"])
          
          se.parts <- pstr.se(n, sd, est.rho)
          df <- unname(se.parts[,"rate.nprime"]) - 2
          se <- unname(se.parts[,"pstr.se"])
        } else{
          stop("test not found!")
        }
        
        #---------------------------------- Serial - Theoretical power ----------------------------------
        if(type=="TH"){
          
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
                                    lower.tail=FALSE), digits=.sigfig)
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
                                     lower.tail=FALSE), digits=.sigfig)
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
          
          #----------------------------------- Serial - Estimated power -----------------------------------
        } else if(type=="S"){
          
          # Determine effect size to achieve given power by searching from 0:.dmax
          notFound <- TRUE  # changes to false when correct effect size for given power is found
          upper <- .dmax
          lower <- 0
          increment <- 1
          while(notFound){
            upperNotFound <- TRUE  # changes to false when upper bound for effect size is found
            while(upperNotFound){
              estpow.up <- round(est.power(upper, df, tmp.dat1$DIFF, se, SIGNIF=alpha), digits=.sigfig)
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
              estpow.low <- round(est.power(lower, df, tmp.dat1$DIFF, se, SIGNIF=alpha), digits=.sigfig)
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
              increment <- increment / 10  # If effect size not found, check
              # smaller range with greater precision
            }
          }
          
          #----------------------------------- Usual - Estimated power ------------------------------------
        } else if(type=="U"){
          
          if(test=="st") {
            test2 <- "t"
            c <- sqrt( 2/n )
            usu.se <- sd*c
            usu.df <- 2*n-2
          } else if(test=="str") {
            test2 <- "tr"
            usu.se <- unname(str.se(n, sd, 0)[,"str.se"])
            usu.df <- 2*n-4
          } else if(test=="pst") {
            test2 <- "pt"
            c <- sqrt( 1/n )
            usu.se <- sd*c
            usu.df <- n-1
          } else if(test=="pstr") {
            test2 <- "ptr"
            usu.se <- unname(pstr.se(n, sd, 0)[,"pstr.se"])
            usu.df <- n-2
          } else{
            stop("Sorry, test not found")
          }
          
          # Determine effect size to achieve given power by searching from 0:.dmax
          notFound <- TRUE
          upper <- .dmax
          lower <- 0
          increment <- 1
          while(notFound){
            upperNotFound <- TRUE  # changes to false when upper bound for effect size is found
            while(upperNotFound){
              estpow.up <- round(est.power(upper, usu.df, tmp.dat1$DIFF,
                                           usu.se, SIGNIF=alpha), digits=.sigfig)
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
              estpow.low <- round(est.power(lower, usu.df, tmp.dat1$DIFF, 
                                            usu.se, SIGNIF=alpha), digits=.sigfig)
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
              increment <- increment / 10  # If effect size not found, check
              # smaller range with greater precision
            }
          }
          
        }
        
        # Append calculated effect sizes to dataset
        tmp <- cbind(type, n, test, rho, effectsize)
        EFFECTSIZES <- rbind(EFFECTSIZES, tmp)
      }
    }
  }
}

# Writes effect sizes dataset
fwrite(EFFECTSIZES, "EFFECTSIZES for allSUMSTAT.csv")

#==========================================================================================================
#**********************************************************************************************************
# N-of-1 t-tests: Effect sizes reformat ----
# Run time: less than 1 minute
#           on a Dell Precision M6800 with an 
#           Intel Core i7-8650U CPU @ 1.90 GHz
#           with 32.0 GB of RAM
#           R version 3.5.3
#**********************************************************************************************************
#==========================================================================================================


#==========================================================================================================
#------------------------------------------------- Set-up -------------------------------------------------
#==========================================================================================================
# Read in dataset
dat <- fread("EFFECTSIZES for allSUMSTAT.csv")

# Parameters
#  .tests is the four serial t-tests. Can take values {'st', 'str', 'pst', 'pstr'}, as long
#     as they exist in the generated dataset. Default has all 4.
#       -st: 2-sample serial t-test for level-change
#       -str: 2-sample serial t-test for rate-change
#       -pst: paired serial t-test for level-change
#       -pstr: paired serial t-test for rate-change
.tests <- c("st", "str", "pst", "pstr")

# Dataset to hold reformatted effect sizes
EFFECTSIZESreformat <- data.table()


#==========================================================================================================
#---------------------------------------- Effect size reformatting ----------------------------------------
#==========================================================================================================
# Loop through tests
for(TEST in .tests){
  # Effect sizes for usual t-tests
  du <- subset(dat, test==TEST & type=="U", select=c(n, test, rho, effectsize))
  colnames(du) <- c("n", "test", "rho", "du")
  
  # Effect sizes for serial t-tests
  ds <- subset(dat, test==TEST & type=="S", select=c(n, test, rho, effectsize))
  colnames(ds) <- c("n", "test", "rho", "ds")
  
  # Theoretical effect sizes
  dth <- subset(dat, test==TEST & type=="TH", select=c(n, test, rho, effectsize))
  colnames(dth) <- c("n", "test", "rho", "dth")
  
  # Merge and append
  tmp.dat1 <- merge(merge(du, ds), dth)
  EFFECTSIZESreformat <- rbind(EFFECTSIZESreformat, tmp.dat1)
}

# Writes reformatted effect sizes dataset
fwrite(EFFECTSIZESreformat, "EFFECTSIZES-reformat for allSUMSTAT.csv")

#==========================================================================================================
#**********************************************************************************************************
# N-of-1 t-tests: Figures 3-6 ----
# Run time: less than 1 minute
#           on a Dell OptiPlex 990 with an 
#           Intel Core i7-2600 CPU @ 3.40 GHz
#           with 8.0 GB of RAM
#           R version 3.4.1
#
# Note: Please make sure files are available from Effect sizes reformat and Type I errors
#**********************************************************************************************************
#==========================================================================================================


#==========================================================================================================
#------------------------------------------------- Set-up -------------------------------------------------
#==========================================================================================================

# Type I Errors
typeIdat <- fread("TYPEIERRORS for allSUMSTAT.csv")
typeIdat.tab <- as.data.table(typeIdat)
setkey(typeIdat.tab, test, rho, n)
n.original <- c(4:12, 30, 50, 100)
n.new <- c(4:15)  # for plotting purposes
typeIdat.tab$n.new <- n.new[match(typeIdat.tab$n, n.original)]

# Effect sizes
ESdat <- fread("EFFECTSIZES-reformat for allSUMSTAT.csv")
ESdat.tab <- as.data.table(ESdat)
setkey(ESdat.tab, test, rho, n)
n.original <- c(4:12, 30, 50, 100)
n.new <- c(4:15)  # for plotting purposes
ESdat.tab$n.new <- n.new[match(ESdat.tab$n, n.original)]

# Change to correct folder
# setwd("<insert location of datasets outputted from simulation>")
setwd("C:/Temp for simulation data")

# Parameters
#  .rho is the serial correlation. Can take values between (-1,1), as long as they exist in
#     the generated dataset. Setting it to the same as the simulation is recommended.
#     Default is {-.33, 0, .33, .67}.
#  .color is the colors for each value of serial correlation.
#     Default is {"orangered", "orange", "orange3", "orange4"}. As rho increases, color will
#     become less red and more brown.
.rho <- c(-.33, 0, .33, .67)
.color <- c("orangered", "orange", "orange3", "orange4")


#==========================================================================================================
#------------------------------------------------ Figures -------------------------------------------------
#==========================================================================================================

#--------------------------------------------- st - Figure 4 ----------------------------------------------
# ***Set up ----
tiff("N-of-1 t-tests Figure 4.tiff", bg="white", width=3200, height=3200, compression="lzw")
layout(matrix(c(1,2,3), 3, 1), heights=c(1,.7,1))
resize <- 4.5
par(cex=resize, cex.axis=1.25, cex.main=1.25, cex.sub=1.25)

# ***First panel - type I errors ----
for(ITER in 1:length(.rho)){
  RHO <- .rho[ITER]
  COLOR <- .color[ITER]
  
  dat.tmp <- typeIdat.tab[.("st", RHO)]
  
  if(ITER==1){
    par(mar=c(4.1, 6.1, .6, 2.1))
    plot(serialtypeI ~ n.new, data=dat.tmp, type="l", col=COLOR, lwd=resize*3,
         xlab="", ylab="", yaxt="n", xaxt="n", xlim=c(4,15), ylim=c(0,.35), frame.plot=FALSE)
    axis(side=1, at=seq(4,15), lwd=0, lwd.ticks=resize, labels=2*c(4:12, 30, 50, 100))
    axis(side=2, at=seq(0,.35,.05), las=2, lwd=resize)
    mtext("Type I Error", side=2, line=4, cex=resize*1.5)
    mtext(expression(italic(m[A]+m[B])), side=1, line=2.5, cex=resize*1.5)
    leg0 <- legend("topright", legend=c("Serial", "Usual"), lty=1:2, lwd=resize*3, bty="n", seg.len=4)
    rect(xleft=leg0$rect$left, ybottom=leg0$rect$top-leg0$rect$h,
         xright=leg0$rect$left+leg0$rect$w, ytop=leg0$rect$top, lwd=resize)
    leg <- legend(x=leg0$rect$left, y=leg0$rect$top, legend=.rho, col=.color,
                  lwd=resize*3, ncol=2, bty="n", xjust=1)
    text(labels=expression(paste("True ", rho, " =")), x=leg$rect$left, y=leg$rect$top-.5*leg$rect$h,
         adj=c(1,0.5), cex=1.25)
    par(cex=resize*1.5)
    txtwidth <- strwidth(expression(paste("True ", rho, " =")))
    par(cex=resize)
    rect(xleft=leg$rect$left-txtwidth, ybottom=leg$rect$top-leg$rect$h,
         xright=leg$rect$left+leg$rect$w, ytop=leg$rect$top, lwd=resize)
    abline(h=.05, col="red", lwd=resize)
    lines(usualtypeI ~ n.new, data=dat.tmp, lwd=resize*3, lty=2, col=COLOR) 
    par(xpd=NA)
    lines(y = c(par("usr")[3],par("usr")[4]), x = c(par("usr")[1],par("usr")[1]), lwd=resize) # Left y axis
    lines(y = c(par("usr")[3],par("usr")[4]), x = c(par("usr")[2],par("usr")[2]), lwd=resize) # Right y axis
    lines(y = c(par("usr")[4],par("usr")[4]), x = c(par("usr")[1],par("usr")[2]), lwd=resize) # Top x axis
    lines(y = c(par("usr")[3],par("usr")[3]), x = c(par("usr")[1],12.25), lwd=resize) # Bottom x axis 1
    lines(y = c(par("usr")[3],par("usr")[3]), x = c(12.75,par("usr")[2]), lwd=resize) # Bottom x axis 2
    text("l", x = 12.25, y = par("usr")[3], srt = -45, font = 2)
    text("l", x = 12.75, y = par("usr")[3], srt = -45, font = 2)
    par(xpd=FALSE)
  }
  else{
    lines(serialtypeI ~ n.new, data=dat.tmp, col=COLOR, lwd=resize*3) 
    lines(usualtypeI ~ n.new, data=dat.tmp, col=COLOR, lwd=resize*3, lty=2) 
  }
}

# ***Second panel - theoretical effect sizes ----
for(ITER in 1:length(.rho)){
  RHO <- .rho[ITER]
  COLOR <- .color[ITER]
  
  dat.tmp <- ESdat.tab[.("st", RHO)] 
  
  if(ITER==1){
    par(mar=c(4.1, 6.1, .1, 2.1))
    plot(dth ~ n.new, data=dat.tmp, type="l", col=COLOR, lwd=resize*3,
         xlab="", ylab="", yaxt="n", xaxt="n", xlim=c(4,15), ylim=c(0,7), frame.plot=FALSE)
    axis(side=1, at=seq(4,15), lwd=0, lwd.ticks=resize, labels=2*c(4:12, 30, 50, 100))
    axis(side=2, at=seq(0,7,1), las=2, lwd=resize)
    mtext(expression(paste("Theoretical ", delta)), side=2, line=4, cex=resize*1.5)
    mtext(expression(italic(m[A]+m[B])), side=1, line=2.5, cex=resize*1.5)
    leg <- legend("topright", legend=.rho, col=.color, lwd=resize*3, ncol=2, bty="n")
    text(labels=expression(paste("True ", rho, " =")), x=leg$rect$left, y=leg$rect$top-.5*leg$rect$h,
         adj=c(1,0.5), cex=1.25)
    par(cex=resize*1.5)
    txtwidth <- strwidth(expression(paste("True ", rho, " =")))
    par(cex=resize)
    rect(xleft=leg$rect$left-txtwidth, ybottom=leg$rect$top-leg$rect$h,
         xright=leg$rect$left+leg$rect$w, ytop=leg$rect$top, lwd=resize)
    par(xpd=NA)
    lines(y = c(par("usr")[3],par("usr")[4]), x = c(par("usr")[1],par("usr")[1]), lwd=resize) # Left y axis
    lines(y = c(par("usr")[3],par("usr")[4]), x = c(par("usr")[2],par("usr")[2]), lwd=resize) # Right y axis
    lines(y = c(par("usr")[4],par("usr")[4]), x = c(par("usr")[1],par("usr")[2]), lwd=resize) # Top x axis
    lines(y = c(par("usr")[3],par("usr")[3]), x = c(par("usr")[1],12.25), lwd=resize) # Bottom x axis 1
    lines(y = c(par("usr")[3],par("usr")[3]), x = c(12.75,par("usr")[2]), lwd=resize) # Bottom x axis 2
    text("l", x = 12.25, y = par("usr")[3], srt = -45, font = 2)
    text("l", x = 12.75, y = par("usr")[3], srt = -45, font = 2)
    par(xpd=FALSE)
  }
  else{
    lines(dth ~ n.new, data=dat.tmp, col=COLOR, lwd=resize*3)
  }
}

# ***Third panel - effect size ratios ----
for(ITER in 1:length(.rho)){
  RHO <- .rho[ITER]
  COLOR <- .color[ITER]
  
  dat.tmp <- ESdat.tab[.("st", RHO)]
  
  ratio.s <- dat.tmp$ds / dat.tmp$dth
  ratio.u <- dat.tmp$du / dat.tmp$dth
  
  if(ITER==1){
    par(mar=c(3.6, 6.1, .1, 2.1))
    plot(ratio.s ~ dat.tmp$n.new, type="l", lty=1, lwd=resize*3, yaxt="n", xaxt="n", col=COLOR,
         ylab="", xlab="", xlim=c(4,15), ylim=c(.25, 2.25), frame.plot=FALSE)
    axis(side=1, at=seq(4,15), lwd=0, lwd.ticks=resize, labels=2*c(4:12, 30, 50, 100))
    axis(side=2, at=seq(.25, 2.25, .25), las=2, lwd=resize)
    mtext(expression(paste("Ratio to Theoretical ", delta)), side=2, line=4, cex=resize*1.5)
    mtext(expression(italic(m[A]+m[B])), side=1, line=2.5, cex=resize*1.5)
    abline(h=1, col="red", lwd=resize)
    lines(ratio.u ~ dat.tmp$n.new, lwd=resize*3, lty=2, col=COLOR)
    leg0 <- legend("topright", legend=c("Serial", "Usual"), lty=1:2, lwd=resize*3, bty="n", seg.len=4)
    rect(xleft=leg0$rect$left, ybottom=leg0$rect$top-leg0$rect$h,
         xright=leg0$rect$left+leg0$rect$w, ytop=leg0$rect$top, lwd=resize)
    leg <- legend(x=leg0$rect$left, y=leg0$rect$top, legend=.rho, col=.color,
                  lwd=resize*3, ncol=2, bty="n", xjust=1)
    text(labels=expression(paste("True ", rho, " =")), x=leg$rect$left, y=leg$rect$top-.5*leg$rect$h,
         adj=c(1,0.5), cex=1.25)
    par(cex=resize*1.5)
    txtwidth <- strwidth(expression(paste("True ", rho, " =")))
    par(cex=resize)
    rect(xleft=leg$rect$left-txtwidth, ybottom=leg$rect$top-leg$rect$h,
         xright=leg$rect$left+leg$rect$w, ytop=leg$rect$top, lwd=resize)
    par(xpd=NA)
    lines(y = c(par("usr")[3],par("usr")[4]), x = c(par("usr")[1],par("usr")[1]), lwd=resize) # Left y axis
    lines(y = c(par("usr")[3],par("usr")[4]), x = c(par("usr")[2],par("usr")[2]), lwd=resize) # Right y axis
    lines(y = c(par("usr")[4],par("usr")[4]), x = c(par("usr")[1],par("usr")[2]), lwd=resize) # Top x axis
    lines(y = c(par("usr")[3],par("usr")[3]), x = c(par("usr")[1],12.25), lwd=resize) # Bottom x axis 1
    lines(y = c(par("usr")[3],par("usr")[3]), x = c(12.75,par("usr")[2]), lwd=resize) # Bottom x axis 2
    text("l", x = 12.25, y = par("usr")[3], srt = -45, font = 2)
    text("l", x = 12.75, y = par("usr")[3], srt = -45, font = 2)
    par(xpd=FALSE)
  }
  else{
    lines(ratio.s ~ dat.tmp$n.new, lwd=resize*3, lty=1, col=COLOR) 
    lines(ratio.u ~ dat.tmp$n.new, lwd=resize*3, lty=2, col=COLOR)
  }
}

dev.off()


#-------------------------------------------- str - Figure 6 ----------------------------------------------
# ***Set up ----
tiff("N-of-1 t-tests Figure 6.tiff", bg="white", width=3200, height=3200, compression="lzw")
layout(matrix(c(1,2,3), 3, 1), heights=c(1,.7,1))
resize <- 4.5
par(cex=resize, cex.axis=1.25, cex.main=1.25, cex.sub=1.25)

# ***First panel - Type I errors ----
for(ITER in 1:length(.rho)){
  RHO <- .rho[ITER]
  COLOR <- .color[ITER]
  
  dat.tmp <- typeIdat.tab[.("str", RHO, c(5:12, 30, 50, 100))] 
  
  if(ITER==1){
    par(mar=c(4.1, 6.1, .6, 2.1))
    plot(serialtypeI ~ n.new, data=dat.tmp, type="l", col=COLOR, lwd=resize*3,
         xlab="", ylab="", yaxt="n", xaxt="n", xlim=c(5,15), ylim=c(0,.3), frame.plot=FALSE)
    axis(side=1, at=seq(5,15), lwd=0, lwd.ticks=resize, labels=2*c(5:12, 30, 50, 100))
    axis(side=2, at=seq(0,.3,.05), las=2, lwd=resize)
    mtext("Type I Error", side=2, line=4, cex=resize*1.5)
    mtext(expression(italic(m[A]+m[B])), side=1, line=2.5, cex=resize*1.5)
    leg0 <- legend("topright", legend=c("Serial", "Usual"), lty=1:2, lwd=resize*3, bty="n", seg.len=4)
    rect(xleft=leg0$rect$left, ybottom=leg0$rect$top-leg0$rect$h,
         xright=leg0$rect$left+leg0$rect$w, ytop=leg0$rect$top, lwd=resize)
    leg <- legend(x=leg0$rect$left, y=leg0$rect$top, legend=.rho, col=.color,
                  lwd=resize*3, ncol=2, bty="n", xjust=1)
    text(labels=expression(paste("True ", rho, " =")), x=leg$rect$left, y=leg$rect$top-.5*leg$rect$h,
         adj=c(1,0.5), cex=1.25)
    par(cex=resize*1.5)
    txtwidth <- strwidth(expression(paste("True ", rho, " =")))
    par(cex=resize)
    rect(xleft=leg$rect$left-txtwidth, ybottom=leg$rect$top-leg$rect$h,
         xright=leg$rect$left+leg$rect$w, ytop=leg$rect$top, lwd=resize)
    abline(h=.05, col="red", lwd=resize)
    lines(usualtypeI ~ n.new, data=dat.tmp, lwd=resize*3, lty=2, col=COLOR)
    par(xpd=NA)
    lines(y = c(par("usr")[3],par("usr")[4]), x = c(par("usr")[1],par("usr")[1]), lwd=resize) # Left y axis
    lines(y = c(par("usr")[3],par("usr")[4]), x = c(par("usr")[2],par("usr")[2]), lwd=resize) # Right y axis
    lines(y = c(par("usr")[4],par("usr")[4]), x = c(par("usr")[1],par("usr")[2]), lwd=resize) # Top x axis
    lines(y = c(par("usr")[3],par("usr")[3]), x = c(par("usr")[1],12.25), lwd=resize) # Bottom x axis 1
    lines(y = c(par("usr")[3],par("usr")[3]), x = c(12.75,par("usr")[2]), lwd=resize) # Bottom x axis 2
    text("l", x = 12.25, y = par("usr")[3], srt = -45, font = 2)
    text("l", x = 12.75, y = par("usr")[3], srt = -45, font = 2)
    par(xpd=FALSE)
  }
  else{
    lines(serialtypeI ~ n.new, data=dat.tmp, col=COLOR, lwd=resize*3) 
    lines(usualtypeI ~ n.new, data=dat.tmp, col=COLOR, lwd=resize*3, lty=2) 
  }
}

# ***Second panel - theoretical effect sizes ----
for(ITER in 1:length(.rho)){
  RHO <- .rho[ITER]
  COLOR <- .color[ITER]
  
  dat.tmp <- ESdat.tab[.("str", RHO, c(5:12, 30, 50, 100))]
  
  if(ITER==1){
    par(mar=c(4.1, 6.1, .1, 2.1))
    plot(dth ~ n.new, data=dat.tmp, type="l", col=COLOR, lwd=resize*3,
         xlab="", ylab="", yaxt="n", xaxt="n", xlim=c(5,15), ylim=c(0,3), frame.plot=FALSE)
    axis(side=1, at=seq(5,15), lwd=0, lwd.ticks=resize, labels=2*c(5:12, 30, 50, 100)) 
    axis(side=2, at=seq(0,3,1), las=2, lwd=resize)
    mtext(expression(paste("Theoretical ", delta)), side=2, line=4, cex=resize*1.5)
    mtext(expression(italic(m[A]+m[B])), side=1, line=2.5, cex=resize*1.5)
    leg <- legend("topright", legend=.rho, col=.color, lwd=resize*3, ncol=2, bty="n")
    text(labels=expression(paste("True ", rho, " =")), x=leg$rect$left, y=leg$rect$top-.5*leg$rect$h,
         adj=c(1,0.5), cex=1.25)
    par(cex=resize*1.5)
    txtwidth <- strwidth(expression(paste("True ", rho, " =")))
    par(cex=resize)
    rect(xleft=leg$rect$left-txtwidth, ybottom=leg$rect$top-leg$rect$h,
         xright=leg$rect$left+leg$rect$w, ytop=leg$rect$top, lwd=resize)
    par(xpd=NA)
    lines(y = c(par("usr")[3],par("usr")[4]), x = c(par("usr")[1],par("usr")[1]), lwd=resize) # Left y axis
    lines(y = c(par("usr")[3],par("usr")[4]), x = c(par("usr")[2],par("usr")[2]), lwd=resize) # Right y axis
    lines(y = c(par("usr")[4],par("usr")[4]), x = c(par("usr")[1],par("usr")[2]), lwd=resize) # Top x axis
    lines(y = c(par("usr")[3],par("usr")[3]), x = c(par("usr")[1],12.25), lwd=resize) # Bottom x axis 1
    lines(y = c(par("usr")[3],par("usr")[3]), x = c(12.75,par("usr")[2]), lwd=resize) # Bottom x axis 2
    text("l", x = 12.25, y = par("usr")[3], srt = -45, font = 2)
    text("l", x = 12.75, y = par("usr")[3], srt = -45, font = 2)
    par(xpd=FALSE)
  }
  else{
    lines(dth ~ n.new, data=dat.tmp, col=COLOR, lwd=resize*3) 
  }
}

# ***Third panel - effect size ratios ----
for(ITER in 1:length(.rho)){
  RHO <- .rho[ITER]
  COLOR <- .color[ITER]
  
  dat.tmp <- ESdat.tab[.("str", RHO, c(5:12, 30, 50, 100))] 
  
  ratio.s <- dat.tmp$ds / dat.tmp$dth
  ratio.u <- dat.tmp$du / dat.tmp$dth
  
  if(ITER==1){
    par(mar=c(3.6, 6.1, .1, 2.1))
    plot(ratio.s ~ dat.tmp$n.new, type="l", lty=1, lwd=resize*3, yaxt="n", xaxt="n", col=COLOR,
         ylab="", xlab="", xlim=c(5,15), ylim=c(.25, 1.5), frame.plot=FALSE)
    axis(side=1, at=seq(5,15), lwd=0, lwd.ticks=resize, labels=2*c(5:12, 30, 50, 100)) 
    axis(side=2, at=seq(.25, 1.5, .25), las=2, lwd=resize)
    mtext(expression(paste("Ratio to Theoretical ", delta)), side=2, line=4, cex=resize*1.5)
    mtext(expression(italic(m[A]+m[B])), side=1, line=2.5, cex=resize*1.5)
    abline(h=1, col="red", lwd=resize)
    lines(ratio.u ~ dat.tmp$n.new, lwd=resize*3, lty=2, col=COLOR)
    leg0 <- legend("bottomright", legend=c("Serial", "Usual"), lty=1:2, lwd=resize*3, bty="n", seg.len=4)
    rect(xleft=leg0$rect$left, ybottom=leg0$rect$top-leg0$rect$h,
         xright=leg0$rect$left+leg0$rect$w, ytop=leg0$rect$top, lwd=resize)
    leg <- legend(x=leg0$rect$left, y=leg0$rect$top, legend=.rho, col=.color,
                  lwd=resize*3, ncol=2, bty="n", xjust=1)
    text(labels=expression(paste("True ", rho, " =")), x=leg$rect$left, y=leg$rect$top-.5*leg$rect$h,
         adj=c(1,0.5), cex=1.25)
    par(cex=resize*1.5)
    txtwidth <- strwidth(expression(paste("True ", rho, " =")))
    par(cex=resize)
    lines(y = c(leg$rect$top-leg$rect$h,leg$rect$top),
          x = c(leg$rect$left-txtwidth,leg$rect$left-txtwidth), lwd=resize) # Left y box
    lines(y = c(leg$rect$top-leg$rect$h,leg$rect$top),
          x = c(leg$rect$left+leg$rect$w,leg$rect$left+leg$rect$w), lwd=resize) # Right y box
    lines(y = c(leg$rect$top,leg$rect$top),
          x = c(leg$rect$left-txtwidth,leg$rect$left+leg$rect$w), lwd=resize) # Top x box
    lines(y = c(leg$rect$top-leg$rect$h,leg$rect$top-leg$rect$h),
          x = c(leg$rect$left-txtwidth,12.25), lwd=resize) # Bottom x box 1
    lines(y = c(leg$rect$top-leg$rect$h,leg$rect$top-leg$rect$h),
          x = c(12.75,leg$rect$left+leg$rect$w), lwd=resize) # Bottom x box 2
    par(xpd=NA)
    lines(y = c(par("usr")[3],par("usr")[4]), x = c(par("usr")[1],par("usr")[1]), lwd=resize) # Left y axis
    lines(y = c(par("usr")[3],par("usr")[4]), x = c(par("usr")[2],par("usr")[2]), lwd=resize) # Right y axis
    lines(y = c(par("usr")[4],par("usr")[4]), x = c(par("usr")[1],par("usr")[2]), lwd=resize) # Top x axis
    lines(y = c(par("usr")[3],par("usr")[3]), x = c(par("usr")[1],12.25), lwd=resize) # Bottom x axis 1
    lines(y = c(par("usr")[3],par("usr")[3]), x = c(12.75,par("usr")[2]), lwd=resize) # Bottom x axis 2
    text("l", x = 12.25, y = par("usr")[3], srt = -45, font = 2)
    text("l", x = 12.75, y = par("usr")[3], srt = -45, font = 2)
    par(xpd=FALSE)
  }
  else{
    lines(ratio.s ~ dat.tmp$n.new, lwd=resize*3, lty=1, col=COLOR) 
    lines(ratio.u ~ dat.tmp$n.new, lwd=resize*3, lty=2, col=COLOR)
  }
}

dev.off()


#-------------------------------------------- pst - Figure 3 ----------------------------------------------
# ***Set up ----
tiff("N-of-1 t-tests Figure 3.tiff", bg="white", width=3200, height=3200, compression="lzw")
layout(matrix(c(1,2,3), 3, 1), heights=c(1,.7,1))
resize <- 4.5
par(cex=resize, cex.axis=1.25, cex.main=1.25, cex.sub=1.25)

# ***First panel - Type I errors ----
for(ITER in 1:length(.rho)){
  RHO <- .rho[ITER]
  COLOR <- .color[ITER]
  
  dat.tmp <- typeIdat.tab[.("pst", RHO)]
  
  if(ITER==1){
    par(mar=c(4.1, 6.1, .6, 2.1))
    plot(serialtypeI ~ n.new, data=dat.tmp, type="l", col=COLOR, lwd=resize*3,
         xlab="", ylab="", yaxt="n", xaxt="n", xlim=c(4,15), ylim=c(0,.35), frame.plot=FALSE)
    axis(side=1, at=seq(4,15), lwd=0, lwd.ticks=resize, labels=c(4:12, 30, 50, 100))
    axis(side=2, at=seq(0,.35,.05), las=2, lwd=resize)
    mtext("Type I Error", side=2, line=4, cex=resize*1.5)
    mtext(expression(italic("m")), side=1, line=2.5, cex=resize*1.5)
    leg0 <- legend("topright", legend=c("Serial", "Usual"), lty=1:2, lwd=resize*3, bty="n", seg.len=4)
    rect(xleft=leg0$rect$left, ybottom=leg0$rect$top-leg0$rect$h,
         xright=leg0$rect$left+leg0$rect$w, ytop=leg0$rect$top, lwd=resize)
    leg <- legend(x=leg0$rect$left, y=leg0$rect$top, legend=.rho, col=.color,
                  lwd=resize*3, ncol=2, bty="n", xjust=1)
    text(labels=expression(paste("True ", rho, " =")), x=leg$rect$left, y=leg$rect$top-.5*leg$rect$h,
         adj=c(1,0.5), cex=1.25)
    par(cex=resize*1.5)
    txtwidth <- strwidth(expression(paste("True ", rho, " =")))
    par(cex=resize)
    rect(xleft=leg$rect$left-txtwidth, ybottom=leg$rect$top-leg$rect$h,
         xright=leg$rect$left+leg$rect$w, ytop=leg$rect$top, lwd=resize)
    abline(h=.05, col="red", lwd=resize)
    lines(usualtypeI ~ n.new, data=dat.tmp, lwd=resize*3, lty=2, col=COLOR) 
    par(xpd=NA)
    lines(y = c(par("usr")[3],par("usr")[4]), x = c(par("usr")[1],par("usr")[1]), lwd=resize) # Left y axis
    lines(y = c(par("usr")[3],par("usr")[4]), x = c(par("usr")[2],par("usr")[2]), lwd=resize) # Right y axis
    lines(y = c(par("usr")[4],par("usr")[4]), x = c(par("usr")[1],par("usr")[2]), lwd=resize) # Top x axis
    lines(y = c(par("usr")[3],par("usr")[3]), x = c(par("usr")[1],12.25), lwd=resize) # Bottom x axis 1
    lines(y = c(par("usr")[3],par("usr")[3]), x = c(12.75,par("usr")[2]), lwd=resize) # Bottom x axis 2
    text("l", x = 12.25, y = par("usr")[3], srt = -45, font = 2)
    text("l", x = 12.75, y = par("usr")[3], srt = -45, font = 2)
    par(xpd=FALSE)
  }
  else{
    lines(serialtypeI ~ n.new, data=dat.tmp, col=COLOR, lwd=resize*3) 
    lines(usualtypeI ~ n.new, data=dat.tmp, col=COLOR, lwd=resize*3, lty=2) 
  }
}

# ***Second panel - theoretical effect sizes ----
for(ITER in 1:length(.rho)){
  RHO <- .rho[ITER]
  COLOR <- .color[ITER]
  
  dat.tmp <- ESdat.tab[.("pst", RHO)] 
  
  if(ITER==1){
    par(mar=c(4.1, 6.1, .1, 2.1))
    plot(dth ~ n.new, data=dat.tmp, type="l", col=COLOR, lwd=resize*3,
         xlab="", ylab="", yaxt="n", xaxt="n", xlim=c(4,15), ylim=c(0,10), frame.plot=FALSE)
    axis(side=1, at=seq(4,15), lwd=0, lwd.ticks=resize, labels=c(4:12, 30, 50, 100))
    axis(side=2, at=seq(0,10,2), las=2, lwd=resize)
    mtext(expression(paste("Theoretical ", delta)), side=2, line=4, cex=resize*1.5)
    mtext(expression(italic("m")), side=1, line=2.5, cex=resize*1.5)
    leg <- legend("topright", legend=.rho, col=.color, lwd=resize*3, ncol=2, bty="n")
    text(labels=expression(paste("True ", rho, " =")), x=leg$rect$left, y=leg$rect$top-.5*leg$rect$h,
         adj=c(1,0.5), cex=1.25)
    par(cex=resize*1.5)
    txtwidth <- strwidth(expression(paste("True ", rho, " =")))
    par(cex=resize)
    rect(xleft=leg$rect$left-txtwidth, ybottom=leg$rect$top-leg$rect$h,
         xright=leg$rect$left+leg$rect$w, ytop=leg$rect$top, lwd=resize)
    par(xpd=NA)
    lines(y = c(par("usr")[3],par("usr")[4]), x = c(par("usr")[1],par("usr")[1]), lwd=resize) # Left y axis
    lines(y = c(par("usr")[3],par("usr")[4]), x = c(par("usr")[2],par("usr")[2]), lwd=resize) # Right y axis
    lines(y = c(par("usr")[4],par("usr")[4]), x = c(par("usr")[1],par("usr")[2]), lwd=resize) # Top x axis
    lines(y = c(par("usr")[3],par("usr")[3]), x = c(par("usr")[1],12.25), lwd=resize) # Bottom x axis 1
    lines(y = c(par("usr")[3],par("usr")[3]), x = c(12.75,par("usr")[2]), lwd=resize) # Bottom x axis 2
    text("l", x = 12.25, y = par("usr")[3], srt = -45, font = 2)
    text("l", x = 12.75, y = par("usr")[3], srt = -45, font = 2)
    par(xpd=FALSE)
  }
  else{
    lines(dth ~ n.new, data=dat.tmp, col=COLOR, lwd=resize*3)
  }
}

# ***Third panel - effect size ratios ----
for(ITER in 1:length(.rho)){
  RHO <- .rho[ITER]
  COLOR <- .color[ITER]
  
  dat.tmp <- ESdat.tab[.("pst", RHO)]
  
  ratio.s <- dat.tmp$ds / dat.tmp$dth
  ratio.u <- dat.tmp$du / dat.tmp$dth
  
  if(ITER==1){
    par(mar=c(3.6, 6.1, .1, 2.1))
    plot(ratio.s ~ dat.tmp$n.new, type="l", lty=1, lwd=resize*3, yaxt="n", xaxt="n", col=COLOR,
         ylab="", xlab="", xlim=c(4,15), ylim=c(0, 2), frame.plot=FALSE)
    axis(side=1, at=seq(4,15), lwd=0, lwd.ticks=resize, labels=c(4:12, 30, 50, 100))
    axis(side=2, at=seq(0, 2, .25), las=2, lwd=resize)
    mtext(expression(paste("Ratio to Theoretical ", delta)), side=2, line=4, cex=resize*1.5)
    mtext(expression(italic("m")), side=1, line=2.5, cex=resize*1.5)
    abline(h=1, col="red", lwd=resize)
    lines(ratio.u ~ dat.tmp$n.new, lwd=resize*3, lty=2, col=COLOR)
    leg0 <- legend("topright", legend=c("Serial", "Usual"), lty=1:2, lwd=resize*3, bty="n", seg.len=4)
    rect(xleft=leg0$rect$left, ybottom=leg0$rect$top-leg0$rect$h,
         xright=leg0$rect$left+leg0$rect$w, ytop=leg0$rect$top, lwd=resize)
    leg <- legend(x=leg0$rect$left, y=leg0$rect$top, legend=.rho, col=.color,
                  lwd=resize*3, ncol=2, bty="n", xjust=1)
    text(labels=expression(paste("True ", rho, " =")), x=leg$rect$left, y=leg$rect$top-.5*leg$rect$h,
         adj=c(1,0.5), cex=1.25)
    par(cex=resize*1.5)
    txtwidth <- strwidth(expression(paste("True ", rho, " =")))
    par(cex=resize)
    rect(xleft=leg$rect$left-txtwidth, ybottom=leg$rect$top-leg$rect$h,
         xright=leg$rect$left+leg$rect$w, ytop=leg$rect$top, lwd=resize)
    par(xpd=NA)
    lines(y = c(par("usr")[3],par("usr")[4]), x = c(par("usr")[1],par("usr")[1]), lwd=resize) # Left y axis
    lines(y = c(par("usr")[3],par("usr")[4]), x = c(par("usr")[2],par("usr")[2]), lwd=resize) # Right y axis
    lines(y = c(par("usr")[4],par("usr")[4]), x = c(par("usr")[1],par("usr")[2]), lwd=resize) # Top x axis
    lines(y = c(par("usr")[3],par("usr")[3]), x = c(par("usr")[1],12.25), lwd=resize) # Bottom x axis 1
    lines(y = c(par("usr")[3],par("usr")[3]), x = c(12.75,par("usr")[2]), lwd=resize) # Bottom x axis 2
    text("l", x = 12.25, y = par("usr")[3], srt = -45, font = 2)
    text("l", x = 12.75, y = par("usr")[3], srt = -45, font = 2)
    par(xpd=FALSE)
  }
  else{
    lines(ratio.s ~ dat.tmp$n.new, lwd=resize*3, lty=1, col=COLOR) 
    lines(ratio.u ~ dat.tmp$n.new, lwd=resize*3, lty=2, col=COLOR)
  }
}

dev.off()


#-------------------------------------------- pstr - Figure 5 ---------------------------------------------
# ***Set up ----
tiff("N-of-1 t-tests Figure 5.tiff", bg="white", width=3200, height=3200, compression="lzw")
layout(matrix(c(1,2,3), 3, 1), heights=c(1,.7,1))
resize <- 4.5
par(cex=resize, cex.axis=1.25, cex.main=1.25, cex.sub=1.25)

# ***First panel - Type I errors ----
for(ITER in 1:length(.rho)){
  RHO <- .rho[ITER]
  COLOR <- .color[ITER]
  
  dat.tmp <- typeIdat.tab[.("pstr", RHO, c(5:12, 30, 50, 100))] 
  
  if(ITER==1){
    par(mar=c(4.1, 6.1, .6, 2.1))
    plot(serialtypeI ~ n.new, data=dat.tmp, type="l", col=COLOR, lwd=resize*3,
         xlab="", ylab="", yaxt="n", xaxt="n", xlim=c(5,15), ylim=c(0,.35), frame.plot=FALSE)
    axis(side=1, at=seq(5,15), lwd=0, lwd.ticks=resize, labels=c(5:12, 30, 50, 100))
    axis(side=2, at=seq(0,.35,.05), las=2, lwd=resize)
    mtext("Type I Error", side=2, line=4, cex=resize*1.5)
    mtext(expression(italic("m")), side=1, line=2.5, cex=resize*1.5)
    leg0 <- legend("topright", legend=c("Serial", "Usual"), lty=1:2, lwd=resize*3, bty="n", seg.len=4)
    rect(xleft=leg0$rect$left, ybottom=leg0$rect$top-leg0$rect$h,
         xright=leg0$rect$left+leg0$rect$w, ytop=leg0$rect$top, lwd=resize)
    leg <- legend(x=leg0$rect$left, y=leg0$rect$top, legend=.rho, col=.color,
                  lwd=resize*3, ncol=2, bty="n", xjust=1)
    text(labels=expression(paste("True ", rho, " =")), x=leg$rect$left, y=leg$rect$top-.5*leg$rect$h,
         adj=c(1,0.5), cex=1.25)
    par(cex=resize*1.5)
    txtwidth <- strwidth(expression(paste("True ", rho, " =")))
    par(cex=resize)
    rect(xleft=leg$rect$left-txtwidth, ybottom=leg$rect$top-leg$rect$h,
         xright=leg$rect$left+leg$rect$w, ytop=leg$rect$top, lwd=resize)
    abline(h=.05, col="red", lwd=resize)
    lines(usualtypeI ~ n.new, data=dat.tmp, lwd=resize*3, lty=2, col=COLOR)
    par(xpd=NA)
    lines(y = c(par("usr")[3],par("usr")[4]), x = c(par("usr")[1],par("usr")[1]), lwd=resize) # Left y axis
    lines(y = c(par("usr")[3],par("usr")[4]), x = c(par("usr")[2],par("usr")[2]), lwd=resize) # Right y axis
    lines(y = c(par("usr")[4],par("usr")[4]), x = c(par("usr")[1],par("usr")[2]), lwd=resize) # Top x axis
    lines(y = c(par("usr")[3],par("usr")[3]), x = c(par("usr")[1],12.25), lwd=resize) # Bottom x axis 1
    lines(y = c(par("usr")[3],par("usr")[3]), x = c(12.75,par("usr")[2]), lwd=resize) # Bottom x axis 2
    text("l", x = 12.25, y = par("usr")[3], srt = -45, font = 2)
    text("l", x = 12.75, y = par("usr")[3], srt = -45, font = 2)
    par(xpd=FALSE)
  }
  else{
    lines(serialtypeI ~ n.new, data=dat.tmp, col=COLOR, lwd=resize*3) 
    lines(usualtypeI ~ n.new, data=dat.tmp, col=COLOR, lwd=resize*3, lty=2) 
  }
}

# ***Second panel - theoretical effect sizes ----
for(ITER in 1:length(.rho)){
  RHO <- .rho[ITER]
  COLOR <- .color[ITER]
  
  dat.tmp <- ESdat.tab[.("pstr", RHO, c(5:12, 30, 50, 100))]
  
  if(ITER==1){
    par(mar=c(4.1, 6.1, .1, 2.1))
    plot(dth ~ n.new, data=dat.tmp, type="l", col=COLOR, lwd=resize*3,
         xlab="", ylab="", yaxt="n", xaxt="n", xlim=c(5,15), ylim=c(0,10), frame.plot=FALSE)
    axis(side=1, at=seq(5,15), lwd=0, lwd.ticks=resize, labels=c(5:12, 30, 50, 100)) 
    axis(side=2, at=seq(0,10,2), las=2, lwd=resize)
    mtext(expression(paste("Theoretical ", delta)), side=2, line=4, cex=resize*1.5)
    mtext(expression(italic("m")), side=1, line=2.5, cex=resize*1.5)
    leg <- legend("topright", legend=.rho, col=.color, lwd=resize*3, ncol=2, bty="n")
    text(labels=expression(paste("True ", rho, " =")), x=leg$rect$left, y=leg$rect$top-.5*leg$rect$h,
         adj=c(1,0.5), cex=1.25)
    par(cex=resize*1.5)
    txtwidth <- strwidth(expression(paste("True ", rho, " =")))
    par(cex=resize)
    rect(xleft=leg$rect$left-txtwidth, ybottom=leg$rect$top-leg$rect$h,
         xright=leg$rect$left+leg$rect$w, ytop=leg$rect$top, lwd=resize)
    par(xpd=NA)
    lines(y = c(par("usr")[3],par("usr")[4]), x = c(par("usr")[1],par("usr")[1]), lwd=resize) # Left y axis
    lines(y = c(par("usr")[3],par("usr")[4]), x = c(par("usr")[2],par("usr")[2]), lwd=resize) # Right y axis
    lines(y = c(par("usr")[4],par("usr")[4]), x = c(par("usr")[1],par("usr")[2]), lwd=resize) # Top x axis
    lines(y = c(par("usr")[3],par("usr")[3]), x = c(par("usr")[1],12.25), lwd=resize) # Bottom x axis 1
    lines(y = c(par("usr")[3],par("usr")[3]), x = c(12.75,par("usr")[2]), lwd=resize) # Bottom x axis 2
    text("l", x = 12.25, y = par("usr")[3], srt = -45, font = 2)
    text("l", x = 12.75, y = par("usr")[3], srt = -45, font = 2)
    par(xpd=FALSE)
  }
  else{
    lines(dth ~ n.new, data=dat.tmp, col=COLOR, lwd=resize*3) 
  }
}

# ***Third panel - effect size ratios ----
for(ITER in 1:length(.rho)){
  RHO <- .rho[ITER]
  COLOR <- .color[ITER]
  
  dat.tmp <- ESdat.tab[.("pstr", RHO, c(5:12, 30, 50, 100))] 
  
  ratio.s <- dat.tmp$ds / dat.tmp$dth
  ratio.u <- dat.tmp$du / dat.tmp$dth
  
  if(ITER==1){
    par(mar=c(3.6, 6.1, .1, 2.1))
    plot(ratio.s ~ dat.tmp$n.new, type="l", lty=1, lwd=resize*3, yaxt="n", xaxt="n", col=COLOR,
         ylab="", xlab="", xlim=c(5,15), ylim=c(0, 1.75), frame.plot=FALSE)
    axis(side=1, at=seq(5,15), lwd=0, lwd.ticks=resize, labels=c(5:12, 30, 50, 100)) 
    axis(side=2, at=seq(0, 1.75, .25), las=2, lwd=resize)
    mtext(expression(paste("Ratio to Theoretical ", delta)), side=2, line=4, cex=resize*1.5)
    mtext(expression(italic("m")), side=1, line=2.5, cex=resize*1.5)
    abline(h=1, col="red", lwd=resize)
    lines(ratio.u ~ dat.tmp$n.new, lwd=resize*3, lty=2, col=COLOR)
    leg0 <- legend("topright", legend=c("Serial", "Usual"), lty=1:2, lwd=resize*3, bty="n", seg.len=4)
    rect(xleft=leg0$rect$left, ybottom=leg0$rect$top-leg0$rect$h,
         xright=leg0$rect$left+leg0$rect$w, ytop=leg0$rect$top, lwd=resize)
    leg <- legend(x=leg0$rect$left, y=leg0$rect$top, legend=.rho, col=.color,
                  lwd=resize*3, ncol=2, bty="n", xjust=1)
    text(labels=expression(paste("True ", rho, " =")), x=leg$rect$left, y=leg$rect$top-.5*leg$rect$h,
         adj=c(1,0.5), cex=1.25)
    par(cex=resize*1.5)
    txtwidth <- strwidth(expression(paste("True ", rho, " =")))
    par(cex=resize)
    rect(xleft=leg$rect$left-txtwidth, ybottom=leg$rect$top-leg$rect$h,
         xright=leg$rect$left+leg$rect$w, ytop=leg$rect$top, lwd=resize)
    par(xpd=NA)
    lines(y = c(par("usr")[3],par("usr")[4]), x = c(par("usr")[1],par("usr")[1]), lwd=resize) # Left y axis
    lines(y = c(par("usr")[3],par("usr")[4]), x = c(par("usr")[2],par("usr")[2]), lwd=resize) # Right y axis
    lines(y = c(par("usr")[4],par("usr")[4]), x = c(par("usr")[1],par("usr")[2]), lwd=resize) # Top x axis
    lines(y = c(par("usr")[3],par("usr")[3]), x = c(par("usr")[1],12.25), lwd=resize) # Bottom x axis 1
    lines(y = c(par("usr")[3],par("usr")[3]), x = c(12.75,par("usr")[2]), lwd=resize) # Bottom x axis 2
    text("l", x = 12.25, y = par("usr")[3], srt = -45, font = 2)
    text("l", x = 12.75, y = par("usr")[3], srt = -45, font = 2)
    par(xpd=FALSE)
  }
  else{
    lines(ratio.s ~ dat.tmp$n.new, lwd=resize*3, lty=1, col=COLOR) 
    lines(ratio.u ~ dat.tmp$n.new, lwd=resize*3, lty=2, col=COLOR)
  }
}

dev.off()
