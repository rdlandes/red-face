#==========================================================================================================
#**********************************************************************************************************
# N-of-1 t-tests: Delay discounting example
# Run time: less than 1 minute
#           on a Dell OptiPlex 990 with an 
#           Intel Core i7-2600 CPU @ 3.40 GHz
#           with 8.0 GB of RAM
#           R version 3.4.1
#
# The data come from 
# [22] Landes RD, Christensen DR, Bickel WK. Delay discounting decreases in those completing
#      treatment for opioid dependence. Exp Clin Psychopharmacol. 2012; 
#      20(4):302-309. 
# Original results are summarized in Table 3 of [22].
#
# The original regression method used refers to 
# [23] Landes, R.D., Pitcock, J.A., Yi, R. and Bickel, W.K. (2010)  Analytical methods to detect 
#      within-individual changes in discounting. Experimental and clinical psychopharmacology, 18(2), 
#      pp. 175-183.
#
#**********************************************************************************************************
#==========================================================================================================

#==========================================================================================================
#------------------------------------------------- Set-up -------------------------------------------------
#==========================================================================================================
# Change to correct folder

#setwd("<insert location of example analysis data file>")
setwd("I:/Tang N-of-1 (8-30-2018 ff)")

# Load library
library(data.table)

# Read in dataset
#dat <- fread("Landes et al. (2012) Table 3 select data.csv")
pairdat <- fread("N-of-1 t-tests Delay discounting raw data.csv")
dat.tab <- as.data.table(pairdat)
setkey(dat.tab, PATIENT)


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
#---------------------------------- Paired serial t-test for level-change ----------------------------------
# Returns degrees of freedom, estimated serial correlation, mean difference,
# standard error of mean difference, t-statistic, 2-sided and 1-sided p-values.
#
# Parameters
#  v1 is the vector of data under condition 1.
#  v0 is the vector of data under condition 0. 
#   If v1 is a vector of differences, v0 is a vector of 0s by default. 
#  rhotype is the type of serial correlation estimate to use. If "WF" (Wayne Fuller's) is not
#    specified, rhow must be inputted and will be used as the serial correlation estimate.
#  rhow is an inputted serial correlation. Only used if rhotype is not specified.
pst.fxn <- function(v1, v0=rep(0, length(v1)), rhotype, rhow){
  # Calculating vector of differences
  diff <- v1 - v0
  
  # Counting length of differences
  n <- length(diff)
  
  # Constructing the variables for the linear model
  one <- rep(1, n)
  y <- diff
  
  # Fitting the linear model & obtaining summary stats
  m <- summary( lm(y ~ one - 1))
  int  <- m$coefficients[1,1]
  sd.y <- m$sigma 
  names(int)  <- NULL
  names(sd.y) <- NULL
  
  # Obtaining residuals & computing correlation
  res <- m$residuals
  if(rhotype=="WF"){
    rho <- WF.r(res)
  } else{
    rho <- rhow
  }
  
  # Test
  c <- (n + 2*rho^(n+1) - n*rho^2 - 2*rho) / (n^2*(rho-1)^2)
  b <- ( n - (n+2*rho^(n+1)-n*rho^2-2*rho)/(n*(rho-1)^2) ) / (n-1)
  nprime <- 1 / (1 - (n-1)*b/n)
  se <- sd.y * sqrt( c/ b )
  DF <- nprime - 1
  sd <- sd.y
  t <- int / se
  pval <- 2*pt(-abs(t), df=DF)
  pval.lo <- pt(t, df=DF)
  pval.up <- pt(-t, df=DF)
  
  # Return summary
  pstSum <- cbind(DF, rho, int, sd, se, t, pval, pval.lo, pval.up)
  return(pstSum)
}

#---------------------------------- Paired serial t-test for rate-change -----------------------------------
# Returns degrees of freedom, estimated serial correlation, slope of difference,
# standard error of slope of difference, t-statistic, 2-sided and 1-sided p-values.
#
# Parameters
#  v1 is the vector of data under condition 1.
#  v0 is the vector of data under condition 0. 
#   If v1 is a vector of differences, v0 is a vector of 0s by default. 
#  rhotype is the type of serial correlation estimate to use. If "WF" (Wayne Fuller's) is not
#    specified, rhow must be inputted and will be used as the serial correlation estimate.
#  rhow is an inputted serial correlation. Only used if rhotype is not specified.
pstr.fxn <- function(v1, v0=rep(0, length(v1)), rhotype, rhow){
  # Calculating vector of differences
  diff <- v1 - v0
  
  # Counting length of differences
  n <- length(diff)
  
  # Constructing the variables for the linear model
  one <- rep(1, n)
  x <- 1:n - mean(1:n)
  y <- diff
  
  # Fitting the linear model & obtaining summary stats
  m <- summary( lm(y ~ one + x - 1))
  int  <- m$coefficients[1,1]
  slp    <- m$coefficients[2,1]
  sd.y <- m$sigma 
  names(int)  <- NULL
  names(slp)    <- NULL
  names(sd.y) <- NULL
  
  # Obtaining residuals & computing correlation
  res <- m$residuals
  if(rhotype=="WF"){
    rho <- WF.r(res)
  } else{
    rho <- rhow
  }
  
  # Test
  c <- 12/(n^2-1)^2 * ( - 6*rho*(rho+1)^2*(rho^n-1)/(n^2*(rho-1)^4) +
                               2*rho*(6*rho^(n+1)+6*rho^n+rho^2-2*rho+1)/(n*(rho-1)^3) -
                               6*rho*(rho^n+1)/(rho-1)^2 - 2*n*rho/(rho-1) + (n^2-1)/n )
  b <- ( ( n^4*(rho-1)^4 + 2*n^3*(rho+1)*(rho-1)^3 + n^2*(rho-1)^2*(4*rho^(n+1)-rho^2+10*rho-1)
                - 2*n*(rho^2-1)*(6*rho^(n+1)+rho^2-2*rho+1) + 2*(3*rho^3+10*rho^2+7*rho+4)*(rho^(n+1)-rho) ) / 
                (n*(n-2)*(n^2-1)*(rho-1)^4) )
  nprime <- 2 / (1-b*(n-2)/n)
  se <- sd.y * sqrt( c/b )
  DF <- nprime - 2
  sd <- sd.y
  t <- slp / se
  pval <- 2*pt(-abs(t), df=DF)
  pval.lo <- pt(t, df=DF)
  pval.up <- pt(-t, df=DF)
  
  # Return summary
  pstrSum <- cbind(DF, rho, slp, sd, se, t, pval, pval.lo, pval.up)
  return(pstrSum)
}

#--------------------------------- 2-sample serial t-test for level-change ---------------------------------
# Returns degrees of freedom, estimated serial correlation, difference in means,
# standard error of difference in means, t-statistic, 2-sided p-value.
#
# Parameters
#  v0 is the vector of data under condition 0.
#  v1 is the vector of data under condition 1.
#  rhotype is the type of serial correlation estimate to use. If "WF" (Wayne Fuller's) is not
#    specified, rhow must be inputted and will be used as the serial correlation estimate.
#  rhow is an inputted serial correlation. Only used if rhotype is not specified.
st.fxn <- function(v0, v1, rhotype, rhow){
  # Counting length of sequences from 2 "independent" conditions
  n0 <- length(v0)
  n1 <- length(v1)
  
  # Constructing the variables for the linear model
  one0 <- c(rep(1, n0), rep(0, n1))
  one1 <- c(rep(0, n0), rep(1, n1))
  y <- c(v0, v1)
  
  # Fitting the linear model & obtaining summary stats
  m <- summary( lm(y ~ one0 + one1 -1))
  int0 <- m$coefficients[1,1]
  int1 <- m$coefficients[2,1]
  sd.y <- m$sigma 
  names(int0) <- NULL
  names(int1) <- NULL
  names(sd.y) <- NULL
  int <- int1-int0
  
  # Obtaining residuals & computing correlation
  res0 <- m$residuals[1:n0]
  res1 <- m$residuals[-(1:n0)]
  if(rhotype=="WF") {
    rho1 <- WF.r(res1)
    rho0 <- WF.r(res0)
    rho <- (n0*rho0 + n1*rho1) / (n0 + n1)
  } else{
    rho <- rhow
  }
  
  # Test
  c0 <- (n0 + 2*rho^(n0+1) - n0*rho^2 - 2*rho) / (n0^2*(rho-1)^2)
  b0 <- ( n0 - (n0+2*rho^(n0+1)-n0*rho^2-2*rho)/(n0*(rho-1)^2) ) / (n0-1)
  nprime0 <- 1 / (1 - (n0-1)*b0/n0)
  
  c1 <- (n1 + 2*rho^(n1+1) - n1*rho^2 - 2*rho) / (n1^2*(rho-1)^2)
  b1 <- ( n1 - (n1+2*rho^(n1+1)-n1*rho^2-2*rho)/(n1*(rho-1)^2) ) / (n1-1)
  nprime1 <- 1 / (1 - (n1-1)*b1/n1)
  se <- sd.y * sqrt( (c0/b0) + (c1/b1) )
  sd <- sd.y
  DF <- nprime0 + nprime1 - 2
  t <- int / se
  pval <- 2*pt(-abs(t), df=DF)
  pval.lo <- pt(t, df=DF)
  pval.up <- pt(-t, df=DF)
  
  # Return summary
  stSum <- cbind(DF, rho, int, sd, se, t, pval, pval.lo, pval.up)
  return(stSum)
}

#--------------------------------- 2-sample serial t-test for rate-change ----------------------------------
# Returns degrees of freedom, estimated serial correlation, difference in slopes,
# standard error of difference in slopes, t-statistic, 2-sided p-value.
#
# Parameters
#  v0 is the vector of data under condition 0.
#  v1 is the vector of data under condition 1.
#  rhotype is the type of serial correlation estimate to use. If "WF" (Wayne Fuller's) is not
#    specified, rhow must be inputted and will be used as the serial correlation estimate.
#  rhow is an inputted serial correlation. Only used if rhotype is not specified.
str.fxn <- function(v0, v1, rhotype, rhow){
  # Counting length of sequences from 2 "independent" conditions
  n0 <- length(v0)
  n1 <- length(v1)
  
  # Constructing the variables for the linear model
  one0 <- c(rep(1, n0), rep(0, n1))
  one1 <- c(rep(0, n0), rep(1, n1))
  x0 <- c( 1:n0 - mean(1:n0), 
           rep(0, n1) )
  x1 <- c( rep(0, n0),
           1:n1 - mean(1:n1) )
  y <- c(v0, v1)
  
  # Fitting the linear model & obtaining summary stats
  m <- summary( lm(y ~ one0 + one1 + x0 + x1 -1))
  slp0   <- m$coefficients[3,1]
  slp1   <- m$coefficients[4,1]
  sd.y <- m$sigma 
  names(slp0)   <- NULL
  names(slp1)   <- NULL
  names(sd.y) <- NULL
  slp.d <- slp1-slp0
  
  # Obtaining residuals & computing correlation
  res0 <- m$residuals[1:n0]
  res1 <- m$residuals[-(1:n0)]
  if(rhotype=="WF") {
    rho1 <- WF.r(res1)
    rho0 <- WF.r(res0)
    rho <- (n0*rho0 + n1*rho1) / (n0 + n1)
  } else{
    rho <- rhow
  }
  
  # Test
  c0 <- 12/(n0^2-1)^2 * ( - 6*rho*(rho+1)^2*(rho^n0-1)/(n0^2*(rho-1)^4) +
                          2*rho*(6*rho^(n0+1)+6*rho^n0+rho^2-2*rho+1)/(n0*(rho-1)^3) -
                          6*rho*(rho^n0+1)/(rho-1)^2 - 2*n0*rho/(rho-1) + (n0^2-1)/n0 )
  b0 <- ( ( n0^4*(rho-1)^4 + 2*n0^3*(rho+1)*(rho-1)^3 + n0^2*(rho-1)^2*(4*rho^(n0+1)-rho^2+10*rho-1)
           - 2*n0*(rho^2-1)*(6*rho^(n0+1)+rho^2-2*rho+1) + 2*(3*rho^3+10*rho^2+7*rho+4)*(rho^(n0+1)-rho) ) / 
           (n0*(n0-2)*(n0^2-1)*(rho-1)^4) )
  nprime0 <- 2 / (1-b0*(n0-2)/n0)
  c1 <- 12/(n1^2-1)^2 * ( - 6*rho*(rho+1)^2*(rho^n1-1)/(n1^2*(rho-1)^4) +
                          2*rho*(6*rho^(n1+1)+6*rho^n1+rho^2-2*rho+1)/(n1*(rho-1)^3) -
                          6*rho*(rho^n1+1)/(rho-1)^2 - 2*n1*rho/(rho-1) + (n1^2-1)/n1 )
  b1 <- ( ( n1^4*(rho-1)^4 + 2*n1^3*(rho+1)*(rho-1)^3 + n1^2*(rho-1)^2*(4*rho^(n1+1)-rho^2+10*rho-1)
           - 2*n1*(rho^2-1)*(6*rho^(n1+1)+rho^2-2*rho+1) + 2*(3*rho^3+10*rho^2+7*rho+4)*(rho^(n1+1)-rho) ) / 
           (n1*(n1-2)*(n1^2-1)*(rho-1)^4) )
  nprime1 <- 2 / (1-b1*(n1-2)/n1)
  
  se <- sd.y * sqrt((c0/b0) + (c1/b1)) 
  DF <- nprime0 + nprime1 - 4
  sd <- sd.y
  t <- slp.d / se
  pval <- 2*pt(-abs(t), df=DF)
  pval.lo <- pt(t, df=DF)
  pval.up <- pt(-t, df=DF)
  
  # Return summary
  strSum <- cbind(DF, rho, slp.d,  sd, se, t, pval, pval.lo, pval.up)
  return(strSum)
}

#==========================================================================================================
#----------------------------------------------- Table 2 ------------------------------------------------
#==========================================================================================================


#--------------------------------- Table 2 data from Patient #1390  -----------------------------------
# These are the raw data in Table 2 

Delay <- unique(dat.tab$DELAY)
Pre.Y <- dat.tab[PATIENT==1390,]$Y0
Post.Y <- dat.tab[PATIENT==1390,]$Y1
Diff <- Pre.Y - Post.Y

Table2.data <- rbind(Delay, Pre.Y, Post.Y, Diff)

#--------------------------------- Table 2 data from Patient #1390  -----------------------------------
# We conduct tests for level- and rate-change with the st and str tests, each at alpha = .05/2.
# Based on Landes & Tang (in preparation), we may also consider pairing the discounting data
# from the 2 time periods. Hence, we perform level- and rate-change tests with pst and pstr. 
# 2-sample tests for level- and rate-change
st.1390 <- st.fxn(Pre.Y, Post.Y, rhotype="WF", rhow=0)
str.1390 <- str.fxn(Pre.Y, Post.Y, rhotype="WF", rhow=0)

# paired tests for level- and rate-change
pst.1390 <- pst.fxn(Diff, rhotype="WF", rhow=0)
pstr.1390 <- pstr.fxn(Diff,  rhotype="WF", rhow=0)


Table2.sumstat <- round( rbind( c( st.1390[c(4,2)], str.1390[c(4,2)] ),
                                c( pst.1390[c(4,2)], pstr.1390[c(4,2)] )), 3)
colnames(Table2.sumstat) <- c("Level.s", "Level.r", "Rate.s", "Rate.r")

Table2.data
Table2.sumstat



#==========================================================================================================
#----------------------------------------------- Table 3 ------------------------------------------------
#==========================================================================================================


#------------------------ Reproducing Table 1 of Landes et al. (2012) -------------------------------------

# The delays on which the indifference points are (nonlinearly) regressed.
X <- rep(unique(dat.tab$DELAY),2)
# Indicator of post-treatment indifference points.
POST <- sort( rep(c(0,1),8) )
# Model Formula  - from Equation 4 in Landes et al. (2010).
mazur1 <- formula( Y ~ 1/( 1 + X*exp(POST*ln.GAMMA + ln.K0) ) )

patients <- unique( dat.tab[dat.tab$BAD_DATA == 0,]$PATIENT )

regr.results <- data.frame()
for( patient in patients){ 
  Y <- c( dat.tab[PATIENT==patient,]$Y0, dat.tab[PATIENT==patient,]$Y1 )
  tmp.dat <- as.data.frame(cbind(X, POST, Y)) 
  colnames(tmp.dat) <- c("X","POST","Y")
  
  start1 <- list(ln.K0 = log( 1/60), ln.GAMMA = 0 )
  fit1 <- summary(nls(mazur1, data=tmp.dat, start= start1)) 
  pval <- fit1$coefficients[2,4]
  signif <- (pval<=.05)
  regr.results <- rbind(regr.results, c(patient, pval, signif))
}
colnames(regr.results) <- c("PATIENT", "p.value","signif")

Table3.orig.results <- table(regr.results$signif)


#--------- Re-analysis of 119 patients using serial t-tests ---------------------
# Dataset to hold results
RESULTS <- data.table()

# Loop through each individual person
for(SUBJ in patients){
  tmp.dat1 <- dat.tab[.(SUBJ)]  # Data 
  
  st.Sum <- as.data.table(st.fxn(v1=tmp.dat1$Y1, v0=tmp.dat1$Y0, rhotype="WF"))  # 2-sample level
  tmp.st <- c( st.Sum[[2]], st.Sum[[3]], st.Sum[[4]], st.Sum[[7]])
  
  str.Sum <- as.data.table(str.fxn(v1=tmp.dat1$Y1, v0=tmp.dat1$Y0, rhotype="WF"))  # 2-sample rate
  tmp.str <- c( str.Sum[[2]], str.Sum[[3]], str.Sum[[4]], str.Sum[[7]])
  
  #tmp.2samp <- 10*(tmp.str[4]<=.025) + (tmp.st[4]<=.025)
  
  pst.Sum <- as.data.table(pst.fxn(v1=tmp.dat1$Y1, v0=tmp.dat1$Y0, rhotype="WF"))  # Paired level
  tmp.pst <- c( pst.Sum[[2]], pst.Sum[[3]], pst.Sum[[4]], pst.Sum[[7]])
  
  pstr.Sum <- as.data.table(pstr.fxn(v1=tmp.dat1$Y1, v0=tmp.dat1$Y0, rhotype="WF"))  # Paired rate
  tmp.pstr <- c( pstr.Sum[[2]], pstr.Sum[[3]], pstr.Sum[[4]], pstr.Sum[[7]])
  
  #tmp.pair <- 10*(tmp.pstr[4]<=.025) + (tmp.pst[4]<=.025)
  
  # Append results to dataset
  tmp <- cbind(SUBJ, t(tmp.st), t(tmp.str), t(tmp.pst), t(tmp.pstr))
  colnames(tmp) <- c("SUBJ", "st.rho", "st.diff", "st.sd", "st.p", 
                             "str.rho", "str.diff", "str.sd", "str.p" , 
                             "pst.rho", "pst.diff", "pst.sd", "pst.p" , 
                             "pstr.rho", "pstr.diff", "pstr.sd", "pstr.p")
  RESULTS <- rbind(RESULTS, tmp)
}

table(RESULTS$signif.2samp)

# Summarizing serial t-test results for Table 3
row1 <- c(sum(RESULTS$st.p<=0.025) +  sum(RESULTS$str.p<=0.025) - sum(RESULTS$st.p<=0.025 & RESULTS$str.p<=0.025 ), NA, NA, NA)
row2 <- c(sum(RESULTS$st.p<=0.025), round( quantile(RESULTS$st.rho, probs=c(0.5, 0.25, 0.75)), 2))
row3 <- c(sum(RESULTS$str.p<=0.025), round( quantile(RESULTS$str.rho, probs=c(0.5, 0.25, 0.75)), 2))
row4 <- c(sum(RESULTS$pst.p<=0.025) +  sum(RESULTS$pstr.p<=0.025) - sum(RESULTS$pst.p<=0.025 & RESULTS$pstr.p<=0.025 ), NA, NA, NA)
row5 <- c(sum(RESULTS$pst.p<=0.025), round( quantile(RESULTS$pst.rho, probs=c(0.5, 0.25, 0.75)), 2))
row6 <- c(sum(RESULTS$pstr.p<=0.025), round( quantile(RESULTS$pstr.rho, probs=c(0.5, 0.25, 0.75)), 2))

# Building the serial t-test results in Table 3
Table3.serial.t <- rbind(row1, row2, row3, row4, row5, row6)
colnames(Table3.serial.t) <- c("Number significant", "Median r", "25th percentile r", "75th percentile r")
rownames(Table3.serial.t) <- c("2-sample serial combined", "2-sample serial t for level-change", "2-sample serial t for rate-change",
                      "Paired serial combined", "Paired serial t for level-change", "Paired serial t for rate-change")

Table3.orig.results
Table3.serial.t




#------------------------ Returning to Table 2 -------------------------------------

#--- Original result for Patient 1390
Y <- c( dat.tab[PATIENT==1390,]$Y0, dat.tab[PATIENT==1390,]$Y1 )
tmp.dat <- as.data.frame(cbind(X, POST, Y)) 
colnames(tmp.dat) <- c("X","POST","Y")

start1 <- list(ln.K0 = log( 1/60), ln.GAMMA = 0 )
fit1 <- summary(nls(mazur1, data=tmp.dat, start= start1)) 
orig.1390 <- fit1$coefficients[2,]
orig.1390

#--- Results from serial t-tests for Patient 1390
st.1390 # 2-sample level-change
str.1390 # 2-sample rate-change
pst.1390 # paired level-change
pstr.1390 # paired rate-change
