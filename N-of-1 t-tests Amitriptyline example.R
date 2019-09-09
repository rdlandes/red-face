#==========================================================================================================
#**********************************************************************************************************
# N-of-1 t-tests: Amitriptyline example
# Run time: less than 1 minute
#           on a Dell OptiPlex 990 with an 
#           Intel Core i7-2600 CPU @ 3.40 GHz
#           with 8.0 GB of RAM
#           R version 3.4.1
#
# The data come from Table 1 in 
# [20] Jaeschke R, Adachi J, Guyatt G, Keller J, Wong B. Clinical usefulness of amitriptyline in 
#      fibromyalgia: The results of 23 N-of-1 randomized controlled trials. The Journal of Rheumatology. 
#      1991;18(3):447-451. 

# A re-analysis of these data using a Bayesian hierarchical model is found in Table 1 of 
# [21] Zucker DR, Schmid CH, McIntosh MW, D'Agostino RB, Selker HP, Lau J. Combining single patient 
#      (N-of-1) trials to estimate population treatment effects and to evaluate individual patient 
#      responses to treatment. Journal of Clinical Epidemiology. 1997;50(4):401-410.  
#**********************************************************************************************************
#==========================================================================================================

#------------------------------- Needed functions -------------------------------

# Computes Fuller's rho estimate for input residuals
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

# Performs paired serial t-test for level-change on input DIFF = Y_A - Y_B vector
pstfxn <- function(diff) {
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
  rho <- WF.r(res)
  # In case serial correlation fails to compute, set all to NA
  if(is.na(rho)){
    d <- NA
    sd.d <- NA
  }
  # c, b, m', and se calculations
  level.c <- (n + 2*rho^(n+1) - n*rho^2 - 2*rho) / (n^2*(rho-1)^2)
  level.b <- ( n - (n+2*rho^(n+1)-n*rho^2-2*rho)/(n*(rho-1)^2) ) / (n-1)
  level.nprime <- 1 / (1 - (n-1)*level.b/n)
  pst.se <- sd.d * sqrt( level.c/level.b )
  # Serial t statistics and upper-tail p-value
  pst.t <- d / pst.se
  pst.pval <- pt(pst.t, df=level.nprime-1, lower.tail=FALSE)
  
  # Gathering summary statistics
  pstSum <- cbind(d, sd.d, rho, pst.t, pst.pval)
  return(pstSum)
}

#----------------- Usual & Serial t-tests for 6 N-of-1 trials  ------------------
# Results in Table 1 of manuscript

Patient9 <- c(0.05,-0.22,0.57,0.36)
t.test(Patient9, alternative=("greater"))
pstfxn(Patient9)

Patient18 <- c(0.64,1.08,-0.36,0.79,-0.64,1.50)
t.test(Patient18, alternative=("greater"))
pstfxn(Patient18)

Patient23 <- c(1.22,1.07,-0.08,0.50)
t.test(Patient23, alternative=("greater"))
pstfxn(Patient23)

Patient17 <- c(-0.08,0.86,1.07,1.15)
t.test(Patient17, alternative=("greater"))
pstfxn(Patient17)

Patient15 <- c(0.86,1.43,0.65,1.86)
t.test(Patient15, alternative=("greater"))
pstfxn(Patient15)

Patient12 <- c(4.29,3.15,0.78,4.49)
t.test(Patient12, alternative=("greater"))
pstfxn(Patient12)
