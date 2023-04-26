############################################################
# Functions for frequentist lasso estimate,
# original Bayesian lasso, fast Bayesian lasso,
# original spike-slab sampler, and fast spike-slab sampler
############################################################

#--------- Load, otherwise install and load packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(lars, tidyverse, foreach, coda, glmnet, robustHD)


#--------- Check for OS and install/load appropriate parallel package
ifelse(Sys.info()["sysname"] == "Windows",
       pacman::p_load(doParallel), # for Windows
       pacman::p_load(doMC))       # for Mac/Linux


#--------- Check for OS and setup appropriate parallel system
cores = detectCores()
if(Sys.info()["sysname"] == "Windows"){
  cl <- makeCluster(cores-2)
  registerDoParallel(cl) # for Windows
}else{
  registerDoMC(cores-2)
}
getDoParWorkers()




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Compute lasso estimator
lasso.est <- function(x,y,lambda){
  fit.lars <- lars(x,y,type="lasso",normalize=F,intercept=F)
  predict(fit.lars,type="coefficients",mode="lambda",s=lambda)$coefficients
}

# Get a value of lambda that yields a specified number of nonzero coefficient estimates
get.lambda <- function(x,y,number.nz){
  fit.lars <- lars(x,y,type="lasso",normalize=F,intercept=F)
  path.segment <- rev(which(fit.lars$df==number.nz))[1]
  lambda.knot.seq <- c(fit.lars$lambda,0)
  if(is.na(path.segment)){
    stop("target number of nonzero coefficient estimates cannot be attained")
  }
  if(path.segment==1){
    return(2*lambda.knot.seq[1])
  }else{
    return((lambda.knot.seq[path.segment-1]+lambda.knot.seq[path.segment])/2)
  }
}

# Get all values of lambda that yield a specified number of nonzero coefficient estimates
get.lambda.all <- function(x,y,number.nz){
  fit.lars <- lars(x,y,type="lasso",normalize=F,intercept=F)
  path.segments <- which(fit.lars$df==number.nz)
  lambda.knot.seq <- c(Inf,fit.lars$lambda,0)
  num.segments <- length(path.segments)
  result <- matrix(NA,nrow=num.segments,ncol=2)
  for(i in 1:num.segments){
    segment <- path.segments[i]
    result[i,] <- lambda.knot.seq[c(segment+1,segment)]
  }
  result[,1] <- rev(result[,1])
  result[,2] <- rev(result[,2])
  colnames(result) <- c("lwr","upr")
  return(result)
}









#%%%%%%%%%%%%%%%%%%%%%%%% (2) BAYESIAN LASSO %%%%%%%%%%%%%%%%%%%%%%

############################################################
# Generate MCMC output for Figure 2, which is also included
# (with additional output) in left panel of Figure 3
############################################################

#%%%%%%%%%%%%%%%%%%%%% BAYESIAN LASSO %%%%%%%%%%%%%%%%%%%%%%%




# Gibbs iteration functions for both Bayesian lassos
# Note: The versions have separate functions, as opposed to being different
#       options of the same function, since the latter would require checking any such
#       options every time the function is called, i.e., in every MCMC iteration.
# Note: The values of XTY, n, and p can obviously be calculated from X and Y, but they
#       are included as inputs to avoid recalculating them every time the function is
#       called, i.e., in every MCMC iteration.
iter.bl.original <- function(beta,sigma2,X,Y,XTY,n,p,lambda){
  d.tau.inv <- rinvgaussian(p,sqrt(lambda^2*sigma2/beta^2),lambda^2)
  d.tau <- as.vector(1/d.tau.inv)
  
  ### Efficient algorithm to generate beta when n < p
  u <- rnorm(p)*sqrt(d.tau)
  v <- X%*%u + rnorm(n)
  U <- chol(X%*%diag(d.tau)%*%t(X)+diag(rep(1,n)))
  w <- backsolve(U,backsolve(t(U),((Y/sqrt(sigma2))-v),upper.tri=F))
  beta.new <- sqrt(sigma2)*(u + d.tau*(t(X)%*%w))
  
  sigma2.new <- (sum((Y-drop(X%*%beta.new))^2)+
                   sum(beta.new^2*d.tau.inv))/rchisq(1,n+p-1)
  return(list(beta=beta.new,sigma2=sigma2.new))
}

iter.bl.fast <- function(beta,sigma2,X,Y,XTY,n,p,lambda){
  d.tau.inv <- rinvgaussian(p,sqrt(lambda^2*sigma2/beta^2),lambda^2)
  d.tau <- as.vector(1/d.tau.inv)
  
  u <- rnorm(p)*sqrt(d.tau)
  v <- X%*%u + rnorm(n)
  U <- chol(X%*%diag(d.tau)%*%t(X)+diag(rep(1,n)))
  w <- backsolve(U,backsolve(t(U),-v,upper.tri=F))
  
  ### Efficient algorithm to generate sigma when n < p
  beta.tilde <- t(X)%*%backsolve(U,backsolve(t(U),Y,upper.tri=F))
  beta.tilde <- d.tau*beta.tilde
  sigma2.new <- (sum(Y^2)-sum(XTY*beta.tilde))/rchisq(1,n-1)
  
  ### Efficient algorithm to generate beta when n < p
  beta.new <- beta.tilde+sqrt(sigma2.new)*(u + d.tau*(t(X)%*%w))
  
  return(list(beta=beta.new,sigma2=sigma2.new))
} 

# Run original and fast Bayesian lassos
run.bl <- function(X,Y,lambda,K,M,fast=F){
  XTY <- drop(t(X)%*%Y)
  n <- dim(X)[1]
  #print(n)
  p <- dim(X)[2]
  #print(p)
  iter.bl <- get(paste("iter.bl.",ifelse(fast,"fast","original"),sep=""))
  foreach(chain = 0:(M-1), .packages = c('coda'), .combine = rbind) %dopar%{
    source("utils.R") # fix for helper function
    beta <- rep(1,p) #init beta
    sigma2 <- 1      #init sigma
    beta.chain <- matrix(NA,nrow=K,ncol=p)
    sigma2.chain <- rep(NA,K)
    
    ## monitor chains
    start_time <- proc.time()
    for(k in 1:K){
      iter.result <- iter.bl(beta,sigma2,X,Y,XTY,n,p,lambda)
      beta <- iter.result$beta
      sigma2 <- iter.result$sigma2
      beta.chain[k,] <- beta
      sigma2.chain[k] <- sigma2
    }
    
    ## get values to return
    time = as.double(proc.time() - start_time)[3]                #time elapsed for chain
    lag1ac = autocorr(as.mcmc(sigma2.chain[-c(1:(K * 0.1))]))[2] #chop of 10% as burn-in
    eff = effectiveSize(as.mcmc(sigma2.chain[-c(1: (K * 0.1))])) #get effective sample size post burnin
    print(paste("chain",chain+1,"of",M,"complete at",date()))
    flush.console()
    
    return(list(chn    = sigma2.chain,
                time   = time,
                lag1ac = lag1ac,
                eff    = eff))
  }
}





# n and p settings
n. <- 75
p. <- n.*c(seq(0.5, 0.9, 0.1), 1:5)

# Approximate pairwise correlation of covariates
r. <- 0.2

# Sparsity of "true" coefficients for generating data
s. <- 0.2

# Grid of settings for runs
runs <- cbind(n.,p.,r.,s.)
n.runs <- dim(runs)[1]

# Regularization parameter
lambda <- 1

# Number of MCMC iterations per chain
K <- 1.5e4

# Number of separate chains to run
M <- 6

# "True" beta and noise (for generating Y)
beta.true <- function(){
  temp <- rep(0,p)
  p.nz <- ceiling(p*s)
  if(p.nz>0){
    temp[1:p.nz] <- rt(p.nz,2)
  }
  return(temp)
}
noise.true <- function(){rt(n,4)}

#--------- Check for OS and setup appropriate parallel system
getDoParWorkers()

#------- Storage
time_10_bl_slow <- matrix(NA, nrow = nrow(runs), ncol = M) #nrow is no. of datasets
time_10_bl_fast <- matrix(NA, nrow = nrow(runs), ncol = M) 

eff_10_bl_slow  <- matrix(NA, nrow = nrow(runs), ncol = M) #nrow is no. of datasets
eff_10_bl_fast  <- matrix(NA, nrow = nrow(runs), ncol = M) 

lag1ac_bl_slow  <- matrix(NA, nrow = nrow(runs), ncol = M) #nrow is no. of datasets
lag1ac_bl_fast  <- matrix(NA, nrow = nrow(runs), ncol = M) 

#----Spike-and-Slab
set.seed(511)
for(i in 0:(n.runs-1)){
  n <- runs[i+1,1]
  p <- runs[i+1,2]
  r <- runs[i+1,3]
  s <- runs[i+1,4]
  
  Xvarhalf <- diag(sqrt(1-r),p)+matrix((sqrt(1+(p-1)*r)-sqrt(1-r))/p,nrow=p,ncol=p)
  X.raw <- matrix(rnorm(n*p),nrow=n,ncol=p)%*%Xvarhalf
  X <- scale(X.raw)*sqrt(n/(n-1))
  X <- matrix(as.vector(X),n,p)
  XTX <- t(X)%*%X
  Y.raw <- drop(X%*%beta.true()+noise.true())
  Y <- Y.raw-mean(Y.raw)
  XTY <- drop(t(X)%*%Y)
  YTY <- sum(Y^2)
  
  
  #--- original blasso
  ## ress_slow is 3BGs and ress_fast is 2BGs
  ress_slow <- run.bl(X,Y,lambda,K,M,fast=F)
  time_10_bl_slow[i+1, ] <- unlist(ress_slow[, "time"])
  eff_10_bl_slow[i+1, ]  <- unlist(ress_slow[, "eff"])
  lag1ac_bl_slow[i+1, ]  <- unlist(ress_slow[, "lag1ac"])
  
  #--- original blasso
  ress_fast <- run.bl(X,Y,lambda,K,M,fast=T)
  time_10_bl_fast[i+1, ] <- unlist(ress_fast[, "time"])
  eff_10_bl_fast[i+1, ]  <- unlist(ress_fast[, "eff"])
  lag1ac_bl_fast[i+1, ]  <- unlist(ress_fast[, "lag1ac"])
  
  flush.console()
}

if(Sys.info()["sysname"] == "Windows") stopCluster(cl)



############################################################
# Generate plots for Blasso
############################################################

#-------------- AutoCorrplot
npratio <- c(seq(0.5, 1, by = 0.1), seq(2, 5, by = 1))
bl_fast_ac_10 <- rowMeans(lag1ac_bl_fast)
bl_slow_ac_10 <- rowMeans(lag1ac_bl_slow)

par(mfcol=c(1, 2))
plot.new(); grid(); par(new = TRUE)
plot(1:10,bl_fast_ac_10, pch=8,type="b", ylim=c(0,1),  col = "darkorange",
     cex=.9, cex.axis=.9, cex.lab=1, xaxt="n",
     xlab=expression(frac(italic(p), italic(n))),ylab="Autocorrelation")
axis(1, at=1:10, labels=npratio)
points(1:10,bl_slow_ac_10,pch=2,type="b", col = "dodgerblue")
legend(1.032*max(3.4),(1+.75)/2,xjust=1,yjust=0.2,pch=c(8, 2),
       c("2BG","3BG"),cex=.9)




#--------- Effective Size plot
npratio <- c(seq(0.5, 1, by = 0.1), seq(2, 5, by = 1))
bl_fast_eff_10 <- eff_10_bl_fast / (time_10_bl_fast * 0.9) #accounting for burn-in
bl_slow_eff_10 <- eff_10_bl_slow / (time_10_bl_slow * 0.9)

eff_fast <- log10(rowMeans(bl_fast_eff_10))
eff_slow <- log10(rowMeans(bl_slow_eff_10))

plot.new(); grid(); par(new = TRUE)
plot(1:10,eff_fast,pch=8,type="b",ylim=c(-3,5), yaxt = "n", col = "darkorange",
     cex=.9,cex.axis=.9,cex.lab=1,xaxt="n",
     xlab=expression(frac(italic(p), italic(n))),ylab="Efficiency (in log scale)")
axis(1, at=1:10, labels=npratio)
axis(2, at=-3:5, labels=c(-3:5))
points(1:10,eff_slow,pch=2,type="b", col = "dodgerblue")
legend(1.032*max(10),(.6 + 3.4),xjust=1,yjust=0.2,pch=c(8, 2),
       c("2BG","3BG"),cex=.9) 


#-------- Save
# Save multiple objects
#save(eff_fast, eff_slow, file = "data_50.RData")
save.image(file = "data_75.RData")

rm(list = ls())
# Load the workspace again
#load("data_50.RData")








#<!-------- STOP ------------>



#%%%%%%%%%%%%%%%%% PROTEIN-GENE DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%
#---------- load package and data
library("robustHD")
data("nci60")  # contains matrices 'protein' and 'gene'

#----------- define response variable
y <- protein[, 92]

#-------- screen most correlated predictor variables
correlations <- apply(gene, 2, corHuber, y)
keep <- partialOrder(abs(correlations), 100, decreasing = TRUE)
x <- gene[, keep]
n <- dim(x)[1]
p <- dim(x)[2]
X <- scale(x)*sqrt(n/(n-1))
Y <- y-mean(y)

# Desired sparsity of lasso estimate (for choosing lambda)
sparsity <- 0.5

# Number of MCMC iterations per chain
K <- 1.8e4

# Number of separate chains to run
M <- 1


set.seed(511)
n.runs <- length(sparsity)
for(i in 0:(n.runs-1)){
  s <- sparsity[i+1]
  p.nz <- ceiling(min(n,p)*s)
  lambda <- get.lambda(X,Y,p.nz)
  
  res_bl_slow <- run.bl(X,Y,lambda,K,M,fast=F)
  res_bl_fast <- run.bl(X,Y,lambda,K,M,fast=T)
  
  flush.console()
}

############################################################
# Plotting stuff
############################################################

#-------------- AutoCorrplot
mcmc_chain_slow_bl <- as.mcmc(res_bl_slow$chn)
mcmc_chain_fast_bl <- as.mcmc(res_bl_fast$chn)

par(mfcol=c(1, 2))
plot.new(); grid(); par(new = TRUE)
traceplot(mcmc_chain_fast_bl, col = alpha("darkorange", 0.5))

d <- density(mcmc_chain_fast_bl)
plot.new(); grid(); par(new = TRUE)
plot(d, main="")
polygon(d, col = alpha("darkorange", 0.5), border="black")


plot.new(); grid(); par(new = TRUE)
traceplot(mcmc_chain_slow_bl, col = alpha("dodgerblue", 0.5))

d <- density(mcmc_chain_slow_bl)
plot.new(); grid(); par(new = TRUE)
plot(d, main="")
polygon(d, col = alpha("dodgerblue", 0.5), border="black")


#-------- Time
time_fast_bl <- res_bl_fast$time*0.9
time_slow_bl <- res_bl_slow$time*0.9

#--------- Effective Size plot
eff_fast_bl <- res_bl_fast$eff
eff_slow_bl <- res_bl_slow$eff

ratio_fast_bl <- eff_fast_bl/time_fast_bl
ratio_slow_bl <- eff_slow_bl/time_slow_bl

#-------- Lag1ac
lag1ac_fast_bl <- res_bl_fast$lag1ac
lag1ac_slow_bl <- res_bl_slow$lag1ac



