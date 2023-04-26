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
  registerdoMC(cores-2)
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

# Slightly modified version of rinvgauss function from statmod package
rinvgauss. <- function(n,mean=1,shape=NULL,dispersion=1) {
	if(!is.null(shape)){dispersion <- 1/shape}
	mu <- rep_len(mean,n)
	phi <- rep_len(dispersion,n)
	r <- rep_len(0,n)
	i <- (mu>0 & phi>0)
	if(!all(i)){
		r[!i] <- NA
		n <- sum(i)
	}
	phi[i] <- phi[i]*mu[i]
	Y <- rchisq(n,df=1)
	X1 <- 1+phi[i]/2*(Y-sqrt(4*Y/phi[i]+Y^2))
	# Note: The line above should yield all X1>0, but it occasionally doesn't due to
	#		numerical precision issues.  The line below detects this and recomputes
	#		the relevant elements of X1 using a 2nd-order Taylor expansion of the 
	#		sqrt function, which is a good approximation whenever the problem occurs.
	if(any(X1<=0)){X1[X1<=0] <- (1/(Y*phi[i]))[X1<=0]}
	firstroot <- as.logical(rbinom(n,size=1L,prob=1/(1+X1)))
	r[i][firstroot] <- X1[firstroot]
	r[i][!firstroot] <- 1/X1[!firstroot]
	mu*r
}

# Draw from inverse-Gaussian distribution while avoiding potential numerical problems
rinvgaussian <- function(n,m,l){
	m. <- m/sqrt(m*l)
	l. <- l/sqrt(m*l)
	sqrt(m*l)*rinvgauss.(n,m.,l.)
}











#%%%%%%%%%%%%%%%%%%%%% (1) SPIKE-AND-SLAB %%%%%%%%%%%%%%%%%%%%%%%


# Gibbs iteration functions for both spike-slab samplers
# Note: The versions have separate functions, as opposed to being different
#       options of the same function, since the latter would require checking any such
#       options every time the function is called, i.e., in every MCMC iteration.
# Note: The values of XTY, n, and p can obviously be calculated from X and Y, but they
#       are included as inputs to avoid recalculating them every time the function is
#       called, i.e., in every MCMC iteration.
iter.ss.original <- function(beta,sigma2,X,Y,XTY,n,p,w,kappa,zeta){
	w.t <- 1/(1+(1/w-1)*sqrt(kappa)*exp(beta^2*(1/kappa-1)/(2*sigma2*zeta)))
	d.tau <- zeta*(1+(kappa-1)*rbinom(p,1,w.t))
	
	### Efficient algorithm to generate beta when n < p
	u <- rnorm(p)*sqrt(d.tau)
	v <- X%*%u + rnorm(n)
	U <- chol(X%*%diag(d.tau)%*%t(X)+diag(rep(1,n)))
	ww <- backsolve(U,backsolve(t(U),((Y/sqrt(sigma2))-v),upper.tri=F))
	beta.new <- sqrt(sigma2)*(u + d.tau*(t(X)%*%ww))
  
  sigma2.new <- (sum((Y-drop(X%*%beta.new))^2)+sum(beta.new^2/d.tau))/rchisq(1,n+p-1)
	return(list(beta=beta.new,sigma2=sigma2.new))
}
iter.ss.fast     <- function(beta,sigma2,X,Y,XTY,n,p,w,kappa,zeta){
	w.t <- 1/(1+(1/w-1)*sqrt(kappa)*exp(beta^2*(1/kappa-1)/(2*sigma2*zeta)))
	d.tau <- zeta*(1+(kappa-1)*rbinom(p,1,w.t))
	
	u <- rnorm(p)*sqrt(d.tau)
	v <- X%*%u + rnorm(n)
	U <- chol(X%*%diag(d.tau)%*%t(X)+diag(rep(1,n)))
	ww <- backsolve(U,backsolve(t(U),-v,upper.tri=F))
	
	### Efficient algorithm to generate sigma when n < p
	beta.tilde <- t(X)%*%backsolve(U,backsolve(t(U),Y,upper.tri=F))
	beta.tilde <- d.tau*beta.tilde
	sigma2.new <- (sum(Y^2)-sum(XTY*beta.tilde))/rchisq(1,n-1)
	
	### Efficient algorithm to generate beta when n < p
	beta.new <- beta.tilde+sqrt(sigma2.new)*(u + d.tau*(t(X)%*%ww))
	
  return(list(beta=beta.new,sigma2=sigma2.new))
} 

# Run original and fast spike-slab samplers
run.ss <- function(X,Y,w,kappa,zeta,K,M,fast=F){	
	XTY <- drop(t(X)%*%Y)
	n <- dim(X)[1]
	p <- dim(X)[2]
	iter.ss <- get(paste("iter.ss.",ifelse(fast,"fast","original"),sep=""))
	foreach(chain = 0:(M-1), .packages = c('coda'), .combine = rbind) %dopar%{
		beta <- rep(1,p)
		sigma2 <- 1
		beta.chain <- matrix(NA,nrow=K,ncol=p)
		sigma2.chain <- rep(NA,K)
		
		## monitor chains
		start_time <- proc.time()
		for(k in 1:K){
				iter.result <- iter.ss(beta,sigma2,X,Y,XTY,n,p,w,kappa,zeta)
				beta <- iter.result$beta
				sigma2 <- iter.result$sigma2
				beta.chain[k,] <- beta
				sigma2.chain[k] <- sigma2
		}
		
		## get valus to return
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

#------------------------------------------------------------
############################################################
# Generate MCMC for n=50
############################################################

# n and p settings
n. <- 100 # repeat for n={50, 100, 150}
p. <- n.*c(seq(0.5, 0.9, 0.1), 1:5)


# Approximate pairwise correlation of covariates
r. <- 0.2

# Sparsity of "true" coefficients for generating data
s. <- 0.2

# Grid of settings for runs
runs <- cbind(n.,p.,r.,s.)
n.runs <- dim(runs)[1]

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
cores = detectCores()
if(Sys.info()["sysname"] == "Windows"){
  cl <- makeCluster(cores-1)
  registerDoParallel(cl) # for Windows
}else{
  registerdoMC(cores-2)
}
getDoParWorkers()

#------- Storage
time_10_ss_slow <- matrix(NA, nrow = nrow(runs), ncol = M) #nrow is no. of datasets
time_10_ss_fast <- matrix(NA, nrow = nrow(runs), ncol = M) 

eff_10_ss_slow  <- matrix(NA, nrow = nrow(runs), ncol = M) #nrow is no. of datasets
eff_10_ss_fast  <- matrix(NA, nrow = nrow(runs), ncol = M) 

lag1ac_ss_slow  <- matrix(NA, nrow = nrow(runs), ncol = M) #nrow is no. of datasets
lag1ac_ss_fast  <- matrix(NA, nrow = nrow(runs), ncol = M) #nrow is no. of datasets

#----Spike-and-Slab
set.seed(111)
for(i in 0:(n.runs-1)){
  n <- runs[i+1,1]
  p <- runs[i+1,2]
  r <- runs[i+1,3]
  s <- runs[i+1,4]
  
  # Hyperparameters for spike-and-slab
  w <- rep(1/2,p)
  kappa <- rep(100,p)
  zeta <- rep(0.01,p)
  
  Xvarhalf <- diag(sqrt(1-r),p)+matrix((sqrt(1+(p-1)*r)-sqrt(1-r))/p,nrow=p,ncol=p)
  X.raw <- matrix(rnorm(n*p),nrow=n,ncol=p)%*%Xvarhalf
  X <- scale(X.raw)*sqrt(n/(n-1))
  XTX <- t(X)%*%X
  Y.raw <- drop(X%*%beta.true()+noise.true())
  Y <- Y.raw-mean(Y.raw)
  XTY <- drop(t(X)%*%Y)
  YTY <- sum(Y^2)
  
  #--- original spike-and-slab
  ## ress_slow is 3BGs and ress_fast is 2BGs
  ress_slow <- run.ss(X,Y,w,kappa,zeta,K,M,fast=F)
  time_10_ss_slow[i+1, ] <- unlist(ress_slow[, "time"])
  eff_10_ss_slow[i+1, ]  <- unlist(ress_slow[, "eff"])
  lag1ac_ss_slow[i+1, ]  <- unlist(ress_slow[, "lag1ac"])
  
  
  #--- new spike-andslab
  ress_fast <- run.ss(X,Y,w,kappa,zeta,K,M,fast=T)
  time_10_ss_fast[i+1, ] <- unlist(ress_fast[, "time"])
  eff_10_ss_fast[i+1, ]  <- unlist(ress_fast[, "eff"])
  lag1ac_ss_fast[i+1, ]  <- unlist(ress_fast[, "lag1ac"])
  
  flush.console()
}

if(Sys.info()["sysname"] == "Windows") stopCluster(cl)
registerDoSEQ()






############################################################
# Generate plots 
############################################################

#-------------- AutoCorrplot


par(mfcol=c(1, 2))
plot.new(); grid(); par(new = TRUE)
plot(1:10,ss_fast_ac_10, pch=8,type="b", ylim=c(0,1),  col = "darkorange",
     cex=.9, cex.axis=.9, cex.lab=1, xaxt="n",
     xlab=expression(frac(italic(p), italic(n))),ylab="Autocorrelation")
axis(1, at=1:10, labels=npratio)
points(1:10,ss_slow_ac_10,pch=2,type="b", col = "dodgerblue")
legend(1.032*max(3.4),(1+.75)/2,xjust=1,yjust=0.2,pch=c(8, 2),
       c("2BG","3BG"),cex=.9)




#--------- Effective Size plot
npratio <- c(seq(0.5, 1, by = 0.1), seq(2, 5, by = 1))
ss_fast_eff_10 <- eff_10_ss_fast / (time_10_ss_fast * 0.9) #accounting for burn-in
ss_slow_eff_10 <- eff_10_ss_slow / (time_10_ss_slow * 0.9)

eff_fast <- log10(rowMeans(ss_fast_eff_10))
eff_slow <- log10(rowMeans(ss_slow_eff_10))

plot.new(); grid(); par(new = TRUE)
plot(1:10,eff_fast,pch=8,type="b",ylim=c(-3,5), yaxt = "n", col = "darkorange",
     cex=.9,cex.axis=.9,cex.lab=1,xaxt="n",
     xlab=expression(frac(italic(p), italic(n))),ylab="Efficiency (in log scale)")
axis(1, at=1:10, labels=npratio)
axis(2, at=-3:5, labels=c(-3:5))
points(1:10,eff_slow,pch=2,type="b", col = "dodgerblue")
legend(1.032*max(10),(.6 + 3.4),xjust=1,yjust=0.2,pch=c(8, 2),
       c("2BG","3BG"),cex=.9) 

dev.off()
#-------- Save
# Save multiple objects
#save(eff_fast, eff_slow, file = "data_50.RData")
save.image(file = "data_100.RData")

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

# Hyperparameters for spike-and-slab
w <- rep(1/2,p)
kappa <- rep(100,p)
zeta <- rep(0.005,p)

# Number of MCMC iterations per chain
K <- 1.8e4

# Number of separate chains to run
M <- 1


set.seed(111)
n.runs <- 1
for(i in 0:(n.runs-1)){

	res_ss_slow <- run.ss(X,Y,w,kappa,zeta,K,M,fast=F)
	res_ss_fast <- run.ss(X,Y,w,kappa,zeta,K,M,fast=T)
	
	flush.console()
}




############################################################
# Plotting stuff
############################################################

#-------------- AutoCorrplot
mcmc_chain_slow_ss <- as.mcmc(res_ss_slow$chn)
mcmc_chain_fast_ss <- as.mcmc(res_ss_fast$chn)

par(mfcol=c(1, 2))
plot.new(); grid(); par(new = TRUE)
traceplot(mcmc_chain_fast_ss, col = alpha("darkorange", 0.5))

d <- density(mcmc_chain_fast_ss)
plot.new(); grid(); par(new = TRUE)
plot(d, main="")
polygon(d, col = alpha("darkorange", 0.5), border="black")


plot.new(); grid(); par(new = TRUE)
traceplot(mcmc_chain_slow_ss, col = alpha("dodgerblue", 0.5))

d <- density(mcmc_chain_slow_ss)
plot.new(); grid(); par(new = TRUE)
plot(d, main="")
polygon(d, col = alpha("dodgerblue", 0.5), border="black")


#-------- Time
time_fast <- res_ss_fast$time*0.9
time_slow <- res_ss_slow$time*0.9

#--------- Effective Size plot
eff_fast <- res_ss_fast$eff
eff_slow <- res_ss_slow$eff

ratio_fast <- eff_fast/time_fast
ratio_slow <- eff_slow/time_slow

#-------- Lag1ac
lag1ac_fast <- res_ss_fast$lag1ac
lag1ac_slow <- res_ss_slow$lag1ac







