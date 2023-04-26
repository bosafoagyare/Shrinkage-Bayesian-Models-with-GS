############################################################
# Functions for frequentist lasso estimate,
# original Bayesian lasso, fast Bayesian lasso,
# original spike-slab sampler, and fast spike-slab sampler
############################################################

#--------- Load, otherwise install and load packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(lars, tidyverse, foreach, coda, glmnet)


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







############################################################
# Generate plots 
############################################################

#-------------- AutoCorrplot
npratio <- c(seq(0.5, 1, by = 0.1), seq(2, 5, by = 1))
ss_fast_ac_10 <- rowMeans(lag1ac_ss_fast)
ss_slow_ac_10 <- rowMeans(lag1ac_ss_slow)

par(mfcol=c(1, 2))
plot.new(); grid(); par(new = TRUE)
plot(1:10,ss_fast_ac_10, pch=8,type="b", ylim=c(0,1),  col = "darkorange",
     cex=.9, cex.axis=.9, cex.lab=1, xaxt="n",
     xlab=expression(frac(italic(p), italic(n))),ylab="Autocorrelation")
axis(1, at=1:10, labels=npratio)
points(1:10,ss_slow_ac_10,pch=2,type="b", col = "dodgerblue")
legend(1.032*max(2.6),(1+.75)/2,xjust=1,yjust=0.2,pch=c(8, 2),
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


#-------- Save
# Save multiple objects
#save(eff_fast, eff_slow, file = "data_50.RData")
save.image(file = "data_100.RData")

rm(list = ls())
# Load the workspace again
#load("data_50.RData")


































































#<!-------- STOP ------------>













############################################################
# Generate MCMC output for Subsection 5.1
# (gene expression data)
############################################################

library(flare)
data(eyedata)
n <- dim(x)[1]
p <- dim(x)[2]
X <- scale(x)*sqrt(n/(n-1))
Y <- y-mean(y)

# Desired sparsity of lasso estimate (for choosing lambda)
sparsity <- 0.5

# Number of MCMC iterations per chain
K <- 20000

# Number of separate chains to run
M <- 1

# Random seed (digits 10-12 of pi)
set.seed(589)

timestamp <- format(Sys.time(),"%Y%m%d%H%M%S")
n.runs <- length(sparsity)
dir.create(timestamp)
setwd(timestamp)
print(paste("Writing output to",timestamp))
flush.console()
rundigits <- max(1,ceiling(log(n.runs,10)))  # digits needed for run label strings
for(i in 0:(n.runs-1)){
	s <- sparsity[i+1]
	p.nz <- ceiling(min(n,p)*s)
	lambda <- get.lambda(X,Y,p.nz)
	runtext <- substring(format(i/(10^rundigits),nsmall=rundigits),3)
	outfile.stem.orig <- paste("eye-orig-",runtext,sep="")
	run.bl(X,Y,lambda,K,M,outfile.stem.orig,fast=F,keep.beta=F,write.each=T)
	outfile.stem.fast <- paste("eye-fast-",runtext,sep="")
	run.bl(X,Y,lambda,K,M,outfile.stem.fast,fast=T,keep.beta=F,write.each=T)
	print(paste("Run",i+1,"of",n.runs,"complete at",date()))
	flush.console()
}
setwd("..")



############################################################
# Calculate results for Subsection 5.1
# (gene expression data)
# using MCMC output generated by code above
############################################################

outputdir <- timestamp

setwd(outputdir)
orig.sigma2.chain <- read.table("eye-orig-0-0-s.txt")
fast.sigma2.chain <- read.table("eye-fast-0-0-s.txt")
setwd("..")

# Fraction of output to discard as burn-in
burn <- 1/11

K.burn <- round(K*burn)
orig.sigma2.chain <- orig.sigma2.chain[-(1:K.burn),]
fast.sigma2.chain <- fast.sigma2.chain[-(1:K.burn),]

library(coda)
orig.sigma2.mcmc <- mcmc(orig.sigma2.chain)
fast.sigma2.mcmc <- mcmc(fast.sigma2.chain)
print(acf(orig.sigma2.chain,1,plot=F)$acf[2,1,1])
print(acf(fast.sigma2.chain,1,plot=F)$acf[2,1,1])
print(effectiveSize(orig.sigma2.mcmc))
print(effectiveSize(fast.sigma2.mcmc))
#print(gelman.diag(orig.sigma2.mcmc))
#print(gelman.diag(fast.sigma2.mcmc))
#gelman.plot(orig.sigma2.mcmc)
#gelman.plot(fast.sigma2.mcmc)
plot.new()
grid()
par(new = TRUE)
traceplot(orig.sigma2.mcmc)
traceplot(fast.sigma2.mcmc)
data.frame(chain = orig.sigma2.chain) %>%
ggplot(aes(x = chain)) +
  geom_histogram(fill = "dodgerblue",
                 color = "black",
                 alpha = 0.3) + 
  theme_bw()

data.frame(chain = orig.sigma2.chain) %>%
ggplot(aes(x = chain)) +
  geom_density(fill = "dodgerblue",
               alpha = 0.3) + 
  theme_bw()

data.frame(chain = orig.sigma2.chain) %>%
ggplot(aes(x = chain)) +
  geom_boxplot(fill = "dodgerblue", alpha = 0.3) +
  theme_bw()




############################################################
# Generate MCMC output for Subsection 5.2
# (infrared spectroscopy data)
############################################################

library(ppls)
data(cookie)
x <- cookie[1:40,1:700]
y <- cookie[1:40,701]
n <- dim(x)[1]
p <- dim(x)[2]
X <- scale(x)*sqrt(n/(n-1))
Y <- y-mean(y)

# Desired sparsity of lasso estimate (for choosing lambda)
sparsity <- 0.5

# Number of MCMC iterations per chain
K <- 20000

# Number of separate chains to run
M <- 1

# Random seed (digits 13-15 of pi)
set.seed(793)

timestamp <- format(Sys.time(),"%Y%m%d%H%M%S")
n.runs <- length(sparsity)
dir.create(timestamp)
setwd(timestamp)
print(paste("Writing output to",timestamp))
flush.console()
rundigits <- max(1,ceiling(log(n.runs,10)))  # digits needed for run label strings
for(i in 0:(n.runs-1)){
	s <- sparsity[i+1]
	p.nz <- ceiling(min(n,p)*s)
	lambda.intervals <- get.lambda.all(X,Y,p.nz)
	lambda <- mean(lambda.intervals[dim(lambda.intervals)[1],])
	runtext <- substring(format(i/(10^rundigits),nsmall=rundigits),3)
	outfile.stem.orig <- paste("cookie-orig-",runtext,sep="")
	run.bl(X,Y,lambda,K,M,outfile.stem.orig,fast=F,keep.beta=F,write.each=T)
	outfile.stem.fast <- paste("cookie-fast-",runtext,sep="")
	run.bl(X,Y,lambda,K,M,outfile.stem.fast,fast=T,keep.beta=F,write.each=T)
	print(paste("Run",i+1,"of",n.runs,"complete at",date()))
	flush.console()
}
setwd("..")



############################################################
# Calculate results for Subsection 5.2
# (infrared spectroscopy data)
# using MCMC output generated by code above
############################################################

outputdir <- timestamp

setwd(outputdir)
orig.sigma2.chain <- read.table("cookie-orig-0-0-s.txt")
fast.sigma2.chain <- read.table("cookie-fast-0-0-s.txt")
setwd("..")

# Fraction of output to discard as burn-in
burn <- 1/11

K.burn <- round(K*burn)
orig.sigma2.chain <- orig.sigma2.chain[-(1:K.burn),]
fast.sigma2.chain <- fast.sigma2.chain[-(1:K.burn),]

library(coda)
orig.sigma2.mcmc <- mcmc(orig.sigma2.chain)
fast.sigma2.mcmc <- mcmc(fast.sigma2.chain)
print(acf(orig.sigma2.chain,1,plot=F)$acf[2,1,1])
print(acf(fast.sigma2.chain,1,plot=F)$acf[2,1,1])
print(effectiveSize(orig.sigma2.mcmc))
print(effectiveSize(fast.sigma2.mcmc))



############################################################
# Generate MCMC output for Subsection 5.3
# (communities and crime data)
############################################################

communities.url1 <- "http://archive.ics.uci.edu/ml/machine-learning-databases/"
communities.url2 <- "communities/communities.data"
communities.url <- paste(communities.url1,communities.url2,sep="")
communities <- read.csv(communities.url,header=F)

p.comm <- dim(communities)[2]
num.missing <- rep(NA,p.comm)
for(j in 1:p.comm){
	num.missing[j] <- sum(communities[,j]=="?")
}

# Columns 1:5 are non-predictive
# Columns 6:30, 32:101, 119:121, 126 are predictive with no missing values
# Column 31 is predictive with 1 missing value
# Columns 102:118, 122:125, 127 are predictive with 1675 missing values
# Column 128 is the response with no missing values

n.all <- dim(communities)[1]
to.include.1 <- c(6:30,32:101,119:121,126)
to.include.2 <- c()
for(j in 1:p.comm){
	if(is.element(j,to.include.1)){
		if(max(table(communities[,j]))<n.all/2){
			to.include.2 <- append(to.include.2,j)
		}
	}
}
cor.X.Y <- rep(0,p.comm)
for(j in 1:p.comm){
	if(is.element(j,to.include.2)){
		cor.X.Y[j] <- cor(communities[,j],communities[,p.comm])
	}
}
to.include.3 <- sort(order(-abs(cor.X.Y))[1:50])

X.1.all <- communities[,to.include.3]
p.1 <- dim(X.1.all)[2]

# Random seed (digits 16-18 of pi)
set.seed(238)

# Choose random subsample of 10 observations for use
n <- 10
obs.include <- sort(sample(n.all,n))

X.1.raw <- X.1.all[obs.include,]
X.1 <- scale(X.1.raw)*sqrt(n/(n-1))
X.2.raw <- matrix(NA,nrow=n,ncol=p.1*(p.1+1)/2)
j.ctr <- 1
for(j1 in 1:p.1){
	for(j2 in j1:p.1){
		X.2.raw[,j.ctr] <- X.1[,j1]*X.1[,j2]
		j.ctr <- j.ctr+1
	}
}
X.2 <- scale(X.2.raw)*sqrt(n/(n-1))
X <- cbind(X.1,X.2)
p <- dim(X)[2]
Y.all <- communities[,p.comm]
Y.raw <- Y.all[obs.include]
Y <- Y.raw-mean(Y.raw)

# Desired sparsity of lasso estimate (for choosing lambda)
sparsity <- 0.5

# Number of MCMC iterations per chain
K <- 20000

# Number of separate chains to run
M <- 1

# Random seed (again) (digits 19-21 of pi)
set.seed(462)

timestamp <- format(Sys.time(),"%Y%m%d%H%M%S")
n.runs <- length(sparsity)
dir.create(timestamp)
setwd(timestamp)
print(paste("Writing output to",timestamp))
flush.console()
rundigits <- max(1,ceiling(log(n.runs,10)))  # digits needed for run label strings
for(i in 0:(n.runs-1)){
	s <- sparsity[i+1]
	p.nz <- ceiling(min(n,p)*s)
	lambda <- get.lambda(X,Y,p.nz)
	runtext <- substring(format(i/(10^rundigits),nsmall=rundigits),3)
	outfile.stem.orig <- paste("communities-orig-",runtext,sep="")
	run.bl(X,Y,lambda,K,M,outfile.stem.orig,fast=F,keep.beta=F,write.each=T)
	outfile.stem.fast <- paste("communities-fast-",runtext,sep="")
	run.bl(X,Y,lambda,K,M,outfile.stem.fast,fast=T,keep.beta=F,write.each=T)
	print(paste("Run",i+1,"of",n.runs,"complete at",date()))
	flush.console()
}
setwd("..")



############################################################
# Calculate results for Subsection 5.3
# (communities and crime data)
# using MCMC output generated by code above
############################################################

outputdir <- timestamp

setwd(outputdir)
orig.sigma2.chain <- read.table("communities-orig-0-0-s.txt")
fast.sigma2.chain <- read.table("communities-fast-0-0-s.txt")
setwd("..")

# Fraction of output to discard as burn-in
burn <- 1/11

K.burn <- round(K*burn)
orig.sigma2.chain <- orig.sigma2.chain[-(1:K.burn),]
fast.sigma2.chain <- fast.sigma2.chain[-(1:K.burn),]

library(coda)
orig.sigma2.mcmc <- mcmc(orig.sigma2.chain)
fast.sigma2.mcmc <- mcmc(fast.sigma2.chain)
print(acf(orig.sigma2.chain,1,plot=F)$acf[2,1,1])
print(acf(fast.sigma2.chain,1,plot=F)$acf[2,1,1])
print(effectiveSize(orig.sigma2.mcmc))
print(effectiveSize(fast.sigma2.mcmc))



############################################################
# Generate MCMC output for Subsection 5.1
# (gene expression data) for spike-and-slab
############################################################

library(flare)
data(eyedata)
n <- dim(x)[1]
p <- dim(x)[2]
X <- scale(x)*sqrt(n/(n-1))
Y <- y-mean(y)

# Hyperparameters for spike-and-slab
w <- rep(1/2,p)
kappa <- rep(100,p)
zeta <- rep(0.00002,p)

# Number of MCMC iterations per chain
K <- 20000

# Number of separate chains to run
M <- 1

# Random seed (digits 22-24 of pi)
set.seed(643)

timestamp <- format(Sys.time(),"%Y%m%d%H%M%S")
n.runs <- 1
dir.create(timestamp)
setwd(timestamp)
print(paste("Writing output to",timestamp))
flush.console()
rundigits <- max(1,ceiling(log(n.runs,10)))  # digits needed for run label strings
for(i in 0:(n.runs-1)){
	runtext <- substring(format(i/(10^rundigits),nsmall=rundigits),3)
	outfile.stem.orig <- paste("ss-eye-orig-",runtext,sep="")
	run.ss(X,Y,w,kappa,zeta,K,M,outfile.stem.orig,fast=F,keep.beta=F,write.each=T)
	outfile.stem.fast <- paste("ss-eye-fast-",runtext,sep="")
	run.ss(X,Y,w,kappa,zeta,K,M,outfile.stem.fast,fast=T,keep.beta=F,write.each=T)
	print(paste("Run",i+1,"of",n.runs,"complete at",date()))
	flush.console()
}
setwd("..")



############################################################
# Calculate results for Subsection 5.1 
# (gene expression data) for spike-and-slab
# using MCMC output generated by code above
############################################################

outputdir <- timestamp

setwd(outputdir)
orig.sigma2.chain <- read.table("ss-eye-orig-0-0-s.txt")
fast.sigma2.chain <- read.table("ss-eye-fast-0-0-s.txt")
setwd("..")

library(coda)
orig.sigma2.mcmc <- mcmc(orig.sigma2.chain)
fast.sigma2.mcmc <- mcmc(fast.sigma2.chain)
print(acf(orig.sigma2.chain,1,plot=F)$acf[2,1,1])
print(acf(fast.sigma2.chain,1,plot=F)$acf[2,1,1])
print(effectiveSize(orig.sigma2.mcmc))
print(effectiveSize(fast.sigma2.mcmc))



############################################################
# Generate MCMC output for left panel of Figure S1
############################################################

# n and p settings
n. <- rep(c(50,100,200),each=5)
p. <- n.*2*(1:5)

# Approximate pairwise correlation of covariates
r. <- 0.8

# Sparsity of "true" coefficients for generating data
s. <- 0.2

# Grid of settings for runs
runs <- cbind(n.,p.,r.,s.)
n.runs <- dim(runs)[1]

# Regularization parameter
lambda <- 1

# Number of MCMC iterations per chain
K <- 20000

# Number of separate chains to run
M <- 10

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

# Random seed (digits 28-30 of pi)
set.seed(279)

timestamp1 <- format(Sys.time(),"%Y%m%d%H%M%S")
dir.create(timestamp1)
setwd(timestamp1)
print(paste("Writing output to",timestamp))
flush.console()
rundigits <- max(1,ceiling(log(n.runs,10)))  # digits needed for run label strings
for(i in 0:(n.runs-1)){
	n <- runs[i+1,1]
	p <- runs[i+1,2]
	r <- runs[i+1,3]
	s <- runs[i+1,4]
	Xvarhalf <- diag(sqrt(1-r),p)+matrix((sqrt(1+(p-1)*r)-sqrt(1-r))/p,nrow=p,ncol=p)
	X.raw <- matrix(rnorm(n*p),nrow=n,ncol=p)%*%Xvarhalf
	X <- scale(X.raw)*sqrt(n/(n-1))
	XTX <- t(X)%*%X
	Y.raw <- drop(X%*%beta.true()+noise.true())
	Y <- Y.raw-mean(Y.raw)
	XTY <- drop(t(X)%*%Y)
	YTY <- sum(Y^2)
	runtext <- substring(format(i/(10^rundigits),nsmall=rundigits),3)
	outfile.stem.orig <- paste("acvp-orig-",runtext,sep="")
	run.bl(X,Y,lambda,K,M,outfile.stem.orig,fast=F,keep.beta=F,write.each=T)
	outfile.stem.fast <- paste("acvp-fast-",runtext,sep="")
	run.bl(X,Y,lambda,K,M,outfile.stem.fast,fast=T,keep.beta=F,write.each=T)
	print(paste("Run",i+1,"of",n.runs,"complete at",date()))
	flush.console()
}
setwd("..")

# n and p settings
n. <- rep(c(50,100,200),each=5)
p. <- n.*c(2,3,4,7,14)/5

# Approximate pairwise correlation of covariates
r. <- 0.8

# Sparsity of "true" coefficients for generating data
s. <- 0.2

# Grid of settings for runs
runs <- cbind(n.,p.,r.,s.)
n.runs <- dim(runs)[1]

# Regularization parameter
lambda <- 1

# Number of MCMC iterations per chain
K <- 20000

# Number of separate chains to run
M <- 10

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

# Random seed (digits 31-33 of pi)
set.seed(502)

timestamp2 <- format(Sys.time(),"%Y%m%d%H%M%S")
dir.create(timestamp2)
setwd(timestamp2)
print(paste("Writing output to",timestamp))
flush.console()
rundigits <- max(1,ceiling(log(n.runs,10)))  # digits needed for run label strings
for(i in 0:(n.runs-1)){
	n <- runs[i+1,1]
	p <- runs[i+1,2]
	r <- runs[i+1,3]
	s <- runs[i+1,4]
	Xvarhalf <- diag(sqrt(1-r),p)+matrix((sqrt(1+(p-1)*r)-sqrt(1-r))/p,nrow=p,ncol=p)
	X.raw <- matrix(rnorm(n*p),nrow=n,ncol=p)%*%Xvarhalf
	X <- scale(X.raw)*sqrt(n/(n-1))
	XTX <- t(X)%*%X
	Y.raw <- drop(X%*%beta.true()+noise.true())
	Y <- Y.raw-mean(Y.raw)
	XTY <- drop(t(X)%*%Y)
	YTY <- sum(Y^2)
	runtext <- substring(format(i/(10^rundigits),nsmall=rundigits),3)
	outfile.stem.orig <- paste("acvp-orig-",runtext,sep="")
	run.bl(X,Y,lambda,K,M,outfile.stem.orig,fast=F,keep.beta=F,write.each=T)
	outfile.stem.fast <- paste("acvp-fast-",runtext,sep="")
	run.bl(X,Y,lambda,K,M,outfile.stem.fast,fast=T,keep.beta=F,write.each=T)
	print(paste("Run",i+1,"of",n.runs,"complete at",date()))
	flush.console()
}
setwd("..")



############################################################
# Generate plot for left panel of Figure S1
# using MCMC output generated by code above
############################################################

outputdir1 <- timestamp1
outputdir2 <- timestamp2

# n and p settings
n.1 <- rep(c(50,100,200),each=5)
p.1 <- n.1*2*(1:5)
n.2 <- rep(c(50,100,200),each=5)
p.2 <- n.2*c(2,3,4,7,14)/5

# Approximate pairwise correlation of covariates
r. <- 0.8

# Sparsity of "true" coefficients for generating data
s. <- 0.2

# Grid of settings for runs
runs1 <- cbind(n.1,p.1,r.,s.)
runs2 <- cbind(n.2,p.2,r.,s.)
n.runs1 <- dim(runs1)[1]
n.runs2 <- dim(runs2)[1]

# Regularization parameter
lambda <- 1

# Number of MCMC iterations per chain
K <- 20000

# Number of separate chains to run
M <- 10

# Fraction of output to discard as burn-in
burn <- 1/11

K.burn <- round(K*burn)
ac.orig.single1 <- rep(NA,n.runs1*M)
ac.fast.single1 <- rep(NA,n.runs1*M)
ac.orig.single2 <- rep(NA,n.runs2*M)
ac.fast.single2 <- rep(NA,n.runs2*M)
chain.run.id1 <- rep(1:n.runs1,each=M)
chain.run.id2 <- rep(1:n.runs2,each=M)
run1digits <- max(1,ceiling(log(n.runs1,10)))	# digits needed for run label strings
run2digits <- max(1,ceiling(log(n.runs2,10)))	# digits needed for run label strings
chaindigits <- max(1,ceiling(log(M,10)))		# digits needed for chain label strings
setwd(outputdir1)
for(i in 0:(n.runs1-1)){
	runtext <- substring(format(i/(10^run1digits),nsmall=run1digits),3)
	for(j in 0:(M-1)){
		chaintext <- substring(format(j/(10^chaindigits),nsmall=chaindigits),3)
		filename.orig.ij <- paste("acvp-orig-",runtext,"-",chaintext,"-s.txt",sep="")
		sigma2.chain.orig.ij <- read.table(filename.orig.ij)[-(1:K.burn),1]
		ac.orig.single1[i*M+j+1] <- acf(sigma2.chain.orig.ij,lag.max=1,plot=F)$acf[2]
		filename.fast.ij <- paste("acvp-fast-",runtext,"-",chaintext,"-s.txt",sep="")
		sigma2.chain.fast.ij <- read.table(filename.fast.ij)[-(1:K.burn),1]
		ac.fast.single1[i*M+j+1] <- acf(sigma2.chain.fast.ij,lag.max=1,plot=F)$acf[2]
	}
}
setwd("..")
setwd(outputdir2)
for(i in 0:(n.runs2-1)){
	runtext <- substring(format(i/(10^run2digits),nsmall=run2digits),3)
	for(j in 0:(M-1)){
		chaintext <- substring(format(j/(10^chaindigits),nsmall=chaindigits),3)
		filename.orig.ij <- paste("acvp-orig-",runtext,"-",chaintext,"-s.txt",sep="")
		sigma2.chain.orig.ij <- read.table(filename.orig.ij)[-(1:K.burn),1]
		ac.orig.single2[i*M+j+1] <- acf(sigma2.chain.orig.ij,lag.max=1,plot=F)$acf[2]
		filename.fast.ij <- paste("acvp-fast-",runtext,"-",chaintext,"-s.txt",sep="")
		sigma2.chain.fast.ij <- read.table(filename.fast.ij)[-(1:K.burn),1]
		ac.fast.single2[i*M+j+1] <- acf(sigma2.chain.fast.ij,lag.max=1,plot=F)$acf[2]
	}
}
setwd("..")
ac.orig1 <- tapply(ac.orig.single1,chain.run.id1,mean)
ac.fast1 <- tapply(ac.fast.single1,chain.run.id1,mean)
ac.orig2 <- tapply(ac.orig.single2,chain.run.id2,mean)
ac.fast2 <- tapply(ac.fast.single2,chain.run.id2,mean)
n. <- c(n.1,n.2)
p. <- c(p.1,p.2)
ac.orig <- c(ac.orig1,ac.orig2)
ac.fast <- c(ac.fast1,ac.fast2)

plotfilename <- paste("acvpn-",timestamp1,"-",timestamp2,".pdf",sep="")
min.ac <- min(ac.orig,ac.fast,0)
pdf(plotfilename,width=5,height=5)
par(mar=c(4,4,1,1),mgp=c(2.75,1.00,0),pty="s")
plot(log(p./n.),ac.orig,pch=1,ylim=c(min.ac,1),cex=1.4,cex.axis=1.4,cex.lab=1.4,
	xlab=expression(log(italic(p)/italic(n))),ylab="Autocorrelation",xaxt="n")
axis(1,seq(-1,2,by=0.5),c("-1.0","","0.0","","1.0","","2.0"),cex.axis=1.4)
points(log(p./n.),ac.fast,pch=2)
legend("topleft",c("Original","Fast"),pch=c(1,2),cex=1.4)
dev.off()



############################################################
# Generate MCMC output for center panel of Figure S1
############################################################

# n and p settings
n. <- rep(c(50,100,200),each=5)
p. <- n.*2*(1:5)

# Approximate pairwise correlation of covariates
r. <- 0.2

# Sparsity of "true" coefficients for generating data
s. <- 0.8

# Grid of settings for runs
runs <- cbind(n.,p.,r.,s.)
n.runs <- dim(runs)[1]

# Regularization parameter
lambda <- 1

# Number of MCMC iterations per chain
K <- 20000

# Number of separate chains to run
M <- 10

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

# Random seed (digits 34-36 of pi)
set.seed(884)

timestamp1 <- format(Sys.time(),"%Y%m%d%H%M%S")
dir.create(timestamp1)
setwd(timestamp1)
print(paste("Writing output to",timestamp))
flush.console()
rundigits <- max(1,ceiling(log(n.runs,10)))  # digits needed for run label strings
for(i in 0:(n.runs-1)){
	n <- runs[i+1,1]
	p <- runs[i+1,2]
	r <- runs[i+1,3]
	s <- runs[i+1,4]
	Xvarhalf <- diag(sqrt(1-r),p)+matrix((sqrt(1+(p-1)*r)-sqrt(1-r))/p,nrow=p,ncol=p)
	X.raw <- matrix(rnorm(n*p),nrow=n,ncol=p)%*%Xvarhalf
	X <- scale(X.raw)*sqrt(n/(n-1))
	XTX <- t(X)%*%X
	Y.raw <- drop(X%*%beta.true()+noise.true())
	Y <- Y.raw-mean(Y.raw)
	XTY <- drop(t(X)%*%Y)
	YTY <- sum(Y^2)
	runtext <- substring(format(i/(10^rundigits),nsmall=rundigits),3)
	outfile.stem.orig <- paste("acvp-orig-",runtext,sep="")
	run.bl(X,Y,lambda,K,M,outfile.stem.orig,fast=F,keep.beta=F,write.each=T)
	outfile.stem.fast <- paste("acvp-fast-",runtext,sep="")
	run.bl(X,Y,lambda,K,M,outfile.stem.fast,fast=T,keep.beta=F,write.each=T)
	print(paste("Run",i+1,"of",n.runs,"complete at",date()))
	flush.console()
}
setwd("..")

# n and p settings
n. <- rep(c(50,100,200),each=5)
p. <- n.*c(2,3,4,7,14)/5

# Approximate pairwise correlation of covariates
r. <- 0.2

# Sparsity of "true" coefficients for generating data
s. <- 0.8

# Grid of settings for runs
runs <- cbind(n.,p.,r.,s.)
n.runs <- dim(runs)[1]

# Regularization parameter
lambda <- 1

# Number of MCMC iterations per chain
K <- 20000

# Number of separate chains to run
M <- 10

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

# Random seed (digits 37-39 of pi)
set.seed(197)

timestamp2 <- format(Sys.time(),"%Y%m%d%H%M%S")
dir.create(timestamp2)
setwd(timestamp2)
print(paste("Writing output to",timestamp))
flush.console()
rundigits <- max(1,ceiling(log(n.runs,10)))  # digits needed for run label strings
for(i in 0:(n.runs-1)){
	n <- runs[i+1,1]
	p <- runs[i+1,2]
	r <- runs[i+1,3]
	s <- runs[i+1,4]
	Xvarhalf <- diag(sqrt(1-r),p)+matrix((sqrt(1+(p-1)*r)-sqrt(1-r))/p,nrow=p,ncol=p)
	X.raw <- matrix(rnorm(n*p),nrow=n,ncol=p)%*%Xvarhalf
	X <- scale(X.raw)*sqrt(n/(n-1))
	XTX <- t(X)%*%X
	Y.raw <- drop(X%*%beta.true()+noise.true())
	Y <- Y.raw-mean(Y.raw)
	XTY <- drop(t(X)%*%Y)
	YTY <- sum(Y^2)
	runtext <- substring(format(i/(10^rundigits),nsmall=rundigits),3)
	outfile.stem.orig <- paste("acvp-orig-",runtext,sep="")
	run.bl(X,Y,lambda,K,M,outfile.stem.orig,fast=F,keep.beta=F,write.each=T)
	outfile.stem.fast <- paste("acvp-fast-",runtext,sep="")
	run.bl(X,Y,lambda,K,M,outfile.stem.fast,fast=T,keep.beta=F,write.each=T)
	print(paste("Run",i+1,"of",n.runs,"complete at",date()))
	flush.console()
}
setwd("..")



############################################################
# Generate plot for center panel of Figure S1
# using MCMC output generated by code above
############################################################

outputdir1 <- timestamp1
outputdir2 <- timestamp2

# n and p settings
n.1 <- rep(c(50,100,200),each=5)
p.1 <- n.1*2*(1:5)
n.2 <- rep(c(50,100,200),each=5)
p.2 <- n.2*c(2,3,4,7,14)/5

# Approximate pairwise correlation of covariates
r. <- 0.2

# Sparsity of "true" coefficients for generating data
s. <- 0.8

# Grid of settings for runs
runs1 <- cbind(n.1,p.1,r.,s.)
runs2 <- cbind(n.2,p.2,r.,s.)
n.runs1 <- dim(runs1)[1]
n.runs2 <- dim(runs2)[1]

# Regularization parameter
lambda <- 1

# Number of MCMC iterations per chain
K <- 11000

# Number of separate chains to run
M <- 10

# Fraction of output to discard as burn-in
burn <- 1/11

K.burn <- round(K*burn)
ac.orig.single1 <- rep(NA,n.runs1*M)
ac.fast.single1 <- rep(NA,n.runs1*M)
ac.orig.single2 <- rep(NA,n.runs2*M)
ac.fast.single2 <- rep(NA,n.runs2*M)
chain.run.id1 <- rep(1:n.runs1,each=M)
chain.run.id2 <- rep(1:n.runs2,each=M)
run1digits <- max(1,ceiling(log(n.runs1,10)))	# digits needed for run label strings
run2digits <- max(1,ceiling(log(n.runs2,10)))	# digits needed for run label strings
chaindigits <- max(1,ceiling(log(M,10)))		# digits needed for chain label strings
setwd(outputdir1)
for(i in 0:(n.runs1-1)){
	runtext <- substring(format(i/(10^run1digits),nsmall=run1digits),3)
	for(j in 0:(M-1)){
		chaintext <- substring(format(j/(10^chaindigits),nsmall=chaindigits),3)
		filename.orig.ij <- paste("acvp-orig-",runtext,"-",chaintext,"-s.txt",sep="")
		sigma2.chain.orig.ij <- read.table(filename.orig.ij)[-(1:K.burn),1]
		ac.orig.single1[i*M+j+1] <- acf(sigma2.chain.orig.ij,lag.max=1,plot=F)$acf[2]
		filename.fast.ij <- paste("acvp-fast-",runtext,"-",chaintext,"-s.txt",sep="")
		sigma2.chain.fast.ij <- read.table(filename.fast.ij)[-(1:K.burn),1]
		ac.fast.single1[i*M+j+1] <- acf(sigma2.chain.fast.ij,lag.max=1,plot=F)$acf[2]
	}
}
setwd("..")
setwd(outputdir2)
for(i in 0:(n.runs2-1)){
	runtext <- substring(format(i/(10^run2digits),nsmall=run2digits),3)
	for(j in 0:(M-1)){
		chaintext <- substring(format(j/(10^chaindigits),nsmall=chaindigits),3)
		filename.orig.ij <- paste("acvp-orig-",runtext,"-",chaintext,"-s.txt",sep="")
		sigma2.chain.orig.ij <- read.table(filename.orig.ij)[-(1:K.burn),1]
		ac.orig.single2[i*M+j+1] <- acf(sigma2.chain.orig.ij,lag.max=1,plot=F)$acf[2]
		filename.fast.ij <- paste("acvp-fast-",runtext,"-",chaintext,"-s.txt",sep="")
		sigma2.chain.fast.ij <- read.table(filename.fast.ij)[-(1:K.burn),1]
		ac.fast.single2[i*M+j+1] <- acf(sigma2.chain.fast.ij,lag.max=1,plot=F)$acf[2]
	}
}
setwd("..")
ac.orig1 <- tapply(ac.orig.single1,chain.run.id1,mean)
ac.fast1 <- tapply(ac.fast.single1,chain.run.id1,mean)
ac.orig2 <- tapply(ac.orig.single2,chain.run.id2,mean)
ac.fast2 <- tapply(ac.fast.single2,chain.run.id2,mean)
n. <- c(n.1,n.2)
p. <- c(p.1,p.2)
ac.orig <- c(ac.orig1,ac.orig2)
ac.fast <- c(ac.fast1,ac.fast2)

plotfilename <- paste("acvpn-",timestamp1,"-",timestamp2,".pdf",sep="")
min.ac <- min(ac.orig,ac.fast,0)
pdf(plotfilename,width=5,height=5)
par(mar=c(4,4,1,1),mgp=c(2.75,1.00,0),pty="s")
plot(log(p./n.),ac.orig,pch=1,ylim=c(min.ac,1),cex=1.4,cex.axis=1.4,cex.lab=1.4,
	xlab=expression(log(italic(p)/italic(n))),ylab="Autocorrelation",xaxt="n")
axis(1,seq(-1,2,by=0.5),c("-1.0","","0.0","","1.0","","2.0"),cex.axis=1.4)
points(log(p./n.),ac.fast,pch=2)
legend("topleft",c("Original","Fast"),pch=c(1,2),cex=1.4)
dev.off()



############################################################
# Generate MCMC output for right panel of Figure S1
############################################################

# n and p settings
n. <- rep(c(50,100,200),each=5)
p. <- n.*2*(1:5)

# Approximate pairwise correlation of covariates
r. <- 0.8

# Sparsity of "true" coefficients for generating data
s. <- 0.8

# Grid of settings for runs
runs <- cbind(n.,p.,r.,s.)
n.runs <- dim(runs)[1]

# Regularization parameter
lambda <- 1

# Number of MCMC iterations per chain
K <- 20000

# Number of separate chains to run
M <- 10

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

# Random seed (digits 40-42 of pi)
set.seed(169)

timestamp1 <- format(Sys.time(),"%Y%m%d%H%M%S")
dir.create(timestamp1)
setwd(timestamp1)
print(paste("Writing output to",timestamp))
flush.console()
rundigits <- max(1,ceiling(log(n.runs,10)))  # digits needed for run label strings
for(i in 0:(n.runs-1)){
	n <- runs[i+1,1]
	p <- runs[i+1,2]
	r <- runs[i+1,3]
	s <- runs[i+1,4]
	Xvarhalf <- diag(sqrt(1-r),p)+matrix((sqrt(1+(p-1)*r)-sqrt(1-r))/p,nrow=p,ncol=p)
	X.raw <- matrix(rnorm(n*p),nrow=n,ncol=p)%*%Xvarhalf
	X <- scale(X.raw)*sqrt(n/(n-1))
	XTX <- t(X)%*%X
	Y.raw <- drop(X%*%beta.true()+noise.true())
	Y <- Y.raw-mean(Y.raw)
	XTY <- drop(t(X)%*%Y)
	YTY <- sum(Y^2)
	runtext <- substring(format(i/(10^rundigits),nsmall=rundigits),3)
	outfile.stem.orig <- paste("acvp-orig-",runtext,sep="")
	run.bl(X,Y,lambda,K,M,outfile.stem.orig,fast=F,keep.beta=F,write.each=T)
	outfile.stem.fast <- paste("acvp-fast-",runtext,sep="")
	run.bl(X,Y,lambda,K,M,outfile.stem.fast,fast=T,keep.beta=F,write.each=T)
	print(paste("Run",i+1,"of",n.runs,"complete at",date()))
	flush.console()
}
setwd("..")

# n and p settings
n. <- rep(c(50,100,200),each=5)
p. <- n.*c(2,3,4,7,14)/5

# Approximate pairwise correlation of covariates
r. <- 0.8

# Sparsity of "true" coefficients for generating data
s. <- 0.8

# Grid of settings for runs
runs <- cbind(n.,p.,r.,s.)
n.runs <- dim(runs)[1]

# Regularization parameter
lambda <- 1

# Number of MCMC iterations per chain
K <- 11000

# Number of separate chains to run
M <- 10

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

# Random seed (digits 43-45 of pi)
set.seed(399)

timestamp2 <- format(Sys.time(),"%Y%m%d%H%M%S")
dir.create(timestamp2)
setwd(timestamp2)
print(paste("Writing output to",timestamp))
flush.console()
rundigits <- max(1,ceiling(log(n.runs,10)))  # digits needed for run label strings
for(i in 0:(n.runs-1)){
	n <- runs[i+1,1]
	p <- runs[i+1,2]
	r <- runs[i+1,3]
	s <- runs[i+1,4]
	Xvarhalf <- diag(sqrt(1-r),p)+matrix((sqrt(1+(p-1)*r)-sqrt(1-r))/p,nrow=p,ncol=p)
	X.raw <- matrix(rnorm(n*p),nrow=n,ncol=p)%*%Xvarhalf
	X <- scale(X.raw)*sqrt(n/(n-1))
	XTX <- t(X)%*%X
	Y.raw <- drop(X%*%beta.true()+noise.true())
	Y <- Y.raw-mean(Y.raw)
	XTY <- drop(t(X)%*%Y)
	YTY <- sum(Y^2)
	runtext <- substring(format(i/(10^rundigits),nsmall=rundigits),3)
	outfile.stem.orig <- paste("acvp-orig-",runtext,sep="")
	run.bl(X,Y,lambda,K,M,outfile.stem.orig,fast=F,keep.beta=F,write.each=T)
	outfile.stem.fast <- paste("acvp-fast-",runtext,sep="")
	run.bl(X,Y,lambda,K,M,outfile.stem.fast,fast=T,keep.beta=F,write.each=T)
	print(paste("Run",i+1,"of",n.runs,"complete at",date()))
	flush.console()
}
setwd("..")



############################################################
# Generate plot for right panel of Figure S1
# using MCMC output generated by code above
############################################################

outputdir1 <- timestamp1
outputdir2 <- timestamp2

# n and p settings
n.1 <- rep(c(50,100,200),each=5)
p.1 <- n.1*2*(1:5)
n.2 <- rep(c(50,100,200),each=5)
p.2 <- n.2*c(2,3,4,7,14)/5

# Approximate pairwise correlation of covariates
r. <- 0.8

# Sparsity of "true" coefficients for generating data
s. <- 0.8

# Grid of settings for runs
runs1 <- cbind(n.1,p.1,r.,s.)
runs2 <- cbind(n.2,p.2,r.,s.)
n.runs1 <- dim(runs1)[1]
n.runs2 <- dim(runs2)[1]

# Regularization parameter
lambda <- 1

# Number of MCMC iterations per chain
K <- 11000

# Number of separate chains to run
M <- 10

# Fraction of output to discard as burn-in
burn <- 1/11

K.burn <- round(K*burn)
ac.orig.single1 <- rep(NA,n.runs1*M)
ac.fast.single1 <- rep(NA,n.runs1*M)
ac.orig.single2 <- rep(NA,n.runs2*M)
ac.fast.single2 <- rep(NA,n.runs2*M)
chain.run.id1 <- rep(1:n.runs1,each=M)
chain.run.id2 <- rep(1:n.runs2,each=M)
run1digits <- max(1,ceiling(log(n.runs1,10)))	# digits needed for run label strings
run2digits <- max(1,ceiling(log(n.runs2,10)))	# digits needed for run label strings
chaindigits <- max(1,ceiling(log(M,10)))		# digits needed for chain label strings
setwd(outputdir1)
for(i in 0:(n.runs1-1)){
	runtext <- substring(format(i/(10^run1digits),nsmall=run1digits),3)
	for(j in 0:(M-1)){
		chaintext <- substring(format(j/(10^chaindigits),nsmall=chaindigits),3)
		filename.orig.ij <- paste("acvp-orig-",runtext,"-",chaintext,"-s.txt",sep="")
		sigma2.chain.orig.ij <- read.table(filename.orig.ij)[-(1:K.burn),1]
		ac.orig.single1[i*M+j+1] <- acf(sigma2.chain.orig.ij,lag.max=1,plot=F)$acf[2]
		filename.fast.ij <- paste("acvp-fast-",runtext,"-",chaintext,"-s.txt",sep="")
		sigma2.chain.fast.ij <- read.table(filename.fast.ij)[-(1:K.burn),1]
		ac.fast.single1[i*M+j+1] <- acf(sigma2.chain.fast.ij,lag.max=1,plot=F)$acf[2]
	}
}
setwd("..")
setwd(outputdir2)
for(i in 0:(n.runs2-1)){
	runtext <- substring(format(i/(10^run2digits),nsmall=run2digits),3)
	for(j in 0:(M-1)){
		chaintext <- substring(format(j/(10^chaindigits),nsmall=chaindigits),3)
		filename.orig.ij <- paste("acvp-orig-",runtext,"-",chaintext,"-s.txt",sep="")
		sigma2.chain.orig.ij <- read.table(filename.orig.ij)[-(1:K.burn),1]
		ac.orig.single2[i*M+j+1] <- acf(sigma2.chain.orig.ij,lag.max=1,plot=F)$acf[2]
		filename.fast.ij <- paste("acvp-fast-",runtext,"-",chaintext,"-s.txt",sep="")
		sigma2.chain.fast.ij <- read.table(filename.fast.ij)[-(1:K.burn),1]
		ac.fast.single2[i*M+j+1] <- acf(sigma2.chain.fast.ij,lag.max=1,plot=F)$acf[2]
	}
}
setwd("..")
ac.orig1 <- tapply(ac.orig.single1,chain.run.id1,mean)
ac.fast1 <- tapply(ac.fast.single1,chain.run.id1,mean)
ac.orig2 <- tapply(ac.orig.single2,chain.run.id2,mean)
ac.fast2 <- tapply(ac.fast.single2,chain.run.id2,mean)
n. <- c(n.1,n.2)
p. <- c(p.1,p.2)
ac.orig <- c(ac.orig1,ac.orig2)
ac.fast <- c(ac.fast1,ac.fast2)

plotfilename <- paste("acvpn-",timestamp1,"-",timestamp2,".pdf",sep="")
min.ac <- min(ac.orig,ac.fast,0)
pdf(plotfilename,width=5,height=5)
par(mar=c(4,4,1,1),mgp=c(2.75,1.00,0),pty="s")
plot(log(p./n.),ac.orig,pch=1,ylim=c(min.ac,1),cex=1.4,cex.axis=1.4,cex.lab=1.4,
	xlab=expression(log(italic(p)/italic(n))),ylab="Autocorrelation",xaxt="n")
axis(1,seq(-1,2,by=0.5),c("-1.0","","0.0","","1.0","","2.0"),cex.axis=1.4)
points(log(p./n.),ac.fast,pch=2)
legend("topleft",c("Original","Fast"),pch=c(1,2),cex=1.4)
dev.off()

