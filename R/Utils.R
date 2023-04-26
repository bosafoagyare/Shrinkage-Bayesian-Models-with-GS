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