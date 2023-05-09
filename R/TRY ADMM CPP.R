#-------------------------------------
## Distributed Longitudinal ADMM
library(Rcpp)
library(RcppArmadillo)

# The soft-thresholding solution of r
shrinkcpp <- function(u, v) {
  w <- (1 + sign(u - v)) / 2 * (u - v) - (1 + sign(-u - v)) / 2 * (-u - v)
  return(w)
}

# The first order derivative of checkloss function 
lossweight <- function(u, tau) {
  num <- length(u)
  loss <- rep(0, num)
  for (i in seq_along(u)) {
    if (u[i] > 0) {
      loss[i] <- -tau
    } else {
      loss[i] <- 1 - tau
    }
  }
  return(loss)
}



# The WQR-ADMM algorithm
WQRCPP <- function(x, y, tau, rho, eps, maxstep, intercept = FALSE, warmstart = FALSE) {
  
  nsum <- nrow(x)
  if(intercept) {
    x <- cbind(rep(1, nsum), x)
  }
  p <- ncol(x)
  betaini <- rep(0, p)
  zini <- rep(0, p)
  uini <- rep(0, p)
  vini <- rep(0, nsum)
  beta <- rep(0, p)
  r <- rep(0, nsum)
  z <- rep(0, p)
  u <- rep(0, p)
  v <- rep(0, nsum)
  xbeta <- rep(0, nsum)
  time <- 0
  time_warm <- 0
  time_prep <- 0
  time_iter <- 0
  final <- rep(0, p + 2)
  
  if(warmstart) {
    start_warm <- Sys.time()
    zini <- solve(t(x) %*% x) %*% t(x) %*% y
    xbeta <- x %*% zini
    finish_warm <- Sys.time()
    elapsed_warm <- finish_warm - start_warm
    time_warm <- as.numeric(elapsed_warm, units = "secs")
    time <- time + time_warm
  }
  
  tmp <- matrix(0, p, p)
  start_prep <- Sys.time()
  if(nsum > p) {
    tmp <- solve(t(x) %*% x + diag(p))
  } else {
    tmp <- diag(p) - t(x) %*% solve(x %*% t(x) + diag(nsum))
  }
  finish_prep <- Sys.time()
  elapsed_prep <- finish_prep - start_prep
  time_prep <- as.numeric(elapsed_prep, units = "secs")
  time <- time + time_prep
  
  iteration <- 0
  distance <- 1
  
  while((distance > eps) && (iteration < maxstep)) {
    
    start_iter <- Sys.time()
    #update beta
    beta <- zini + uini / rho
    #update r
    r <- shrinkcpp(vini / rho + y - xbeta - 0.5 * (2 * tau - 1) / rho, 0.5 * rep(1, nsum) / rho)
    #update z
    z <- tmp %*% (t(x) %*% (y - r + vini / rho) + beta - uini / rho)
    #update u and v
    xbeta <- x %*% z
    u <- uini + rho * (z - beta)
    v <- vini + rho * (y - xbeta - r)
    finish_iter <- Sys.time()
    elapsed_iter <- finish_iter - start_iter
    time_iter <- as.numeric(elapsed_iter, units = "secs")
    time <- time + time_iter
    
    if(sum(abs(betaini)) == 0) {
      distance <- 1
    } else {
      distance <- sum(abs(beta - betaini)) / sum(abs(betaini))
    }
    
    betaini <- beta
    zini <- z
    uini <- u
    vini <- v
    iteration <- iteration + 1
    
  }
  
  final[1:p] <- beta
  final[p + 1] <- time
  return(final)
  
}



WeightCPP <- function(x, y, rep, tau, betahat, eps, maxstep, intercept, type) {
  
  nsum <- nrow(x)
  if(intercept) {
    x <- cbind(rep(1, nsum), x)
  }
  p <- ncol(x)
  n <- length(rep)
  time_prep1 <- 0
  time_prep2 <- 0
  time_iter <- 0
  time_wei <- 0
  time <- 0
  
  # calculate D1 and D2
  e <- rep(0, nsum)
  xM <- matrix(0, nrow = p, ncol = nsum)
  D1 <- matrix(0, nrow = p, ncol = p)
  D2 <- matrix(0, nrow = p, ncol = p)
  start_prep1 <- Sys.time()
  e <- y - x %*% betahat
  sigma <- 1 / sqrt(rowSums(x^2) / n)
  psi <- sigma * dnorm(e * sigma)
  xw <- sweep(x, 1, sqrt(psi), "*")
  D1 <- t(xw) %*% xw / n
  counter <- 1
  for (i in 1:n) {
    xi <- x[counter:(counter + rep[i] - 1), ]
    M <- matrix(0, nrow = rep[i], ncol = rep[i])
    if(type == "exchangeable") {
      M <- matrix(1, nrow = rep[i], ncol = rep[i])
      diag(M) <- rep(0, rep[i])
    } else {
      diag(M) <- rep(1, rep[i] - 1)
      diag(M[-1, -1]) <- 1
    }
    xM[, counter:(counter + rep[i] - 1)] <- t(xi) %*% M
    counter <- counter + rep[i]
  }
  xMw <- sweep(xM, 2, sqrt(psi), "*")
  D2 <- xMw %*% xw / n
  D <- D2 %*% solve(D1)
  finish_prep1 <- Sys.time()
  elapsed_prep1 <- finish_prep1 - start_prep1
  time_prep1 <- as.numeric(elapsed_prep1)
  time <- time + time_prep1
  
  # calculate g
  g <- matrix(0, nrow = n, ncol = p)
  counter <- 1
  start_prep2 <- Sys.time()
  for (i in 1:n) {
    xi <- x[counter:(counter + rep[i] - 1), ]
    lossi <- lossweight(e[counter:(counter + rep[i] - 1)], tau)
    g[i, ] <- t((xM[, counter:(counter + rep[i] - 1)] - D %*% t(xi)) %*% lossi)
    counter <- counter + rep[i]
  }
  finish_prep2 <- Sys.time()
  elapsed_prep2 <- finish_prep2 - start_prep2
  time_prep2 <- as.numeric(elapsed_prep2)
  time <- time + time_prep2
  
  # update lambda
  lambdaini <- rep(0, p)
  lambda <- rep(0, p)
  iteration <- 0
  distance <- 1
  
  while (iteration < maxstep & distance > eps) {
    S <- rep(0, p)
    F <- matrix(0, nrow = p, ncol = p)
    
    start_iter <- Sys.time()
    tmp <- 1 / (g * lambdaini + rep(1, n))
    S <- t(g) %*% tmp
    gw <- sweep(g, 1, tmp, FUN = "*")
    F <- t(gw) %*% gw
    lambda <- lambdaini + solve(F) %*% S
    finish_iter <- Sys.time()
    
    elapsed_iter <- finish_iter - start_iter
    time_iter <- as.numeric(elapsed_iter, units = "secs")
    time <- time + time_iter
    distance <- sum(abs(lambda - lambdaini))
    lambdaini <- lambda
    iteration <- iteration + 1
  }
  
  # calculate weight
  weight <- rep(0, n)
  
  start_wei <- Sys.time()
  weight <- 1 / (n * (rep(1, n) + g %*% lambda))
  finish_wei <- Sys.time()
  
  elapsed_wei <- finish_wei - start_wei
  time_wei <- as.numeric(elapsed_wei, units = "secs")
  time <- time + time_wei
  
  # return final value
  final <- c(rep(0, n + 2))
  final[1:n] <- weight
  final[n + 1] <- time
  final[n + 2] <- iteration
  return(final)
  
}




# Compute the weights
WS <- function(x, y, rep, tau, betahat, intercept, corrtype = "ar1", epsw = 1e-06, maxstepw = 100) {
  n <- length(rep)
  nsum <- nrow(x)
  weight <- WeightCPP(x, y, rep, tau, betahat, epsw, maxstepw, intercept, corrtype)
  weightnew <- rep(0, nsum)
  counter <- 1
  
  for (j in seq_along(rep)) {
    weightnew[counter:(counter + rep[j] - 1)] <- weight[j] * rep(1, rep[j])
    counter <- counter + rep[j]
  }
  
  final <- list(Weight = weightnew, Time_W = weight[nsum + 1])
  return(final)
}

#---------------- PARALLEL --------




library(Matrix)
shrinkcpp1 <- function(u, v){
  w <- (1 + sign(u - v)) / 2 * (u - v) - (1 + sign(-u - v)) / 2 * (-u - v)
  return(w)
}

deriv <- function(beta, a, lambda, penalty){
  p <- length(beta)
  df <- rep(0, p)
  if(penalty == "scad"){
    for(j in 1:p){
      if(abs(beta[j]) <= lambda) df[j] <- lambda
      else if(abs(beta[j]) <= a * lambda) df[j] <- (a * lambda - abs(beta[j])) / (a - 1)
    }
  } else if(penalty == "mcp"){
    for(j in 1:p){
      if(abs(beta[j]) <= a * lambda) df[j] <- lambda - abs(beta[j]) / a
    }
  } else{
    df <- rep(lambda, p)
  }
  return(df)
}

QRADM <- function(xr, yr, ta, rhor, lambdar, iter, intercept, M, penalty, a){
  maxit <- iter
  tau <- ta
  rho <- rhor
  lambda <- lambdar
  alpha <- 1.7
  x <- as.matrix(xr)
  y <- as.vector(yr)
  n <- nrow(x)
  p <- ncol(x)
  ni <- n / M
  lambda <- lambda / n
  rho <- rho / n
  r <- rep(0, n)
  df <- rep(0, p)
  u <- rep(0, n)
  beta <- rep(0, p)
  xbetai <- rep(0, ni)
  betaold <- beta
  if(intercept){
    p <- p + 1
    x <- cbind(rep(1, n), x)
    comparev <- rep(0, 3)
    beta <- rep(0, p)
  } else{
    comparev <- rep(0, 2)
    beta <- rep(0, p)
  }

  r <- y - x %*% beta
  ABSTOL <- 1e-7
  RELTOL <- 1e-4
  rnorm <- NULL
  epspri <- NULL
  snorm <- NULL
  epsdual <- NULL
  betai <- matrix(0, p, M)
  etai <- matrix(0, p, M)
  dat <- array(0, dim = c(p, p, M))
  for (i in 1:M) {
    tmp <- NULL
    xi <- x[(ni * i - ni + 1):(ni * i), ]
    if (ni > p && intercept) {
      tmp <- solve(t(xi) %*% xi + diag(p))
    } else {
      tmp <- diag(p) - xi %*% solve(xi %*% t(xi) + diag(ni))
    }
    dat[, , i] <- tmp
  }
  iteration <- 0
  while (iteration < maxit) {
    beta_avg <- apply(betai, 1, mean)
    eta_avg <- apply(etai, 1, mean)
    betaold <- beta
    if (intercept) {
      df <- deriv(beta[-1], a, lambda, penalty)
      beta[1] <- beta_avg[1] + eta_avg[1] / rho
      beta[-1] <- shrinkcpp1(beta_avg[-1] + eta_avg[-1] / rho, df / (rho * M))
    } else {
      df <- deriv(beta, a, lambda, penalty)
      beta <- shrinkcpp1(beta_avg + eta_avg / rho, df / (rho * M))
    }
    for (i in 1:M) {
      yi <- y[(ni * i - ni + 1):(ni * i)]
      ui <- u[(ni * i - ni + 1):(ni * i)]
      beta_i <- betai[, i]
      ri <- r[(ni * i - ni + 1):(ni * i)]
      xi <- x[(ni * i - ni + 1):(ni * i), ]
      xbetai <- alpha * xi %*% beta_i + (1 - alpha) * (yi - ri)
      r[(ni * i - ni + 1):(ni * i)] <- shrinkcpp1(ui / rho + yi - xbetai - 0.5 * (2 * tau - 1) / (n * rho), 0.5 * rep(1, ni) / (n * rho))
      betai[, i] <- dat[, , i] %*% (t(xi) %*% (yi - r[(ni * i - ni + 1):(ni * i)] + ui / rho) - eta_avg)
      etai[, i] <- etai[, i] + rho * (betai[, i] - beta)
    }
    iteration <- iteration + 1
  }

  if(intercept){
    rnorm = sqrt(sum((y-x%*%beta-r)^2))
    snorm = sqrt(sum((rho*x[,1:(p-1)]*(beta[1:(p-1)]-betaold[1:(p-1)]))^2))
    comparev=c(sqrt(sum((x[,1:(p-1)]*beta[1:(p-1)])^2)),sqrt(sum((-r)^2)),sqrt(sum((y-beta[1])^2)))
    epspri = sqrt(n)*ABSTOL + RELTOL*max(comparev)
    epsdual = sqrt(n)*ABSTOL + RELTOL*sqrt(sum(u^2))
  }
  else{
    rnorm = sqrt(sum((y-x%*%beta-r)^2))
    snorm = sqrt(sum((rho*x*(beta-betaold))^2))
    comparev=c(sqrt(sum((x*beta)^2)),sqrt(sum((-r)^2)))
    epspri = sqrt(n)*ABSTOL + RELTOL*max(comparev)
    epsdual = sqrt(n)*ABSTOL + RELTOL*sqrt(sum(u^2))
  }
  if (rnorm < epspri && snorm < epsdual)
    iteration = maxit+1
  else
    iteration = iteration + 1
}
betanew=beta
return(betanew)
}












