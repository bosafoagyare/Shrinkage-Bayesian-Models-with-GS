#(1) %%%%%%%%%% BAYESIAN LASSO %%%%%%%%%%%%%%
#----------- 3-Block Algorithm
iter_bl <- function(beta, sigma2, X, Y, XTY, n, p, lambda){
  d_tau_inv <- rinvgauss(p, mu = sqrt(lambda^2*sigma2/beta^2), shape = lambda^2)
  d_tau <- 1/d_tau_inv
  
  # Efficient algorithm to generate beta when n < p
  u <- rnorm(p) * sqrt(d_tau)
  v <- X %*% u + rnorm(n)
  U <- chol(X %*% diag(d_tau) %*% t(X) + diag(n))
  w <- backsolve(U, backsolve(t(U), (Y/sqrt(sigma2)) - v, upper.tri = FALSE))
  beta_new <- sqrt(sigma2) * (u + d_tau * (t(X) %*% w))
  
  sigma2_new <- (sum((Y - X %*% beta_new)^2) + sum(beta_new^2 * d_tau_inv)) / rchisq(1, n + p - 1)
  
  return(list(beta = beta_new, sigma2 = sigma2_new))
}




#---------- 2-Block Algorithm
iter_bl_fast <- function(beta, sigma2, X, Y, XTY, n, p, lambda){
  d_tau_inv <- rinvgaussian(p, sqrt(lambda^2*sigma2/beta^2), lambda^2)
  d_tau <- as.vector(1/d_tau_inv)
  
  u <- rnorm(p)*sqrt(d_tau)
  v <- X %*% u + rnorm(n)
  U <- chol(X %*% diag(d_tau) %*% t(X) + diag(rep(1, n)))
  w <- backsolve(U, backsolve(t(U), -v, upper.tri = FALSE))
  
  # Efficient algorithm to generate sigma when n < p
  beta_tilde <- crossprod(X, w)
  beta_tilde <- d_tau * beta_tilde
  sigma2_new <- (sum(Y^2) - sum(XTY * beta_tilde)) / rchisq(1, n - 1)
  
  # Efficient algorithm to generate beta when n < p
  beta_new <- beta_tilde + sqrt(sigma2_new) * (u + d_tau * crossprod(t(X), w))
  
  return(list(beta = beta_new, sigma2 = sigma2_new))
}
