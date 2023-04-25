
#(1) %%%%%%%% Data generators %%%%%%%%%%%%%%%%%
## (i) true beta
## s: sparsity proportion
generate_beta_true <- function(p, s){
  beta_true <- numeric(p)
  n_nonzeros <- ceiling(p*s)
  if(n_nonzeros > 0){
    beta_true[1:n_nonzeros] <- rt(n_nonzeros, 2)
  }
  return(beta_true)
}

## (ii) error terms
generate_noise_true <- function(n){
  rnorm(n)
}


#------ Inverse Gaussian distribution
rinvgauss_mod <- function(n, mean = 1, shape = NULL, dispersion = 1) {
  if (!is.null(shape)) {
    dispersion <- 1 / shape
  }
  mu <- rep(mean, length.out = n)
  phi <- rep(dispersion, length.out = n)
  r <- numeric(n)
  i <- (mu > 0 & phi > 0)
  if (!all(i)) {
    r[!i] <- NA
    n <- sum(i)
  }
  phi[i] <- phi[i] * mu[i]
  Y <- rchisq(n, df = 1)
  X1 <- 1 + phi[i] / 2 * (Y - sqrt(4 * Y / phi[i] + Y^2))
  X1[!is.finite(X1)] <- (1 / (Y[!is.finite(X1)] * phi[i][!is.finite(X1)]))
  first_root <- as.logical(rbinom(n, size = 1, prob = 1 / (1 + X1)))
  r[first_root] <- X1[first_root]
  r[!first_root] <- 1 / X1[!first_root]
  r[i] <- r[i] * mu[i]
  return(r)
}

#------ Robust Inverse Gaussian Draws
rinvgaussian <- function(n, mean, lambda){
  m <- mean / sqrt(mean * lambda)
  l <- lambda / sqrt(mean * lambda)
  return(sqrt(mean * lambda) * rinvgauss_mod(n, m, l))
}





# (2) %%%%%%%%%%%%%% LASSO ESTIMATTOR %%%%%%%%%%%%%%%
# (i) Compute lasso estimator
compute_lasso_estimator <- function(x, y, lambda) {
  lasso_fit <- lars(x, y, type = "lasso", normalize = FALSE, intercept = FALSE)
  predict(lasso_fit, type = "coefficients", mode = "lambda", s = lambda)$coefficients
}

# (ii) Obtain the value of lambda that results in a specific number of non-zero coefficient estimates.
get_lambda_for_nonzero <- function(x, y, num_nonzero) {
  lasso_fit <- lars(x, y, type = "lasso", normalize = FALSE, intercept = FALSE)
  nonzero_path <- rev(which(lasso_fit$df == num_nonzero))[1]
  lambda_sequence <- c(lasso_fit$lambda, 0)
  if(is.na(nonzero_path)) {
    stop("The target number of nonzero coefficients cannot be attained.")
  }
  if(nonzero_path == 1) {
    return(2 * lambda_sequence[1])
  } else {
    return((lambda_sequence[nonzero_path-1] + lambda_sequence[nonzero_path]) / 2)
  }
}

# (iii) Obtain all value of lambda that results in a specific number of non-zero coefficient estimates.
# get_lambda_all <- function(x, y, number_nonzero){
#   fit <- glmnet(x, y, family="gaussian", alpha=1, intercept=FALSE, standardize=FALSE)
#   nonzero_indices <- max.col(abs(coef(fit)), ties.method="first")[-1]
#   lambda_values <- fit$lambda[nonzero_indices]
#   lambda_values_sorted <- sort(lambda_values, decreasing=TRUE)
#   lambda_pairs <- matrix(NA, nrow=length(lambda_values_sorted)-1, ncol=2)
#   for(i in 1:(length(lambda_values_sorted)-1)){
#     lambda_pairs[i,] <- c(lambda_values_sorted[i+1], lambda_values_sorted[i])
#   }
#   colnames(lambda_pairs) <- c("lwr", "upr")
#   return(lambda_pairs)
# }

get_lambda_all <- function(x, y, number_nonzero){
  lars_fit <- lars(x, y, type = "lasso", normalize = FALSE, intercept = FALSE)
  nonzero_indices <- which.max(lars_fit$df >= number_nonzero)
  if (length(nonzero_indices) == 0){
    return(matrix(numeric(0), nrow = 0, ncol = 2))
  }
  path_lambdas <- rev(c(Inf, lars_fit$lambda, 0))
  result <- matrix(numeric(0), nrow = length(nonzero_indices), ncol = 2)
  for(i in 1:length(nonzero_indices)){
    idx <- nonzero_indices[i]
    result[i,] <- path_lambdas[c(idx+1, idx)]
  }
  colnames(result) <- c("lwr", "upr")
  return(result)
}
