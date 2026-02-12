# Note: np::npreg is used for eta; the optimizer and bound-finder are implemented in C++.

library(Rcpp)
library(np)
library(data.table)
library(readr)
library(dplyr)

# ------------------------------
# C++ code compiled via Rcpp
# ------------------------------
Rcpp::sourceCpp("./mmd_cpp.cpp")

# ------------------------------
# R wrapper: MMD_bounds_cpp
# ------------------------------
# v0_col, v1_col: numeric vectors length n. the minimum and the maximum of unknown v.

MMD_bounds_cpp <- function(formula, data, v0_col, v1_col,
                           box_radius = 1.0,
                           bfgs_maxiter = 1000,
                           bfgs_gtol = 1e-6,
                           bfgs_ftol = 1e-12,
                           bisect_iter = 60,
                           verbose = FALSE) {
  
  # 1. Process Formula and Data
  cl <- match.call()
  mf <- model.frame(formula, data)
  y <- model.response(mf)
  
  # We extract the model matrix but FORCE removal of intercept because
  # the C++ code handles the intercept separately as pptr[0] (gamma0).
  full_x_mat <- model.matrix(formula, mf)
  intercept_idx <- which(colnames(full_x_mat) == "(Intercept)")
  
  if (length(intercept_idx) > 0) {
    x_mat <- full_x_mat[, -intercept_idx, drop = FALSE]
  } else {
    x_mat <- full_x_mat
  }
  
  v0 <- data[[v0_col]]
  v1 <- data[[v1_col]]
  
  n <- as.integer(length(y))
  d <- as.integer(ncol(x_mat))
  p <- 2 + d # [gamma0, gammaV, beta_x1, ..., beta_xd]

  # 2. Estimate eta (Nuisance Parameter)
  # For eta = E[Y | v0, v1, X], we use all observed information.
  # We use the clean x_mat (no constant intercept) to avoid the np warning.
  df_np <- data.frame(yn = as.numeric(y), v0 = as.numeric(v0), v1 = as.numeric(v1))
  if (d > 0) {
    df_np <- cbind(df_np, as.data.frame(x_mat))
  }
  
  # Construct formula for np: yn ~ v0 + v1 + x1 + ...
  rhs_np <- c("v0", "v1", colnames(x_mat))
  fml_np <- reformulate(rhs_np, response = "yn")

  if(verbose) cat("Estimating eta via nonparametric regression...\n")
  bw <- np::npregbw(fml_np, data = df_np, regtype = "lc")
  model_eta <- np::npreg(bws = bw)
  eta <- as.numeric(predict(model_eta, newdata = df_np[, -1, drop = FALSE]))

  # 3. Optimization (C++)
  par0 <- rep(0, p)
  if(verbose) cat("Starting BFGS optimization...\n")
  opt <- bfgs_opt_cpp(as.matrix(x_mat), as.numeric(v0), as.numeric(v1), as.numeric(eta),
                      par0,
                      maxiter = as.integer(bfgs_maxiter),
                      gtol = bfgs_gtol,
                      ftol = bfgs_ftol)

  par_opt <- as.numeric(opt$par)
  Q_min <- as.numeric(opt$value)
  threshold <- Q_min + (log(n) / n)

  # 4. Bounds via Bisection (C++)
  if(verbose) cat("Computing coordinate-wise bounds...\n")
  bounds_res <- coord_bounds_bisect_cpp(as.matrix(x_mat),
                                        as.numeric(v0), as.numeric(v1), as.numeric(eta),
                                        par_opt,
                                        threshold,
                                        box_radius = as.numeric(box_radius),
                                        max_bisect_iter = as.integer(bisect_iter),
                                        tol = 1e-8)

  # 5. Result Construction
  # Mapping names: [gamma0, gammaV, x_covariates]
  coeff_names <- c("(Intercept)", paste0("latent_", v0_col, "_", v1_col), colnames(x_mat))
  names(par_opt) <- coeff_names
  
  res <- list(
    coefficients = par_opt,
    bounds = data.frame(
      lower = as.numeric(bounds_res$lower),
      upper = as.numeric(bounds_res$upper),
      row.names = coeff_names
    ),
    Q_min = Q_min,
    threshold = threshold,
    n = n,
    call = cl,
    convergence = opt$conv
  )
  
  class(res) <- "mmd_bounds"
  return(res)
}

# ------------------------------
# S3 Methods
# ------------------------------

print.mmd_bounds <- function(x, ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("Coefficients:\n")
  print(x$coefficients)
  cat("\nMinimum Q:", round(x$Q_min, 6), "\n")
}

summary.mmd_bounds <- function(object, ...) {
  ans <- list(
    call = object$call,
    stats = data.frame(
      Estimate = object$coefficients,
      Lower_Bound = object$bounds$lower,
      Upper_Bound = object$bounds$upper
    ),
    Q_min = object$Q_min,
    n = object$n,
    conv = object$convergence
  )
  class(ans) <- "summary.mmd_bounds"
  return(ans)
}

print.summary.mmd_bounds <- function(x, ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("Sample Size: ", x$n, "\n")
  cat("Optimization converged:", x$conv, "\n")
  cat("Minimum Q:   ", round(x$Q_min, 6), "\n\n")
  cat("Point Estimates and Identification Bounds:\n")
  print(x$stats)
  cat("\n")
}

# ------------------------------
# Small demo
# ------------------------------
if (FALSE) {
  # example simulation
  set.seed(63130)
  N <- 10**4 # Population size

  # Population {{{
  ## Truth
  #gamma <- c(1, 1, -1)
  gamma <- c(1, 2, rep(-1, 5))
  v <- rnorm(N, 0, 2)
  # v <- runif(N, -2, 3) # try uniform distribution
  
  ## Observed
  v1 = ceiling(v)
  v0 <- v1 - 1
  # x <- rnorm(n, 1, 4)
  # x <- runif(N, 0, 5) # try uniform distribution
  # multivariate uniform distribution
  d <- 5
  x <- matrix(runif(N * d), nrow = N, ncol = d)
  ## Error
  eps <- rnorm(N, 0, 1)
  
  ## y
  y <- gamma[1] + gamma[2]*v + x %*% gamma[3:length(gamma)] + eps
  
  ## pop
  pop <- data.frame(y, x, v0, v1)
  # }}}

  fit <- MMD_bounds_cpp(y ~ X1 + X2 + X3 + X4 + X5, 
    data = pop, v0_col = "v0", v1_col = "v1")
  summary(fit)
}
