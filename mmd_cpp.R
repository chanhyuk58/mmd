library(Rcpp)
library(np)         
library(nloptr)     
library(stats)

# 1. Load C++ Backend
Rcpp::sourceCpp("mmd_cpp.cpp")

# 2. Main Estimation Function
MMD_bounds <- function(formula, data, v0_col, v1_col,  # {{{
                       method = c("projection", "grid"),
                       grid_list = NULL,
                       alpha = 0.05,           # Critical value
                       B = 200,                # Subsample rep for CI 
                       b_exponent = 0.8,       # Subsample proportion for CI
                       verbose = TRUE) {
  
  method <- match.arg(method)
  cl <- match.call()
  
  # --- 1. Data Preparation ---
  if(verbose) cat(">> [1/5] Preparing data...\n")
  mf <- model.frame(formula, data)
  y  <- as.numeric(model.response(mf)) 
  full_x_mat <- model.matrix(formula, mf)
  intercept_idx <- which(colnames(full_x_mat) == "(Intercept)")
  x_mat <- if(length(intercept_idx) > 0) full_x_mat[, -intercept_idx, drop=FALSE] else full_x_mat
  
  v0 <- data[[v0_col]]; v1 <- data[[v1_col]]
  n <- nrow(x_mat); d <- ncol(x_mat)
  param_names <- c("(Intercept)", paste0("latent_", v0_col), colnames(x_mat))
  p <- length(param_names)
  
  # --- 2. Nuisance Parameter (eta) ---
  # nonparametric kernel regression using np package
  if(verbose) cat(">> [2/5] Estimating eta (E[Y|X,V])...\n")
  df_np <- data.frame(yn = y, v0 = v0, v1 = v1)
  if (d > 0) df_np <- cbind(df_np, as.data.frame(x_mat))
  fml_np <- as.formula(paste("yn ~", paste(names(df_np)[-1], collapse = " + ")))
  bw <- np::npregbw(fml_np, data = df_np, regtype = "lc", bwmethod = "cv.aic") 
  eta <- fitted(np::npreg(bws = bw))
  
  # --- 3. Whole Sample Minimum for Q ---
  if(verbose) cat(">> [3/5] Finding unconstrained sample minimum...\n")
  obj_fun_R <- function(par) Q_obj_cpp(as.matrix(x_mat), v0, v1, par, eta)
  opt_unconstr <- optim(rep(0, p), obj_fun_R, method = "Nelder-Mead", control = list(maxit = 5000))
  Q_min_hat <- opt_unconstr$value
  theta_hat <- opt_unconstr$par 
  
  # --- 4. Subsampling for CI ---
  tau_CI <- NA
  if (B > 0) {
    if(verbose) cat(">> [4/5] Subsampling for critical values...\n")
    b_size <- floor(n^b_exponent); W_stats <- numeric(B)
    for(i in 1:B) {
      idx <- sample(0:(n-1), b_size, replace = FALSE)
      Q_b_min <- optim(theta_hat, function(p) Q_obj_subsample_cpp(as.matrix(x_mat), v0, v1, p, eta, idx), 
                       method = "Nelder-Mead", control = list(maxit = 500))$value
      Q_b_theta <- Q_obj_subsample_cpp(as.matrix(x_mat), v0, v1, theta_hat, eta, idx)
      W_stats[i] <- b_size * (Q_b_theta - Q_b_min)
    }
    tau_CI <- Q_min_hat + (quantile(W_stats, 1 - alpha, names = FALSE) / n)
  }
  
  tau_ID <- Q_min_hat + (log(n) / n)
  
  # --- 5. Estimation Branching ---
  bounds_ID <- matrix(NA, p, 2, dimnames = list(param_names, c("Lower", "Upper")))
  bounds_CI <- matrix(NA, p, 2, dimnames = list(param_names, c("Lower", "Upper")))
  
  if (method == "projection") {
    if(verbose) cat(">> [5/5] Running Projection Method (nloptr)...\n")
    opts <- list("algorithm" = "NLOPT_LN_COBYLA", "xtol_rel" = 1e-4, "maxeval" = 1000)
    for(k in 1:p) {
      # ID Set
      bounds_ID[k, 1] <- nloptr(theta_hat, function(x) x[k], eval_g_ineq = function(x) obj_fun_R(x)-tau_ID, opts=opts)$objective
      bounds_ID[k, 2] <- -nloptr(theta_hat, function(x) -x[k], eval_g_ineq = function(x) obj_fun_R(x)-tau_ID, opts=opts)$objective
      # CI
      if (B > 0) {
        bounds_CI[k, 1] <- nloptr(theta_hat, function(x) x[k], eval_g_ineq = function(x) obj_fun_R(x)-tau_CI, opts=opts)$objective
        bounds_CI[k, 2] <- -nloptr(theta_hat, function(x) -x[k], eval_g_ineq = function(x) obj_fun_R(x)-tau_CI, opts=opts)$objective
      }
    }
  } else {
    if(verbose) cat(">> [5/5] Running Grid Search (Baseline)...\n")
    if(is.null(grid_list)) stop("grid_list required for grid search.")
    
    # Create Grid
    full_grid <- as.matrix(expand.grid(grid_list))
    if(verbose) cat("       Grid size:", nrow(full_grid), "points.\n")
    
    # Evaluate Grid in C++
    q_values <- Q_eval_grid_cpp(as.matrix(x_mat), v0, v1, full_grid, eta)
    
    # Filter and find bounds
    id_points <- full_grid[q_values <= tau_ID, , drop=FALSE]
    if(nrow(id_points) > 0) {
      bounds_ID[, 1] <- apply(id_points, 2, min)
      bounds_ID[, 2] <- apply(id_points, 2, max)
    }
    
    if (B > 0) {
      ci_points <- full_grid[q_values <= tau_CI, , drop=FALSE]
      if(nrow(ci_points) > 0) {
        bounds_CI[, 1] <- apply(ci_points, 2, min)
        bounds_CI[, 2] <- apply(ci_points, 2, max)
      }
    }
  }
  
  res <- list(param_names=param_names, point_est=theta_hat, bounds_ID=as.data.frame(bounds_ID), 
              bounds_CI=as.data.frame(bounds_CI), Q_min=Q_min_hat, thresholds=c(ID=tau_ID, CI=tau_CI), 
              n=n, method=method, B=B, call=cl)
  class(res) <- "mmd_results"
  return(res)
}
# }}}

# 3. S3 Methods
print.mmd_results <- function(x, ...) {
  cat("\nMMD Estimator (Method:", x$method, ")\n")
  cat("Min Q:", round(x$Q_min, 6), "| n:", x$n, "\n\n")
  cat("Estimated Identified Set:\n"); print(round(x$bounds_ID, 4))
  if (x$B > 0) { cat("\n95% Confidence Intervals:\n"); print(round(x$bounds_CI, 4)) }
}

summary.mmd_results <- function(object, ...) {
  tab <- data.frame(Point_Est = object$point_est, ID_L = object$bounds_ID$Lower, ID_U = object$bounds_ID$Upper)
  if (object$B > 0) { tab$CI_L = object$bounds_CI$Lower; tab$CI_U = object$bounds_CI$Upper }
  ans <- list(table = tab, method = object$method, n = object$n, Q_min = object$Q_min)
  class(ans) <- "summary.mmd_results"; return(ans)
}

print.summary.mmd_results <- function(x, ...) {
  cat("\nMMD SUMMARY (", x$method, ")\n")
  cat("n:", x$n, "| Min Q:", round(x$Q_min, 6), "\n\n")
  print(round(x$table, 4)); cat("\n")
}
