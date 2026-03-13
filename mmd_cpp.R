library(Rcpp)
library(np)
library(dplyr)
library(nloptr)

# ------------------------------
# C++ code compiled via Rcpp
# ------------------------------
Rcpp::sourceCpp("./mmd_cpp.cpp")

# ------------------------------
# R wrapper: MMD_bounds_cpp
# ------------------------------
# v0_col, v1_col: numeric vectors length n. the minimum and the maximum of unknown v.

MMD_bounds <- function(formula, data, v0_col, v1_col,
                       alpha = 0.05,            # Critical value
                       B = 200,                 # Bootstrap reps for CI calculation
                       b_exponent = 0.8,        # Proportion of sample used for bootstrap CI
                       verbose = TRUE) {
  
  # --- 1. Data Preparation ---
  cl <- match.call()
  if(verbose) cat(">> [1/5] Preparing data...\n")
  
  mf <- model.frame(formula, data)
  y  <- as.numeric(model.response(mf)) 
  
  full_x_mat <- model.matrix(formula, mf)
  intercept_idx <- which(colnames(full_x_mat) == "(Intercept)")
  
  if (length(intercept_idx) > 0) {
    x_mat <- full_x_mat[, -intercept_idx, drop = FALSE]
  } else {
    x_mat <- full_x_mat
  }
  
  v0 <- data[[v0_col]]
  v1 <- data[[v1_col]]
  n <- nrow(x_mat)
  d <- ncol(x_mat)
  
  param_names <- c("(Intercept)", paste0("latent_", v0_col), colnames(x_mat))
  p <- length(param_names)
  
  # --- 2. Nuisance Parameter Estimation (eta) ---
  if(verbose) cat(">> [2/5] Estimating nuisance parameter eta...\n")
  
  df_np <- data.frame(yn = y, v0 = v0, v1 = v1)
  if (d > 0) df_np <- cbind(df_np, as.data.frame(x_mat))
  fml_np <- as.formula(paste("yn ~", paste(names(df_np)[-1], collapse = " + ")))
  
  # May tweak bwmethod for speed. Currently get it from CV.
  bw <- np::npregbw(fml_np, data = df_np, regtype = "lc", bwmethod = "cv.aic") 
  model_eta <- np::npreg(bws = bw)
  eta <- fitted(model_eta)
  
  # --- 3. Unconstrained Minimization (Point Estimates) ---
  if(verbose) cat(">> [3/5] Finding unconstrained sample minimum...\n")
  
  obj_fun_R <- function(par) {
    Q_obj_cpp(as.matrix(x_mat), as.numeric(v0), as.numeric(v1), par, as.numeric(eta))
  }
  
  init_par <- rep(0, p)
  opt_unconstr <- optim(init_par, obj_fun_R, method = "Nelder-Mead", control = list(maxit = 5000))
  
  Q_min_hat <- opt_unconstr$value
  theta_hat <- opt_unconstr$par 
  
  if(verbose) cat(sprintf("       Sample Minimum Q_n = %.6f\n", Q_min_hat))
  
  # --- 4. Subsampling (CONDITIONAL ON B > 0) ---
  # Bootstrap subsampling is needed to use Kaidd, Molinari, and Stoyse (2019) framework.
  # c_alpha is the thresholds.
  c_alpha <- NA
  tau_CI <- NA
  
  if (B > 0) {
    if(verbose) cat(">> [4/5] Subsampling for critical values (B =", B, ")...\n")
    
    b_size <- floor(n^b_exponent)
    W_stats <- numeric(B)
    
    obj_fun_sub_R <- function(par, indices_0_based) {
      Q_obj_subsample_cpp(as.matrix(x_mat), as.numeric(v0), as.numeric(v1), 
                          par, as.numeric(eta), indices_0_based)
    }
    
    for(i in 1:B) {
      idx <- sample(0:(n-1), b_size, replace = FALSE)
      
      # We restart optim from theta_hat for speed
      opt_sub <- optim(theta_hat, obj_fun_sub_R, indices_0_based = idx, 
                       method = "Nelder-Mead", control = list(maxit = 500))
      
      Q_b_min <- opt_sub$value
      Q_b_theta_hat <- obj_fun_sub_R(theta_hat, idx)
      
      W_stats[i] <- b_size * (Q_b_theta_hat - Q_b_min)
    }
    
    # Thresholds obtained
    c_alpha <- quantile(W_stats, probs = 1 - alpha, names = FALSE)
    tau_CI <- Q_min_hat + (c_alpha / n)
    
    if(verbose) cat(sprintf("       c_alpha (95%%) = %.4f\n", c_alpha))
  } else {
    if(verbose) cat(">> [4/5] Skipping subsampling (B=0). CI will not be computed.\n")
  }
  
  # --- 5. Projection Method ---
  if(verbose) cat(">> [5/5] Running projections...\n")
  
  # Thresholds for Identified Set: log(n) / n -> From Chernozhukov et al. (2017)
  epsilon_n <- log(n) / n
  tau_ID <- Q_min_hat + epsilon_n
  
  bounds_ID <- matrix(NA, nrow = p, ncol = 2, dimnames = list(param_names, c("Lower", "Upper")))
  bounds_CI <- matrix(NA, nrow = p, ncol = 2, dimnames = list(param_names, c("Lower", "Upper")))
  
  opts <- list("algorithm" = "NLOPT_LN_COBYLA", "xtol_rel" = 1.0e-4, "maxeval" = 1000)
  
  # Loop over parameters
  # Using constrained optimization, we find bounds for each parameter. 
  # Thus we have to loop k times. 
  for(k in 1:p) {
    eval_f_min <- function(x) { x[k] }
    eval_f_max <- function(x) { -x[k] }
    
    # --- A. Estimated Identified Set ---
    eval_g_ID <- function(x) { obj_fun_R(x) - tau_ID }
    
    res_L_ID <- nloptr(x0 = theta_hat, eval_f = eval_f_min, eval_g_ineq = eval_g_ID, opts = opts)
    res_U_ID <- nloptr(x0 = theta_hat, eval_f = eval_f_max, eval_g_ineq = eval_g_ID, opts = opts)
    
    bounds_ID[k, 1] <- res_L_ID$objective
    bounds_ID[k, 2] <- -res_U_ID$objective 
    
    # --- B. Confidence Interval (Only if B > 0) ---
    if (B > 0 && !is.na(tau_CI)) {
      eval_g_CI <- function(x) { obj_fun_R(x) - tau_CI }
      
      res_L_CI <- nloptr(x0 = theta_hat, eval_f = eval_f_min, eval_g_ineq = eval_g_CI, opts = opts)
      res_U_CI <- nloptr(x0 = theta_hat, eval_f = eval_f_max, eval_g_ineq = eval_g_CI, opts = opts)
      
      bounds_CI[k, 1] <- res_L_CI$objective
      bounds_CI[k, 2] <- -res_U_CI$objective
    }
  }
  
  # --- 6. Compile Results ---
  res <- list(
    param_names = param_names,
    point_est = theta_hat, 
    bounds_ID = as.data.frame(bounds_ID),
    bounds_CI = as.data.frame(bounds_CI), # Might contain NAs if B=0
    Q_min = Q_min_hat,
    thresholds = c(ID = tau_ID, CI = tau_CI),
    n = n,
    alpha = alpha,
    B = B, # Store B to check later in summary
    call = cl
  )
  
  class(res) <- "mmd_results"
  return(res)
}

# ------------------------------------------------------------------------------
# 3. S3 Methods
# ------------------------------------------------------------------------------

print.mmd_results <- function(x, ...) {
  cat("\nModified Minimum Distance Estimator (Manski-Tamer)\n")
  cat("--------------------------------------------------\n")
  cat("Call:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("Sample Size: ", x$n, "\n")
  cat("Min Q(theta):", round(x$Q_min, 6), "\n\n")
  
  cat("Estimated Identified Set:\n")
  print(round(x$bounds_ID, 4))
  
  if (x$B > 0) {
    cat("\n95% Confidence Intervals:\n")
    print(round(x$bounds_CI, 4))
  } else {
    cat("\n(Confidence Intervals not computed because B=0)\n")
  }
  cat("\n")
}

summary.mmd_results <- function(object, ...) {
  # May only include identified set bounds 
  if (object$B > 0) {
    tab <- data.frame(
      Point_Est = object$point_est,
      ID_Lower = object$bounds_ID$Lower,
      ID_Upper = object$bounds_ID$Upper,
      CI_Lower = object$bounds_CI$Lower,
      CI_Upper = object$bounds_CI$Upper
    )
  } else {
    tab <- data.frame(
      Point_Est = object$point_est,
      ID_Lower = object$bounds_ID$Lower,
      ID_Upper = object$bounds_ID$Upper
    )
  }
  
  ans <- list(
    call = object$call,
    n = object$n,
    Q_min = object$Q_min,
    thresholds = object$thresholds,
    table = tab,
    alpha = object$alpha,
    B = object$B
  )
  class(ans) <- "summary.mmd_results"
  return(ans)
}

print.summary.mmd_results <- function(x, ...) {
  cat("\n=================================================================\n")
  cat(" MMD ESTIMATION RESULTS (Projection Method)\n")
  cat("=================================================================\n")
  cat("Sample Size:   ", x$n, "\n")
  cat("Min Objective: ", sprintf("%.6f", x$Q_min), "\n")
  
  if (x$B > 0) {
    cat("Thresholds:    ID =", sprintf("%.6f", x$thresholds["ID"]), 
        " | CI =", sprintf("%.6f", x$thresholds["CI"]), "\n")
  } else {
    cat("Thresholds:    ID =", sprintf("%.6f", x$thresholds["ID"]), 
        " | CI = NA (B=0)\n")
  }
  
  cat("Confidence Lvl:", (1 - x$alpha) * 100, "%\n\n")
  
  cat("Parameter Bounds:\n")
  print(round(x$table, 4))
  cat("\n")
  cat("=================================================================\n")
}
