library(Rcpp)
library(np)         
library(nloptr)     
library(stats)

if(file.exists("mmd_cpp.cpp")) {
  Rcpp::sourceCpp("mmd_cpp.cpp")
} else {
  stop("Error: 'mmd_cpp.cpp' not found.")
}

MMD_bounds <- function(formula, data, v0_col, v1_col,
                       method = c("projection", "grid"),
                       grid_radius = 0.05,   
                       alpha = 0.05,
                       B = 200, 
                       b_exponent = 0.8,
                       n_starts = 5,         
                       param_bounds = 100,   # Projection bounds
                       K_folds = 5,          # eta Cross-Fitting
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
  
  # --- 2. Nuisance Parameter (eta) via CROSS-FITTING ---
  if(verbose) cat(sprintf(">> [2/5] Estimating eta via %d-Fold Cross-Fitting...\n", K_folds))
  df_np <- data.frame(yn = y, v0 = v0, v1 = v1)
  if (d > 0) df_np <- cbind(df_np, as.data.frame(x_mat))
  fml_np <- as.formula(paste("yn ~", paste(names(df_np)[-1], collapse = " + ")))
  
  # eta calculation
  bw <- np::npregbw(fml_np, data = df_np, regtype = "lc", bwmethod = "cv.aic") 
  
  folds <- sample(rep(1:K_folds, length.out = n))
  eta <- numeric(n)
  
  for(k in 1:K_folds) {
    train_df <- df_np[folds != k, , drop = FALSE]
    test_df  <- df_np[folds == k, , drop = FALSE]
    
    mod_k <- np::npreg(bws = bw, txdat = train_df[,-1], tydat = train_df[,1])
    eta[folds == k] <- predict(mod_k, newdata = test_df[,-1])
  }
  
  # --- 3. Sample Minimum (Multi-Start Subplex) ---
  if(verbose) cat(sprintf(">> [3/5] Finding global sample minimum (%d starts)...\n", n_starts))
  obj_fun_R <- function(par) Q_obj_cpp(as.matrix(x_mat), v0, v1, par, eta)
  
  opts_unconstr <- list("algorithm" = "NLOPT_LN_SBPLX", "xtol_rel" = 1e-7, "maxeval" = 10000)
  
  best_Q <- Inf
  best_theta <- rep(0, p)
  
  for(i in 1:n_starts) {
    start_val <- if(i == 1) rep(0, p) else rnorm(p, 0, 0.5)
    
    opt <- nloptr(x0 = start_val, eval_f = obj_fun_R, opts = opts_unconstr)
    
    if(opt$objective < best_Q) {
      best_Q <- opt$objective
      best_theta <- opt$solution
    }
  }
  Q_min_hat <- best_Q
  theta_hat <- best_theta 
  names(theta_hat) <- param_names
  
  if(verbose) cat(sprintf("       Global Sample Minimum (Q_min_hat) = %.6f\n", Q_min_hat))
  
  # --- 4. Subsampling for CI (BOBYQA) ---
  # This follows the projection method proposed by Kaido, Molinari, and Stoye 2019
  tau_CI <- NA
  if (B > 0) {
    if(verbose) cat(">> [4/5] Subsampling for critical values...\n")
    b_size <- floor(n^b_exponent); W_stats <- numeric(B)
    
    opts_sub <- list("algorithm" = "NLOPT_LN_BOBYQA", "xtol_rel" = 1e-5, "maxeval" = 1000)
    
    for(i in 1:B) {
      idx <- sample(0:(n-1), b_size, replace = FALSE)
      
      # Find subsample minimum starting from theta_hat
      # This replaces global Q_min. Standard practice.
      sub_obj <- function(p) Q_obj_subsample_cpp(as.matrix(x_mat), v0, v1, p, eta, idx)
      opt_sub <- nloptr(x0 = theta_hat, eval_f = sub_obj, opts = opts_sub)
      Q_b_min <- opt_sub$objective
      
      # Evaluate at full-sample estimate
      Q_b_theta <- sub_obj(theta_hat)
      W_stats[i] <- b_size * (Q_b_theta - Q_b_min)
    }
    tau_CI <- Q_min_hat + (quantile(W_stats, 1 - alpha, names = FALSE) / n)
  }
  
  tau_ID <- Q_min_hat + (log(n) / n)
  if(verbose) cat(sprintf("       Threshold (tau_ID) = %.6f\n", tau_ID))
  
  # --- 5. Estimation Branching ---
  bounds_ID <- matrix(NA, p, 2, dimnames = list(param_names, c("Lower", "Upper")))
  bounds_CI <- matrix(NA, p, 2, dimnames = list(param_names, c("Lower", "Upper")))
  
  if (method == "projection") {
    if(verbose) cat(">> [5/5] Running Projection Method (nloptr)...\n")
    opts_proj <- list("algorithm" = "NLOPT_LN_COBYLA", "xtol_rel" = 1e-6, "maxeval" = 5000)
    
    # Bounded projections to prevent infinite sets
    lb <- rep(-param_bounds, p)
    ub <- rep(param_bounds, p)
    
    for(k in 1:p) {
      bounds_ID[k, 1] <- nloptr(theta_hat, function(x) x[k], eval_g_ineq = function(x) obj_fun_R(x)-tau_ID, lb=lb, ub=ub, opts=opts_proj)$objective
      bounds_ID[k, 2] <- -nloptr(theta_hat, function(x) -x[k], eval_g_ineq = function(x) obj_fun_R(x)-tau_ID, lb=lb, ub=ub, opts=opts_proj)$objective
      if (B > 0) {
        bounds_CI[k, 1] <- nloptr(theta_hat, function(x) x[k], eval_g_ineq = function(x) obj_fun_R(x)-tau_CI, lb=lb, ub=ub, opts=opts_proj)$objective
        bounds_CI[k, 2] <- -nloptr(theta_hat, function(x) -x[k], eval_g_ineq = function(x) obj_fun_R(x)-tau_CI, lb=lb, ub=ub, opts=opts_proj)$objective
      }
    }
  } else {
    if(verbose) cat(">> [5/5] Running 3-Step Adaptive Grid Search...\n")
    
    generate_local_grids <- function(centers_mat, radius, pts_per_dim) {
      offsets <- as.matrix(expand.grid(rep(list(seq(-radius, radius, length.out = pts_per_dim)), p)))
      C_rep <- centers_mat[rep(1:nrow(centers_mat), each = nrow(offsets)), , drop=FALSE]
      O_rep <- offsets[rep(1:nrow(offsets), times = nrow(centers_mat)), , drop=FALSE]
      return(unique(C_rep + O_rep))
    }
    
    # --- STEP 1: Coarse Grid  ---
    r1 <- grid_radius; pts1 <- 5
    grid1 <- generate_local_grids(matrix(theta_hat, nrow=1), r1, pts1)
    q1 <- Q_eval_grid_cpp(as.matrix(x_mat), v0, v1, grid1, eta)
    
    # Prune based on quantiles
    q_range1 <- max(q1) - min(q1)
    thresh1 <- max(tau_ID, min(q1) + 0.05 * q_range1)
    top_50_val <- sort(q1)[min(length(q1), 50)]
    survivors1 <- grid1[q1 <= max(thresh1, top_50_val), , drop=FALSE]
    
    # --- STEP 2: Medium Grid ---
    r2 <- r1 / (pts1 - 1); pts2 <- 3
    grid2 <- generate_local_grids(survivors1, r2, pts2)
    q2 <- Q_eval_grid_cpp(as.matrix(x_mat), v0, v1, grid2, eta)
    
    # prune again
    q_range2 <- max(q2) - min(q2)
    thresh2 <- max(tau_ID, min(q2) + 0.01 * q_range2)
    top_50_val2 <- sort(q2)[min(length(q2), 50)]
    survivors2 <- grid2[q2 <= max(thresh2, top_50_val2), , drop=FALSE]
    
    # --- STEP 3: Fine Grid ---
    r3 <- r2 / (pts2 - 1); pts3 <- 3
    grid3 <- generate_local_grids(survivors2, r3, pts3)
    q3 <- Q_eval_grid_cpp(as.matrix(x_mat), v0, v1, grid3, eta)
    
    id_points <- grid3[q3 <= tau_ID, , drop=FALSE]
    if(nrow(id_points) > 0) {
      bounds_ID[, 1] <- apply(id_points, 2, min)
      bounds_ID[, 2] <- apply(id_points, 2, max)
    } else {
      warning("No points met the strict tau_ID threshold in Step 3. Returning NAs.")
    }
    
    if (B > 0) {
      ci_points <- grid3[q3 <= tau_CI, , drop=FALSE]
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

# ------------------------------------------------------------------------------
# 3. S3 Methods
# ------------------------------------------------------------------------------

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
