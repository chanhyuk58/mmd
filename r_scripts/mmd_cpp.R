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
                       grid_radius = 1.0,       # Decides the size of candidates
                       grid_points = 10000,    # Eval Points for each Parameters
                       alpha = 0.05,
                       B = 200, 
                       b_exponent = 0.8,
                       n_starts = 5,         
                       param_bounds = 100,   
                       K_folds = 5,          
                       verbose = TRUE) {
  
  method <- match.arg(method)
  cl <- match.call()
  
  # --- 1. Data Preparation  ---
  if(verbose) cat(">> [1/5] Preparing data...\n")
  mf <- model.frame(formula, data)
  y  <- as.numeric(model.response(mf)) 
  full_x_mat <- model.matrix(formula, mf)
  intercept_idx <- which(colnames(full_x_mat) == "(Intercept)")
  x_mat <- if(length(intercept_idx) > 0) full_x_mat[, -intercept_idx, drop=FALSE] else full_x_mat
  
  valid_idx <- as.integer(rownames(mf))
  v0 <- as.numeric(data[[v0_col]][valid_idx])
  v1 <- as.numeric(data[[v1_col]][valid_idx])
  
  n <- nrow(x_mat); d <- ncol(x_mat)
  param_names <- c("(Intercept)", paste0("latent_", v0_col), colnames(x_mat))
  p <- length(param_names)
  
  # --- 2. Nuisance Parameter (eta) ---
  if(verbose) cat(sprintf(">> [2/5] Estimating eta via %d-Fold Cross-Fitting...\n", K_folds))
  df_np <- data.frame(yn = y, v0 = v0, v1 = v1)
  if (d > 0) df_np <- cbind(df_np, as.data.frame(x_mat))
  fml_np <- as.formula(paste("yn ~", paste(names(df_np)[-1], collapse = " + ")))
  
  bw_obj <- np::npregbw(fml_np, data = df_np, regtype = "lc", bwmethod = "cv.aic") 
  raw_bws <- bw_obj$bw 
  
  folds <- sample(rep(1:K_folds, length.out = n))
  eta <- numeric(n)
  
  for(k in 1:K_folds) {
    train_df <- df_np[folds != k, , drop = FALSE]
    test_df  <- df_np[folds == k, , drop = FALSE]
    
    # Train model
    mod_k <- np::npreg(bws = raw_bws, txdat = train_df[, -1, drop = FALSE], tydat = train_df[, 1])
    
    # Validate
    pred_vals <- as.numeric(predict(mod_k, newdata = test_df[, -1, drop = FALSE]))
    eta[folds == k] <- pred_vals
  }

  # --- 3. Sample Minimum (Multi-Start Subplex) ---
  if(verbose) cat(sprintf(">> [3/5] Finding global sample minimum (%d starts)...\n", n_starts))
  obj_fun_R <- function(par) {
    Q_obj_cpp(as.matrix(x_mat), as.numeric(v0), as.numeric(v1), as.numeric(par), as.numeric(eta))
  }
  
  opts_unconstr <- list("algorithm" = "NLOPT_LN_SBPLX", "xtol_rel" = 1e-7, "maxeval" = 10000)
  best_Q <- Inf; best_theta <- rep(0, p)
  
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
  
  # --- 4. Subsampling for CI ---
  # Following Kaido, Molinari, and Stoye (2019)
  tau_CI <- NA
  if (B > 0) {
    if(verbose) cat(">> [4/5] Subsampling for critical values...\n")
    b_size <- floor(n^b_exponent); W_stats <- numeric(B)
    opts_sub <- list("algorithm" = "NLOPT_LN_BOBYQA", "xtol_rel" = 1e-5, "maxeval" = 1000)
    
    for(i in 1:B) {
      idx <- sample(0:(n-1), b_size, replace = FALSE)
      sub_obj <- function(p) Q_obj_subsample_cpp(as.matrix(x_mat), v0, v1, p, eta, idx)
      opt_sub <- nloptr(x0 = theta_hat, eval_f = sub_obj, opts = opts_sub)
      W_stats[i] <- b_size * (sub_obj(theta_hat) - opt_sub$objective)
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
    lb <- rep(-param_bounds, p); ub <- rep(param_bounds, p)
    
    run_proj <- function(obj_dir, tau) {
      opt <- nloptr(theta_hat, eval_f = obj_dir, eval_g_ineq = function(x) obj_fun_R(x)-tau, lb=lb, ub=ub, opts=opts_proj)
      if(opt$status < 0 || opt$status == 5) warning("Projection optimizer failed or hit maxeval. Bounds may be inaccurate.")
      return(opt$objective)
    }
    
    for(k in 1:p) {
      bounds_ID[k, 1] <- run_proj(function(x) x[k], tau_ID)
      bounds_ID[k, 2] <- -run_proj(function(x) -x[k], tau_ID)
      if (B > 0) {
        bounds_CI[k, 1] <- run_proj(function(x) x[k], tau_CI)
        bounds_CI[k, 2] <- -run_proj(function(x) -x[k], tau_CI)
      }
    }
  } else {
    if(verbose) cat(sprintf(">> [5/5] Running Center-Out Profile Grid Search...\n"))
    opts_prof <- list("algorithm" = "NLOPT_LN_BOBYQA", "xtol_rel" = 1e-4, "maxeval" = 5000)
    
    for(k in 1:p) {
      if(verbose) cat(sprintf("       Profiling %s...\n", param_names[k]))
      
      eval_prof <- function(par_minus_k, v) {
        full_par <- numeric(p); full_par[-k] <- par_minus_k; full_par[k] <- v
        return(Q_obj_cpp(as.matrix(x_mat), v0, v1, full_par, eta))
      }
      
      grid_right <- seq(theta_hat[k], theta_hat[k] + grid_radius, length.out = grid_points)
      grid_left  <- seq(theta_hat[k], theta_hat[k] - grid_radius, length.out = grid_points)
      
      q_prof_right <- numeric(grid_points)
      q_prof_left  <- numeric(grid_points)
      
      # Sweep Right
      current_guess <- theta_hat[-k]
      for(i in 1:grid_points) {
        opt <- nloptr(x0 = current_guess, eval_f = function(x) eval_prof(x, grid_right[i]), opts = opts_prof)
        q_prof_right[i] <- opt$objective
        current_guess <- opt$solution
      }
      
      # Sweep Left
      current_guess <- theta_hat[-k]
      for(i in 1:grid_points) {
        opt <- nloptr(x0 = current_guess, eval_f = function(x) eval_prof(x, grid_left[i]), opts = opts_prof)
        q_prof_left[i] <- opt$objective
        current_guess <- opt$solution
      }
      
      # Combine results
      full_grid_k <- c(rev(grid_left), grid_right[-1])
      full_q_prof <- c(rev(q_prof_left), q_prof_right[-1])
      
      # Find Bounds
      valid_ID <- full_grid_k[full_q_prof <= tau_ID]
      if(length(valid_ID) > 0) {
        bounds_ID[k, ] <- c(min(valid_ID), max(valid_ID))
      }
      
      if (B > 0) {
        valid_CI <- full_grid_k[full_q_prof <= tau_CI]
        if(length(valid_CI) > 0) {
          bounds_CI[k, ] <- c(min(valid_CI), max(valid_CI))
        }
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
