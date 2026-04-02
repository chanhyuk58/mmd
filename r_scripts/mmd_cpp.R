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
                       grid_radius = 0.5,    # Radius in standardized space
                       grid_points = 100,    
                       alpha = 0.05,
                       B = 200, 
                       b_exponent = 0.8,
                       n_starts = 5,         
                       param_bounds = 20,    
                       K_folds = 5,          
                       strict_bw = FALSE,    # CV for npreg Bandwidth setting
                       verbose = TRUE) {
  
  method <- match.arg(method)
  cl <- match.call()
  
  # --- 1. Data Preparation ---
  if(verbose) cat(">> [1/5] Preparing data...\n")
  
  fml_full <- update(formula, paste("~ . +", v0_col, "+", v1_col))
  mf <- model.frame(fml_full, data)
  
  # Extract Y
  y <- as.numeric(model.response(mf)) 
  
  # Extract X matrix
  full_x_mat <- model.matrix(formula, mf)
  intercept_idx <- which(colnames(full_x_mat) == "(Intercept)")
  x_mat <- if(length(intercept_idx) > 0) full_x_mat[, -intercept_idx, drop=FALSE] else full_x_mat
  
  # Extract v0 and v1
  v0 <- as.numeric(mf[[v0_col]])
  v1 <- as.numeric(mf[[v1_col]])
  
  n <- nrow(x_mat); d <- ncol(x_mat)
  param_names <- c("(Intercept)", paste0("latent_", v0_col), colnames(x_mat))
  p <- length(param_names)
  
  # --- 1b. Invisible Full Standardization ---
  if(verbose) cat(">> [1b/5] Invisible standardization...\n")
  mu_v <- mean(c(v0, v1), na.rm = TRUE)
  sd_v <- sd(c(v0, v1), na.rm = TRUE)
  if(sd_v == 0 || is.na(sd_v)) sd_v <- 1.0
  
  v0_std <- (v0 - mu_v) / sd_v
  v1_std <- (v1 - mu_v) / sd_v
  
  mu_x <- apply(x_mat, 2, mean, na.rm = TRUE)
  sd_x <- apply(x_mat, 2, sd, na.rm = TRUE)
  sd_x[sd_x == 0 | is.na(sd_x)] <- 1.0 
  
  x_mat_std <- scale(x_mat, center = mu_x, scale = sd_x)
  scale_vec <- c(1.0, sd_v, sd_x)
  
  # --- 2. Nuisance Parameter (eta) via cross-fitting ---
  if(verbose) cat(sprintf(">> [2/5] Estimating eta via %d-Fold Cross-Fitting...\n", K_folds))
  df_np <- data.frame(yn = y, v0 = v0_std, v1 = v1_std)
  if (d > 0) df_np <- cbind(df_np, as.data.frame(x_mat_std))
  fml_np <- as.formula(paste("yn ~", paste(names(df_np)[-1], collapse = " + ")))
  
  folds <- sample(rep(1:K_folds, length.out = n))
  eta <- numeric(n)
  
  if (!strict_bw) {
    if(verbose) cat("       (Using full-sample bandwidth for speed)\n")
    bw_obj <- np::npregbw(fml_np, data = df_np, regtype = "lc", bwmethod = "cv.aic") 
    raw_bws <- bw_obj$bw 
  }
  
  for(k in 1:K_folds) {
    train_df <- df_np[folds != k, , drop = FALSE]
    test_df  <- df_np[folds == k, , drop = FALSE]
    
    if (strict_bw) {
      bw_k <- np::npregbw(fml_np, data = train_df, regtype = "lc", bwmethod = "cv.aic")
      mod_k <- np::npreg(bws = bw_k$bw, txdat = train_df[,-1, drop=FALSE], tydat = train_df[,1])
    } else {
      mod_k <- np::npreg(bws = raw_bws, txdat = train_df[,-1, drop=FALSE], tydat = train_df[,1])
    }
    eta[folds == k] <- as.numeric(predict(mod_k, newdata = test_df[,-1, drop=FALSE]))
  }
  
  # --- 3. Sample Minimum (Multi-Start BOBYQA) ---
  if(verbose) cat(sprintf(">> [3/5] Finding global sample minimum (%d starts)...\n", n_starts))
  
  obj_fun_R <- function(par_std) {
    Q_obj_cpp(as.matrix(x_mat_std), as.numeric(v0_std), as.numeric(v1_std), as.numeric(par_std), as.numeric(eta))
  }
  
  opts_unconstr <- list("algorithm" = "NLOPT_LN_BOBYQA", "xtol_rel" = 1e-7, "maxeval" = 10000)
  best_Q <- Inf; theta_hat_std <- rep(0, p)
  
  for(i in 1:n_starts) {
    start_val <- if(i == 1) rep(0, p) else rnorm(p, 0, 0.5)
    opt <- nloptr(x0 = start_val, eval_f = obj_fun_R, opts = opts_unconstr)
    if(opt$objective < best_Q) {
      best_Q <- opt$objective
      theta_hat_std <- opt$solution
    }
  }
  Q_min_hat <- best_Q
  if(verbose) cat(sprintf("       Global Sample Minimum (Q_min_hat) = %.6f\n", Q_min_hat))
  
  # Calculate Raw Point Estimate
  theta_hat_raw <- numeric(p)
  theta_hat_raw[2] <- theta_hat_std[2] / sd_v
  if(d > 0) theta_hat_raw[3:p] <- theta_hat_std[3:p] / sd_x
  theta_hat_raw[1] <- theta_hat_std[1] - (theta_hat_raw[2] * mu_v) - sum(theta_hat_raw[3:p] * mu_x)
  names(theta_hat_raw) <- param_names
  
  # --- 4. Subsampling for CI ---
  # Chernozhucov et al 2018
  tau_CI <- NA
  if (B > 0) {
    if(verbose) cat(">> [4/5] Subsampling for critical values...\n")
    b_size <- floor(n^b_exponent); W_stats <- numeric(B)
    opts_sub <- list("algorithm" = "NLOPT_LN_BOBYQA", "xtol_rel" = 1e-5, "maxeval" = 1000)
    
    for(i in 1:B) {
      idx <- sample(0:(n-1), b_size, replace = FALSE)
      sub_obj <- function(p_std) Q_obj_subsample_cpp(as.matrix(x_mat_std), as.numeric(v0_std), as.numeric(v1_std), as.numeric(p_std), as.numeric(eta), idx)
      opt_sub <- nloptr(x0 = theta_hat_std, eval_f = sub_obj, opts = opts_sub)
      W_stats[i] <- b_size * (sub_obj(theta_hat_std) - opt_sub$objective)
    }
    tau_CI <- Q_min_hat + (quantile(W_stats, 1 - alpha, names = FALSE) / n)
  }
  
  tau_ID <- Q_min_hat + (log(n) / n)
  if(verbose) cat(sprintf("       Threshold (tau_ID) = %.6f\n", tau_ID))
  
  # --- 5. Estimation Branching ---
  bounds_ID_raw <- matrix(NA, p, 2, dimnames = list(param_names, c("Lower", "Upper")))
  bounds_CI_raw <- matrix(NA, p, 2, dimnames = list(param_names, c("Lower", "Upper")))
  
  if (method == "projection") {
    if(verbose) cat(">> [5/5] Running Directional Projection Method (AUGLAG)...\n")
    
    local_opts <- list("algorithm" = "NLOPT_LN_BOBYQA", "xtol_rel" = 1e-5)
    opts_proj <- list("algorithm" = "NLOPT_LN_AUGLAG", "local_opts" = local_opts, "xtol_rel" = 1e-5, "maxeval" = 5000)
    lb <- rep(-param_bounds, p); ub <- rep(param_bounds, p)
    
    run_proj <- function(obj_dir, tau) {
      opt <- nloptr(theta_hat_std, eval_f = obj_dir, eval_g_ineq = function(x) obj_fun_R(x)-tau, lb=lb, ub=ub, opts=opts_proj)
      if (obj_fun_R(opt$solution) > tau + 1e-4) return(NA) 
      return(opt$objective)
    }
    
    for(k in 1:p) {
      if (k == 1) {
        eval_dir <- function(x) {
          raw_slopes_x <- if(d > 0) (x[3:p] / sd_x) else 0
          x[1] - (x[2] / sd_v * mu_v) - sum(raw_slopes_x * mu_x)
        }
      } else if (k == 2) {
        eval_dir <- function(x) x[2] / sd_v
      } else {
        eval_dir <- function(x) x[k] / sd_x[k - 2]
      }
      
      bounds_ID_raw[k, 1] <- run_proj(eval_dir, tau_ID)
      bounds_ID_raw[k, 2] <- -run_proj(function(x) -eval_dir(x), tau_ID)
      
      if (B > 0) {
        bounds_CI_raw[k, 1] <- run_proj(eval_dir, tau_CI)
        bounds_CI_raw[k, 2] <- -run_proj(function(x) -eval_dir(x), tau_CI)
      }
    }
    
  } else {
    if(verbose) cat(sprintf(">> [5/5] Running Directional Profile Grid Search...\n"))
    opts_prof <- list("algorithm" = "NLOPT_LN_BOBYQA", "xtol_rel" = 1e-4, "maxeval" = 500)
    
    for(k in 1:p) {
      if(verbose) cat(sprintf("       Profiling %s...\n", param_names[k]))
      
      radius_raw <- grid_radius / scale_vec[k] 
      grid_right <- seq(theta_hat_raw[k], theta_hat_raw[k] + radius_raw, length.out = grid_points)
      grid_left  <- seq(theta_hat_raw[k], theta_hat_raw[k] - radius_raw, length.out = grid_points)
      
      q_prof_right <- numeric(grid_points); q_prof_left  <- numeric(grid_points)
      
      if (k == 1) {
        eval_prof <- function(slopes_std, v_raw) {
          raw_slopes_x <- if(d > 0) (slopes_std[-1] / sd_x) else 0
          beta0_std <- v_raw + (slopes_std[1] / sd_v * mu_v) + sum(raw_slopes_x * mu_x)
          return(obj_fun_R(c(beta0_std, slopes_std)))
        }
      } else {
        eval_prof <- function(other_params_std, v_raw) {
          v_std <- v_raw * scale_vec[k] 
          full_par_std <- numeric(p)
          full_par_std[-k] <- other_params_std
          full_par_std[k] <- v_std
          return(obj_fun_R(full_par_std))
        }
      }
      
      current_guess <- theta_hat_std[-k]
      for(i in 1:grid_points) {
        opt <- nloptr(x0 = current_guess, eval_f = function(x) eval_prof(x, grid_right[i]), opts = opts_prof)
        q_prof_right[i] <- opt$objective; current_guess <- opt$solution
      }
      
      current_guess <- theta_hat_std[-k]
      for(i in 1:grid_points) {
        opt <- nloptr(x0 = current_guess, eval_f = function(x) eval_prof(x, grid_left[i]), opts = opts_prof)
        q_prof_left[i] <- opt$objective; current_guess <- opt$solution
      }
      
      full_grid_k_raw <- c(rev(grid_left), grid_right[-1])
      full_q_prof <- c(rev(q_prof_left), q_prof_right[-1])
      
      valid_ID_raw <- full_grid_k_raw[full_q_prof <= tau_ID]
      if(length(valid_ID_raw) > 0) {
        bounds_ID_raw[k, ] <- c(min(valid_ID_raw), max(valid_ID_raw))
        if(min(valid_ID_raw) == grid_left[grid_points] || max(valid_ID_raw) == grid_right[grid_points]) {
          warning(sprintf("Identified set for %s hit the grid boundary. Increase grid_radius.", param_names[k]))
        }
      }
      
      if (B > 0) {
        valid_CI_raw <- full_grid_k_raw[full_q_prof <= tau_CI]
        if(length(valid_CI_raw) > 0) {
          bounds_CI_raw[k, ] <- c(min(valid_CI_raw), max(valid_CI_raw))
        }
      }
    }
  }
  
  res <- list(param_names=param_names, point_est=theta_hat_raw, bounds_ID=as.data.frame(bounds_ID_raw), 
              bounds_CI=as.data.frame(bounds_CI_raw), Q_min=Q_min_hat, thresholds=c(ID=tau_ID, CI=tau_CI), 
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
