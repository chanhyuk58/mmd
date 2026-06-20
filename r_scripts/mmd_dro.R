library(lpSolve)
library(nloptr)
library(stats)
library(graphics)
library(grDevices)

# @param formula A formula object.
# @param data A data.frame containing the variables.
# @param v0_col String. Name of lower bound column.
# @param v1_col String. Name of upper bound column.
# @param family String. "gaussian" (LPM/Linear) or "binomial" (Logit).
# @param scale String. "absolute" (years/units cap) or "proportion" (0 to 1 scaling).
# @param delta Numeric vector of sensitivity values. Defaults to most conservative if NULL.
# @param B Integer. Number of bootstrap replicates for CI (0 to skip).
# @param b_exponent Exponent for subsample size.
# @param alpha Significance level.
# @param n_starts Number of random starts for global min.
# @param grid_radius Radius in standardized space.
# @param grid_points Points per parameter for profiling.
# @param param_bounds Absolute bounds for parameters to prevent infinity.
# @param verbose Logical. Print progress?

mmd_dro <- function(formula, data, v0_col, v1_col,
                    family = c("gaussian", "binomial"),
                    method = c("projection", "profile"),
                    scale = c("absolute", "proportion"),
                    delta = NULL,
                    B = 0,
                    b_exponent = 0.8,
                    alpha = 0.05,
                    n_starts = 5,
                    grid_radius = 0.5,
                    grid_points = 100,
                    param_bounds = 20,
                    verbose = TRUE) {
  
  family <- match.arg(family)
  method <- match.arg(method)
  scale_type <- match.arg(scale)
  cl <- match.call()
  
  # --- 1. Data Preparation (NA-Safe) ---
  if(verbose) cat(">> [1/5] Preparing and validating data...\n")
  fml_full <- update(formula, paste("~ . +", v0_col, "+", v1_col))
  mf <- model.frame(fml_full, data)
  
  y <- as.numeric(model.response(mf))
  full_x_mat <- model.matrix(formula, mf)
  intercept_idx <- which(colnames(full_x_mat) == "(Intercept)")
  x_mat <- if(length(intercept_idx) > 0) full_x_mat[, -intercept_idx, drop=FALSE] else full_x_mat
  
  valid_idx <- as.integer(rownames(mf))
  v0 <- as.numeric(data[[v0_col]][valid_idx])
  v1_max <- as.numeric(data[[v1_col]][valid_idx]) 
  
  n <- nrow(x_mat); d <- ncol(x_mat)
  param_names <- c("(Intercept)", paste0("latent_", v0_col), colnames(x_mat))
  p <- length(param_names)
  
  # --- 1b. INVISIBLE PURE STANDARDIZATION ---
  uncensored_v0 <- v0[v1_max == v0]
  
  if (length(uncensored_v0) > 1 && !is.na(sd(uncensored_v0, na.rm = TRUE)) && sd(uncensored_v0, na.rm = TRUE) > 0) {
    mu_v <- mean(uncensored_v0, na.rm = TRUE)
    sd_v <- sd(uncensored_v0, na.rm = TRUE)
  } else {
    mu_v <- mean(v0, na.rm = TRUE)
    sd_v <- sd(v0, na.rm = TRUE)
  }
  if(sd_v == 0 || is.na(sd_v)) sd_v <- 1.0
  
  v0_std <- (v0 - mu_v) / sd_v
  v1_max_std <- (v1_max - mu_v) / sd_v
  
  if (d > 0) {
    mu_x <- colMeans(x_mat, na.rm = TRUE)
    sd_x <- sapply(as.data.frame(x_mat), sd, na.rm = TRUE)
    sd_x[sd_x == 0 | is.na(sd_x)] <- 1.0
    x_mat_std <- scale(x_mat, center = mu_x, scale = sd_x)
  } else {
    mu_x <- numeric(0)
    sd_x <- numeric(0)
    x_mat_std <- x_mat
  }
  scale_vec <- c(1.0, sd_v, sd_x)
  
  # Set up directional projection vectors
  c_proj_list <- lapply(1:p, function(k) {
    c_proj <- rep(0, p)
    if (k == 1) {
      raw_slopes_x <- if(d > 0) (1 / sd_x) else 0
      c_proj[1] <- 1
      c_proj[2] <- -mu_v / sd_v
      if (d > 0) c_proj[3:p] <- -mu_x / sd_x
    } else if (k == 2) {
      c_proj[2] <- 1 / sd_v
    } else {
      c_proj[k] <- 1 / sd_x[k - 2]
    }
    return(c_proj)
  })
  
  if (is.null(delta)) {
    delta <- if (scale_type == "absolute") Inf else 1.0
  }
  
  # --- 1c. Vectorized L2 DRO Loss Function (For Binomial/Logit) ---
  Q_dro_R <- function(par_std, v0_s, v1_max_s, x_s, y_val, d_val) {
    g0 <- par_std[1]
    gV <- par_std[2]
    fx <- if (ncol(x_s) > 0) as.numeric(x_s %*% par_std[3:length(par_std)]) else 0
    
    unobs <- v1_max - v0
    v1_act <- if (scale_type == "absolute") v0 + pmin(unobs, d_val) else v0 + d_val * unobs
    v1_act_s <- (v1_act - mu_v) / sd_v
    
    f0 <- g0 + gV * v0_s + fx
    f1 <- g0 + gV * v1_act_s + fx
    
    L0 <- (y_val - f0)^2
    L1 <- (y_val - f1)^2
    return(sum(pmax(L0, L1)) / length(y_val))
  }
  
  # --- 2. Core Estimation Solver ---
  solve_dro_bounds <- function(v0_s, v1_s, x_s, y_val, d_val) {
    
    # ==========================================================================
    # --- GAUSSIAN FAMILY (LINEAR PROGRAMMING) ---
    # ==========================================================================
    if (family == "gaussian") {
      # Decision variables: x = (theta_pos [p], theta_neg [p], t [n])
      # Objective: Minimize sum(t_i)/n + 1e-6 * sum(theta_pos + theta_neg)
      obj_coeff <- c(rep(1e-6, 2*p), rep(1/n, n))
      
      # Construct LP Constraint Matrix (4n x (2p+n))
      A_theta <- matrix(0, nrow = 4*n, ncol = p)
      unobs <- v1_max - v0
      v1_act <- if (scale_type == "absolute") v0 + pmin(unobs, d_val) else v0 + d_val * unobs
      v1_act_s <- (v1_act - mu_v) / sd_v
      
      Z0 <- cbind(1, v0_s, x_s)
      Z1 <- cbind(1, v1_act_s, x_s)
      
      A_theta[seq(1, 4*n, 4), ] <- Z0
      A_theta[seq(2, 4*n, 4), ] <- -Z0
      A_theta[seq(3, 4*n, 4), ] <- Z1
      A_theta[seq(4, 4*n, 4), ] <- -Z1
      
      A_t <- matrix(0, nrow = 4*n, ncol = n)
      A_t[matrix(c(1:(4*n), rep(1:n, each = 4)), ncol = 2)] <- 1
      
      A_full <- cbind(A_theta, -A_theta, A_t)
      
      rhs <- numeric(4*n)
      rhs[seq(1, 4*n, 4)] <- y_val
      rhs[seq(2, 4*n, 4)] <- -y_val
      rhs[seq(3, 4*n, 4)] <- y_val
      rhs[seq(4, 4*n, 4)] <- -y_val
      
      # All 4n constraints are strictly ">=" in the epigraph formulation
      const_dir <- rep(">=", 4*n)
      
      # Solve Global Minimum
      opt_min <- lpSolve::lp("min", obj_coeff, A_full, const_dir, rhs)
      Q_min <- opt_min$objval
      theta_hat_std <- opt_min$solution[1:p] - opt_min$solution[(p+1):(2*p)]
      
      # Calculate Threshold
      tau_ID <- Q_min + (log(n) / n)
      
      bounds <- matrix(NA, p, 2, dimnames = list(param_names, c("Lower", "Upper")))
      
      if (method == "projection") {
        # --- PROJECTION METHOD ---
        A_proj <- rbind(A_full, c(rep(0, 2*p), rep(1, n)))
        rhs_proj <- c(rhs, tau_ID * n)
        dir_proj <- c(const_dir, "<=")
        
        for (k in 1:p) {
          c_proj <- c(c_proj_list[[k]], -c_proj_list[[k]], rep(0, n))
          
          opt_L <- lpSolve::lp("min", c_proj, A_proj, dir_proj, rhs_proj)
          bounds[k, 1] <- if (opt_L$status == 0) opt_L$objval else NA
          
          opt_U <- lpSolve::lp("max", c_proj, A_proj, dir_proj, rhs_proj)
          bounds[k, 2] <- if (opt_U$status == 0) opt_U$objval else NA
        }
      } else {
        # --- PROFILE METHOD ---
        theta_hat_raw <- numeric(p)
        theta_hat_raw[2] <- theta_hat_std[2] / sd_v
        if(d > 0) theta_hat_raw[3:p] <- theta_hat_std[3:p] / sd_x
        theta_hat_raw[1] <- theta_hat_std[1] - (theta_hat_raw[2] * mu_v) - sum(theta_hat_raw[3:p] * mu_x)
        
        # Row 1 is the threshold constraint. Row 2 is the parameter equality.
        A_prof <- rbind(A_full, c(rep(0, 2*p), rep(1, n)), rep(0, 2*p + n)) 
        dir_prof <- c(const_dir, "<=", "=")
        
        for (k in 1:p) {
          radius_raw <- grid_radius / scale_vec[k]
          grid_k <- seq(theta_hat_raw[k] - radius_raw, theta_hat_raw[k] + radius_raw, length.out = grid_points)
          valid_pts <- c()
          
          A_prof[4*n + 2, ] <- c(c_proj_list[[k]], -c_proj_list[[k]], rep(0, n))
          
          for (g in grid_k) {
            rhs_prof <- c(rhs, tau_ID * n, g)
            opt_prof <- lpSolve::lp("min", obj_coeff, A_prof, dir_prof, rhs_prof)
            if (opt_prof$status == 0) valid_pts <- c(valid_pts, g)
          }
          
          bounds[k, 1] <- if (length(valid_pts) > 0) min(valid_pts) else NA
          bounds[k, 2] <- if (length(valid_pts) > 0) max(valid_pts) else NA
        }
      }
      
      # Convert point estimate back to raw scale
      theta_hat_raw <- numeric(p)
      theta_hat_raw[2] <- theta_hat_std[2] / sd_v
      if(d > 0) theta_hat_raw[3:p] <- theta_hat_std[3:p] / sd_x
      theta_hat_raw[1] <- theta_hat_std[1] - (theta_hat_raw[2] * mu_v) - sum(theta_hat_raw[3:p] * mu_x)
      names(theta_hat_raw) <- param_names
      
      return(list(bounds = bounds, theta_raw = theta_hat_raw, Q_min = Q_min, tau_ID = tau_ID))
    }
    
    # ==========================================================================
    # --- BINOMIAL FAMILY (NON-LINEAR PROGRAMMING) ---
    # ==========================================================================
    if (family == "binomial") {
      # Non-linear optimizer logic placeholder
    }
  }
  
  # --- 3. Estimation over Delta ---
  if(verbose) cat(">> [2/4] Solving direct DRO bounds...\n")
  results_list <- list()
  
  for (d_val in delta) {
    fit <- solve_dro_bounds(v0_std, v1_max_std, x_mat_std, y, d_val)
    results_list[[as.character(d_val)]] <- list(
      bounds_ID = as.data.frame(fit$bounds), 
      point_est = fit$theta_raw,
      Q_min = fit$Q_min,
      tau_ID = fit$tau_ID,
      delta = d_val
    )
  }
  
  # --- 4. Parallel-Safe Bootstrap for Confidence Intervals ---
  if (B > 0) {
    if(verbose) cat(sprintf(">> [3/4] Bootstrapping bounds (B = %d)...\n", B))
    boot_bounds <- array(NA, dim = c(p, 2, B, length(delta)))
    
    for (b in 1:B) {
      boot_idx <- sample(1:n, n, replace = TRUE)
      y_b <- y[boot_idx]
      v0_s_b <- v0_std[boot_idx]
      v1_s_b <- v1_max_std[boot_idx]
      x_s_b <- x_mat_std[boot_idx, , drop=FALSE]
      
      for (d_idx in seq_along(delta)) {
        fit_b <- solve_dro_bounds(v0_s_b, v1_s_b, x_s_b, y_b, delta[d_idx])
        boot_bounds[, , b, d_idx] = fit_b$bounds
      }
    }
    
    for (d_idx in seq_along(delta)) {
      bounds_CI <- matrix(NA, p, 2, dimnames = list(param_names, c("Lower", "Upper")))
      for (k in 1:p) {
        bounds_CI[k, 1] <- quantile(boot_bounds[k, 1, , d_idx], probs = alpha / 2, na.rm = TRUE)
        bounds_CI[k, 2] <- quantile(boot_bounds[k, 2, , d_idx], probs = 1 - alpha / 2, na.rm = TRUE)
      }
      results_list[[as.character(delta[d_idx])]]$bounds_CI <- as.data.frame(bounds_CI)
    }
  }
  
  res <- list(results=results_list, formula=formula, family=family, scale=scale_type, delta=delta, n=n, B=B, call=cl)
  class(res) <- "mmd_results"
  return(res)
}


# ------------------------------------------------------------------------------
# S3 Print and Summary Methods
# ------------------------------------------------------------------------------

print.mmd_results <- function(x, ...) {
  cat("\nMMD-DRO Robust Estimator (Linear Programming)\n")
  cat("--------------------------------------------\n")
  cat("Scale Type:  ", x$scale, "\n")
  cat("Sample Size: ", x$n, "\n\n")
  
  max_delta <- as.character(max(x$delta))
  cat("Estimated Bounds (at Delta =", max_delta, "):\n")
  print(round(x$results[[max_delta]]$bounds_ID, 4))
  cat("\n")
}

summary.mmd_results <- function(object, ...) {
  summary_list <- lapply(object$results, function(res) {
    df <- res$bounds_ID
    df$Parameter <- rownames(df)
    df$Delta <- res$delta
    if (!is.null(res$bounds_CI)) {
      df$CI_Lower <- res$bounds_CI$Lower
      df$CI_Upper <- res$bounds_CI$Upper
    }
    return(df)
  })
  
  summary_table <- do.call(rbind, summary_list)
  rownames(summary_table) <- NULL
  
  ans <- list(
    table = summary_table,
    n = object$n,
    scale = object$scale,
    B = object$B,
    call = object$call
  )
  class(ans) <- "summary.mmd_results"
  return(ans)
}

print.summary.mmd_results <- function(x, ...) {
  cat("\nMMD-DRO SUMMARY TABLE\n")
  cat("=================================================================\n")
  cat("Scale: ", x$scale, " | n: ", x$n, "\n\n")
  if (x$B > 0) {
    print(round(x$table[, c("Parameter", "Delta", "Lower", "Upper", "CI_Lower", "CI_Upper")], 4))
  } else {
    print(round(x$table[, c("Parameter", "Delta", "Lower", "Upper")], 4))
  }
  cat("=================================================================\n\n")
}

# ------------------------------------------------------------------------------
# S3 Plot Method (The Censoring Sensitivity Plot)
# ------------------------------------------------------------------------------

plot.mmd_results <- function(x, parameter, ...) {
  if (length(x$delta) < 2) {
    stop("Error: Plotting requires a vector of Delta values (sensitivity sweep).")
  }
  
  sum_tab <- summary(x)$table %>% filter(Parameter == parameter)
  
  if (nrow(sum_tab) == 0) {
    stop("Error: Parameter not found in results. Available: ", paste(unique(summary(x)$table$Parameter), collapse=", "))
  }
  
  # Standardize label depending on scale type
  x_label <- if (x$scale == "absolute") {
    expression(paste("Max Allowed Unobserved Peace (", delta, " years)"))
  } else {
    expression(paste("Proportion of Lifespan Unobserved (", delta, ")"))
  }
  
  plot(NULL, xlim = range(x$delta), ylim = range(c(sum_tab$Lower, sum_tab$Upper), na.rm=TRUE),
       xlab = x_label,
       ylab = paste("Causal Bounds for", parameter),
       main = paste("Censoring Sensitivity Plot:", parameter),
       bty = "n", las = 1, cex.lab = 1.1)
  
  # Shaded identified set polygon
  polygon(c(sum_tab$Delta, rev(sum_tab$Delta)), 
          c(sum_tab$Lower, rev(sum_tab$Upper)), 
          col = rgb(0.2, 0.4, 0.8, 0.2), border = NA)
  
  # Bound lines
  lines(sum_tab$Delta, sum_tab$Lower, col = "royalblue4", lwd = 2)
  lines(sum_tab$Delta, sum_tab$Upper, col = "royalblue4", lwd = 2)
  
  # If CI is available, plot it as dotted outer lines
  if (!is.null(sum_tab$CI_Lower)) {
    lines(sum_tab$Delta, sum_tab$CI_Lower, col = "red", lty = 3, lwd = 1.5)
    lines(sum_tab$Delta, sum_tab$CI_Upper, col = "red", lty = 3, lwd = 1.5)
  }
  
  # Dashed line at 0 (The significance threshold)
  abline(h = 0, col = "grey50", lty = 3, lwd = 1.5)
  
  # Legend
  legend_labels <- c("Identified Set", "Significance Limit (0)")
  legend_cols <- c("royalblue4", "grey50")
  legend_ltys <- c(1, 3)
  legend_fills <- c(rgb(0.2, 0.4, 0.8, 0.2), NA)
  
  if (!is.null(sum_tab$CI_Lower)) {
    legend_labels <- c(legend_labels, "95% Confidence Interval")
    legend_cols <- c(legend_cols, "red")
    legend_ltys <- c(legend_ltys, 3)
    legend_fills <- c(legend_fills, NA)
  }
  
  legend("topright", legend = legend_labels, fill = legend_fills, border = NA,
         col = legend_cols, lty = legend_ltys, lwd = 1.5, bty = "n", cex = 0.9)
}
