library(lpSolve)
library(nloptr)
library(stats)

mmd_dro <- function(formula, data, v0_col, v1_col,
                    family = c("gaussian", "binomial"),
                    method = c("projection", "profile"),
                    scale = c("absolute", "proportion"),
                    delta = NULL,
                    B = 0,
                    b_exponent = 0.8,
                    alpha = 0.05,
                    grid_radius = 0.5,
                    grid_points = 100,
                    verbose = TRUE) {
  
  family <- match.arg(family)
  method <- match.arg(method)
  scale_type <- match.arg(scale)
  cl <- match.call()
  
  # --- 1. Data Preparation ---
  if(verbose) cat(">> [1/5] Preparing and validating data...\n")
  fml_full <- update(formula, paste("~ . +", v0_col, "+", v1_col))
  mf <- model.frame(fml_full, data)
  
  y <- as.numeric(model.response(mf))
  full_x_mat <- model.matrix(formula, mf)
  intercept_idx <- which(colnames(full_x_mat) == "(Intercept)")
  x_mat <- if(length(intercept_idx) > 0) full_x_mat[, -intercept_idx, drop=FALSE] else full_x_mat
  
  valid_rows <- rownames(mf)
  v0 <- as.numeric(data[valid_rows, v0_col])
  v1_max <- as.numeric(data[valid_rows, v1_col]) 
  
  n <- nrow(x_mat); d <- ncol(x_mat)
  param_names <- c("(Intercept)", paste0("latent_", v0_col), colnames(x_mat))
  p <- length(param_names)
  
  # --- 1b. INVISIBLE PURE STANDARDIZATION ---
  uncensored_v0 <- v0[v1_max == v0]
  mu_v <- mean(uncensored_v0, na.rm = TRUE)
  sd_v <- sd(uncensored_v0, na.rm = TRUE)
  if(sd_v == 0 || is.na(sd_v)) sd_v <- 1.0
  
  v0_std <- (v0 - mu_v) / sd_v
  v1_max_std <- (v1_max - mu_v) / sd_v
  
  mu_x <- apply(x_mat, 2, mean, na.rm = TRUE)
  sd_x <- apply(x_mat, 2, sd, na.rm = TRUE)
  sd_x[sd_x == 0 | is.na(sd_x)] <- 1.0
  x_mat_std <- scale(x_mat, center = mu_x, scale = sd_x)
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
  
  # --- 2. Core Estimation Loop over Delta ---
  results_list <- list()
  
  for (d_val in delta) {
    unobs <- v1_max - v0
    v1_act <- if (scale_type == "absolute") v0 + pmin(unobs, d_val) else v0 + d_val * unobs
    v1_act_s <- (v1_act - mu_v) / sd_v
    
    Z0 <- cbind(1, v0_std, x_mat_std)
    Z1 <- cbind(1, v1_act_s, x_mat_std)
    
    if (family == "gaussian") {
      # LP Setup: Decision variables = (theta_std [p], t [n])
      obj_coeff <- c(rep(0, p), rep(1/n, n))
      
      A_theta <- matrix(0, nrow = 4*n, ncol = p)
      A_theta[seq(1, 4*n, 4), ] <- Z0
      A_theta[seq(2, 4*n, 4), ] <- -Z0
      A_theta[seq(3, 4*n, 4), ] <- Z1
      A_theta[seq(4, 4*n, 4), ] <- -Z1
      
      A_t <- matrix(0, nrow = 4*n, ncol = n)
      A_t[matrix(c(1:(4*n), rep(1:n, each = 4)), ncol = 2)] <- 1
      A_full <- cbind(A_theta, A_t)
      
      rhs <- numeric(4*n)
      rhs[seq(1, 4*n, 4)] <- y
      rhs[seq(2, 4*n, 4)] <- -y
      rhs[seq(3, 4*n, 4)] <- y
      rhs[seq(4, 4*n, 4)] <- -y
      
      # Solve Global Minimum
      opt_min <- lpSolve::lp("min", obj_coeff, A_full, ">=", rhs, unconst.cols = 1:p)
      Q_min <- opt_min$objval
      theta_hat_std <- opt_min$solution[1:p]
      
      # Calculate Threshold
      tau_ID <- Q_min + (log(n) / n)
      
      bounds_ID <- matrix(NA, p, 2, dimnames = list(param_names, c("Lower", "Upper")))
      
      if (method == "projection") {
        # --- PROJECTION METHOD (DIRECT LP) ---
        A_proj <- rbind(A_full, c(rep(0, p), rep(1, n)))
        rhs_proj <- c(rhs, tau_ID * n)
        dir_proj <- c(rep(">=", 4*n), "<=")
        
        for (k in 1:p) {
          c_proj <- c(c_proj_list[[k]], rep(0, n))
          opt_L <- lpSolve::lp("min", c_proj, A_proj, dir_proj, rhs_proj, unconst.cols = 1:p)
          bounds_ID[k, 1] <- opt_L$objval
          opt_U <- lpSolve::lp("max", c_proj, A_proj, dir_proj, rhs_proj, unconst.cols = 1:p)
          bounds_ID[k, 2] <- opt_U$objval
        }
      } else {
        # --- PROFILE GRID METHOD (FAST LP SWEEPS) ---
        # Get raw estimates for centering the grid
        theta_hat_raw <- numeric(p)
        theta_hat_raw[2] <- theta_hat_std[2] / sd_v
        if(d > 0) theta_hat_raw[3:p] <- theta_hat_std[3:p] / sd_x
        theta_hat_raw[1] <- theta_hat_std[1] - (theta_hat_raw[2] * mu_v) - sum(theta_hat_raw[3:p] * mu_x)
        
        # We append: sum(t_i)/n <= tau_ID AND theta_k = g (the grid point)
        # Note: theta_k_std = c_proj_std * x
        A_prof <- rbind(A_full, c(rep(0, p), rep(1, n)), rep(0, p + n)) # placeholder last row
        dir_prof <- c(rep(">=", 4*n), "<=", "=")
        
        for (k in 1:p) {
          radius_raw <- grid_radius / scale_vec[k]
          grid_k <- seq(theta_hat_raw[k] - radius_raw, theta_hat_raw[k] + radius_raw, length.out = grid_points)
          valid_pts <- c()
          
          # Set up the scaling mapping for the constraint row
          c_map_std <- rep(0, p + n)
          if (k == 1) {
            c_map_std[1] <- 1
            c_map_std[2] <- -mu_v / sd_v
            if (d > 0) c_map_std[3:p] <- -mu_x / sd_x
          } else if (k == 2) {
            c_map_std[2] <- 1 / sd_v
          } else {
            c_map_std[k] <- 1 / sd_x[k - 2]
          }
          A_prof[4*n + 2, ] <- c_map_std
          
          for (g in grid_k) {
            rhs_prof <- c(rhs, tau_ID * n, g)
            opt_prof <- lpSolve::lp("min", obj_coeff, A_prof, dir_prof, rhs_prof, unconst.cols = 1:p)
            if (opt_prof$status == 0) valid_pts <- c(valid_pts, g)
          }
          
          bounds_ID[k, 1] <- if (length(valid_pts) > 0) min(valid_pts) else NA
          bounds_ID[k, 2] <- if (length(valid_pts) > 0) max(valid_pts) else NA
        }
      }
      
      theta_hat_raw <- numeric(p)
      theta_hat_raw[2] <- theta_hat_std[2] / sd_v
      if(d > 0) theta_hat_raw[3:p] <- theta_hat_std[3:p] / sd_x
      theta_hat_raw[1] <- theta_hat_std[1] - (theta_hat_raw[2] * mu_v) - sum(theta_hat_raw[3:p] * mu_x)
      names(theta_hat_raw) <- param_names
      
      results_list[[as.character(d_val)]] <- list(
        bounds_ID = as.data.frame(bounds_ID),
        point_est = theta_hat_raw,
        Q_min = Q_min,
        delta = d_val
      )
    }
  }
  
  res <- list(results=results_list, formula=formula, family=family, scale=scale_type, delta=delta, n=n, B=B, call=cl)
  class(res) <- "mmd_results"
  return(res)
}
