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

  generate_data <- function(J, T_full, T_obs, beta_0, beta_1, beta_2, # {{{
                            gamma_1, gamma_2, gamma_3, rho_x1_t) {
  
    T_start <- T_full - T_obs + 1
  
    # ── Country-level latent variable that creates X1-duration correlation ──
    Z <- rnorm(J, 0, 1)
  
    # ── Generate X1 (correlated with duration via Z) and X2 (independent) ──
    # X1 depends on Z: countries with higher Z have higher X1
    # Z also affects the hazard (through X1), creating the correlation
    X1 <- rho_x1_t * Z + sqrt(1 - rho_x1_t^2) * rnorm(J, 0, 1)
    X2 <- rnorm(J, 0, 1)
  
    # ── Simulate event histories for each country ──
    # We build the panel year by year: for each country-year, compute Pr(event),
    # draw the outcome, and update the duration counter.
  
    records <- list()
  
    for (j in 1:J) {
      t_star <- 0  # True time since last event (or since country "birth")
  
      for (year in 1:T_full) {
        # Linear predictor
        linpred <- beta_0 + beta_1 * X1[j] + beta_2 * X2[j] +
                   gamma_1 * t_star + gamma_2 * t_star^2 + gamma_3 * t_star^3
  
        # Probability of event (logistic CDF)
        prob <- 1 / (1 + exp(-linpred))
  
        # Draw outcome
        y <- rbinom(1, 1, prob)
  
        # Record this observation
        records[[length(records) + 1]] <- data.frame(
          country = j, year = year, y = y,
          X1 = X1[j], X2 = X2[j],
          t_star = t_star,   # True duration-to-date
          Z = Z[j]
        )
  
        # Update duration counter
        if (y == 1) {
          t_star <- 0  # Reset after event
        } else {
          t_star <- t_star + 1
        }
      }
    }
  
    full_data <- bind_rows(records)
  
    # ── Now create the "observed" data: only years T_start to T_full ──
    obs_data <- full_data %>% filter(year >= T_start)
  
    # ── Code the duration as the researcher would ──
    # The researcher starts t at 0 for each country in the first observed year,
    # and resets to 0 after each observed event.
    obs_data <- obs_data %>%
      group_by(country) %>%
      mutate(
        # Identify events in the observed window
        event_cumsum = cumsum(lag(y, default = 0)),
        # For the first spell (before any observed event), t starts at 0
        # For subsequent spells, t resets after each event
        spell_id = event_cumsum,
        t_coded = ave(seq_along(y), spell_id, FUN = function(x) seq_along(x) - 1)
      ) %>%
      ungroup()
  
    # ── Identify left-censored spells ──
    # A spell is left-censored if it's the first spell for a country AND
  
    # the country existed before T_start AND no event happened in year T_start-1
    # (i.e., we don't know when the previous event was).
    # In our setup, every country exists from year 1, so the first spell is
    # censored unless an event happened in year T_start - 1 or the country
    # had its time counter at 0 in year T_start.
  
    obs_data <- obs_data %>%
      group_by(country) %>%
      mutate(
        is_first_spell = (spell_id == 0),
        # The spell is censored if t* > t in the first observed year
        first_year_t_star = t_star[1],
        first_year_t_coded = t_coded[1],
        is_left_censored = is_first_spell & (first_year_t_star > first_year_t_coded)
      ) %>%
      ungroup()
  
    # ── Compute spell-level information for the auxiliary model ──
    # For uncensored spells: we know the full spell length
    # For left-censored spells: we only know the coded portion
  
    spell_info <- obs_data %>%
      group_by(country, spell_id) %>%
      summarize(
        X1 = first(X1),
        X2 = first(X2),
        Z  = first(Z),
        is_left_censored = first(is_left_censored),
        spell_length_coded = max(t_coded) + 1,  # Coded spell length
        spell_length_true  = max(t_star) - min(t_star) + 1, # True spell length within window
        # True total spell length (including pre-window portion)
        t_star_at_start = min(t_star),
        t_star_at_end   = max(t_star),
        ended_with_event = last(y) == 1,
        # The true measurement error for this spell
        tau = ifelse(is_left_censored, first(t_star) - first(t_coded), 0),
        .groups = "drop"
      )
  
    # For uncensored complete spells, the true spell length = t* at event + 1
    # We need complete uncensored spells for the auxiliary model
    spell_info <- spell_info %>%
      mutate(
        # Total true duration of the spell (for uncensored spells)
        true_total_duration = t_star_at_end + 1,
        # Is this spell "complete" (ended with an event)?
        is_complete = ended_with_event,
        # For the survival model: observed duration and censoring indicator
        surv_time = ifelse(!is_left_censored,
                           spell_length_coded,
                           spell_length_coded),
        surv_event = as.integer(ended_with_event)
      )
  
    # Country age (years since "birth" = year 1 in our simulation)
    obs_data <- obs_data %>%
      mutate(country_age = T_full)  # All countries exist from year 1
  
    return(list(
      obs_data   = obs_data,
      full_data  = full_data,
      spell_info = spell_info,
      X1 = X1, X2 = X2
    ))
  }
  # }}}


  # Get a Sample
  set.seed(63130)
  dat <- generate_data(J, T_full, T_obs, beta_0, beta_1, beta_2,
                       gamma_1, gamma_2, gamma_3, rho_x1_t)
  
  obs  <- dat$obs_data
  spells <- dat$spell_info
  
  cat("── Dataset Summary ──\n")
  cat("Total observations:", nrow(obs), "\n")
  cat("Events (y=1):", sum(obs$y), "(", round(100*mean(obs$y), 1), "%)\n")
  cat("Left-censored observations:", sum(obs$is_left_censored),
      "(", round(100*mean(obs$is_left_censored), 1), "%)\n")
  cat("\n── Spell Summary ──\n")
  cat("Total spells:", nrow(spells), "\n")
  cat("Left-censored spells:", sum(spells$is_left_censored), "\n")
  cat("Uncensored complete spells:", sum(!spells$is_left_censored & spells$is_complete), "\n")
  cat("Uncensored right-censored spells:", sum(!spells$is_left_censored & !spells$is_complete), "\n")


  fit <- MMD_bounds_cpp(y ~ X1 + X2 + X3 + X4 + X5, 
    data = pop, v0_col = "v0", v1_col = "v1")
  summary(fit)
}
