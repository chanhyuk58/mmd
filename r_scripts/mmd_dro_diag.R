library(dplyr)
library(stats)

#' MMD Pre-Flight Diagnostics
#'
#' @param formula A formula object.
#' @param data A data.frame containing the variables.
#' @param v0_col String. Name of lower bound column.
#' @param v1_col String. Name of upper bound column.
#'
#' @return An S3 object of class 'mmd_diag'
#' @export

mmd_diagnose <- function(formula, data, v0_col, v1_col) {
  
  cl <- match.call()
  
  # 1. Data Preparation
  fml_full <- update(formula, paste("~ . +", v0_col, "+", v1_col))
  mf <- model.frame(fml_full, data)
  
  full_x_mat <- model.matrix(formula, mf)
  intercept_idx <- which(colnames(full_x_mat) == "(Intercept)")
  x_mat <- if(length(intercept_idx) > 0) full_x_mat[, -intercept_idx, drop=FALSE] else full_x_mat
  
  v0 <- as.numeric(mf[[v0_col]])
  v1 <- as.numeric(mf[[v1_col]])
  n <- length(v0)
  
  width <- v1 - v0
  if (any(width < 0)) {
    stop("Data Error: Lower bound v0 exceeds upper bound v1 for some observations.")
  }
  
  # Manski & Tamer (2002, Econometrica, Eq 16): Interval width directly governs set identification.
  cens_rate <- mean(width > 0)
  avg_width_full <- mean(width)
  avg_width_cens <- if (sum(width > 0) > 0) mean(width[width > 0]) else 0
  max_width <- max(width)
  
  # Conditional Information Slack Index (ISI): Noise-to-Signal ratio.
  # Normalizes expected censoring interval by the standard deviation 
  # of the uncorrupted, exact latent observations.
  sd_v0_uncens <- if (sum(width == 0) > 1) sd(v0[width == 0]) else 0
  isi <- if (sd_v0_uncens > 0) avg_width_cens / sd_v0_uncens else Inf
  
  # Multicollinearity (Kappa): Bontemps, Magnac, and Maurin (2012, Econometrica) prove 
  # that the covariance structure of covariates directly inflates the identified set bounds.
  X_v0 <- cbind(v0 = v0, x_mat)
  R <- cor(X_v0, use = "pairwise.complete.obs")
  eig_vals <- eigen(R)$values
  kappa <- if (min(eig_vals) > 0) sqrt(max(eig_vals) / min(eig_vals)) else Inf
  
  res <- list(
    n = n,
    cens_rate = cens_rate,
    avg_width_full = avg_width_full,
    avg_width_cens = avg_width_cens,
    max_width = max_width,
    isi = isi,
    kappa = kappa,
    corr_matrix = R,
    call = cl
  )
  
  class(res) <- "mmd_diag"
  return(res)
}

# ------------------------------------------------------------------------------
# S3 Print Method
# ------------------------------------------------------------------------------

print.mmd_diag <- function(x, ...) {
  cat("\n=================================================================\n")
  cat(" MMD PRE-FLIGHT DATA DIAGNOSTICS\n")
  cat("=================================================================\n")
  cat("Sample Size (Complete Cases): ", x$n, "\n\n")
  
  cat("1. CENSORING STRUCTURE:\n")
  cat("  Censoring Rate (v1 > v0):   ", sprintf("%.2f%%", x$cens_rate * 100), "\n")
  cat("  Average Censoring Width:    ", sprintf("%.2f years", x$avg_width_full), "\n")
  cat("  Avg Width (Censored Only):  ", sprintf("%.2f years", x$avg_width_cens), "\n")
  cat("  Maximum Censoring Width:    ", sprintf("%.2f years", x$max_width), "\n\n")
  
  cat("2. IDENTIFICATION STRENGTH:\n")
  cat("  Information Slack Index (ISI):", sprintf("%.4f", x$isi), "\n")
  cat("  System Collinearity (Kappa):  ", sprintf("%.4f", x$kappa), "\n\n")
  
  cat("3. METHODOLOGICAL PRE-DIAGNOSIS & ADVICE:\n")
  
  if (x$cens_rate > 0.5) {
    cat("  [!] WARNING: Very high censoring (>50%). Expect wide parameter bounds.\n")
  } else {
    cat("  [+] GOOD: Censoring rate is under 50%. Uncensored cases will help anchor the bounds.\n")
  }
  
  if (x$isi > 1.0) {
    cat("  [!] WARNING: ISI > 1.0. Unobserved history is larger than observed variation.\n")
    cat("      Bounds for the latent variable and intercept are likely to cross zero.\n")
  } else {
    cat("  [+] GOOD: ISI < 1.0. Observed variation dominates. Expect reasonably tight bounds.\n")
  }
  
  if (x$kappa >= 30) {
    cat("  [!] WARNING: Severe multicollinearity (Kappa >= 30). This will dramatically\n")
    cat("      inflate the bounds because the parameters can easily compensate for each other.\n")
    cat("      Consider dropping redundant covariates.\n")
  } else if (x$kappa >= 15) {
    cat("  [*] NOTE: Moderate multicollinearity (15 <= Kappa < 30). Bounds may be slightly inflated.\n")
  } else {
    cat("  [+] GOOD: Low multicollinearity (Kappa < 15). Optimizer will be highly stable.\n")
  }
  
  cat("=================================================================\n\n")
}
