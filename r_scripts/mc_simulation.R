library(doFuture)
library(foreach)
library(dplyr)
library(Rcpp)
library(splines)
library(nloptr)
library(ranger)

cat("Packages loaded ...\n")

mc_reps <- 100

cat(">> Master: Compiling C++ backend...\n")
Rcpp::sourceCpp("./mmd_cpp.cpp", cacheDir = tempdir())
cat(">> Master: Compilation complete.\n")

lsb_hosts <- Sys.getenv("LSB_HOSTS")
n_cores <- if (lsb_hosts != "") {
  length(strsplit(lsb_hosts, " ")[[1]])
} else {
  max(1, parallel::detectCores() - 1)
}

registerDoFuture()
plan(multisession, workers = n_cores)
cat(sprintf(">> Registered %d workers.\n", n_cores))

cat(">> Starting Monte Carlo...\n")

set.seed(4875995)
results_list <- foreach(i = 1:mc_reps,
                        .options.future = list(seed = TRUE),
                        .errorhandling = "stop") %dofuture% {

  Rcpp::sourceCpp("./mmd_cpp.cpp", cacheDir = tempdir())
  source("./mmd_cpp.R")
  source("./generate_pop.R")

  sim <- gen_pop(
    J = 80,
    T_full = 80,
    birth_range = c(1, 40),
    obs_start_range = c(40, 70),
    mean_gdp = 8.5,
    mean_pop = 16.0,
    gdp_shock = 0.05
  )
  pop <- sim$data
  true_params <- sim$true_params

  # --- 1. Estimation with Cross-fitted Eta (ranger) ---
  fit_estimated <- MMD_bounds(
    onset ~ log_gdp + democ + eth_het + log_pop + log_ref,
    data = pop, v0_col = "v0", v1_col = "v1",
    method = "projection",
    B = 0,
    eta_method = "ranger",
    family = "lpm",
    verbose = FALSE
  )

  # --- 2. Estimation with Oracle (Analytical) Eta ---
  fit_oracle <- MMD_bounds(
    onset ~ log_gdp + democ + eth_het + log_pop + log_ref,
    data = pop, v0_col = "v0", v1_col = "v1",
    method = "projection",
    B = 0,
    eta_oracle = pop$eta_oracle,
    family = "lpm",
    verbose = FALSE
  )

  rep_res <- data.frame(
    iteration = i,
    param     = names(true_params),
    truth     = as.numeric(true_params),
    
    est_point  = fit_estimated$point_est,
    est_low    = fit_estimated$bounds_ID$Lower,
    est_upp    = fit_estimated$bounds_ID$Upper,
    
    orc_point  = fit_oracle$point_est,
    orc_low    = fit_oracle$bounds_ID$Lower,
    orc_upp    = fit_oracle$bounds_ID$Upper
  )

  rep_res <- rep_res %>%
    mutate(
      est_covered = (truth >= est_low & truth <= est_upp),
      orc_covered = (truth >= orc_low & truth <= orc_upp),
      est_width   = est_upp - est_low,
      orc_width   = orc_upp - orc_low
    )

  return(rep_res)
}

all_results <- bind_rows(results_list)

summary_stats <- all_results %>%
  group_by(param) %>%
  summarize(
    True_Value    = first(truth),
    Avg_Est_Point = mean(est_point, na.rm = TRUE),
    Avg_Est_Low   = mean(est_low, na.rm = TRUE),
    Avg_Est_Upp   = mean(est_upp, na.rm = TRUE),
    Est_Coverage  = mean(est_covered, na.rm = TRUE),
    Avg_Est_Width = mean(est_width, na.rm = TRUE),
    
    Avg_Orc_Point = mean(orc_point, na.rm = TRUE),
    Avg_Orc_Low   = mean(orc_low, na.rm = TRUE),
    Avg_Orc_Upp   = mean(orc_upp, na.rm = TRUE),
    Orc_Coverage  = mean(orc_covered, na.rm = TRUE),
    Avg_Orc_Width = mean(orc_width, na.rm = TRUE),
    mc_reps       = n()
  )

print(summary_stats)

save(all_results, summary_stats, file = paste0("mc_comparison_", Sys.Date(), ".rda"))
