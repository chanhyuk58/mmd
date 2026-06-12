library(doFuture)
library(foreach)
library(dplyr)
library(stats)

cat(">> Packages loaded successfully...\n")

# Source our pure R modular files
source("./generate_pop.R")
source("./mmd_dro_diag.R")
source("./mmd_dro.R")

cat(">> Functions loaded successfully...\n")

# --- 1. Setup Parallel Cluster ---
lsb_hosts <- Sys.getenv("LSB_HOSTS")
n_cores <- if (lsb_hosts != "") {
  length(strsplit(lsb_hosts, " ")[[1]])
} else {
  max(1, parallel::detectCores() - 1)
}

registerDoFuture()
plan(multisession, workers = n_cores)
cat(sprintf(">> Registered %d local socket workers (No SSH or compiler needed).\n", n_cores))

# --- 2. Monte Carlo Parameters ---
mc_reps <- 100
set.seed(4875995)

cat(">> Starting Monte Carlo Simulation...\n")

# --- 3. Parallel Execution ---
results_list <- foreach(i = 1:mc_reps,
                        .options.future = list(seed = TRUE),
                        .errorhandling = "stop") %dofuture% {
  
  source("./generate_pop.R")
  source("./mmd_dro.R")

  # Generate Population
  pop <- gen_pop_simple(n = 1000)
  colnames(pop) <- c("y", "x", "v0", "v1")
  
  # True Parameters from Manski and Tamer
  # Intercept = 1, beta_x = -1, gamma_v = 1
  true_params <- c("(Intercept)" = 1, "x" = -1, "latent_v0" = 1)

  # Estimation A: Direct LP Projection (Exact & Global)
  fit_proj <- mmd_dro(
    y ~ x, data = pop, v0_col = "v0", v1_col = "v1",
    family = "gaussian", method = "projection",
    delta = Inf, B = 0, verbose = FALSE
  )

  # Estimation B: LP-Profile Grid (Centered on Global Min)
  fit_prof <- mmd_dro(
    y ~ x, data = pop, v0_col = "v0", v1_col = "v1",
    family = "gaussian", method = "profile",
    delta = Inf, B = 0,
    grid_radius = 0.5, grid_points = 100,
    verbose = FALSE
  )

  # Compile Results
  rep_res <- data.frame(
    iteration = i,
    param     = names(true_params),
    truth     = as.numeric(true_params),
    # Projection Results
    proj_low  = fit_proj$results[["Inf"]]$bounds_ID$Lower,
    proj_upp  = fit_proj$results[["Inf"]]$bounds_ID$Upper,
    # Profile Results
    prof_low  = fit_prof$results[["Inf"]]$bounds_ID$Lower,
    prof_upp  = fit_prof$results[["Inf"]]$bounds_ID$Upper
  )

  # Calculate Coverage
  rep_res <- rep_res %>%
    mutate(
      proj_covered = (truth >= proj_low & truth <= proj_upp),
      prof_covered = (truth >= prof_low & truth <= prof_upp),
      proj_width   = proj_upp - proj_low
    )

  return(rep_res)
}

# --- 4. Aggregate Results ---
all_results <- bind_rows(results_list)

summary_stats <- all_results %>%
  group_by(param) %>%
  summarize(
    True_Value    = first(truth),
    Avg_Proj_Low  = mean(proj_low, na.rm = TRUE),
    Avg_Proj_Upp  = mean(proj_upp, na.rm = TRUE),
    Proj_Coverage = mean(proj_covered, na.rm = TRUE),
    Avg_Prof_Low  = mean(prof_low, na.rm = TRUE),
    Avg_Prof_Upp  = mean(prof_upp, na.rm = TRUE),
    Prof_Coverage = mean(prof_covered, na.rm = TRUE),
    mc_reps       = n()
  )

print(summary_stats)

# Save
save(all_results, summary_stats, file = paste0("mc_final_simple_", Sys.Date(), ".rda"))
