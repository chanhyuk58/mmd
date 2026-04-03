library(doParallel)
library(foreach)
library(Rcpp)
library(dplyr)
library(tidyr)

# Parse IBM LSF Environment
lsb_hosts <- Sys.getenv("LSB_HOSTS")

if (lsb_hosts != "") {
  node_list <- strsplit(lsb_hosts, " ")[[1]]
  cl <- makePSOCKcluster(node_list)
  cat(sprintf(">> HPC Mode: Registered %d slots on LSF.\n", length(node_list)))
} else {
  # Local Fallback
  n_cores <- parallel::detectCores() - 1
  cl <- makePSOCKcluster(node_list, rshcmd = "blaunch")
  cat(sprintf(">> Local Mode: Registered %d cores.\n", n_cores))
}

registerDoParallel(cl)

# Pre-Compile C++
cat(">> Compiling C++ backend on master node...\n")
cache_path <- paste0(getwd(), "/cpp_cache")
if(!dir.exists(cache_path)) dir.create(cache_path)
Rcpp::sourceCpp("./mmd_cpp.cpp", cacheDir = cache_path)

# Initialize Distributed Workers
cat(">> Initializing remote workers...\n")
clusterEvalQ(cl, {
  library(Rcpp)
  library(splines)
  library(nloptr)
  library(dplyr)
  
  # Load pre-compiled binary
  Rcpp::sourceCpp("./mmd_cpp.cpp", cacheDir = paste0(getwd(), "/cpp_cache"))
  
  # Source logic
  source("./mmd_cpp.R")
  source("./generate_pop.R")
})

# Parallel Monte Carlo Loop
mc_reps <- 100  # Total iterations
cat(sprintf(">> Starting MC loop (%d iterations)...\n", mc_reps))

results_list <- foreach(i = 1:mc_reps, .errorhandling = 'pass') %dopar% {
  
  # A. Generate Data
  sim <- gen_pop(
    J = 500, 
    T_full = 100, 
    birth_range = c(1, 50),
    obs_start_range = c(50, 90),
    beta_0 = 0.2, 
    beta_gdp = -0.10, 
    beta_democ = -0.05, 
    beta_eth = 0.02,
    beta_pop = 0.05,
    beta_ref = 0.05,
    beta_v = -0.02, 
    mean_gdp = 8.5,
    mean_pop = 16.0,
    gdp_shock = 0.05,
    seed = 72938
  )
  pop <- sim$data
  true_params <- sim$true_params
  
  # B. Estimate via Projection Method
  fit_proj <- MMD_bounds(
    onset ~ log_gdp + democ + eth_het + log_pop + log_ref, 
    data = pop, v0_col = "v0", v1_col = "v1",
    method = "projection", B = 0, verbose = FALSE
  )
  
  # C. Estimate via Profile Grid Method
  fit_prof <- MMD_bounds(
    onset ~ log_gdp + democ + eth_het + log_pop + log_ref, 
    data = pop, v0_col = "v0", v1_col = "v1",
    method = "profile", 
    grid_radius = 0.5, 
    grid_points = 100,
    B = 0, verbose = FALSE
  )
  
  # D. Combine Results for this iteration
  res_proj <- fit_proj$bounds_ID %>% 
    mutate(param = rownames(.), method = "projection")
  
  res_prof <- fit_prof$bounds_ID %>% 
    mutate(param = rownames(.), method = "profile")
  
  combined <- bind_rows(res_proj, res_prof) %>%
    mutate(
      rep = i,
      truth = true_params[param],
      covered = (truth >= Lower & truth <= Upper),
      width = Upper - Lower
    )
  
  return(combined)
}

stopCluster(cl)

# --- 5. Final Aggregation and Summary ---
cat(">> Aggregating results...\n")
all_reps_df <- bind_rows(results_list)

# Calculate Summary Statistics
summary_table <- all_reps_df %>%
  group_by(param, method) %>%
  summarize(
    True_Value   = first(truth),
    Mean_Lower   = mean(Lower, na.rm = TRUE),
    Mean_Upper   = mean(Upper, na.rm = TRUE),
    Mean_Width   = mean(width, na.rm = TRUE),
    ID_Coverage  = mean(covered, na.rm = TRUE),
    Fail_Rate    = mean(is.na(Lower)),
    .groups = "drop"
  ) %>%
  arrange(param, method)

# Print Summary to Console
print(as.data.frame(summary_table))

# Save everything
saveRDS(all_reps_df, file = "mc_full_data.rds")
saveRDS(summary_table, file = "mc_summary_table.rds")
write.csv(summary_table, "mmd_comparison_results.csv", row.names = FALSE)

cat(">> Simulation Complete. Files saved to disk.\n")
