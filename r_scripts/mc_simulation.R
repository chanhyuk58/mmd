library(doFuture)
library(foreach)
library(dplyr)
library(Rcpp)
library(splines)
library(nloptr)

cat("Packages loaded ...\n")
  
# Set Up
mc_reps <- 100                 

# --- 1. Master Compilation ---
cache_path <- "./cpp_cache"
if(dir.exists(cache_path)) unlink(cache_path, recursive = TRUE)
dir.create(cache_path)

cat(">> Master: Compiling C++ backend...\n")
Rcpp::sourceCpp("./mmd_cpp.cpp", cacheDir = cache_path)
cat(">> Master: Compilation complete.\n")

# --- 2. Setup Cluster ---
lsb_hosts <- Sys.getenv("LSB_HOSTS")
n_cores <- if (lsb_hosts != "") length(strsplit(lsb_hosts, " ")[[1]]) else (parallel::detectCores() - 1)

registerDoFuture()
cat(sprintf(">> Registered %d local sockets.\n", n_cores))

# --- 3. Worker Initialization ---
cat(">> Initializing workers one-by-one to prevent filesystem congestion...\n")

worker_ids <- 1:n_cores

for (id in worker_ids) {
  f <- future({
    library(Rcpp)
    library(splines)
    library(nloptr)
    library(dplyr)
    
    Rcpp::sourceCpp("./mmd_cpp.cpp", cacheDir = "./cpp_cache")
    
    source("./mmd_cpp.R")
    source("./generate_pop.R")
    
    return(TRUE)
  })
  value(f)
  
  if (id %% 10 == 0) cat(sprintf("   [%d/%d] workers ready...\n", id, n_cores))
  Sys.sleep(0.1) # Breathe
}

cat(">> All workers initialized successfully. Starting Monte Carlo...\n")

# Monte Carlo Simulation
results_list <- foreach(i = 1:mc_reps, .errorhandling = "stop") %dopar% {
  
  # Generate Data
  sim <- gen_pop(
    J = 100, 
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

  cat("Population Generated ...\n")
  
  # Projection Method
  fit_proj <- MMD_bounds(
    onset ~ log_gdp + democ + eth_het + log_pop + log_ref, 
    data = pop, v0_col = "v0", v1_col = "v1",
    method = "projection",
    B = 0,
    verbose = FALSE
  )

  cat("Projection method is done ...\n")
  
  # Profile Method
  fit_prof <- MMD_bounds(
    onset ~ log_gdp + democ + eth_het + log_pop + log_ref, 
    data = pop, v0_col = "v0", v1_col = "v1",
    method = "profile",
    grid_radius = 0.5,
    grid_points = 100, 
    B = 0,
    verbose = FALSE
  )

  cat("Profile method is done ...\n")
  
  # Results
  rep_res <- data.frame(
    iteration = i,
    param     = names(true_params),
    truth     = as.numeric(true_params),
    # Projection Results
    proj_est  = fit_proj$point_est,
    proj_low  = fit_proj$bounds_ID$Lower,
    proj_upp  = fit_proj$bounds_ID$Upper,
    # Profile Results
    prof_low  = fit_prof$bounds_ID$Lower,
    prof_upp  = fit_prof$bounds_ID$Upper
  )
  
  # Coverage
  rep_res <- rep_res %>%
    mutate(
      proj_covered = (truth >= proj_low & truth <= proj_upp),
      prof_covered = (truth >= prof_low & truth <= prof_upp),
      proj_width   = proj_upp - proj_low
    )
  
  return(rep_res)
}

stopCluster(cl)

# MC results
all_results <- bind_rows(results_list)

summary_stats <- all_results %>%
  group_by(param) %>%
  summarize(
    True_Value    = first(truth),
    Avg_Proj_Est  = mean(proj_est),
    Avg_Proj_Low  = mean(proj_low),
    Avg_Proj_Upp  = mean(proj_upp),
    Proj_Coverage   = mean(proj_covered),
    Avg_Prof_Est  = mean(prof_est),
    Avg_Prof_Low  = mean(prof_low),
    Avg_Prof_Upp  = mean(prof_upp),
    Prof_Coverage   = mean(prof_covered),
    n = n()
  )

print(summary_stats)

save(all_results, summary_stats, file = paste0("mc_final_", Sys.Date(), ".rda"))
