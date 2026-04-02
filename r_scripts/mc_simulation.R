library("txtplot")
library(parallel)
options(np.cores = detectCores() - 1)

source("./generate_pop.R")
source("./mmd_cpp.R")

# The bounds is still very large. 

mc <- 1

for (i in 1:mc) {
  # Create a population
    pop <- gen_pop(
                   J = 100, # Count of countries
                   T_full = 100, # Max true duration / history
                   birth_range = c(1, 50),
                   obs_start_range = c(50, 60),
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
                   )$data
  # head(pop)
  cat("Head of population:\n")
  print(head(pop))
  cat("Number of Total Civil War Onset:\n")
  print(sum(pop$onset))
  cat("Proportion of Total Civil War Onset:\n")
  print(mean(pop$onset))
  cat("Summary of the length of V and the Density plot:\n")
  print(summary(pop$v1 - pop$v0))
  txtdensity((pop$v1 - pop$v0))

  # Estimate bounds
  bounds_projection <- MMD_bounds(
    onset ~ log_gdp + democ + eth_het + log_pop + log_ref, 
    data = pop, v0_col = "v0", v1_col = "v1",
    method = "projection",
    # grid_list = grid_list,
    alpha = 0.05,            # Critical value
    B = 0,                 # Bootstrap reps for CI calculation
    # b_exponent = 0.8,        # Proportion of sample used for bootstrap CI
    verbose = TRUE
  )
  bounds_grid <- MMD_bounds(
    onset ~ log_gdp + democ + eth_het + log_pop + log_ref, 
    data = pop, v0_col = "v0", v1_col = "v1",
    method = "grid",
    grid_radius = 1.0,      # This is more than enough for LPM model, since beta =< 1
    grid_points = 10000,    # Eval Points for each Parameters
    alpha = 0.05,            # Critical value
    B = 0,                 # Bootstrap reps for CI calculation
    # b_exponent = 0.8,        # Proportion of sample used for bootstrap CI
    verbose = TRUE
  )
  cat("Summary of Bounds from Pojection Method:\n")
  print(bounds_projection)
  cat("Summary of Bounds from Grid-search Method:\n")
  print(bounds_grid)
  save(bounds_projection, file = paste("./projection", "_", Sys.Date(), ".rda"))
  save(bounds_grid, file = paste("./grid", "_", Sys.Date(), ".rda"))
}

proc.time()
gc()
