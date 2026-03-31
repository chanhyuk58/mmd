library("txtplot")
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
  )$data
  # head(pop)
  print(sum(pop$onset))
  print(mean(pop$onset))
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
    alpha = 0.05,            # Critical value
    B = 0,                 # Bootstrap reps for CI calculation
    # b_exponent = 0.8,        # Proportion of sample used for bootstrap CI
    verbose = TRUE
  )
  print(bounds_projection)
  print(bounds_grid)
  save(bounds_projection, file = paste("./projection", "_", Sys.Date(), ".rda"))
  save(bounds_grid, file = paste("./grid", "_", Sys.Date(), ".rda"))
}

proc.time()
gc()
