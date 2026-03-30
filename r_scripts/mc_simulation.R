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
    beta_0 = 0.2, 
    gamma_v = -0.05, 
    beta_gdp = 0.00, 
    beta_democ = 0.00, 
    beta_eth = 0.00,
    beta_pop = 0.05,
    beta_ref = 0.00,
    # beta_democ = -0.05, 
    # beta_eth = 0.02,
    # beta_pop = 0.05,
    # beta_ref = 0.05,
  )$data
  # head(pop)
  print(sum(pop$onset))
  print(mean(pop$onset))
  print(summary(pop$v1 - pop$v0))
  txtdensity((pop$v1 - pop$v0))

  # Estimate bounds
  stp <- 0.01
  grid_list <- list(
    "(Intercept)" = seq(-0.1, 0.4, stp),
    "latent_v0"  = seq(-0.2, 0.1, stp),
    "log_gdp"   = seq(-0.0, 0.0, stp),
    "democ"     = seq(-0.0, 0.0, stp),
    "eth_het"   = seq(-0.0, 0.0, stp),
    "log_pop"   = seq(-0.1, 0.2, stp),
    "log_ref"   = seq(-0.0, 0.0, stp)
  )
  method <- "projection"
  bounds <- MMD_bounds(
    onset ~ log_gdp + democ + eth_het + log_pop + log_ref, 
    data = pop, v0_col = "v0", v1_col = "v1",
    method = method,
    # grid_list = grid_list,
    alpha = 0.05,            # Critical value
    B = 0,                 # Bootstrap reps for CI calculation
    # b_exponent = 0.8,        # Proportion of sample used for bootstrap CI
    verbose = TRUE
  )
  print(bounds)
  save(bounds, file = paste("./", method, "_", Sys.Date(), ".rda"))
}
