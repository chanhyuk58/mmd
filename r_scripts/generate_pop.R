library(dplyr)
library(tidyr)

gen_pop <- function(
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
  gdp_shock = 0.05
) {

  # 1. Meta-Data & Fixed Effects
  country_meta <- data.frame(country = 1:J) %>%
    mutate(
      birth_year = sample(birth_range[1]:birth_range[2], J, replace = TRUE),
      raw_obs_start = sample(obs_start_range[1]:obs_start_range[2], J, replace = TRUE),
      obs_start_year = pmax(birth_year, raw_obs_start),
      eth_het = runif(J, -1, 1)
    )

  # 2. Initialize Time-Varying Vectors
  N_total <- J * T_full
  ids <- rep(1:J, each = T_full)
  years <- rep(1:T_full, times = J)

  log_gdp <- numeric(N_total)
  log_pop <- numeric(N_total)
  log_ref <- numeric(N_total)
  democ <- integer(N_total)
  onset <- integer(N_total)
  true_duration <- numeric(N_total)

  idx_mat <- matrix(1:N_total, nrow = T_full, ncol = J, byrow = FALSE)

  # 3. Initialization (t=1)
  t1 <- idx_mat[1, ]
  log_gdp[t1] <- rnorm(J, mean_gdp, 1.0)
  log_pop[t1] <- rnorm(J, mean_pop, 1.5)
  log_ref[t1] <- rgamma(J, shape = 2, scale = 0.5)
  democ[t1] <- rbinom(J, 1, 0.5)
  true_duration[t1] <- 0

  # 4. Time Loop
  birth_years_vec <- country_meta$birth_year
  eth_het_vec <- country_meta$eth_het

  for(t in 2:T_full) {
    curr <- idx_mat[t, ]
    prev <- idx_mat[t-1, ]

    y_lag <- onset[prev]
    gdp_lag <- log_gdp[prev]
    pop_lag <- log_pop[prev]
    ref_lag <- log_ref[prev]
    dem_lag <- democ[prev]
    dur_lag <- true_duration[prev]

    # Updates (Same dynamics as before)
    shock <- ifelse(y_lag == 1, -gdp_shock, 0)
    log_gdp[curr] <- 0.95 * gdp_lag + 0.05 * mean_gdp + rnorm(J, 0, 0.1) + shock

    growth <- ifelse(y_lag == 1, 0, 0.015)
    log_pop[curr] <- pop_lag + growth + rnorm(J, 0, 0.005)

    log_ref[curr] <- 0.7 * ref_lag + rgamma(J, shape = 1, scale = 0.5)

    p_flip <- ifelse(y_lag == 1, 0.20, 0.01)
    do_flip <- runif(J) < p_flip
    democ[curr] <- ifelse(do_flip, 1 - dem_lag, dem_lag)

    is_birth <- (birth_years_vec == t)
    reset <- (y_lag == 1) | is_birth
    true_duration[curr] <- ifelse(reset, 0, dur_lag + 1)

    # --- LINEAR PROBABILITY MODEL GENERATION ---
    active <- (t >= birth_years_vec)
    if(any(active)) {
      act <- curr[active]

      lin_pred <- beta_0 +
                  beta_gdp * log_gdp[act] +
                  beta_democ * democ[act] +
                  beta_eth * eth_het_vec[active] +
                  beta_pop * log_pop[act] +
                  beta_ref * log_ref[act] +
                  beta_v * true_duration[act]

      # Clamp probabilities to [0, 1]
      probs <- pmin(pmax(lin_pred, 0), 1)

      onset[act] <- rbinom(length(act), 1, probs)
    }
  }

  # 5. Construct Data Frame
  full_data <- data.frame(
    country = ids,
    year = years,
    onset = onset,
    log_gdp = log_gdp,
    democ = democ,
    log_pop = log_pop,
    log_ref = log_ref,
    true_duration = true_duration
  ) %>%
    left_join(country_meta, by = "country")

  # 6. Apply Observation Window & Bounds
  obs_data <- full_data %>%
    filter(year >= birth_year) %>%
    filter(year >= obs_start_year) %>%
    arrange(country, year) %>%
    group_by(country) %>%
    mutate(
      spell_id = cumsum(lag(onset, default = 0)),
      v0 = sequence(rle(spell_id)$lengths) - 1,
      unobserved_history = obs_start_year - birth_year,
      v1 = ifelse(spell_id > 0, v0, v0 + unobserved_history)
    ) %>%
    ungroup() %>%
    select(country, year, onset, log_gdp, democ, eth_het, log_pop, log_ref,
           v0, v1, true_duration, birth_year, obs_start_year)

  true_params <- c(
    "(Intercept)" = beta_0,
    "latent_v" = beta_v,
    "log_gdp" = beta_gdp,
    "democ" = beta_democ,
    "eth_het" = beta_eth,
    "log_pop" = beta_pop,
    "log_ref" = beta_ref
  )

  return(list(data = obs_data, true_params = true_params, meta = country_meta))
}

gen_pop_simple <- function(
                           n = 1000,
                           p = 5
                           ){
    gamma = c(1, 1, rep(-1, p))
    v <- runif(n, -2, 3) # try uniform distribution

    ## Observed
    v1 = ceiling(v)
    v0 <- v1 - 1
    # x <- runif(N, 0, 5) # try uniform distribution
    # multivariate uniform distribution
    x <- matrix(runif(n * p), nrow = n, ncol = p)
    ## Error
    eps <- rnorm(n, 0, 1)

    ## y
    y <- gamma[1] + gamma[2]*v + x %*% gamma[3:length(gamma)] + eps

    df <- data.frame(y, x, v0, v1)
    return(df)
}

# Sample test
if (FALSE) {
  pop <- gen_pop(
    J = 100, T_full = 100,
    birth_range = c(1, 50),
    obs_start_range = c(50, 90),
    beta_0 = 0.10,
  )

  print(pop$data)
  print(sum(pop$data$onset))
  print(mean(pop$data$onset))
}
