library(dplyr)
library(readr)
library(foreign)
library(here)
library(txtplot)

# Load replication data
sg <- foreign::read.dta(
  here("gleditschsalehyan", "R Replication", "gleditsch_salehyan_correctedtime.dta")
)

# Load state list and find the minimum establishment year for each ccode
states2024 <- read_csv(here("data", "statelist2024.csv"), col_types = cols()) %>%
  filter(!is.na(ccode)) %>%
  group_by(ccode) %>%
  summarize(styear = min(styear, na.rm = TRUE), .groups = "drop")

# Merge datasets
sg <- sg %>%
  left_join(states2024, by = "ccode")

# Identify Censoring Windows
sg <- sg %>%
  group_by(ccode) %>%
  arrange(year, .by_group = TRUE) %>%
  mutate(
    first_war1_idx = which(nonset == 1)[1],
    first_war1_idx = ifelse(is.na(first_war1_idx), n() + 1, first_war1_idx),
    is_censored1   = row_number() <= first_war1_idx,
    
    first_war2_idx = which(bigconset2 == 1)[1],
    first_war2_idx = ifelse(is.na(first_war2_idx), n() + 1, first_war2_idx),
    is_censored2   = row_number() <= first_war2_idx
  ) %>%
  ungroup()

mean(sg$first_war1_idx > 1) # 0.9655% of the data points are not censored
mean(sg$first_war2_idx > 1) # 0.9725% of the data points are not censored
mean(sg$nonset, na.rm=TRUE)
sd(sg$lngdp, na.rm = TRUE)
sd(sg$lnpop, na.rm = TRUE)

# ------------------------------------------------------------------------------
# Most Conservative (Upper bound is state age)
# ------------------------------------------------------------------------------

sg_most_conserv <- sg %>%
  mutate(
    potential_max1 = ifelse(!is.na(styear) & (year - styear) > peace1, year - styear, peace1),
    potential_max2 = ifelse(!is.na(styear) & (year - styear) > peace2, year - styear, peace2),
    
    peace1_min = peace1,
    peace1_max = ifelse(is_censored1, potential_max1, peace1),
    peace1_width = peace1_max - peace1_min,
    
    peace2_min = peace2,
    peace2_max = ifelse(is_censored2, potential_max2, peace2),
    peace2_width = peace2_max - peace2_min,
  ) %>%
  select(-potential_max1, -potential_max2, -is_censored1, -is_censored2, -first_war1_idx, -first_war2_idx)

summary(sg_most_conserv$peace1_width)
txtdensity(na.omit(sg_most_conserv$peace1_width))

write_csv(sg_most_conserv, here("data", "sg1_most_conserv.csv"))

# ------------------------------------------------------------------------------
# Capped at 5 (Upper bound is min(state age, 5))
# ------------------------------------------------------------------------------

sg_at5 <- sg %>%
  mutate(
    potential_max1 = ifelse(!is.na(styear) & (year - styear) > peace1, year - styear, peace1),
    potential_max2 = ifelse(!is.na(styear) & (year - styear) > peace2, year - styear, peace2),
    
    peace1_min = peace1,
    peace1_max = ifelse(is_censored1, pmax(pmin(potential_max1, 5), peace1), peace1),
    peace1_width = peace1_max - peace1_min,
    
    peace2_min = peace2,
    peace2_max = ifelse(is_censored2, pmax(pmin(potential_max2, 5), peace2), peace2),
    peace2_width = peace2_max - peace2_min,
  ) %>%
  select(-potential_max1, -potential_max2, -is_censored1, -is_censored2, -first_war1_idx, -first_war2_idx)

summary(sg_at5$peace1_width)
txtdensity(na.omit(sg_most_conserv$peace1_width))

write_csv(sg_at5, here("data", "sg1_at5.csv"))

# ------------------------------------------------------------------------------
# Unit Max
# ------------------------------------------------------------------------------

sg_unit_max <- sg %>%
  group_by(ccode) %>%
  mutate(
    country_max1 = max(peace1, na.rm = TRUE),
    country_max1 = ifelse(is.infinite(country_max1), peace1, country_max1),
    
    country_max2 = max(peace2, na.rm = TRUE),
    country_max2 = ifelse(is.infinite(country_max2), peace2, country_max2),
  ) %>%
  ungroup() %>%
  mutate(
    peace1_min = peace1,
    peace1_max = ifelse(is_censored1, pmax(country_max1, peace1), peace1),
    peace1_width = peace1_max - peace1_min,
    
    peace2_min = peace2,
    peace2_max = ifelse(is_censored2, pmax(country_max2, peace2), peace2),
    peace2_width = peace2_max - peace2_min,
  ) %>%
  select(-country_max1, -country_max2, -is_censored1, -is_censored2, -first_war1_idx, -first_war2_idx)

summary(sg_unit_max$peace1_width)
txtdensity(na.omit(sg_most_conserv$peace1_width))

write_csv(sg_unit_max, here("data", "sg1_unit_max.csv"))
