library("data.table")
library("xtable")
# library("texreg")
source("mmd_cpp.R")

# Identified the coding error for the peace1 year interval.
# peace1_max > peace1_min does not hold for every sample.
# Probably due to the year of state's establishment is defined 
# differently across the data. Thus, naive calculation based on 
# one data causes errors.

sg <- foreign::read.dta("../../gleditschsalehyan/R Replication/gleditsch_salehyan_correctedtime.dta")
readr::write_excel_csv(sg, file="../data/sg.csv")
sg <- as.data.table(sg)
head(sg)
names(sg)

# outcome: nonset / bigconset civil war onset at the year. based on the definition of "war"
# explanatory: peace1 peace time duration. This is the left-censored var.

## Replication
t4m2=glm(nonset~logref2+nbcwbin+polityb+polityb2+lngdp+lnpop+het+peace1+s1a+s2a+s3a ,family=binomial(link="logit"), data=sg)
# screenreg(t4m2) # Successful replication of the original result.
summary(t4m2) # Successful replication of the original result.

## Get Min Max

# two categories to deal with: 
# 1) At least 1 civil war recorded during the period.
# 2) No civil war recorded on the data set.
#    - Maximum: From the formation of the country to that year.
#    - Minimum: From the starting year of the dataset.

# Merge with country formation year dataset
states2016 <- fread("../data/states2016.csv")
# head(states2016)
states2016 <- states2016[!is.na(ccode), ]
states2016 <- states2016[, .SD[which.min(styear)], ccode][, .(ccode, styear)]
# states2016[, order(.N), ccode]

sg <- merge(sg, states2016, by=c("ccode"), all.x=T)
# head(sg)
# names(sg)

# The Most Conservative {{{
## Based on `nonset`
# two categories
sg1_1 <- sg[!is.na(nonset), .SD[(seq(.N) <= which.max(nonset == 1) | sum(nonset) == 0), ], ccode
  ][, `:=`(
  peace1_min = peace1,
  peace1_max = ifelse(year - styear > 0, year - styear, peace1)
)]
# other data points
sg1_2 <- sg[!is.na(nonset), .SD[seq(.N) > which.max(nonset == 1) & sum(nonset) > 0, ], ccode
  ][, `:=`(
  peace1_min = peace1,
  peace1_max = peace1
)]
# merge
sg1 <- rbindlist(list(sg1_1, sg1_2))
sg1 <- sg1[order(ccode, year), ]

## Based on `bigconset2`
# two categories
sg1_1 <- sg1[!is.na(bigconset2), .SD[(seq(.N) <= which.max(bigconset2 == 1) | sum(bigconset2) == 0), ], ccode
  ][, `:=`(
  peace2_min = peace2,
  peace2_max = ifelse(year - styear > 0, year - styear, peace2)
)]
# other data points
sg1_2 <- sg1[!is.na(bigconset2), .SD[seq(.N) > which.max(bigconset2 == 1) & sum(bigconset2) > 0, ], ccode
  ][, `:=`(
  peace2_min = peace2,
  peace2_max = peace2
)]
# merge
sg1 <- rbindlist(list(sg1_1, sg1_2))
sg1 <- sg1[order(ccode, year), ]
fwrite(sg1, file="../data/sg1_most_conserv.csv", bom=T)
# }}}

# # Histogram of Max Peace Years {{{
# pdf(width = 5, height = 5, file = "./hist_peace_years_most_conserv.pdf")
# hist(sg1$peace1_max)
# dev.off()
# # }}}

# At 5 {{{
## Based on `nonset`
# two categories
sg1_1 <- sg[!is.na(nonset), .SD[(seq(.N) <= which.max(nonset == 1) | sum(nonset) == 0), ], ccode
  ][, `:=`(
  peace1_min = peace1,
  peace1_max = min(ifelse(year - styear > 0, year - styear, peace1), 5)
)]
# other data points
sg1_2 <- sg[!is.na(nonset), .SD[seq(.N) > which.max(nonset == 1) & sum(nonset) > 0, ], ccode
  ][, `:=`(
  peace1_min = peace1,
  peace1_max = peace1
)]
# merge
sg1 <- rbindlist(list(sg1_1, sg1_2))
sg1 <- sg1[order(ccode, year), ]

## Based on `bigconset2`
# two categories
sg1_1 <- sg1[!is.na(bigconset2), .SD[(seq(.N) <= which.max(bigconset2 == 1) | sum(bigconset2) == 0), ], ccode
  ][, `:=`(
  peace2_min = peace2,
  peace2_max = ifelse(year - styear > 0, year - styear, peace2)
)]
# other data points
sg1_2 <- sg1[!is.na(bigconset2), .SD[seq(.N) > which.max(bigconset2 == 1) & sum(bigconset2) > 0, ], ccode
  ][, `:=`(
  peace2_min = peace2,
  peace2_max = peace2
)]
# merge
sg1 <- rbindlist(list(sg1_1, sg1_2))
sg1 <- sg1[order(ccode, year), ]
fwrite(sg1, file="../data/sg1_at5.csv", bom=T)
# }}}

# At Country Level Max {{{
## Based on `nonset`
# two categories
sg1_2 <- sg[!is.na(nonset), .SD[seq(.N) > which.max(nonset == 1) & sum(nonset) > 0, ], ccode
  ][, `:=`(
  peace1_min = peace1,
  peace1_max = peace1
)]
sg1_1 <- sg[!is.na(nonset), .SD[(seq(.N) <= which.max(nonset == 1) | sum(nonset) == 0), ], ccode
  ][, `:=`(
  peace1_min = peace1,
  peace1_max = max(sg1_2$peace1_max) 
)]
# other data points
# merge
sg1 <- rbindlist(list(sg1_1, sg1_2))
sg1 <- sg1[order(ccode, year), ]

## Based on `bigconset2`
# two categories
sg1_1 <- sg1[!is.na(bigconset2), .SD[(seq(.N) <= which.max(bigconset2 == 1) | sum(bigconset2) == 0), ], ccode
  ][, `:=`(
  peace2_min = peace2,
  peace2_max = ifelse(year - styear > 0, year - styear, peace2)
)]
# other data points
sg1_2 <- sg1[!is.na(bigconset2), .SD[seq(.N) > which.max(bigconset2 == 1) & sum(bigconset2) > 0, ], ccode
  ][, `:=`(
  peace2_min = peace2,
  peace2_max = peace2
)]
# merge
sg1 <- rbindlist(list(sg1_1, sg1_2))
sg1 <- sg1[order(ccode, year), ]
fwrite(sg1, file="../data/sg1_unit_max.csv", bom=T)
# }}}

## Run MMD

Xs <- c("nonset", "logref2", "nbcwbin", "polityb", "polityb2", "lngdp", "lnpop", "het", "peace1_min", "peace1_max")
colnames <- c("Intercept", "Peace Years", "Log Refugees", "Neighbor Civil War", "Polity", "Polity Squared", "Log GDP", "Ethnic Heterogeneity")

# # The Most Conservative Bounds
sg1 <- fread("../data/sg1_most_conserv.csv")
print(mean(sg1$nonset))
length(unique(sg1$stateid))
length(unique(sg1$year))

summary(sg1$peace1_max - sg1$peace1_min)
sg1[sg1$peace1_max < sg1$peace1_min, c("stateid", "year", "peace1", "peace1_max", "peace1_min")]

# gleditschsalehyan1 <- MMD_bounds_cpp(
#   nonset ~ logref2 + nbcwbin + polityb + polityb2 + lngdp + het,
#   v0 = "peace1_min",
#   v1 = "peace1_max",
#   data = na.omit(sg1[ ,..Xs])
# )
#
# print(summary(gleditschsalehyan1))
#
# tab <- cbind(colnames, summary(gleditschsalehyan1)$stats)
# print(xtable(tab), include.rownames = FALSE, file = "../data/Salehyan_mmd_the_most_conservative.csv")

# # Ceiling at 30?
# sg1 <- fread("../data/sg1.csv")
# sg2 <- copy(sg1)
# sg2 <- sg2[, peace1_max := ifelse(peace1_max >= 30, 30, peace1_max)]
#
# gleditschsalehyan2 <- MMD_bounds_cpp(
#   nonset ~ logref2 + nbcwbin + polityb + polityb2 + lngdp + het,
#   v0 = "peace1_min",
#   v1 = "peace1_max",
#   data = na.omit(sg2[ ,..Xs])
# )
#
# print(summary(gleditschsalehyan2))
#
# tab <- cbind(colnames, summary(gleditschsalehyan2)$stats)
# print(xtable(tab), include.rownames = FALSE, 
#   file = "../data/Salehyan_mmd_ceiling_at30.csv")

# Ceiling at 5
sg1 <- fread("../data/sg1_at5.csv")
sg2 <- copy(sg1)
sg2 <- sg2[, peace1_max := ifelse(peace1_max >= 5, 5, peace1_max)]

gleditschsalehyan2 <- MMD_bounds_cpp(
  nonset ~ logref2 + nbcwbin + polityb + polityb2 + lngdp + het,
  v0 = "peace1_min",
  v1 = "peace1_max",
  data = na.omit(sg2[ ,..Xs])
)

print(summary(gleditschsalehyan2))

tab <- cbind(colnames, summary(gleditschsalehyan2)$stats)
print(xtable(tab), include.rownames = FALSE, 
  file = "../data/Salehyan_mmd_ceiling_at5.csv")

# Unit Max
sg1 <- fread("../data/sg1_unit_max.csv")
sg2 <- copy(sg1)
sg2 <- sg2[, peace1_max := ifelse(peace1_max >= 5, 5, peace1_max)]

gleditschsalehyan2 <- MMD_bounds_cpp(
  nonset ~ logref2 + nbcwbin + polityb + polityb2 + lngdp + het,
  v0 = "peace1_min",
  v1 = "peace1_max",
  data = na.omit(sg2[ ,..Xs])
)

print(summary(gleditschsalehyan2))

tab <- cbind(colnames, summary(gleditschsalehyan2)$stats)
print(xtable(tab), include.rownames = FALSE, 
  file = "../data/Salehyan_mmd_ceiling_unit_max.csv")
