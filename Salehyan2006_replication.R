library('data.table')
# library('texreg')

sg <- foreign::read.dta('../../gleditschsalehyan/R Replication/gleditsch_salehyan_correctedtime.dta')
readr::write_excel_csv(sg, file='../data/sg.csv')
sg <- as.data.table(sg)
head(sg)
names(sg)

# outcome: nonset / bigconset civil war onset at the year. based on the definition of 'war'
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
states2016 <- fread('../data/states2016.csv')
# head(states2016)
states2016 <- states2016[!is.na(ccode), ]
states2016 <- states2016[, .SD[which.min(styear)], ccode][, .(ccode, styear)]
states2016[, order(.N), ccode]

sg <- merge(sg, states2016, by=c('ccode'), all.x=T)
# head(sg)
# names(sg)

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
fwrite(sg1, file='../data/sg1.csv', bom=T)

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
fwrite(sg1, file='../data/sg1.csv', bom=T)

# sg1[, which(nonset != bigconset2)]
# (sg1$peace1 != sg1$peace2)
# (sg1$peace1_max != sg1$peace2_max)

names(sg1)

## Run MMD
source("mmd_cpp.R")


gleditschsalehyan1 <- MMD_bounds_cpp(
  nonset ~ logref2 + nbcwbin + polityb + polityb2 + lngdp + het,
  v0 = "peace1_min",
  v1 = "peace1_max",
  data = na.omit(sg1[,c("nonset", "logref2", "nbcwbin", "polityb", "polityb2", "lngdp", "het", "peace1_min", "peace1_max")])
)

print(summary(gleditschsalehyan1))
