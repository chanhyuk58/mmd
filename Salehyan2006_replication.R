library('data.table')
library('texreg')

sg <- foreign::read.dta('../../gleditschsalehyan/R Replication/gleditsch_salehyan_correctedtime.dta')
readr::write_excel_csv(sg, file='../sg.csv')
sg <- as.data.table(sg)
head(sg)
names(sg)

# outcome: nonset / bigconset civil war onset at the year. based on the definition of 'war'
# explanatory: peace1 peace time duration. This is the left-censored var.

## Replication
t4m2=glm(nonset~logref2+nbcwbin+polityb+polityb2+lngdp+lnpop+het+peace1+s1a+s2a+s3a ,family=binomial(link="logit"), data=sg)
screenreg(t4m2) # Successful replication of the original result.

## Get Min Max

# two categories to deal with: 
# 1) At least 1 civil war recorded during the period.
# 2) No civil war recorded on the data set.
#    - Maximum: From the formation of the country to that year.
#    - Minimum: From the starting year of the dataset.

## Merge with country formation year dataset
# cow2iso <- fread('../data/cow2iso.csv')
# cow2iso <- cow2iso[!is.na(cow_id), ]
# cow2iso <- cow2iso[, .SD[which.min(valid_from)], cow_id][, .(cow_id, valid_from)]
# # cow2iso[, order(.N), cow_id]
# names(cow2iso) <- c('ccode', 'form_year')

# sg <- merge(sg, cow2iso, by=c('ccode'), all.x=T)

## Subset data with `nonset` value
sg1 <- sg[!is.na(nonset), ]

# Index to separate two cases.
indx_no_onset <- sg1[, sum(nonset) == 0, ccode][V1 == T, 1][[1]]

### At least 1 civil war
sg1_ex_onset <- sg1[!(ccode %in% indx_no_onset), ]
# filter rows before the first recorded civil war.
sg1_ex_onset <- sg1_ex_onset[sg1_ex_onset[, seq(.N) < which.max(nonset == 1), ccode][, V1], ]

### No civil war
sg1_no_onset <- sg1[(ccode %in% indx_no_onset), ]
sg1_no_onset <- sg1_no_onset[, `:=`(
  peace1_min = peace1,
  peace1_max = 
) ]



sg[, peace1_max := peace1, ccode]
sg[sg[, which(nonset == 1), ccode]$V1, peace1_max := NA, ccode]
sg[ccode== 20,head(.SD, 10), , .SDcols=c('peace1_max', 'peace1', 'nonset')]
