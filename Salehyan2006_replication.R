library('data.table')
library('texreg')

sg <- foreign::read.dta('../../gleditschsalehyan/R Replication/gleditsch_salehyan_correctedtime.dta')
readr::write_excel_csv(sg, file='../sg.csv')
sg <- as.data.table(sg)
head(sg)
names(sg)
# outcome: nonset / bigconset based on the definition of 'war'
# explanatory: peace1

## Replication
t4m2=glm(nonset~logref2+nbcwbin+polityb+polityb2+lngdp+lnpop+het+peace1+s1a+s2a+s3a ,family=binomial(link="logit"), data=sg)
screenreg(t4m2)

## Take the rows before the first onset
sg[, .SD[which(nonset == 1)[1]], ccode, .SDcols=c('stateid', 'nonset', 'peace1', 'year')]
sg[, peace1_max := peace1, ccode]
sg[sg[, which(nonset == 1), ccode]$V1, peace1_max := NA, ccode]
sg[ccode== 20,head(.SD, 10), , .SDcols=c('peace1_max', 'peace1', 'nonset')]
