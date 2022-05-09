# GIMME and scgVAR Analysis of fMRI data
library(gimme);library(graphicalVAR)
list.files('./fMRI Data/')
dat1 = NULL; ind=1
for(i in c(2:5,7,9)){
  temp = read.csv(paste0("./fMRI Data/sub_100",i,"_ses_1_gordon.csv"))
  temp$ID = rep(ind, nrow(temp))
  dat1 = rbind(dat1, temp)
  ind = ind+1
}
dat2 = NULL
for(i in c(1:3,5:7)){
  temp = read.csv(paste0("./fMRI Data/sub_200",i,"_ses_1_gordon.csv"))
  temp$ID = rep(ind, nrow(temp))
  dat2 = rbind(dat2, temp)
  ind = ind+1
}
dat = rbind(dat1, dat2)

## S-GIMME
gimme(data = './fMRI Data D1/', out = './fMRI GIMME/', sep = ',', header = TRUE, subgroup = TRUE)
## scgVAR
gvarclust(data = dat, idvar = 'ID', vars = names(dat)[1:14], 
          gamma = 0.0, dir = paste0('./fMRI GIMME/gamma 0.00/'))
gvarclust(data = dat, idvar = 'ID', vars = names(dat)[1:14], 
          gamma = 0.5, dir = paste0('./fMRI GIMME/gamma 0.50/'))

## Analyzing All Data
## S-GIMME
gimme(data = './fMRI Data/', out = './fMRI GIMME All/', sep = ',', header = TRUE, subgroup = TRUE)
## scgVAR
files = list.files('./fMRI Data/')
dat1 = NULL; ind=1
for(i in files){
  temp = read.csv(paste0("./fMRI Data/",i))
  temp$ID = rep(ind, nrow(temp))
  dat1 = rbind(dat1, temp)
  ind = ind+1
}
gvarclust(data = dat, idvar = 'ID', vars = names(dat)[1:14], 
          gamma = 0.0, dir = paste0('./fMRI GIMME All/gamma 0.00/'))
gvarclust(data = dat, idvar = 'ID', vars = names(dat)[1:14], 
          gamma = 0.5, dir = paste0('./fMRI GIMME All/gamma 0.50/'))


