# Loading Packages
library(gimme);library(graphicalVAR);library(dplyr)
# Loading Data-Generating Function
source('~/work/StateSpaceFunctions.R')
# Setting Directory
thedir = '~/scratch/'
# When running off cluster, l should = x; 
# otherwise, loop over l as many times as you want MC reps
l = x; 
# Sample Size is 50/MC rep
# 10,000 timepoints were simulated at a delta-t = 0.1
# Then, we loop over "theDeltat" to subsample the true data at each sampling interval
SampleSize = 50
maxT <- 10000
Deltat = .1
TimeSeq=seq(0, (maxT-1)*Deltat,by=Deltat)
for (theDeltat in c(0.1,.5,1.0,10)){
  # Reading in data files
  dat1.list = list();dat1.listb = list();dat1 = NULL
  # Subsampling the data by theDeltat
  if(theDeltat == 0.1){
    theIndex = which(TimeSeq %in% tail(TimeSeq, 100))
  }else{theIndex = tail(which(TimeSeq %% theDeltat<.001), 100)}
  # Change these to change how the thresholding goes for s-gimme and scgVAR
  gimme.crit = 0.75; gvar.crit = 0.00
  amt = c((l*25-24):(l*25),(l*25-24):(l*25)+2500)
  # Generating File Storage
  dir.create(paste0(thedir,'/Results/'))
  dir.create(paste0(thedir,'/Results/Last100/'))
  dir.create(paste0(thedir,'/Results/Last100/Large/'))
  dir.create(paste0(thedir,'/Results/Last100/Large/',theDeltat))
  dir.create(paste0(thedir,'/Results/Last100/Large/',theDeltat,'/job ',l,'/'))
  dir.create(paste0(thedir,'/Results/Last100/Large/',theDeltat,'/job ',l,'/',gimme.crit,'N',gvar.crit,'/'))
  filedir = paste0(thedir,'/Results/Last100/Large/',theDeltat,'/job ',l,'/',gimme.crit,'N',gvar.crit,'/')
  ff <- dir(filedir, recursive=TRUE, full.names=TRUE) #Get all file names in that directory
  file.remove(ff)
  # Preparing data for scgVAR and sgimme which use data.frame and lists, respectively
  for(m in 1:SampleSize){
    dat1.list[[m]] = read.table(paste0(thedir,'/JUL22/OUModels/SUBJ',amt[m],'OU_cohen Large.txt'), 
                                sep = ',')[theIndex,]
    dat1 = rbind(dat1, dat1.list[[m]])
    dat1.listb[[m]] = dat1.list[[m]][,3:6] 
  }
  # Running s-gimme and scgVAR
  gimme(dat1.listb[1:SampleSize], out = filedir, 
        subgroup = TRUE, standardize = TRUE, groupcutoff = gimme.crit)
  gvarclust(data = dat1, idvar = 'ID', vars = paste0('y',1:4), gamma = gvar.crit, 
            dir = filedir)
}