#script to compute Inverse-Manhattan-Block distance results for simulated data=
#Eric Schulz and Charley Wu, June 2017

#house keeping
rm(list=ls())

#load packages
packages <- c('plyr', 'ggplot2', 'jsonlite', 'MASS', 'gridExtra')
lapply(packages, require, character.only = TRUE)

#read in data
dat<-read.csv("ExperimentData/simDataBMTLocal.csv")
#dat<-read.csv("ExperimentData/simDataGPLocal.csv")


##############################################################################################################
# COMPUTE MANHATTAN BLOCK DISTANCE
##############################################################################################################
computeIMD <- function(destinationFolder, d){
  cand<-expand.grid('x'=0:10, 'y'=0:10) #candidates)
  k<-1 #counter
  #loop through subjects
  for (h in 1:max(d$id)){
    dsub<-subset(d, id==h)
    #loop through rounds
    for (i in 1:8){
      dblock<-subset(dsub, round==i)
      dist<-matrix(0, nrow=121, ncol=nrow(dblock)-1)
      #loop through trials
      for (j in 1:(nrow(dblock)-1)){
        dist[,j]<-abs(cand$x-dblock$x[j])+abs(cand$y-dblock$y[j])
      }
      dist[dist==0]<-1
      dist<-dist^-1
      write.table(dist, file = paste0(destinationFolder,"/block", k, ".csv"), sep = ",", col.names = FALSE, row.names = FALSE)
      k<-k+1
    }
  }
}

computeIMD('bmtLocalBlockDistance', dat)
#computeIMD('gpLocalBlockDistance', dat)