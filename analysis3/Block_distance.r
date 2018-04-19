#script to Inverse-Manhattan-Block distance results (i.e., inertia model)
#Eric Schulz and Charley Wu July 2018

#house keeping
rm(list=ls())

#load packages
packages <- c('plyr', 'ggplot2', 'jsonlite', 'MASS', 'gridExtra')
lapply(packages, require, character.only = TRUE)

#read in data
dat<-read.csv("ExperimentData/exp3.csv", sep=';')


##############################################################################################################
# IMPORT DATA
##############################################################################################################
#dummy data frame
d<-data.frame(id=numeric(), trial=numeric(), x=numeric(), y=numeric(), 
              z=numeric(), zmax=numeric(), scenario=numeric(), 
              horizon=numeric(), round=numeric())

#loop through data
for (i in 1:nrow(dat)){
  #read JSON of search history
  dfinal <-fromJSON(as.character(dat$searchHistory[i]))
  #sampled x value
  x<-unlist(dfinal$xcollect)
  #sampled y value
  y<-unlist(dfinal$ycollect)
  #sampled z value
  z<-unlist(dfinal$zcollect)
  #length
  len<-sapply(dfinal[1]$xcollect, length)
  #trial 1:20 or 1:40
  trial<-c(1:len[1],1:len[2],1:len[3],1:len[4])
  #round, 1:4
  round<-rep(1:8, len)
  #horizon, 20 or 40
  horizon<-rep(len-1, len)
  #maximum value or cumulative
  scenario<-rep(dat$scenario[i], sum(len))
  #smoothness
  kernel<-rep(dat$kernel[i], sum(len))
  #id number
  id<-rep(i, sum(len))
  #zmax<-unlist(sapply(dfinal$zcollect, maxton))
  #dummy frame
  dummy<-data.frame(id, trial, x, y, z, scenario, horizon, round)
  #bind them together
  d<-rbind(d, dummy)
}

#create a unique listing for each space in the grid
allopts<-expand.grid(0:10, 0:10)
d$chosen<--99
for (i in 1:nrow(d)){ #loop through data and assign the unique listing from allopts to data$chosen based on the chosen x and y values
  d$chosen[i]<-which(d$x[i]==allopts$Var1 & d$y[i]==allopts$Var2)
}

#create an overall id, unique for each subject for each round
d$overall<-1:nrow(d)-with(d, ave(id, paste(id, round), FUN = seq_along))
d<-transform(d,overall=as.numeric(factor(overall)))
d$round <- factor(d$round) #convert round into factor

##############################################################################################################
# COMPUTE MANHATTAN BLOCK DISTANCE
##############################################################################################################
computeIMD <- function(){
  cand<-expand.grid(0:10, 0:10) #candidates
  names(cand)<-c("x","y")
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
      write.table(dist, file = paste0("blockdistance/block", k, ".csv"), sep = ",", col.names = FALSE, row.names = FALSE)
      k<-k+1
    }
  }
}

#computeIMD()

##############################################################################################################
# ESTIMATE MODEL FIT
##############################################################################################################
source('Models.R')

inertiaFit<-function(participant, tau){
  subjD <- subset(d, id==participant)
  #Vector to store negative log likelihods
  nLL <- rep(0,length(8)) 
  for (r in 1:8){ #Begin looping through each round
    #subset of data for round r
    roundD <- subset(subjD, round==r)
    #is this round a short or long horizon?
    horizon = roundD$horizon[1]
    #Observations of subject choice behavior
    chosen <- roundD$chosen
    chosen <- chosen[2:length(chosen)] # trim first observation, since it wasn't a choice but a randomly revealed tile
    #Utilties of each choice
    utilities <- NULL
    #loop through observations
    for (i in 1:horizon){ #skip the last observation, because no choice was made based on that information
      #Inverse Manhatten distance
      out <- inertia(subjD$id[1], i, as.integer(r))
      utilityVec <- ucb(out, c(0)) #use UCB as a proxy, where beta=0
      utilities <- rbind(utilities, utilityVec) # build horizon_length x 121 matrix, where each row holds the utilities of each choice at each decision time in the search horizon
    }
    #Softmax rule
    utilities <- utilities/max(utilities) #scale to max value of one
    p <- exp(utilities/tau)
    p <- p/rowSums(p)
    #avoid underflow by setting a floor and a ceiling
    p <- (pmax(p, 0.00001))
    p <- (pmin(p, 0.99999))
    #Calculate Negative log likelihood
    nLL[r] <- -sum(log(p[cbind(c(1:horizon),chosen)]))
  }#end loop through rounds
  return(sum(nLL))  #Return negative log likelihoods of all observations 
}

inertiaNLLs <- NULL
d$inertia <- 0
for (i in 1:80){
  inertiaNLLs[i] <- inertiaFit(i, 1)
  d$inertia[d$id==i]  <- rep(inertiaNLLs[i], length(d$inertia[d$id==i]))
}
