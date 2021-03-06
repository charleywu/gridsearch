#Generate Simulated Data 1D
#Charley Wu 2017

#############################################################################################################################
# IMPORT DATA AND ENVIRONMENTS
#############################################################################################################################
#house keeping
rm(list=ls())

#load packages
packages <- c('plyr', 'jsonlite', 'ggplot2', 'gridExtra', 'reshape2', "grid", 'matrixcalc')
lapply(packages, require, character.only = TRUE)
source('Models.R')

#clusterid <- as.integer(runif(1, 1, 400))
clusterid <- as.integer(commandArgs(TRUE)[1]) #Cluster id, corresponds to an integer used to indicate which combination of kernel and acquisition function to simulate
set.seed(clusterid)

models <- c("BMT", "localBMT", "GP", "localGP")
rep <- seq(1,100)
opts <- expand.grid(models, rep)
colnames(opts) <- c("model", "rep")

condition <- opts[clusterid,]
acq <- ucb


if (condition$model=='BMT'){
  outputFileName <- paste0('simDataBMT', condition$rep)
  pars <- read.csv('rationalModels/parameters/BMT.csv')
  localize <- FALSE
  k <- bayesianMeanTracker
}else if (condition$model=="localBMT"){
  outputFileName <- paste0('simDataBMTLocal', condition$rep)
  pars <- read.csv('rationalModels/parameters/localBMT.csv')
  localize <- TRUE
  k <- bayesianMeanTracker
}else if (condition$model=='GP'){
  outputFileName <-  paste0('simDataGP', condition$rep)
  pars <- read.csv('rationalModels/parameters/gp.csv')
  localize <- FALSE
  k <- rbf
}else if (condition$model=="localGP"){
  outputFileName <- paste0('simDataGPLocal', condition$rep)
  pars <- read.csv('rationalModels/parameters/localgp.csv')
  localize <- TRUE
  k <- rbf
}


pars$horizon <- factor(pars$horizon, levels=c(5, 10), labels=c("Short", "Long"))

#Participant data
source('dataMunging.R')
d <- dataImport(normalize=FALSE) #participant data
simD <- d #used for simulated data
#add columns for the parameters we used
simD$tau <- NA
simD$kError <- NA
simD$lambda <- NA
simD$beta <- NA


#Environments
setwd("..") #Set into parent folder to read kernel files
setwd("experiment1D")
#load environments from json, unlist, transform into numeric, convert into matrix, and name dimensions
roughEnvironments <- lapply(fromJSON("kernel1.json", flatten=TRUE), FUN=function(x) matrix(as.numeric(unlist(x)), ncol=2, byrow=TRUE, dimnames=list(seq(1,30), c('y', 'x'))))
smoothEnvironments <- lapply(fromJSON("kernel2.json", flatten=TRUE), FUN=function(x) matrix(as.numeric(unlist(x)), ncol=2, byrow=TRUE, dimnames=list(seq(1,30), c('y', 'x'))))
setwd("..")#Step back into analysis folder
setwd("analysis1D")


for (i in 1:81){ #loop thorugh participants
  subjd <- subset(d, id == i) #subset participant
  subjParams <- subset(pars, participant==i) #subset parameters
  if (subjd$kernel[1]=="Rough"){#which class of environments? 
    envClass <-roughEnvironments
  }else{ envClass <-smoothEnvironments }
  envs <- unique(subjd$env) #which specific environments did the subject encounter?
  #counters to relate to leaveoutindex in paramEstimates
  longRound <- 1
  shortRound <- 1
  for (round in 1:16){ #loop through rounds
    envNum <- envs[round] #env number
    horizonLength <- subjd[subjd$env==envNum, ]$Horizon[1] #what was the horizon length
    horizonValue <- subjd[subjd$env==envNum, ]$horizon[1]
    ifelse(horizonLength == "Long", leaveoutindex <- longRound, leaveoutindex <- shortRound)
    ifelse(horizonLength == "Long", longRound <- longRound + 1, shortRound <- shortRound + 1) #increment leaveout index counters
    params <- subset(subjParams, horizon==horizonLength)[leaveoutindex,]
    if(inherits(k, "GP")){lambda <- params$lambda
    }else if(inherits(k, "KalmanFilter")){kError <- params$kError}
    beta <- params$beta
    tau <- params$tau
    location <- sample(1:30,1) #first location is random
    #Observations
    X <- c(location)
    #rewards
    reward<-c()
    reward[1] <- Y <-Ymax<- envClass[[envNum+1]][location,"y"]*100 #add 1 to envNum to go from range 0-19 to 1-20
    prevPost <- NULL #set the previous posterior computation to NULL for the kalman filter
    #beginloop through remaining trials and make decisions based on GP preditions
    for (j in 2:(horizonValue+1)){
      if (inherits(k, 'GP')){#compute GP posterior predictions
        post <- gpr(X.test = 1:30, theta = c(lambda, lambda, 1, 0.0001), X = X, Y = ((Y-50)/100), k = rbf) #scale Y observations to zero mean and variance of 1
      }else if (inherits(k, 'KalmanFilter')){
        post <- bayesianMeanTracker(x = X[j-1] -1, y = (reward[j-1] - 50) / 100, prevPost = prevPost, theta = c(kError))
        #update prevPost for the next round
        prevPost <- post}
      #compute acquisition function evaluation and weight by inverse manhattan distance
      utilityVec <- ucb(post, pars = c(beta))
      #scale to max of prevent overflow by subtracting max
      utilityVec <- utilityVec - max(utilityVec)
      #compute softmax choice probabilities
      p <- exp(utilityVec/tau)
      p <- p/sum(p)
      #Sample next choice
      location <- sample(1:30,1, prob=p, replace=TRUE)
      #update reward, X1, X2, and Y
      reward[j] <- envClass[[envNum+1]][location,"y"] * 100
      X <- c(X, location)
      Y <- c(Y,  reward[j])
      Ymax <- c(Ymax, max(Y))
    }
    #insert data intosimD
    simD[simD$id==i & simD$env==envNum,]$x <- X - 1 #IMPORTANT: rescale back to range of 0-29 from 1-30
    simD[simD$id==i & simD$env==envNum,]$y <- Y
    simD[simD$id==i & simD$env==envNum,]$ymax <- Ymax
    simD[simD$id==i & simD$env==envNum,]$tau <- tau
    if(inherits(k, "GP")){
      simD[simD$id==i & simD$env==envNum,]$lambda <- lambda
    }else if(inherits(k, "KalmanFilter")){
      simD[simD$id==i & simD$env==envNum,]$kError <- kError
    }
    simD[simD$id==i & simD$env==envNum,]$beta <- beta
  }
}
simD$Model <- condition$model

write.csv(simD, paste0('rationalModels/simulatedData/', outputFileName, '.csv'))
