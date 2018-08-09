#Generate Simulated Model Data using Participant parameter estimates 
#Charley Wu 2018

#############################################################################################################################
# IMPORT DATA AND ENVIRONMENTS
#############################################################################################################################
#house keeping
rm(list=ls())

#load packages
packages <- c('plyr', 'jsonlite', 'ggplot2', 'gridExtra', 'reshape2', "grid", 'matrixcalc')
lapply(packages, require, character.only = TRUE)
source("Models.R") #model specifications
source('dataMunging.R')

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

pars$horizon <- factor(pars$horizon, levels=c(20, 40), labels=c("Short", "Long"))

#Participant data
d <- dataImport(normalize=FALSE) #participant data
simD <- d #used for simulated data
#add columns for the parameters used to generate the data
simD$tau <- NA
simD$lambda <- NA
simD$kError <- NA
simD$beta <- NA

#Environments
setwd("..") #Set into parent folder to read kernel files
setwd("experiment3")
#load environments from json, unlist, transform into numeric, convert into matrix, and name dimensions
environments <- fromJSON("environments/agridat.json", flatten=TRUE)
setwd("..")#Step back into analysis folder
setwd("analysis3")


#build manhattan blocks and choice matrix
choices <- manhattan <- expand.grid('x1'=0:10, 'x2'=0:10) 

#############################################################################################################################
# Run simulation
#############################################################################################################################

for (i in 1:80){ #loop thorugh participants
  subjd <- subset(d, id == i) #subset participant
  subjParams <- subset(pars, participant==i) #subset parameters
  envs <- unique(subjd$env) #which specific environments did the subject encounter?
  #counters to relate to leaveoutindex in paramEstimates
  longRound <- 1
  shortRound <- 1
  for (round in 1:8){ #loop through rounds
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
    location <- sample(1:121,1) #first location is random
    #Observations
    X1 <- choices[location,'x1']
    X2<- choices[location,'x2']
    chosen <- c(location)
    #rewards
    reward<-c()
    reward[1] <- Y <- Ymax <- environments[[envNum+1]][location,"y"]*100 #add 1 to envNum to go from range 0-19 to 1-20
    prevPost <- NULL  #set the previous posterior computation to NULL for the kalman filter
    #beginloop through remaining trials and make decisions based on GP preditions
    for (j in 2:(horizonValue+1)){
      if (localize){
        #compute Manhattan Distances
        prev <- unlist(manhattan[location,])
        distance <- sapply(seq(1:121), FUN=function(x) abs(manhattan[x, "x1"] - prev["x1"]) + abs(manhattan[x, "x2"] - prev["x2"]))
        #set values of 0 to 1
        distance[distance==0]<-1
        if (inherits(k, 'GP')){#compute GP posterior predictions
          post <- gpr(X.test = choices, theta = c(lambda, lambda, 1, 0.0001), X = cbind(X1,X2), Y = ((Y-50)/100), k = rbf) #scale Y observations to zero mean and variance of 1
        }else if (inherits(k, 'KalmanFilter')){
          post <- bayesianMeanTracker(x = as.matrix(cbind(X1[j-1],X2[j-1])), y = (reward[j-1]-50)/100, prevPost = prevPost, theta = c(kError))
          #update prevPost for the next round
          prevPost <- post}
        #compute acquisition function evaluation and weight by inverse manhattan distance
        utilityVec <- ucb(post, pars = c(beta)) / distance
      }else{
        if (inherits(k, 'GP')){#compute GP posterior predictions
          post <- gpr(X.test = choices, theta = c(lambda, lambda, 1, 0.0001), X = cbind(X1,X2), Y = ((Y-50)/100), k = rbf) #scale Y observations to zero mean and variance of 1
        }else if (inherits(k, 'KalmanFilter')){
          post <- bayesianMeanTracker(x = as.matrix(cbind(X1[j-1],X2[j-1])), y = (reward[j-1]-50)/100, prevPost = prevPost, theta = c(kError))
          #update prevPost for the next round
          prevPost <- post}
        #compute acquisition function evaluation and weight by inverse manhattan distance
        utilityVec <- ucb(post, pars = c(beta))
      }
      #scale to max of prevent overflow by subtracting max
      utilityVec <- utilityVec - max(utilityVec)
      #compute softmax choice probabilities
      p <- exp(utilityVec/tau)
      p <- p/sum(p)
      #Sample next choice
      location <- sample(1:121,1, prob=p, replace=TRUE)
      #update reward, X1, X2, and Y
      reward[j] <- environments[[envNum+1]][location,"y"] * 100
      X1 <- c(X1, choices[location, 'x1'])
      X2 <- c(X2, choices[location, 'x2'])
      chosen <- c(chosen, location)
      Y <- c(Y,  reward[j])
      Ymax <- c(Ymax, max(Y))
    }
    #insert data intosimD
    simD[simD$id==i & simD$env==envNum,]$x <- X1
    simD[simD$id==i & simD$env==envNum,]$y <- X2
    simD[simD$id==i & simD$env==envNum,]$z <- Y
    simD[simD$id==i & simD$env==envNum,]$zmax <- Ymax
    simD[simD$id==i & simD$env==envNum,]$chosen <- chosen
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

write.csv(simD, paste0('rationalModels/simulatedData/',outputFileName,'.csv'))
