#Mismatch Experiment 3
#Charley Wu 2018

#house keeping
rm(list=ls())

outputfolder <- 'modelResults/mismatch'
reps <- 10000 #replications
lambdaGrid <- seq(0.1,10, 0.1)

#clusterid <- as.integer(runif(1,1,30))
clusterid <- as.integer(commandArgs(TRUE)[1]) 
targetLambda <- lambdaGrid[clusterid] #simulation to try

#load packages
packages <- c('plyr', 'jsonlite', 'ggplot2', 'gridExtra', 'reshape2', 'stargazer', 'coefplot', "grid", 'matrixcalc', 'parallel')
lapply(packages, require, character.only = TRUE)

#Participant data
source('dataMunging.R')
d <- dataImport() #participant data
#GP model specifications
source("Models.R")

#Environments
setwd("..") #Set into parent folder to read kernel files
setwd("experiment3")
#load environments from json, unlist, transform into numeric, convert into matrix, and name dimensions
environments <-fromJSON("environments/agridat.json", flatten=TRUE)
setwd("..")#Step back into analysis folder
setwd("analysis3")

#parameter estimates
gpPars <- read.csv('rationalModels/parameters/gp.csv')
localgpPars <- read.csv('rationalModels/parameters/localgp.csv')
#bmtPars <-  read.csv('rationalModels/parameters/BMT.csv')
#localbmtPars <-  read.csv('rationalModels/parameters/localBMT.csv')

#############################################################################################################################
# GP Model
#############################################################################################################################

gpRationalModel <- function(replications, outputfile, lambda, beta, tau, acq=ucb, kernel=rbf){
  #choice matrix
  choices <- expand.grid(0:10, 0:10) #build choice matrix
  names(choices)<-c("x1","x2")
  #run for smooth environments
  rewards <- mapply(1:replications, FUN=function(x){
    #randomly choose environment
    envNum <- sample(1:20,1) 
    #1st trial is random
    location <- sample(1:121,1)#each location as an integer from 1:121; first value is random, so treat both environments as replications of the same actions
    #Observations
    X1 <- choices[location,'x1']
    X2 <- choices[location,'x2']
    #reward
    reward<- c()
    reward[1] <- Y <- environments[[envNum]][location,"y"]*100
    for (j in 2:41){ #after that, loop through remaining trials and make decisions based on GP preditions
      #compute posterior predictions
      post <- gpr(X.test = choices, theta = c(lambda, lambda, 1, 0.0001), X = cbind(X1,X2), Y = ((Y-50)/100), k = kernel) #scale observed Y to zero mean, variance of 1
      #compute acquisition function evaluation
      utilityVec <- acq(post, pars = c(beta))
      #scale to max of prevent overflow by subtracting max
      utilityVec <- utilityVec - max(utilityVec)
      #compute softmax choice probabilities
      p <- exp(utilityVec/tau)
      p <- p/sum(p)
      #Sample next choice
      location <- sample(1:121,1, prob=p, replace=TRUE)
      #update reward, X1, X2, and Y for both smooth and rough
      reward[j] <- environments[[envNum]][location,"y"] * 100
      X1 <- c(X1, choices[location, 'x1'])
      X2 <- c(X2, choices[location, 'x2'])
      Y <- c(Y,  reward[j])}
    reward})
  #put into dataFrame
  gpDF <- data.frame(trial=seq(1:41),
                     lambda = lambda,
                     meanReward = rowMeans(rewards), #reward averaged over replications
                     meanSE = apply(rewards, 1, FUN = function(x) sd(x)/sqrt(length(x))),
                     maxReward = rowMeans(apply(rewards, 2, FUN=function(x) maxton(x))), #apply maxton function over columns (trials) to compute the max reward at each time t, and then average the max values over repliations
                     maxSE =  apply(apply(rewards, 2, FUN=function(x) maxton(x)), 1, FUN=function(x)  sd(x)/sqrt(length(x)))) 
  #add mmodel label
  gpDF$Model <- rep("GP-UCB", 41)
  #write to csv
  if(!is.null(outputfile)){
    write.csv(gpDF, outputfile)
  }
  return(gpDF)
}


filename <- paste0(outputfolder, '/', clusterid, '.csv')

#start <- proc.time()
gpDF <- gpRationalModel(reps, outputfile = filename , lambda = targetLambda, beta = median(gpPars$beta), tau= median(gpPars$tau))
#proc.time() - start

#############################################################################################################################
# Local-GP Model
#############################################################################################################################

localGPmodel<- function(replications, outputfile, lambda, tau, beta, acq=ucb, kernel=rbf){
  #build manhattan blocks and choice matrix
  manhattan <- expand.grid(0:10, 0:10) 
  names(manhattan)<-c("x1","x2")
  choices <- manhattan
  #choice matrix
  choices <- expand.grid(0:10, 0:10) #build choice matrix
  names(choices)<-c("x1","x2")
  #run for smooth environments
  rewards <- mapply(1:replications, FUN=function(x){
    #randomly choose environment
    envNum <- sample(1:20,1) 
    #1st trial is random
    location <- sample(1:121,1)#each location as an integer from 1:121; first value is random, so treat both environments as replications of the same actions
    #Observations
    X1 <- choices[location,'x1']
    X2 <- choices[location,'x2']
    #reward
    reward<- c()
    reward[1] <- Y <- environments[[envNum]][location,"y"]*100
    for (j in 2:41){ #after that, loop through remaining trials and make decisions based on GP preditions
      #compute Manhattan Distances
      prev <- unlist(manhattan[location,])
      distance <- sapply(seq(1:121), FUN=function(x) abs(manhattan[x, "x1"] - prev["x1"]) + abs(manhattan[x, "x2"] - prev["x2"]))
      #set values of 0 to 1
      distance[distance==0]<-1
      #compute posterior predictions
      post <- gpr(X.test = choices, theta = c(lambda, lambda, 1, 0.0001), X = cbind(X1,X2), Y = ((Y-50)/100), k = kernel) #scale observed Y to zero mean, variance of 1
      #compute acquisition function evaluation and weight by inverse manhattan distance
      utilityVec <- acq(post, pars = c(beta)) /  distance
      #scale to max of prevent overflow by subtracting max
      utilityVec <- utilityVec - max(utilityVec)
      #compute softmax choice probabilities
      p <- exp(utilityVec/tau)
      p <- p/sum(p)
      #Sample next choice
      location <- sample(1:121,1, prob=p, replace=TRUE)
      #update reward, X1, X2, and Y for both smooth and rough
      reward[j] <- environments[[envNum]][location,"y"] * 100
      X1 <- c(X1, choices[location, 'x1'])
      X2 <- c(X2, choices[location, 'x2'])
      Y <- c(Y,  reward[j])}
    reward})
  #put into dataFrame
  lgpDF <- data.frame(trial=seq(1:41),
                      lambda = lambda,
                      meanReward = rowMeans(rewards), #reward averaged over replications
                      meanSE = apply(rewards, 1, FUN = function(x) sd(x)/sqrt(length(x))),
                      maxReward = rowMeans(apply(rewards, 2, FUN=function(x) maxton(x))), #apply maxton function over columns (trials) to compute the max reward at each time t, and then average the max values over repliations
                      maxSE =  apply(apply(rewards, 2, FUN=function(x) maxton(x)), 1, FUN=function(x)  sd(x)/sqrt(length(x)))) 
  #add mmodel label
  #label model
  lgpDF$Model <- rep("LocalGP-UCB", 41)
  #write to csv
  if(!is.null(outputfile)){
    write.csv(lgpDF, outputfile)
  }
  return(lgpDF)
}


filename <- paste0(outputfolder, '/local', clusterid, '.csv')

#start <- proc.time()
#lgpDF <- localGPmodel(reps, outputfile = filename , lambda = targetLambda, beta = median(localgpPars$beta), tau= median(localgpPars$tau))
#proc.time() - start


