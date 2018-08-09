#Charley Wu 2018
#Script to run rational models in comparison to human behavior split by condition

#house keeping
rm(list=ls())

#load packages
packages <- c('plyr', 'jsonlite', 'ggplot2', 'reshape2', "grid", 'matrixcalc')
lapply(packages, require, character.only = TRUE)

reps <- 10000 #replications

#Cluster id from qsub
#DEBUG clusterid <- as.integer(runif(1, 1, 41))
clusterid <- as.integer(commandArgs(TRUE)[1]) #Cluster id, corresponds to an integer used to indicate which combination of kernel and acquisition function to simulate
set.seed(clusterid)

payoffConditions <- c("Cumulative", "Best")
environments <- c("Smooth", "Rough")
horizon <- c(20,40)
models <- c("local", "gp", "localgp", "bmt", "localbmt")
opts <- expand.grid(payoffConditions, environments, horizon, models)
colnames(opts) <- c("payoff", "env", "horizon", "model")

condition <- opts[clusterid,]
#############################################################################################################################
# IMPORT ENVIRONMENTS
#############################################################################################################################

#GP model specifications
source("Models.R")


maxton<-function(x){
  maxn<-rep(0,length(x))
  maxn[1]<-x[1]
  for (i in 2:length(x)){
    if (x[i]>maxn[i-1]){maxn[i]<-x[i]}
    else{maxn[i]<-maxn[i-1]}
  }
  return(maxn)
}

#Environments
setwd("..") #Set into parent folder to read kernel files
setwd("experiment2D")
#load environments from json, unlist, transform into numeric, convert into matrix, and name dimensions
roughEnvironments <- lapply(fromJSON("kernel1.json", flatten=TRUE), FUN=function(x) matrix(as.numeric(unlist(x)), ncol=3, byrow=TRUE, dimnames=list(seq(1,121), c('x2', 'y', 'x1'))))
smoothEnvironments <- lapply(fromJSON("kernel2.json", flatten=TRUE), FUN=function(x) matrix(as.numeric(unlist(x)), ncol=3, byrow=TRUE, dimnames=list(seq(1,121), c('x2', 'y', 'x1'))))
setwd("..")#Step back into analysis folder
setwd("analysis2D")

#parameter estimates
localPars <- read.csv('rationalModels/parameters/local.csv')
gpPars <- read.csv('rationalModels/parameters/gp.csv')
localgpPars <- read.csv('rationalModels/parameters/localgp.csv')
bmtPars <-  read.csv('rationalModels/parameters/BMT.csv')
localbmtPars <-  read.csv('rationalModels/parameters/localBMT.csv')
#############################################################################################################################
# Learning Curve Simulations
#############################################################################################################################

learningCurveSimulation <- function(replications, outputfile, condition){
  if (condition$env=="Smooth"){
    environments <- smoothEnvironments
  }else
    environments <- roughEnvironments
  horizon <- condition$horizon
  model <- condition$model
  #build manhattan blocks
  manhattan <- expand.grid(0:10, 0:10) 
  names(manhattan)<-c("x1","x2")
  #Local Search
  if (model == 'local'){
    #load parameters and subset by condition
    params = subset(localPars, reward==condition$payoff & environment==condition$env & horizon == condition$horizon)
    learningCurve <- sapply(1:replications, FUN=function(x){
      tau <- sample(params$tau,1) #sample a tau value for each replication
      rewards <-c() #store rewards
      envNum <- sample(1:20,1) #randomly choose environment
      #1st trial is random
      location <- sample(1:121,1)
      rewards[1] <- environments[[envNum]][location,"y"]*100 #pull out rewards
      for (j in 2:(condition$horizon+1)){ #after that, loop through remaining trials based on inertia-based transition probabilities
        #given previous location, calculate distances to each other location
        prevLocation <- unlist(manhattan[location,])
        distances <- sapply(seq(1:121), FUN= function(x) abs(manhattan[x,"x1"] - prevLocation["x1"]) + abs(manhattan[x,"x1"] - prevLocation["x1"]))
        #set values of 0 to 1, and then take the inverse
        distances[distances==0]<-1
        distances<-unlist(distances)^-1
        #Calculate softmax transition probabilities, sampling from the vector of given tau values
        p <- exp(distances/tau)
        p <- p/sum(p)
        #select a number between 1:121 based on the probabilities p
        location <- sample(1:121,1, prob=p, replace=TRUE) #update location
        rewards[j] <- environments[[envNum]][location,"y"] * 100}
      rewards}) #return rewards
  }
  #GP-UCB
  else if (model=='gp'){
    #load parameters and subset by condition
    params = subset(gpPars, reward==condition$payoff & environment==condition$env & horizon == condition$horizon)
    #choice matrix
    choices <- expand.grid(0:10, 0:10) #build choice matrix
    names(choices)<-c("x1","x2")
    learningCurve <- sapply(1:replications, FUN=function(x){
      paramSample <- params[sample(1:ncol(params), 1),]
      lambda <- paramSample$lambda
      beta <- paramSample$beta
      tau <- paramSample$tau
      reward <-c() #store rewards
      envNum <- sample(1:20,1) #randomly choose environment
      #1st trial is random
      location <- sample(1:121,1)
      #Observations
      X1 <- choices[location,'x1']
      X2 <- choices[location,'x2']
      reward[1] <- Y <- environments[[envNum]][location,"y"]*100 #pull out rewards
      for (j in 2:(condition$horizon+1)){ #after that, loop through remaining trials based on inertia-based transition probabilities
        #compute posterior predictions
        post <- gpr(X.test = choices, theta = c(lambda, lambda, 1, 0.0001), X =  cbind(X1,X2), Y = ((Y-50)/100), k = rbf) #scale observed Y to zero mean, variance of 1
        #compute acquisition function evaluation
        utilityVec <- ucb(post, pars = c(beta))
        #prevent overflow by subtracting max
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
  }
  #Local GP-UCB
  else if (model=='localgp'){
    #load parameters and subset by condition
    params = subset(localgpPars, reward==condition$payoff & environment==condition$env & horizon == condition$horizon)
    #choice matrix
    choices <- expand.grid(0:10, 0:10) #build choice matrix
    names(choices)<-c("x1","x2")
    learningCurve <- sapply(1:replications, FUN=function(x){
      paramSample <- params[sample(1:ncol(params), 1),]
      lambda <- paramSample$lambda
      beta <- paramSample$beta
      tau <- paramSample$tau
      reward <-c() #store rewards
      envNum <- sample(1:20,1) #randomly choose environment
      #1st trial is random
      location <- sample(1:121,1)
      #Observations
      X1 <- choices[location,'x1']
      X2 <- choices[location,'x2']
      reward[1] <- Y <- environments[[envNum]][location,"y"]*100 #pull out rewards
      for (j in 2:(condition$horizon+1)){ #after that, loop through remaining trials based on inertia-based transition probabilities
        #compute Manhattan Distances
        prev <- unlist(manhattan[location,])
        distance <- sapply(seq(1:121), FUN=function(x) abs(manhattan[x, "x1"] - prev["x1"]) + abs(manhattan[x, "x2"] - prev["x2"]))
        #set values of 0 to 1
        distance[distance==0]<-1
        #compute posterior predictions
        post <- gpr(X.test = choices, theta = c(lambda, lambda, 1, 0.0001), X =  cbind(X1,X2), Y = ((Y-50)/100), k = rbf) #scale observed Y to zero mean, variance of 1
        #compute acquisition function evaluation
        utilityVec <- ucb(post, pars = c(beta))/distance
        #prevent overflow by subtracting max
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
  }
  #BMT
  else if (model=='bmt'){
    #load parameters and subset by condition
    params = subset(bmtPars, reward==condition$payoff & environment==condition$env & horizon == condition$horizon)
    #choice matrix
    choices <- expand.grid(0:10, 0:10) #build choice matrix
    names(choices)<-c("x1","x2")
    learningCurve <- sapply(1:replications, FUN=function(x){
      paramSample <- params[sample(1:ncol(params), 1),]
      kError <- paramSample$kError
      beta <- paramSample$beta
      tau <- paramSample$tau
      reward <-c() #store rewards
      envNum <- sample(1:20,1) #randomly choose environment
      #1st trial is random
      location <- sample(1:121,1)
      #Observations
      X1 <- choices[location,'x1']
      X2 <- choices[location,'x2']
      reward[1] <- Y <- environments[[envNum]][location,"y"]*100 #pull out rewards
      #posterior
      prevPost <- NULL  #initialize previous posterior for sequential updating of BMT posterior
      for (j in 2:(condition$horizon+1)){ #after that, loop through remaining trials based on inertia-based transition probabilities
        #update posterior predictions
        post <- bayesianMeanTracker(x = cbind(X1[j-1],X2[j-1] ), y = (reward[j-1] - 50) / 100, prevPost = prevPost, theta = c(kError))  #IMPORTANT, rescale sampled location between 0:29 because bayesianMeanTracker automatically compensates for the input range by adding +1
        prevPost <- post  #save new posterior as prevPost for next round
        #compute acquisition function evaluation
        utilityVec <- ucb(post, pars = c(beta))
        #prevent overflow by subtracting max
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
  }
  #Local BMT
  else if (model=='localbmt'){
    #load parameters and subset by condition
    params = subset(localbmtPars, reward==condition$payoff & environment==condition$env & horizon == condition$horizon)
    #choice matrix
    choices <- expand.grid(0:10, 0:10) #build choice matrix
    names(choices)<-c("x1","x2")
    learningCurve <- sapply(1:replications, FUN=function(x){
      paramSample <- params[sample(1:ncol(params), 1),]
      kError <- paramSample$kError
      beta <- paramSample$beta
      tau <- paramSample$tau
      reward <-c() #store rewards
      envNum <- sample(1:20,1) #randomly choose environment
      #1st trial is random
      location <- sample(1:121,1)
      #Observations
      X1 <- choices[location,'x1']
      X2 <- choices[location,'x2']
      reward[1] <- Y <- environments[[envNum]][location,"y"]*100 #pull out rewards
      #posterior
      prevPost <- NULL  #initialize previous posterior for sequential updating of BMT posterior
      for (j in 2:(condition$horizon+1)){ #after that, loop through remaining trials based on inertia-based transition probabilities
        #compute Manhattan Distances
        prev <- unlist(manhattan[location,])
        distance <- sapply(seq(1:121), FUN=function(x) abs(manhattan[x, "x1"] - prev["x1"]) + abs(manhattan[x, "x2"] - prev["x2"]))
        #set values of 0 to 1
        distance[distance==0]<-1
        #update posterior predictions
        post <- bayesianMeanTracker(x = cbind(X1[j-1],X2[j-1]), y = (reward[j-1] - 50) / 100, prevPost = prevPost, theta = c(kError))  #IMPORTANT, rescale sampled location between 0:29 because bayesianMeanTracker automatically compensates for the input range by adding +1
        prevPost <- post  #save new posterior as prevPost for next round
        #compute acquisition function evaluation
        utilityVec <- ucb(post, pars = c(beta))/ distance
        #prevent overflow by subtracting max
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
  }
  #Save data
  learningCurveDF <-data.frame(trial=seq(1:(horizon + 1)),
                               Environment=rep(condition$env, horizon+1),
                               PayoffCondition = rep(condition$payoff, horizon+1),
                               Horizon = rep(horizon, horizon+1),
                               meanReward = rowMeans(learningCurve), #reward averaged over replications
                               meanSE = apply(learningCurve, 1, FUN = function(x) sd(x)/sqrt(length(x))),
                               maxReward =  rowMeans(apply(learningCurve, 2, FUN=function(x) maxton(x))), #apply maxton function over columns (trials) to compute the max reward at each time t, and then average the max values over repliations
                               maxSE =  apply(apply(learningCurve, 2, FUN=function(x) maxton(x)), 1, FUN=function(x)  sd(x)/sqrt(length(x))))
  
  #add model 
  learningCurveDF$Model <- rep(condition$model, nrow(learningCurveDF))
  #write to csv
  if(!is.null(outputfile)){
    write.csv(learningCurveDF, outputfile)
  }
  return(learningCurveDF)
}

#run simulations and save data
outputFileName = paste0('rationalModels/learningCurves/',clusterid,'.csv')
learningCurveDF <- learningCurveSimulation(reps, outputFileName, condition)


