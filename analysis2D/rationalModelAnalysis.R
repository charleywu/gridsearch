#Charley Wu 2018
#2D
#Script to run rational models in comparison to human behavior

#############################################################################################################################
# IMPORT DATA AND ENVIRONMENTS
#############################################################################################################################
#house keeping
rm(list=ls())
reps <- 10000 #replications

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
# RANDOM MODEL
# to read the previous simulation from disk:
# randomDF <- read.csv("rationalModels/random.csv")
#############################################################################################################################

randomModel <- function(replications, outputfile){
  #Run for both types of environments
  smoothReward <- mcmapply(1:replications, FUN=function(x) sample(smoothEnvironments[[sample(1:20,1)]][,'y'], 41, replace = TRUE)*100, mc.preschedule=TRUE, mc.cores=7)
  roughReward <- mcmapply(1:replications, FUN=function(x) sample(roughEnvironments[[sample(1:20,1)]][,'y'], 41, replace = TRUE)*100, mc.preschedule=TRUE, mc.cores=7)
  #put into dataFrame
  randomDF <- data.frame(trial=rep(seq(1:41), 2),
                         Environment=c(rep("Smooth", 41), rep("Rough", 41)),
                         meanReward = c(rowMeans(smoothReward), rowMeans(roughReward)), #reward averaged over replications
                         meanSE = c(apply(smoothReward, 1, FUN = function(x) sd(x)/sqrt(length(x))), apply(roughReward, 1, FUN = function(x) sd(x)/sqrt(length(x)))),
                         maxReward =  c(rowMeans(apply(smoothReward, 2, FUN=function(x) maxton(x))),  rowMeans(apply(roughReward, 2, FUN=function(x) maxton(x)))), #apply maxton function over columns (trials) to compute the max reward at each time t, and then average the max values over repliations
                         maxSE =  c(apply(apply(smoothReward, 2, FUN=function(x) maxton(x)), 1, FUN=function(x)  sd(x)/sqrt(length(x))), apply(apply(roughReward, 2, FUN=function(x) maxton(x)), 1, FUN=function(x)  sd(x)/sqrt(length(x)))))
  #add model label
  randomDF$Model <- rep("Random", 82)
  #write output
  if(!is.null(outputfile)){
    write.csv(randomDF, outputfile) 
  }
  return(randomDF)
}

randomDF <- randomModel(reps, "rationalModels/random.csv")

#############################################################################################################################
# Local Search Model
# to read the previous simulation from disk:
# localDF <- read.csv("rationalModels/localSearch.csv")
#############################################################################################################################

localSearch <- function(replications, outputfile, taulist, cores=7){
  #Separate parameters for each enviroment
  smoothTau <- subset(taulist, environment=="Smooth")$tau
  roughTau <-  subset(taulist, environment=="Rough")$tau
  #build manhattan blocks
  manhattan <- expand.grid(0:10, 0:10) 
  names(manhattan)<-c("x1","x2")
  #Run for smooth environments
  smoothReward <- mcmapply(1:replications, FUN=function(x){
    tau <- sample(smoothTau,1) #sample a tau value for each replication
    rewards <-c() #store rewards
    envNum <- sample(1:20,1) #randomly choose environment
    #1st trial is random
    location <- sample(1:121,1)
    rewards[1] <- smoothEnvironments[[envNum]][location,"y"]*100 #pull out rewards
    for (j in 2:41){ #after that, loop through remaining trials based on inertia-based transition probabilities
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
      rewards[j] <- smoothEnvironments[[envNum]][location,"y"] * 100}
    rewards}, mc.preschedule = TRUE, mc.cores=cores) #return rewards
  #run for rough environments
  roughReward <- mcmapply(1:replications, FUN=function(x){
    tau <- sample(roughTau,1) #sample a tau value for each replication
    rewards <-c() #store rewards
    envNum <- sample(1:20,1) #randomly choose environment
    #1st trial is random
    location <- sample(1:121,1)
    rewards[1] <- roughEnvironments[[envNum]][location,"y"]*100 #pull out rewards
    for (j in 2:41){ #after that, loop through remaining trials based on inertia-based transition probabilities
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
      rewards[j] <- roughEnvironments[[envNum]][location,"y"] * 100}
    rewards},  mc.preschedule = TRUE, mc.cores=cores) #return rewards
  #put into dataFrame
  localSearchDF <-data.frame(trial=rep(seq(1:41), 2),
                             Environment=c(rep("Smooth", 41), rep("Rough", 41)),
                             meanReward = c(rowMeans(smoothReward), rowMeans(roughReward)), #reward averaged over replications
                             meanSE = c(apply(smoothReward, 1, FUN = function(x) sd(x)/sqrt(length(x))), apply(roughReward, 1, FUN = function(x) sd(x)/sqrt(length(x)))),
                             maxReward =  c(rowMeans(apply(smoothReward, 2, FUN=function(x) maxton(x))),  rowMeans(apply(roughReward, 2, FUN=function(x) maxton(x)))), #apply maxton function over columns (trials) to compute the max reward at each time t, and then average the max values over repliations
                             maxSE =  c(apply(apply(smoothReward, 2, FUN=function(x) maxton(x)), 1, FUN=function(x)  sd(x)/sqrt(length(x))), apply(apply(roughReward, 2, FUN=function(x) maxton(x)), 1, FUN=function(x)  sd(x)/sqrt(length(x)))))
  
  #add model and environment labels
  localSearchDF$Model <- rep("Local Search", 82)
  localSearchDF$Environment <- factor(localSearchDF$Environment, levels=c("Smooth", "Rough"))
  #write to csv
  if(!is.null(outputfile)){
    write.csv(localSearchDF, outputfile)
  }
  return(localSearchDF)
}

localDF <- localSearch(reps, "rationalModels/localSearch.csv", localPars)

#############################################################################################################################
# GP Model
# to read the previous simulation from disk:
# gpDF <- read.csv("rationalModels/GPUCB.csv")
#############################################################################################################################

gpRationalModel <- function(replications, outputfile, parameters, acq=ucb, kernel=rbf, cores=7){
  smoothPars <- subset(parameters, environment=="Smooth")
  roughPars <- subset(parameters, environment =="Rough")
  #choice matrix
  choices <- expand.grid(0:10, 0:10) #build choice matrix
  names(choices)<-c("x1","x2")
  #run for smooth environments
  smoothReward <- mcmapply(1:replications, FUN=function(x){
    #sample parameters
    params <- parameters[sample(1:ncol(smoothPars), 1),]
    lambda <- params$lambda
    beta <- params$beta
    #tau <- params$tau
    tau <- median(smoothPars$tau)
    #randomly choose environment
    envNum <- sample(1:20,1) 
    #1st trial is random
    location <- sample(1:121,1)#each location as an integer from 1:121; first value is random, so treat both environments as replications of the same actions
    #Observations
    X1 <- choices[location,'x1']
    X2 <- choices[location,'x2']
    #reward
    reward<- c()
    reward[1] <- Y <- smoothEnvironments[[envNum]][location,"y"]*100
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
      reward[j] <- smoothEnvironments[[envNum]][location,"y"] * 100
      X1 <- c(X1, choices[location, 'x1'])
      X2 <- c(X2, choices[location, 'x2'])
      Y <- c(Y,  reward[j])}
  reward}, mc.preschedule = TRUE, mc.cores=cores)
  #run for smooth environments
  roughReward <- mcmapply(1:replications, FUN=function(x){
    #sample parameters
    params <- parameters[sample(1:ncol(roughPars), 1),]
    lambda <- params$lambda
    beta <- params$beta
    #tau <- params$tau
    tau <- median(roughPars$tau)
    #randomly choose environment
    envNum <- sample(1:20,1) 
    #1st trial is random
    location <- sample(1:121,1)#each location as an integer from 1:121; first value is random, so treat both environments as replications of the same actions
    #Observations
    X1 <- choices[location,'x1']
    X2 <- choices[location,'x2']
    #reward
    reward<- c()
    reward[1] <- Y <- roughEnvironments[[envNum]][location,"y"]*100
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
      reward[j] <- roughEnvironments[[envNum]][location,"y"] * 100
      X1 <- c(X1, choices[location, 'x1'])
      X2 <- c(X2, choices[location, 'x2'])
      Y <- c(Y,  reward[j])}
    reward}, mc.preschedule = TRUE, mc.cores=cores)
  #put into dataFrame
  gpDF <-data.frame(trial=rep(seq(1:41), 2),
                       Environment=c(rep("Smooth", 41), rep("Rough", 41)),
                       meanReward = c(rowMeans(smoothReward), rowMeans(roughReward)), #reward averaged over replications
                       meanSE = c(apply(smoothReward, 1, FUN = function(x) sd(x)/sqrt(length(x))), apply(roughReward, 1, FUN = function(x) sd(x)/sqrt(length(x)))),
                       maxReward =  c(rowMeans(apply(smoothReward, 2, FUN=function(x) maxton(x))),  rowMeans(apply(roughReward, 2, FUN=function(x) maxton(x)))), #apply maxton function over columns (trials) to compute the max reward at each time t, and then average the max values over repliations
                       maxSE =  c(apply(apply(smoothReward, 2, FUN=function(x) maxton(x)), 1, FUN=function(x)  sd(x)/sqrt(length(x))), apply(apply(roughReward, 2, FUN=function(x) maxton(x)), 1, FUN=function(x)  sd(x)/sqrt(length(x)))))
  #add mmodel label
  gpDF$Model <- rep("GP-UCB", 82)
  #write to csv
  if(!is.null(outputfile)){
    write.csv(gpDF, outputfile)
  }
  return(gpDF)
}

gpDF <- gpRationalModel(reps, outputfile = "rationalModels/GPUCBmedianTau.csv", parameters = gpPars)

#############################################################################################################################
# Local-GP Model
# to read the previous simulation from disk:
# lgpDF <- read.csv( "rationalModels/localGPUCB.csv")
#############################################################################################################################

localGPmodel<- function(replications, outputfile, parameters, acq=ucb, kernel=rbf, cores=7){
  smoothPars <- subset(parameters, environment=="Smooth")
  roughPars <- subset(parameters, environment =="Rough")
  #build manhattan blocks and choice matrix
  manhattan <- expand.grid(0:10, 0:10) 
  names(manhattan)<-c("x1","x2")
  choices <- manhattan
  #Run for smooth environments
  smoothReward <- mcmapply(1:replications, FUN=function(x){
    #sample parameters
    params <- parameters[sample(1:ncol(smoothPars), 1),]
    lambda <- params$lambda
    beta <- params$beta
    #tau <- params$tau
    tau <- median(smoothPars$tau)
    envNum <- sample(1:20,1) #randomly choose environment
    #1st trial is random
    location <- sample(1:121,1)
    #Observations
    X1 <- choices[location,'x1']
    X2<- choices[location,'x2']
    #rewards
    reward<-c()
    reward[1] <- Y <- smoothEnvironments[[envNum]][location,"y"]*100
    #beginloop through remaining trials and make decisions based on GP preditions
    for (j in 2:41){
      #compute Manhattan Distances
      prev <- unlist(manhattan[location,])
      distance <- sapply(seq(1:121), FUN=function(x) abs(manhattan[x, "x1"] - prev["x1"]) + abs(manhattan[x, "x2"] - prev["x2"]))
      #set values of 0 to 1
      distance[distance==0]<-1
      #compute GP posterior predictions
      post <- gpr(X.test = choices, theta = c(lambda, lambda, 1, 0.0001), X = cbind(X1,X2), Y = ((Y-50)/100), k = kernel) #scale Y observations to zero mean and variance of 1
      #compute acquisition function evaluation and weight by inverse manhattan distance
      utilityVec <- acq(post, pars = c(beta)) / distance
      #scale to max of prevent overflow by subtracting max
      utilityVec <- utilityVec - max(utilityVec)
      #compute softmax choice probabilities
      p <- exp(utilityVec/tau)
      p <- p/sum(p)
      #Sample next choice
      location <- sample(1:121,1, prob=p, replace=TRUE)
      #update reward, X1, X2, and Y
      reward[j] <- smoothEnvironments[[envNum]][location,"y"] * 100
      X1 <- c(X1, choices[location, 'x1'])
      X2 <- c(X2, choices[location, 'x2'])
      Y <- c(Y,  reward[j])}
    reward},mc.preschedule = TRUE, mc.cores=cores)
  #Run for rough environments
  roughReward <- mcmapply(1:replications, FUN=function(x){
    #sample parameters
    params <- parameters[sample(1:ncol(roughPars), 1),]
    lambda <- params$lambda
    beta <- params$beta
    #tau <- params$tau
    tau <- median(roughPars$tau)
    envNum <- sample(1:20,1) #randomly choose environment
    #1st trial is random
    location <- sample(1:121,1)
    #Observations
    X1 <- choices[location,'x1']
    X2<- choices[location,'x2']
    #rewards
    reward<-c()
    reward[1] <- Y <- roughEnvironments[[envNum]][location,"y"]*100
    #beginloop through remaining trials and make decisions based on GP preditions
    for (j in 2:41){
      #compute Manhattan Distances
      prev <- unlist(manhattan[location,])
      distance <- sapply(seq(1:121), FUN=function(x) abs(manhattan[x, "x1"] - prev["x1"]) + abs(manhattan[x, "x2"] - prev["x2"]))
      #set values of 0 to 1
      distance[distance==0]<-1
      #compute GP posterior predictions
      post <- gpr(X.test = choices, theta = c(lambda, lambda, 1, 0.0001), X = cbind(X1,X2), Y = ((Y-50)/100), k = kernel) #scale Y observations to zero mean and variance of 1
      #compute acquisition function evaluation and weight by inverse manhattan distance
      utilityVec <- acq(post, pars = c(beta)) / distance
      #scale to max of prevent overflow by subtracting max
      utilityVec <- utilityVec - max(utilityVec)
      #compute softmax choice probabilities
      p <- exp(utilityVec/tau)
      p <- p/sum(p)
      #Sample next choice
      location <- sample(1:121,1, prob=p, replace=TRUE)
      #update reward, X1, X2, and Y
      reward[j] <- roughEnvironments[[envNum]][location,"y"] * 100
      X1 <- c(X1, choices[location, 'x1'])
      X2 <- c(X2, choices[location, 'x2'])
      Y <- c(Y,  reward[j])}
    reward},mc.preschedule = TRUE, mc.cores=cores)
  #put into dataFrame
  lgpDF  <-data.frame(trial=rep(seq(1:41), 2),
                      Environment=c(rep("Smooth", 41), rep("Rough", 41)),
                      meanReward = c(rowMeans(smoothReward), rowMeans(roughReward)), #reward averaged over replications
                      meanSE = c(apply(smoothReward, 1, FUN = function(x) sd(x)/sqrt(length(x))), apply(roughReward, 1, FUN = function(x) sd(x)/sqrt(length(x)))),
                      maxReward =  c(rowMeans(apply(smoothReward, 2, FUN=function(x) maxton(x))),  rowMeans(apply(roughReward, 2, FUN=function(x) maxton(x)))), #apply maxton function over columns (trials) to compute the max reward at each time t, and then average the max values over repliations
                      maxSE =  c(apply(apply(smoothReward, 2, FUN=function(x) maxton(x)), 1, FUN=function(x)  sd(x)/sqrt(length(x))), apply(apply(roughReward, 2, FUN=function(x) maxton(x)), 1, FUN=function(x)  sd(x)/sqrt(length(x)))))
  #label model
  lgpDF$Model <- rep("LocalGP-UCB", 82)
  #write to csv
  if(!is.null(outputfile)){
    write.csv(lgpDF, outputfile)
  }
  return(lgpDF)
}

lgpDF <- localGPmodel(reps, outputfile = "rationalModels/localGPUCBmedianTau.csv", parameters = localgpPars)

#############################################################################################################################
# BMT-UCB Model
# to read the previous simulation from disk:
# bmtDF <- read.csv("rationalModels/BMTUCB.csv")
#############################################################################################################################

bmtRationalModel <- function(replications, outputfile, parameters, acq=ucb, cores=7){
  smoothPars <- subset(parameters, environment=="Smooth")
  roughPars <- subset(parameters, environment =="Rough")
  #choice matrix
  choices <- expand.grid(0:10, 0:10) #build choice matrix
  names(choices)<-c("x1","x2")
  #run for smooth environments
  smoothReward <- mcmapply(1:replications, FUN=function(x){
    #sample parameters
    params <- parameters[sample(1:ncol(smoothPars), 1),]
    kError <- params$kError
    beta <- params$beta
    #tau <- params$tau
    tau <- median(smoothPars$tau)
    #randomly choose environment
    envNum <- sample(1:20,1) 
    #1st trial is random
    location <- sample(1:121,1)#each location as an integer from 1:121; first value is random, so treat both environments as replications of the same actions
    #Observations
    X1 <- choices[location,'x1']
    X2 <- choices[location,'x2']
    #reward
    reward<- c()
    reward[1] <- Y <- smoothEnvironments[[envNum]][location,"y"]*100
    #posterior
    prevPost <- NULL  #initialize previous posterior for sequential updating of BMT posterior
    for (j in 2:41){ #after that, loop through remaining trials and make decisions based on BMT-UCB preditions
      #update posterior predictions
      post <- bayesianMeanTracker(x = as.matrix(cbind(X1,X2)), y = (reward[j-1]-50)/100, prevPost = prevPost, theta = c(kError))
      prevPost <- post  #save new posterior as prevPost for next round
      #compute acquisition function evaluation
      utilityVec <- acq(post, pars = c(beta))
      #to prevent overflow, subtract max of q(x) vector 
      utilityVec <- utilityVec - max(utilityVec)
      #compute softmax choice probabilities
      p <- exp(utilityVec/tau)
      p <- p/sum(p)
      #Sample next choice
      location <- sample(1:121,1, prob=p, replace=TRUE)
      #update reward, X1, X2, and Y for both smooth and rough
      reward[j] <- smoothEnvironments[[envNum]][location,"y"] * 100
      X1 <- c(X1, choices[location, 'x1'])
      X2 <- c(X2, choices[location, 'x2'])
      Y <- c(Y,  reward[j])}
    reward}, mc.preschedule = TRUE, mc.cores=cores)
  #run for smooth environments
  roughReward <- mcmapply(1:replications, FUN=function(x){
    #sample parameters
    params <- parameters[sample(1:ncol(roughPars), 1),]
    kError <- params$kError
    beta <- params$beta
    #tau <- params$tau
    tau <- median(roughPars$tau)
    #randomly choose environment
    envNum <- sample(1:20,1) 
    #1st trial is random
    location <- sample(1:121,1)#each location as an integer from 1:121; first value is random, so treat both environments as replications of the same actions
    #Observations
    X1 <- choices[location,'x1']
    X2 <- choices[location,'x2']
    #reward
    reward<- c()
    reward[1] <- Y <- roughEnvironments[[envNum]][location,"y"]*100
    #posterior
    prevPost <- NULL  #initialize previous posterior for sequential updating of BMT posterior
    for (j in 2:41){ #after that, loop through remaining trials and make decisions based on GP preditions
      #update posterior predictions
      post <- bayesianMeanTracker(x = as.matrix(cbind(X1,X2)), y = (reward[j-1]-50)/100, prevPost = prevPost, theta = c(kError))
      prevPost <- post  #save new posterior as prevPost for next round
      #compute acquisition function evaluation
      utilityVec <- acq(post, pars = c(beta))
      #prevent overflow by subtracting max
      utilityVec <- utilityVec - max(utilityVec)
      #compute softmax choice probabilities
      p <- exp(utilityVec/tau)
      p <- p/sum(p)
      #Sample next choice
      location <- sample(1:121,1, prob=p, replace=TRUE)
      #update reward, X1, X2, and Y for both smooth and rough
      reward[j] <- roughEnvironments[[envNum]][location,"y"] * 100
      X1 <- c(X1, choices[location, 'x1'])
      X2 <- c(X2, choices[location, 'x2'])
      Y <- c(Y,  reward[j])}
    reward}, mc.preschedule = TRUE, mc.cores=cores)
  #put into dataFrame
  bmtDF <-data.frame(trial=rep(seq(1:41), 2),
                    Environment=c(rep("Smooth", 41), rep("Rough", 41)),
                    meanReward = c(rowMeans(smoothReward), rowMeans(roughReward)), #reward averaged over replications
                    meanSE = c(apply(smoothReward, 1, FUN = function(x) sd(x)/sqrt(length(x))), apply(roughReward, 1, FUN = function(x) sd(x)/sqrt(length(x)))),
                    maxReward =  c(rowMeans(apply(smoothReward, 2, FUN=function(x) maxton(x))),  rowMeans(apply(roughReward, 2, FUN=function(x) maxton(x)))), #apply maxton function over columns (trials) to compute the max reward at each time t, and then average the max values over repliations
                    maxSE =  c(apply(apply(smoothReward, 2, FUN=function(x) maxton(x)), 1, FUN=function(x)  sd(x)/sqrt(length(x))), apply(apply(roughReward, 2, FUN=function(x) maxton(x)), 1, FUN=function(x)  sd(x)/sqrt(length(x)))))
  #add model label
  bmtDF$Model <- rep("BMT-UCB", 82)
  #write to csv
  if(!is.null(outputfile)){
    write.csv(bmtDF, outputfile)
  }
  return(bmtDF)
}

bmtDF <- bmtRationalModel(reps, outputfile = "rationalModels/BMTUCBmedianTau.csv", parameters = bmtPars)

#############################################################################################################################
# Local BMT-UCB Model
# to read the previous simulation from disk:
# localbmtDF <- read.csv("rationalModels/localBMTUCB.csv")
#############################################################################################################################

localbmtRationalModel <- function(replications, outputfile, parameters, acq=ucb, cores=7){
  smoothPars <- subset(parameters, environment=="Smooth")
  roughPars <- subset(parameters, environment =="Rough")
  #build manhattan blocks and choice matrix
  manhattan <- expand.grid(0:10, 0:10) 
  names(manhattan)<-c("x1","x2")
  choices <- manhattan
  #run for smooth environments
  smoothReward <- mcmapply(1:replications, FUN=function(x){
    #sample parameters
    params <- parameters[sample(1:ncol(smoothPars), 1),]
    kError <- params$kError
    beta <- params$beta
    #tau <- params$tau
    tau <- median(smoothPars$tau)
    #randomly choose environment
    envNum <- sample(1:20,1) 
    #1st trial is random
    location <- sample(1:121,1)#each location as an integer from 1:121; first value is random, so treat both environments as replications of the same actions
    #Observations
    X1 <- choices[location,'x1']
    X2 <- choices[location,'x2']
    #reward
    reward<- c()
    reward[1] <- Y <- smoothEnvironments[[envNum]][location,"y"]*100
    #posterior
    prevPost <- NULL  #initialize previous posterior for sequential updating of BMT posterior
    for (j in 2:41){ #after that, loop through remaining trials and make decisions based on BMT-UCB preditions
      #compute Manhattan Distances
      prev <- unlist(manhattan[location,])
      distance <- sapply(seq(1:121), FUN=function(x) abs(manhattan[x, "x1"] - prev["x1"]) + abs(manhattan[x, "x2"] - prev["x2"]))
      #set values of 0 to 1
      distance[distance==0]<-1
      #update posterior predictions
      post <- bayesianMeanTracker(x = as.matrix(cbind(X1,X2)), y = (reward[j-1]-50)/100, prevPost = prevPost, theta = c(kError))
      prevPost <- post  #save new posterior as prevPost for next round
      #compute acquisition function evaluation and weight by inverse manhattan distance
      utilityVec <- acq(post, pars = c(beta)) / distance
      #prevent overflow by subtracting by max
      utilityVec <- utilityVec - max(utilityVec)
      #compute softmax choice probabilities
      p <- exp(utilityVec/tau)
      p <- p/sum(p)
      #Sample next choice
      location <- sample(1:121,1, prob=p, replace=TRUE)
      #update reward, X1, X2, and Y for both smooth and rough
      reward[j] <- smoothEnvironments[[envNum]][location,"y"] * 100
      X1 <- c(X1, choices[location, 'x1'])
      X2 <- c(X2, choices[location, 'x2'])
      Y <- c(Y,  reward[j])}
    reward}, mc.preschedule = TRUE, mc.cores=cores)
  #run for smooth environments
  roughReward <- mcmapply(1:replications, FUN=function(x){
    #sample parameters
    params <- parameters[sample(1:ncol(roughPars), 1),]
    kError <- params$kError
    beta <- params$beta
    #tau <- params$tau
    tau <- median(roughPars$tau)
    #randomly choose environment
    envNum <- sample(1:20,1) 
    #1st trial is random
    location <- sample(1:121,1)#each location as an integer from 1:121; first value is random, so treat both environments as replications of the same actions
    #Observations
    X1 <- choices[location,'x1']
    X2 <- choices[location,'x2']
    #reward
    reward<- c()
    reward[1] <- Y <- roughEnvironments[[envNum]][location,"y"]*100
    #posterior
    prevPost <- NULL  #initialize previous posterior for sequential updating of BMT posterior
    for (j in 2:41){ #after that, loop through remaining trials and make decisions based on GP preditions
      #compute Manhattan Distances
      prev <- unlist(manhattan[location,])
      distance <- sapply(seq(1:121), FUN=function(x) abs(manhattan[x, "x1"] - prev["x1"]) + abs(manhattan[x, "x2"] - prev["x2"]))
      #set values of 0 to 1
      distance[distance==0]<-1
      #update posterior predictions
      post <- bayesianMeanTracker(x = as.matrix(cbind(X1,X2)), y = (reward[j-1]-50)/100, prevPost = prevPost, theta = c(kError))
      prevPost <- post  #save new posterior as prevPost for next round
      #compute acquisition function evaluation and weight by inverse manhattan distance
      utilityVec <- acq(post, pars = c(beta)) / distance
      #prevent overflow by subtracting max
      utilityVec <- utilityVec - max(utilityVec)
      #compute softmax choice probabilities
      p <- exp(utilityVec/tau)
      p <- p/sum(p)
      #Sample next choice
      location <- sample(1:121,1, prob=p, replace=TRUE)
      #update reward, X1, X2, and Y for both smooth and rough
      reward[j] <- roughEnvironments[[envNum]][location,"y"] * 100
      X1 <- c(X1, choices[location, 'x1'])
      X2 <- c(X2, choices[location, 'x2'])
      Y <- c(Y,  reward[j])}
    reward}, mc.preschedule = TRUE, mc.cores=cores)
  #put into dataFrame
  localbmtDF <-data.frame(trial=rep(seq(1:41), 2),
                     Environment=c(rep("Smooth", 41), rep("Rough", 41)),
                     meanReward = c(rowMeans(smoothReward), rowMeans(roughReward)), #reward averaged over replications
                     meanSE = c(apply(smoothReward, 1, FUN = function(x) sd(x)/sqrt(length(x))), apply(roughReward, 1, FUN = function(x) sd(x)/sqrt(length(x)))),
                     maxReward =  c(rowMeans(apply(smoothReward, 2, FUN=function(x) maxton(x))),  rowMeans(apply(roughReward, 2, FUN=function(x) maxton(x)))), #apply maxton function over columns (trials) to compute the max reward at each time t, and then average the max values over repliations
                     maxSE =  c(apply(apply(smoothReward, 2, FUN=function(x) maxton(x)), 1, FUN=function(x)  sd(x)/sqrt(length(x))), apply(apply(roughReward, 2, FUN=function(x) maxton(x)), 1, FUN=function(x)  sd(x)/sqrt(length(x)))))
  #add model label
  localbmtDF$Model <- rep("Local BMT-UCB", 82)
  #write to csv
  if(!is.null(outputfile)){
    write.csv(localbmtDF, outputfile)
  }
  return(localbmtDF)
}

localbmtDF <- localbmtRationalModel(reps, outputfile = "rationalModels/localBMTUCBmedianTau.csv", parameters = localbmtPars)

#############################################################################################################################
# PLOTS
# Load all models
randomDF <- read.csv("rationalModels/random.csv")
gpDF <- read.csv("rationalModels/GPUCB.csv")
lgpDF <- read.csv("rationalModels/localGPUCB.csv")
bmtDF <- read.csv("rationalModels/BMTUCB.csv")
localbmtDF <- read.csv("rationalModels/localBMTUCB.csv")
#Load median tau models
randomDF <- read.csv("rationalModels/random.csv")
gpDF <- read.csv("rationalModels/GPUCBmedianTau.csv")
lgpDF <- read.csv("rationalModels/localGPUCBmedianTau.csv")
bmtDF <- read.csv("rationalModels/BMTUCBmedianTau.csv")
localbmtDF <- read.csv("rationalModels/localBMTUCBmedianTau.csv")
#############################################################################################################################

#Join all DFs together 
rationalDF <- rbind(randomDF, gpDF, lgpDF, bmtDF, localbmtDF)
rationalDF$Environment <- factor(rationalDF$Environment, levels=c("Rough", "Smooth"))


#convert to 5 round average
rationalDF$trial5 <- round((rationalDF$trial+1)/5)*5
rationalDF$trial5 <- ifelse(rationalDF$trial5<5,0,rationalDF$trial5)
rational5DF <- ddply(rationalDF, ~trial5+Model+Environment, summarise, meanReward=mean(meanReward), maxReward=mean(maxReward))

#add human data
source('dataMunging.R') #source data import function
d <- dataImport(normalize=FALSE)
#calculate 5 trial average
d$trial5<-round((d$trial+1)/5)*5
d$trial5<-ifelse(d$trial5<5,0,d$trial5)
dplot5<-ddply(d,~trial5+kernel,summarise,meanReward=mean(z), maxReward=mean(zmax))
colnames(dplot5)[2]<-"Environment" #rename colname
dplot5$Model <- rep("Human", nrow(dplot5))
#Join human with rational models
rational5 <- rbind(rational5DF, dplot5[,c(1,5,2,3,4)]) #join together but reorder columbs of dplot5

#main text plot
rational5 <- subset(rational5, Model = "Random" | Model=="LocalGP-UCB" | Model=="GP-UCB" | Model == "Human" | Model=="BMT-UCB" | Model=="Local BMT-UCB")
#reorder factor levels
rational5$Model <- factor(rational5$Model, levels = c("Random", "BMT-UCB", "Local BMT-UCB", "GP-UCB", "LocalGP-UCB", "Human"))

#Plot of mean Reward
p1<- ggplot(rational5, aes(x=trial5, y=meanReward, col=Model, shape=Model))+
  geom_point(size=2) +
  geom_line(size=.8) +
  #geom_errorbar(aes(ymin=meanReward-meanSE, ymax=meanReward+meanSE), width=0.1, size=1) +
  ylab("Average Reward")+xlab("Trial")+
  theme_classic()+
  facet_wrap(~Environment) +
  #scale_shape_manual(values = c(32, 16,17,15,4,7)) +
  scale_color_manual(values=c("black",  "#F0E442", "#E69F00", "#009E73", "#56B4E9", "#fb9a99"))+
  #scale_color_brewer(palette="Paired", direction=1)+
  coord_cartesian(ylim=c(45,85))+
  theme(text = element_text(size=12,  family="sans"), legend.position="top")+
  theme(legend.position="none", strip.background=element_blank(), legend.key=element_rect(color=NA))
p1
ggsave(filename = "plots/modelPerformanceAvg.pdf", plot = p1,height =2.5, width = 4.2, units = "in", useDingbats=FALSE) 

#ggsave(filename = "plots/modelPerformanceLegend.pdf", plot = p1 + theme(legend.position="right"),height =2.82, width = 8, units = "in", useDingbats=FALSE) 
#Plot of mean Reward
p2<- ggplot(rational5, aes(x=trial5, y=maxReward, col=Model, shape=Model))+
  geom_point(size=3) +
  geom_line(size=1) +
  #geom_errorbar(aes(ymin=meanReward-meanSE, ymax=meanReward+meanSE), width=0.1, size=1) +
  ylab("Maximum Reward")+xlab("Trial")+
  theme_bw()+
  facet_wrap(~Environment) +
  #scale_shape_manual(values = c(32, 16,17,15,4,7)) +
  scale_color_manual(values=c("black", "#E69F00", "#56B4E9", "#009E73","#F0E442", "#fb9a99"))+
  #scale_color_brewer(palette="Paired", direction=1)+
  coord_cartesian(ylim=c(65,100))+
  theme(text = element_text(size=16,  family="serif"), legend.position="top")+
  theme(legend.position="none", strip.background=element_blank(), legend.key=element_rect(color=NA))
p2

#Barplot of Max Reward
maxDF <- subset(rational5, trial5==40)
p3 <- ggplot(maxDF, aes(x=Model, y=maxReward,fill=Model))+
  stat_summary(fun.y = mean, geom = "bar", position = "dodge", color='black') + 
  #geom_errorbar(aes(ymin=meanReward-meanSE, ymax=meanReward+meanSE), width=0.1, size=1) +
  ylab("Max. Reward")+xlab("Trial")+
  theme_bw()+
  scale_fill_manual(values=c("black", "#E69F00", "#56B4E9", "#009E73","#F0E442", "#fb9a99"))+
  coord_cartesian(ylim=c(60,100))+
  facet_wrap(~Environment)+
  theme(text = element_text(size=16,  family="serif"))+
  theme(legend.position="none", strip.background=element_blank(), legend.key=element_rect(color=NA)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
p3
ggsave(filename = "plots/modelPerformanceMax.pdf", plot = p2, height =2.5, width = 6, units = "in", useDingbats=FALSE) 

#############################################################################################################################
# DEBUG: FIGURE OUT WHY GP Model HAS SUBLINEAR REGRET
# to read the previous simulation from disk:
# gpDF <- read.csv("rationalModels/GPUCB-Regret.csv")
#############################################################################################################################

#Dynamic beta
#fixed tau
rationalGPdynamic <- function(replications, outputfile, parameters, acq=ucb, kernel=rbf, cores=7, delta=.1, argmax=TRUE){
  #choice matrix
  choices <- expand.grid(0:10, 0:10) #build choice matrix
  names(choices)<-c("x1","x2")
  #run for smooth environments
  smoothReward <- mcmapply(1:replications, FUN=function(x){
    #sample parameters
    params <- parameters[sample(1:ncol(parameters), 1),]
    lambda <- params[1]
    beta <- params[2]
    tau <-params[3]
    #randomly choose environment
    envNum <- sample(1:20,1) 
    #1st trial is random
    location <- sample(1:121,1)#each location as an integer from 1:121; first value is random, so treat both environments as replications of the same actions
    #Observations
    X1 <- choices[location,'x1']
    X2 <- choices[location,'x2']
    #reward
    reward<- c()
    reward[1] <- Y <- smoothEnvironments[[envNum]][location,"y"]*100
    for (j in 2:101){ #after that, loop through remaining trials and make decisions based on GP preditions
      #compute posterior predictions
      post <- gpr(X.test = choices, theta = c(lambda, lambda, 1, 0.0001), X = cbind(X1,X2), Y = ((Y-50)/100), k = kernel) #scale observed Y to zero mean, variance of 1
      #compute acquisition function evaluation
      beta_t <- sqrt(2 * log(2*(j^2)*(pi^2)*(1/(6*delta))))/5 #Downscaling from Srinivas et al.
      utilityVec <- acq(post, pars = c(beta_t))
      if (argmax==TRUE){
        location <- which.max(utilityVec)
      }else{
        #scale to max of one
        utilityVec <- utilityVec/max(utilityVec)
        #compute softmax choice probabilities
        p <- exp(utilityVec/tau)
        p <- p/sum(p)
        #Sample next choice
        location <- sample(1:121,1, prob=p, replace=TRUE)
      }
      #update reward, X1, X2, and Y 
      reward[j] <- smoothEnvironments[[envNum]][location,"y"] * 100
      X1 <- c(X1, choices[location, 'x1'])
      X2 <- c(X2, choices[location, 'x2'])
      Y <- c(Y,  reward[j])}
    reward}, mc.preschedule = TRUE, mc.cores=cores)
  #run for smooth environments
  roughReward <- mcmapply(1:replications, FUN=function(x){
    #sample parameters
    params <- parameters[sample(1:ncol(parameters), 1),]
    lambda <- params[1]
    beta <- params[2]
    tau <- params[3]
    #randomly choose environment
    envNum <- sample(1:20,1) 
    #1st trial is random
    location <- sample(1:121,1)#each location as an integer from 1:121; first value is random, so treat both environments as replications of the same actions
    #Observations
    X1 <- choices[location,'x1']
    X2 <- choices[location,'x2']
    #reward
    reward<- c()
    reward[1] <- Y <- roughEnvironments[[envNum]][location,"y"]*100
    for (j in 2:101){ #after that, loop through remaining trials and make decisions based on GP preditions
      #compute posterior predictions
      post <- gpr(X.test = choices, theta = c(lambda, lambda, 1, 0.0001), X = cbind(X1,X2), Y = ((Y-50)/100), k = kernel) #scale observed Y to zero mean, variance of 1
      #compute acquisition function evaluation
      beta_t <- 2 * log(2*(j^2)*(pi^2)*(1/(6*delta)))
      utilityVec <- acq(post, pars = c(beta_t))
      if (argmax==TRUE){
        location<- which.max(utilityVec)
      }else{
        #scale to max of one
        utilityVec <- utilityVec/max(utilityVec)
        #compute softmax choice probabilities
        p <- exp(utilityVec/tau)
        p <- p/sum(p)
        #Sample next choice
        location <- sample(1:121,1, prob=p, replace=TRUE)
      }
      #update reward, X1, X2, and Y for both smooth and rough
      reward[j] <- roughEnvironments[[envNum]][location,"y"] * 100
      X1 <- c(X1, choices[location, 'x1'])
      X2 <- c(X2, choices[location, 'x2'])
      Y <- c(Y,  reward[j])}
    reward}, mc.preschedule = TRUE, mc.cores=cores)
  #put into dataFrame
  gpDF <-data.frame(trial=rep(seq(1:101), 2),
                    Environment=c(rep("Smooth", 101), rep("Rough", 101)),
                    meanReward = c(rowMeans(smoothReward), rowMeans(roughReward)), #reward averaged over replications
                    meanSE = c(apply(smoothReward, 1, FUN = function(x) sd(x)/sqrt(length(x))), apply(roughReward, 1, FUN = function(x) sd(x)/sqrt(length(x)))),
                    maxReward =  c(rowMeans(apply(smoothReward, 2, FUN=function(x) maxton(x))),  rowMeans(apply(roughReward, 2, FUN=function(x) maxton(x)))), #apply maxton function over columns (trials) to compute the max reward at each time t, and then average the max values over repliations
                    maxSE =  c(apply(apply(smoothReward, 2, FUN=function(x) maxton(x)), 1, FUN=function(x)  sd(x)/sqrt(length(x))), apply(apply(roughReward, 2, FUN=function(x) maxton(x)), 1, FUN=function(x)  sd(x)/sqrt(length(x)))))
  #add mmodel label
  gpDF$Model <- rep("GP-UCB", 202)
  #write to csv
  if(!is.null(outputfile)){
    write.csv(gpDF, outputfile)
  }
  return(gpDF)
}


gpDF <- rationalGPdynamic(100, outputfile = "rationalModels/GPUCB-Regret3.csv", parameters = gpPars)

ggplot(gpDF, aes(x=trial, y=meanReward, color=Environment)) + geom_line()
