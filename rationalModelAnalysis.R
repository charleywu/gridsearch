#Charley Wu 2016
#Script to run rational models in comparison to human behavior

#############################################################################################################################
# IMPORT DATA AND ENVIRONMENTS
#############################################################################################################################
#house keeping
rm(list=ls())

ptm<-proc.time()
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
setwd("experiment")
#load environments from json, unlist, transform into numeric, convert into matrix, and name dimensions
roughEnvironments <- lapply(fromJSON("kernel1.json", flatten=TRUE), FUN=function(x) matrix(as.numeric(unlist(x)), ncol=3, byrow=TRUE, dimnames=list(seq(1,121), c('x2', 'y', 'x1'))))
smoothEnvironments <- lapply(fromJSON("kernel2.json", flatten=TRUE), FUN=function(x) matrix(as.numeric(unlist(x)), ncol=3, byrow=TRUE, dimnames=list(seq(1,121), c('x2', 'y', 'x1'))))
setwd("..")#Step back into analysis folder
setwd("analysis")

#parameter estimates

localPars <- as.matrix(read.csv('rationalModels/parameters/local.csv'))
gpPars <- matrix(as.matrix(read.csv('rationalModels/parameters/gp.csv')), ncol=3)
localgpPars <- matrix(as.matrix(read.csv('rationalModels/parameters/localgp.csv')), ncol=3)

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
  manhattan <- expand.grid(0:10, 0:10) #build manhattan blocks
  names(manhattan)<-c("x1","x2")
  #Run for smooth environments
  smoothReward <- mcmapply(1:replications, FUN=function(x){
    tau <- sample(taulist,1) #sample a tau value for each replication
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
    tau <- sample(taulist,1) #sample a tau value for each replication
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
  #choice matrix
  choices <- expand.grid(0:10, 0:10) #build choice matrix
  names(choices)<-c("x1","x2")
  #run for smooth environments
  smoothReward <- mcmapply(1:replications, FUN=function(x){
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
    reward[1] <- Y <- smoothEnvironments[[envNum]][location,"y"]*100
    for (j in 2:41){ #after that, loop through remaining trials and make decisions based on GP preditions
      #compute posterior predictions
      post <- gpr(X.test = choices, theta = c(lambda, lambda, 1, 0.0001), X = cbind(X1,X2), Y = ((Y-50)/100), k = kernel) #scale observed Y to zero mean, variance of 1
      #compute acquisition function evaluation
      utilityVec <- acq(post, pars = c(beta))
      #scale to max of one
      utilityVec <- utilityVec/max(utilityVec)
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
    for (j in 2:41){ #after that, loop through remaining trials and make decisions based on GP preditions
      #compute posterior predictions
      post <- gpr(X.test = choices, theta = c(lambda, lambda, 1, 0.0001), X = cbind(X1,X2), Y = ((Y-50)/100), k = kernel) #scale observed Y to zero mean, variance of 1
      #compute acquisition function evaluation
      utilityVec <- acq(post, pars = c(beta))
      #scale to max of one
      utilityVec <- utilityVec/max(utilityVec)
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


gpDF <- gpRationalModel(reps, outputfile = "rationalModels/GPUCB.csv", parameters = gpPars)
#Hyper parameters
#lambda <- 0.5645783
#beta <- 0.4802621
#tau <- 0.02673054 #softmax parameter

#############################################################################################################################
# Inertia-GP Model
# to read the previous simulation from disk:
# lgpDF <- read.csv( "rationalModels/localGPUCB.csv")
#############################################################################################################################

localGPmodel<- function(replications, outputfile, parameters, acq=ucb, kernel=rbf, cores=7){
  #build manhattan blocks and choice matrix
  manhattan <- expand.grid(0:10, 0:10) 
  names(manhattan)<-c("x1","x2")
  choices <- manhattan
  #Run for smooth environments
  smoothReward <- mcmapply(1:replications, FUN=function(x){
    #sample parameters
    params <- parameters[sample(1:ncol(parameters), 1),]
    lambda <- params[1]
    beta <- params[2]
    tau <- params[3]
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
      #scale to max of one
      utilityVec <- utilityVec/max(utilityVec)
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
    params <- parameters[sample(1:ncol(parameters), 1),]
    lambda <- params[1]
    beta <- params[2]
    tau <- params[3]
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
      #scale to max of one
      utilityVec <- utilityVec/max(utilityVec)
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


lgpDF <- localGPmodel(reps, outputfile = "rationalModels/localGPUCB.csv", parameters = localgpPars)
#Hyper parameters
#lambda <- 0.7563726
#beta <- 0.4972964
#tau <- 0.1361702 #softmax parameter


#end timer
proc.time()-ptm
#############################################################################################################################
# PLOTS
#############################################################################################################################

#Join all DFs together 
rationalDF <- rbind(randomDF, localDF, gpDF, lgpDF)
rationalDF$Environment <- factor(rationalDF$Environment, levels=c("Smooth", "Rough"))

#convert to 5 round average
rationalDF$trial5 <- round(rationalDF$trial/5)*5
rationalDF$trial5 < ifelse(rationalDF$trial5<5,0,rationalDF$trial5)
rational5DF <- ddply(rationalDF, ~trial5+Model+Environment, summarise, meanReward=mean(meanReward), maxReward=mean(maxReward))

#add human data
source('dataMunging.R') #source data import function
d <- dataImport(normalize=FALSE)
#calculate 5 trial average
d$trial5<-round(d$trial/5)*5
d$trial5<-ifelse(d$trial5<5,0,d$trial5)
dplot5<-ddply(d,~trial5+kernel,summarise,meanReward=mean(z), maxReward=mean(zmax))
colnames(dplot5)[2]<-"Environment" #rename colname
dplot5$Model <- rep("Human", nrow(dplot5))
#Join human with rational models
rational5 <- rbind(rational5DF, dplot5[,c(1,5,2,3,4)]) #join together but reorder columbs of dplot5

#Plot of mean Reward
p1<- ggplot(rational5, aes(x=trial5, y=meanReward, col=Model, shape=Model, linetype=Environment))+
  geom_point(size=3) +
  geom_line(size=1) +
  #geom_errorbar(aes(ymin=meanReward-meanSE, ymax=meanReward+meanSE), width=0.1, size=1) +
  ylab("Avg. Reward")+xlab("Trial")+
  theme_bw()+
  scale_color_brewer(palette="Paired", direction=1)+
  coord_cartesian(ylim=c(45,75))+
  theme(text = element_text(size=24,  family="serif"), legend.position="top")+
  theme(legend.position="right", strip.background=element_blank(), legend.key=element_rect(color=NA))
p1

#Plot of Max Reward
p2 <- ggplot(rational5, aes(x=trial5, y=maxReward,col=Model, shape=Model, linetype=Environment))+
  geom_point(size=3) +
  geom_line(size=1) +
  #geom_errorbar(aes(ymin=meanReward-meanSE, ymax=meanReward+meanSE), width=0.1, size=1) +
  ylab("Max. Reward")+xlab("Trial")+
  theme_bw()+
  scale_color_brewer(palette="Paired", direction=1)+
  coord_cartesian(ylim=c(75,100))+
  theme(text = element_text(size=24,  family="serif"), legend.position="top")+
  theme(legend.position="none", strip.background=element_blank(), legend.key=element_rect(color=NA))
p2

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

legend1 <- get_legend(p1) #get legend
p1 <- p1 + theme(legend.position="none") #now remove legend from plot

#new composite graph
p3 <- grid.arrange(p1, p2, legend1, ncol=3, widths=c(6, 6, 4), top=textGrob("Model Performance",gp=gpar(fontsize=24, fontfamily="serif")))

ggsave(filename = "plots/rationalModel.pdf", p3, height =4.33, width = 9, units = "in")




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
