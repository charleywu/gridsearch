#Charley Wu 2018
#Experiment 3
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
setwd("experiment3")
#load environments from json, unlist, transform into numeric, convert into matrix, and name dimensions
environments <-fromJSON("environments/agridat.json", flatten=TRUE)
setwd("..")#Step back into analysis folder
setwd("analysis3")

#parameter estimates
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
  rewards <- mcmapply(1:replications, FUN=function(x) sample(environments[[sample(1:20,1)]][,'y'], 41, replace = TRUE)*100, mc.preschedule=TRUE, mc.cores=7)
  #put into dataFrame
  randomDF <- data.frame(trial=seq(1:41),
                         meanReward = rowMeans(rewards), #reward averaged over replications
                         meanSE = apply(rewards, 1, FUN = function(x) sd(x)/sqrt(length(x))),
                         maxReward = rowMeans(apply(rewards, 2, FUN=function(x) maxton(x))), #apply maxton function over columns (trials) to compute the max reward at each time t, and then average the max values over repliations
                         maxSE =  apply(apply(rewards, 2, FUN=function(x) maxton(x)), 1, FUN=function(x)  sd(x)/sqrt(length(x)))) 
  #add model label
  randomDF$Model <- rep("Random", 41)
  #write output
  if(!is.null(outputfile)){
    write.csv(randomDF, outputfile) 
  }
  return(randomDF)
}

randomDF <- randomModel(reps, "rationalModels/random.csv")

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
  rewards <- mcmapply(1:replications, FUN=function(x){
    #sample parameters
    params <- parameters[sample(nrow(parameters), 1), ]
    lambda <- params$lambda
    beta <- params$beta
    #tau <- params$tau
    tau <- median(params$tau)
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
  reward}, mc.preschedule = TRUE, mc.cores=cores)
  #put into dataFrame
  gpDF <- data.frame(trial=seq(1:41),
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

gpDF <- gpRationalModel(reps, outputfile = "rationalModels/GPUCBmedianTau.csv", parameters = gpPars)

#############################################################################################################################
# Local-GP Model
# to read the previous simulation from disk:
# lgpDF <- read.csv( "rationalModels/localGPUCB.csv")
#############################################################################################################################

localGPmodel<- function(replications, outputfile, parameters, acq=ucb, kernel=rbf, cores=7){
  #build manhattan blocks and choice matrix
  manhattan <- expand.grid(0:10, 0:10) 
  names(manhattan)<-c("x1","x2")
  choices <- manhattan
  #choice matrix
  choices <- expand.grid(0:10, 0:10) #build choice matrix
  names(choices)<-c("x1","x2")
  #run for smooth environments
  rewards <- mcmapply(1:replications, FUN=function(x){
    #sample parameters
    params <- parameters[sample(nrow(parameters), 1), ]
    lambda <- params$lambda
    beta <- params$beta
    #tau <- params$tau
    tau <- median(params$tau)
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
    reward}, mc.preschedule = TRUE, mc.cores=cores)
  #put into dataFrame
  lgpDF <- data.frame(trial=seq(1:41),
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

lgpDF <- localGPmodel(reps, outputfile = "rationalModels/localGPUCBmedianTau.csv", parameters = localgpPars)

#############################################################################################################################
# BMT-UCB Model
# to read the previous simulation from disk:
# bmtDF <- read.csv("rationalModels/BMTUCB.csv")
#############################################################################################################################

bmtRationalModel <- function(replications, outputfile, parameters, acq=ucb, cores=7){
  #choice matrix
  choices <- expand.grid(0:10, 0:10) #build choice matrix
  names(choices)<-c("x1","x2")
  #run for smooth environments
  rewards <- mcmapply(1:replications, FUN=function(x){
    #sample parameters
    params <- parameters[sample(nrow(parameters), 1), ]
    kError <- params$kError
    beta <- params$beta
    #tau <- params$tau
    tau <- median(params$tau)
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
      reward[j] <- environments[[envNum]][location,"y"] * 100
      X1 <- c(X1, choices[location, 'x1'])
      X2 <- c(X2, choices[location, 'x2'])
      Y <- c(Y,  reward[j])}
    reward}, mc.preschedule = TRUE, mc.cores=cores)
 
  #put into dataFrame
  bmtDF <- data.frame(trial=seq(1:41),
                      meanReward = rowMeans(rewards), #reward averaged over replications
                      meanSE = apply(rewards, 1, FUN = function(x) sd(x)/sqrt(length(x))),
                      maxReward = rowMeans(apply(rewards, 2, FUN=function(x) maxton(x))), #apply maxton function over columns (trials) to compute the max reward at each time t, and then average the max values over repliations
                      maxSE =  apply(apply(rewards, 2, FUN=function(x) maxton(x)), 1, FUN=function(x)  sd(x)/sqrt(length(x)))) 
  #add model label
  bmtDF$Model <- rep("BMT-UCB", 41)
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
  #build manhattan blocks and choice matrix
  manhattan <- expand.grid(0:10, 0:10) 
  names(manhattan)<-c("x1","x2")
  choices <- manhattan
  #run for smooth environments
  rewards <- mcmapply(1:replications, FUN=function(x){
    #sample parameters
    params <- parameters[sample(nrow(parameters), 1), ]
    kError <- params$kError
    beta <- params$beta
    #tau <- params$tau
    tau <- median(params$tau)
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
    #posterior
    prevPost <- NULL  #initialize previous posterior for sequential updating of BMT posterior
    for (j in 2:41){ #after that, loop through remaining trials and make decisions based on LBMT-UCB preditions
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
      #to prevent overflow, subtract max of q(x) vector 
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
    reward}, mc.preschedule = TRUE, mc.cores=cores)
  
  #put into dataFrame
  localbmtDF <- data.frame(trial=seq(1:41),
                      meanReward = rowMeans(rewards), #reward averaged over replications
                      meanSE = apply(rewards, 1, FUN = function(x) sd(x)/sqrt(length(x))),
                      maxReward = rowMeans(apply(rewards, 2, FUN=function(x) maxton(x))), #apply maxton function over columns (trials) to compute the max reward at each time t, and then average the max values over repliations
                      maxSE =  apply(apply(rewards, 2, FUN=function(x) maxton(x)), 1, FUN=function(x)  sd(x)/sqrt(length(x)))) 
  #add model label
  localbmtDF$Model <- rep("Local BMT-UCB", 41)
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
gpDF <- read.csv("rationalModels/GPUCBmedianTau.csv")
lgpDF <- read.csv("rationalModels/localGPUCBmedianTau.csv")
bmtDF <- read.csv("rationalModels/BMTUCBmedianTau.csv")
localbmtDF <- read.csv("rationalModels/localBMTUCBmedianTau.csv")
#############################################################################################################################

#Join all DFs together 
rationalDF <- rbind(randomDF, gpDF, lgpDF, bmtDF, localbmtDF)


#convert to 5 round average
rationalDF$trial5 <- round((rationalDF$trial+1)/5)*5
rationalDF$trial5 <- ifelse(rationalDF$trial5<5,0,rationalDF$trial5)
rational5DF <- ddply(rationalDF, ~trial5+Model, summarise, meanReward=mean(meanReward), maxReward=mean(maxReward))

#add human data
source('dataMunging.R') #source data import function
d <- dataImport(normalize=FALSE)
#calculate 5 trial average
d$trial5<-round((d$trial+1)/5)*5
d$trial5<-ifelse(d$trial5<5,0,d$trial5)
dplot5<-ddply(d,~trial5,summarise,meanReward=mean(z), maxReward=mean(zmax))
dplot5$Model <- rep("Human", nrow(dplot5))
#Join human with rational models
rational5 <- rbind(rational5DF, dplot5[,c(1,4,2,3)]) #join together but reorder columbs of dplot5

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
  #scale_shape_manual(values = c(32, 16,17,15,4,7)) +
  scale_color_manual(values=c("black",  "#F0E442", "#E69F00", "#009E73", "#56B4E9", "#fb9a99"))+
  #scale_color_brewer(palette="Paired", direction=1)+
  coord_cartesian(ylim=c(45,85))+
  theme(text = element_text(size=12,  family="sans"), legend.position="right")+
  theme(strip.background=element_blank(), legend.key=element_rect(color=NA))
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
  #scale_shape_manual(values = c(32, 16,17,15,4,7)) +
  scale_color_manual(values=c("black", "#E69F00", "#56B4E9", "#009E73","#F0E442", "#fb9a99"))+
  #scale_color_brewer(palette="Paired", direction=1)+
  coord_cartesian(ylim=c(65,100))+
  theme(text = element_text(size=16,  family="serif"), legend.position="top")+
  theme(legend.position="bottom", strip.background=element_blank(), legend.key=element_rect(color=NA))
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
