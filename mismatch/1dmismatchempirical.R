#1D mismatch simulations
#Eric Schulz, 2017

setwd("..")
setwd("analysis1D")
source("Models.R") #Load models
setwd("..")
setwd("mismatch")

packages <- c('plyr', 'MASS', 'matrixcalc')
lapply(packages, require, character.only = TRUE)

clusterid <- as.integer(commandArgs(TRUE)[1]) #retrieve clusterid for running it in parallel across multiple computers

#Median participant estimates
tau<-0.007772736
beta<-0.5038621

Xstar<-matrix(1:30, ncol=1) #Matrix of options
L<-expand.grid(seq(0.1,3,0.1), seq(0.1,3,0.1)) #teacher vs. student lambda combinations

#initalize dataframe
dat<-data.frame(trial=numeric(), round=numeric(), lteach=numeric(), llearn=numeric(), X=numeric(), Y=numeric())

#Assign teacher and student lambdas based on cluster id
lteach<-L[clusterid,1]
llearn<-L[clusterid,2]

tcov<-rbf(Xstar, Xstar, c(lteach,1,0)) #Initialize RBF kernel using teacher lambda

#Simuluate 100 replications
for (j in 1:100){
  f<-mvrnorm(1, rep(0, length(Xstar)), tcov) #sample environment function from GP prior
  f<-f-mean(f) #center mean of function at 0
  X <- sample(1:30,1) #first choice is random
  Y<-f[X] #add to observations
  for (k in 2:11){ #use GP-UCB using student lambda to simulate search behavior for 10 trials
    #compute posterior predictions
    post <- gpr(X.test = Xstar, theta = c(llearn, 1, 0.0001), X = X, Y = Y, k = rbf) #scale observed Y to zero mean, variance of 1
    #compute acquisition function evaluation
    utilityVec <- ucb(post, pars = c(beta))
    #prevent overflow by subtracting max
    utilityVec <- utilityVec - max(utilityVec)
    #compute softmax choice probabilities
    p <- exp(utilityVec/tau)
    p <- p/sum(p)
    #Sample next choice
    location <- sample(1:30,1, prob=p, replace=TRUE)
    #update reward, X1, X2, and Y for both smooth and rough
    reward <- f[location]
    X <- c(X, location)
    Y <- c(Y,  reward)
  }
  dummy<-data.frame(trial=1:11, round=rep(j, 11), lteach=lteach, llearn=llearn, X=X, Y=Y)
  dat<-rbind(dat, dummy)
}

name<-paste0("/home/ucabchu/Scratch/mismatch/data", clusterid) #save data on cluster
write.csv(output, paste0(name, '.csv'))
#End of cluster code



###################################################################
#Run the following on the cluster to compute all 900 teacher-student combinations
for (k in 1:900){
  dummy<-read.csv(paste0("data", k, ".csv"))
  if (k==1){d<-ddply(dummy, ~trial, summarize, m=mean(Y))}
  if (k>1){d<-rbind(d, ddply(dummy, ~trial, summarize, m=mean(Y)))}
}
write.csv(d, "mismatch1d.csv")
