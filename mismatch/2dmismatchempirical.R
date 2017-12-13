#2D mismatch simulations
#Eric Schulz, 2017

setwd("..")
setwd("analysis2D")
source("Models.R") #Load models
setwd("..")
setwd("mismatch")
packages <- c('plyr', 'MASS', 'matrixcalc')
lapply(packages, require, character.only = TRUE)

clusterid <- as.integer(commandArgs(TRUE)[1]) #retrieve clusterid for running it in parallel across multiple computers

#Median participant estimates
tau<-0.0155044
beta<-0.4696316

Xstar<-as.matrix(expand.grid(0:10, 0:10)) #Matrix of options
L<-expand.grid(seq(0.1,3,0.1), seq(0.1,3,0.1)) #teacher vs. student lambda combinations

#initalize dataframe
dat<-data.frame(trial=numeric(), round=numeric(), lteach=numeric(), llearn=numeric(), x1=numeric(),x2=numeric(), y=numeric())

#Assign teacher and student lambdas based on cluster id
lteach<-L[clusterid,1]
llearn<-L[clusterid,2]

tcov<-rbf(Xstar, Xstar, c(lteach,lteach,1,0)) #Initialize RBF kernel using teacher lambda

#Simuluate 100 replications
for (j in 1:100){
  m<-mvrnorm(1, rep(0,nrow(Xstar)), tcov) #sample environment function from GP prior
  f<-data.frame(x1=Xstar[,1], x2=Xstar[,2], y=m) #reshape into dataframe
  f$y<-f$y-mean(f$y) #center mean of function at 0
  samp <- sample(1:121,1) #first choice is random
  Y<-matrix(f$y[samp]) #add to observations
  x1<-f$x1[samp]
  x2<-f$x2[samp]
  X<-matrix(c(x1, x2),ncol=2)
  for (k in 2:40){  #use GP-UCB using student lambda to simulate search behavior for 40 trials
    #compute posterior predictions
    post <- gpr(X.test = Xstar, theta = c(llearn, llearn, 1, 0.0001), X = X, Y = Y, k = rbf) #scale observed Y to zero mean, variance of 1
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
    reward <- f$y[location]
    x1<-f$x1[location]
    x2<-f$x2[location]
    X <- rbind(X, matrix(c(x1, x2), ncol=2))
    Y <- rbind(Y,  matrix(reward))
  }
  dummy<-data.frame(trial=1:40, round=rep(j, 40), lteach=lteach, llearn=llearn, x1=X[,1],x2=X[,1], y=Y)
  dat<-rbind(dat, dummy)
}

name<-paste0("/home/ucabchu/Scratch/mismatch/d2data", clusterid) #save data on cluster
write.csv(dat, paste0(name, '.csv'))
#End of cluster code

###################################################################
#Run the following on the cluster to compute all 900 teacher-student combinations
for (k in 1:900){
  dummy<-read.csv(paste0("d2data", k, ".csv"))
  if (k==1){d<-ddply(dummy, ~trial, summarize, m=mean(y))}
  if (k>1){d<-rbind(d, ddply(dummy, ~trial, summarize, m=mean(y)))}
}
write.csv(d, "mismatch2d.csv")

