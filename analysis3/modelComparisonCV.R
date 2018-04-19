#Modeling results
#Charley Wu March 2018
#Experiment 3

rm(list=ls()) #house keeping

#load packages
packages <- c('plyr', 'jsonlite', 'DEoptim', "matrixcalc", "fields")
lapply(packages, require, character.only = TRUE)
#Source dependencies 
source('Models.R')
source('dataMunging.R')

##############################################################################################################
#Cluster configuration: (1 subject x model combination) per CPU
##############################################################################################################

#IMPORTANT: update batch name
batchName = 'first' #saves output csv files in 'modelResults/batchName/*.csv'

#InertiaWeighting of acquisition functions used in Local-GP models
inertiaWeight <- FALSE

#Cluster id from qsub
# clusterid <- runif(1,1,80)
clusterid <- as.integer(commandArgs(TRUE)[1]) #Cluster id, corresponds to an integer used to indicate which combination of kernel and acquisition function to simulate

#create list of all kernel functions
kernellist<-list(rbf) 
#kernellist<-list(bayesianMeanTracker, rbf) 

#names of all kernel functions
kernelnames<-c("RBF")
#kernelnames<-c("LBMT", "LRBF")

#list of all acquisition functions
#acqlist<-list(ucb, greedyVar, greedyMean, exofimp, poi, pmu)
acqlist<-list(poi)

#names of all acquisition functions
#acqnames<-c("UCB", "GV", "GM", "EXI", "POI", "PMU")
acqnames<-c("POI")

#all combinations of kernels and acquisition functions will be needed
combs<-expand.grid(1:length(kernellist), 1:length(acqlist))

#create a matrix with  combinations of subjectIds and model combinations
subjectComb <- expand.grid(1:80, 1:(length(kernellist) * length(acqlist))) #1:? defines the number of unique models to be analyzed
subjectId <- subjectComb[clusterid,1] #used to identify unique subjects
combId <- subjectComb[clusterid,2] #used to identify a unique model (i.e., combination of GP kernel and acquisition function)

#trim combs if cluster==TRUE
model <- combs[combId,]

set.seed(clusterid) #set seed as the participant id

##############################################################################################################
#Compile Participant Data
##############################################################################################################
data <- dataImport() #sourced from dataMunging.R

if (inertiaWeight){blocks <- readBlockDistance(subjectId)} #if localized models, initialize inertia weights by reading pre-computed csv files
##############################################################################################################
#Model Fitting
##############################################################################################################

#Negative Log Likelihood of model prediction 
#parameters are distinct for each individual model
#subjD = subset of data for a specific subject
#kernel = any of fastKalmanFilter, standard GP (rbf, matern, or oru), or lingp
#acquisition = any of ucb, probofimp, expofimp, or pmu
#horizonLength 20 or 40
#rounds list of numbers between 1:4, relative to the order of rounds with the specified horizon length
#inertiaWeight==TRUE weighs the q(x) of the acquisition function by the inverse manhattan distance
#refactorUCB: whether or not to refactor tau and beta parameters in UCB
modelFit<-function(par, subjD, acquisition, k,  horizonLength, rounds, inertiaWeight=inertiaWeight, refactorUCB=F){
  #Extract and process parameters
  par<-exp(par) #exponentiate parameters to make a non-negative and convex optimization surface
  #last parameter is always inverse temperature for softmax
  tau<-par[length(par)] 
  #Which posterior function to use; therefore, which parameters to use
  if (inherits(k, "KalmanFilter")){ #null kernel indicates kalman filter model
    kNoise <- par[1]
    parVec <- c(kNoise) #Vector of parameters to send to the KF posterior function
  }else if(inherits(k, "GP")){ #lambda
    lambda <- par[1]
    parVec <- c(lambda, lambda, 1, .0001) # Vector of parameters to send to the GP posterior vector, where sF and sN are fixed
  }
  #Additional acquisition function dependent parameters
  if (inherits(acquisition, "UCB")){ #check if UCB is used
    beta <- par[length(par)-1] #If UCB, beta is always 2nd last
    #refactor beta and tau into gamma and beta_star, where gamma = 1/tau and beta_star = beta/tau
    if (refactorUCB==TRUE){
      gamma<-1/tau
      beta_star <- beta/tau
    }
  }
  #which rounds to consider?
  horizondata <- subset(subjD, horizon==horizonLength)
  #which round numbers belong to the selected horizon?
  roundList <- unique(horizondata$round)[as.integer(rounds)]
  #construct training set
  trainingSet <- subset(subjD, round %in% roundList)
  #Vector to store negative log likelihods
  nLL <- rep(0,length(rounds)) 
  for (r in unique(trainingSet$round)){ #Begin looping through each round
    #subset of data for round r
    roundD <- subset(subjD, round==r)
    #is this round a short or long horizon?
    horizon = roundD$horizon[1]
    #Observations of subject choice behavior
    chosen <- roundD$chosen
    chosen <- chosen[2:length(chosen)] # trim first observation, since it wasn't a choice but a randomly revealed tile
    y  <- roundD$z[1:horizon] #trim off the last observation, because it was not used to inform a choice (round already over)
    x1 <- roundD$x[1:horizon]
    x2 <- roundD$y[1:horizon]
    #create observation matrix
    X<-as.matrix(cbind(x1,x2))
    #initialize Xtest
    Xnew<-as.matrix(expand.grid(0:10,0:10))
    #make sure X is a matrix
    X<-as.matrix(X)
    #Utilties of each choice
    utilities <- NULL
    prevPost <- NULL #set the previous posterior computation to NULL for the kalman filter
    #loop through observations
    for (i in 1:horizon){ #skip the last observation, because no choice was made based on that information
      #new observation
      X1<-matrix(X[1:i,], ncol=2)
      y1<-y[1:i]
      #Which posterior function to use
      if (inherits(k, "KalmanFilter")){# kalman filter model
        out<- bayesianMeanTracker(x = X1[i,], y=y[i], prevPost = prevPost, theta = parVec)
        #update prevPost for the next round
        prevPost <- out
      }else if (inherits(k, "linGP")){#linear kernel GP
        out <- lingp(subjD$id[1], i, as.integer(r))
      }else if (inherits(k, "inertia")){#inertia model
        out <- inertia(blocks, i, as.integer(r))
        beta <- 0 #set a dummy beta model to make it work with UCB
      }else if (inherits(k, "WSLS")){ #win stay lose go
        out<- WSLS(y1, X1)
      }else{# GP with length-scale parameterized kernel
        out <- gpr(X.test=Xnew, theta=parVec, X=X1, Y=y1, k=k) #Mu and Sigma predictions for each of the 121 arms; either GP or Kalman filter
      }
      #Slightly different function calls for each acquisition function
      if (inherits(acquisition, "UCB")){ #UCB takes a beta parameter
        if (refactorUCB==TRUE){
          utilityVec <- acquisition(out, c(gamma, beta_star))
        } else{
          utilityVec<-acquisition(out, c(beta))
        }
      }else if(inherits(acquisition, "Imp")){ #ProbImp or ExpectedImp
        y.star <- max(roundD$z[1:i]) #Best revealed solution
        utilityVec <- acquisition(out, y.star)
      }else{ #PMU or any other
        utilityVec <- acquisition(out)
      }
      if (inertiaWeight==TRUE){ #if inertiaWeighting is specified
        #weight the utilityVec by the inverse manhattan distance
        IMD<- inertia(blocks, i, as.integer(r)) 
        utilityVec <- utilityVec * IMD$mu
      }
      utilityVec <- utilityVec - max(utilityVec) #avoid overflow
      utilities <- rbind(utilities, utilityVec) # build horizon_length x 121 matrix, where each row holds the utilities of each choice at each decision time in the search horizon
    }
    #Softmax rule
    if (inherits(acquisition, "UCB") & refactorUCB==TRUE){
      #no need for tau because it is already factored into gamma and beta_star
      p<- exp(utilities)
    }else{#all other acquisition functions still need to be scaled by tau
      p <- exp(utilities/tau)
    }
    p <- p/rowSums(p)
    #avoid underflow by setting a floor and a ceiling
    p <- (pmax(p, 0.00001))
    p <- (pmin(p, 0.99999))
    #Calculate Negative log likelihood
    nLL[which(unique(trainingSet$round)==r)] <- -sum(log(p[cbind(c(1:horizon),chosen)]))
  }#end loop through rounds
  return(sum(nLL))  #Return negative log likelihoods of all observations 
}

##############################################################################################################
#CROSS VALIDATION FUNCTION
##############################################################################################################
#function to plug in to the optimaztion routine
#selector: scalar, indicates a specific participant
#kernel, function, can be "rbf", "oru", "matern", or NULL (Kalman filter)
#acquisition, function, can be "ucb", "probofimp", "expofimp", or "PMU"
#horizonLength is 20 or 40
#leaveoutindex is 1,2,3, or 4
#inertiaWeight: whether nor not acquisition functions are weighted by inertia
cvfun<-function(selector, kernelFun, acquisition, horizonLength, leaveoutindex, inertiaWeight){
  #subselect participant, horizon and rounds not left out
  d1<-subset(data, id==selector)
  #training set
  rounds <- c("1","2","3","4")
  trainingSet <- rounds[! rounds==rounds[leaveoutindex]] #remove round specified by leaveoutindex
  #test set
  testSet <- rounds[leaveoutindex]
  #determine parameters for different model combinations
  if (inherits(kernelFun, "linGP") | inherits(kernelFun, "WSLS")){ #linearGP or WSLS
    nParams <- 0 #no kernel parameters
  }else { #normal GP has lambda and BMT has error variance
    nParams <- 1 
  }
  if (inherits(acquisition, 'UCB')){
    nParams <- nParams + 1 #add beta parameter
  }
  nParams <- nParams + 1 #add one for tau, which is in all models
  if ( inherits(kernelFun, "inertia")){ #special case for the baseline inertia model, which only has tau
    nParams <- 1 
  }
  #Set upper and lower bounds based on nParams
  lbound <- rep(-5, nParams)
  ubound <- rep(5, nParams)
  #Begin cross validation routine
  if (nParams>=2){#if 2 or more parameters
    #TRAINING SET
    fit<-DEoptim(modelFit, lower=lbound, upper=ubound, subjD=d1, k=kernelFun, horizonLength = horizonLength, rounds = trainingSet, acquisition=acquisition,inertiaWeight=inertiaWeight, DEoptim.control(itermax=100))
    paramEstimates <- fit$optim$bestmem #MODEL DEPENDENT PARAMETER ESTIMATES
    #TEST SET
    predict <- modelFit(par=paramEstimates, subjD=d1, acquisition=acquisition, k=kernelFun, horizonLength=horizonLength, rounds=testSet, inertiaWeight=inertiaWeight)
    output <- c(horizonLength, leaveoutindex, predict, fit$optim$bestmem) #horizonLength, leaveoutindex, nLL, parameters....
  } else{
    #TRAINING SET
    fit<-optimize(modelFit, lower=lbound, upper=ubound, subjD=d1, k=kernelFun, horizonLength = horizonLength, rounds = trainingSet, acquisition=acquisition, inertiaWeight=inertiaWeight)
    paramEstimates <- fit$minimum #MODEL DEPENDENT PARAMETER ESTIMATES
    #TEST SET
    predict <- modelFit(par=paramEstimates, subjD=d1, acquisition=acquisition, k=kernelFun, horizonLength=horizonLength, rounds=testSet, inertiaWeight=inertiaWeight)
    output <- c(horizonLength, leaveoutindex, predict, fit$minimum) #horizonLength, leaveoutindex, nLL, parameters....
  }
  return(output) #return optimized value
}

##############################################################################################################
#OPTIMIZATION ROUTINE
##############################################################################################################
output <- c()
#cross-validation routine
for (h in c(20,40)){ #loop through horizons
  for (loo in 1:4){#leave one out
    cv <- cvfun(subjectId, kernelFun=kernellist[[model[[1]]]], acquisition = acqlist[[model[[2]]]], horizonLength = h, leaveoutindex=loo, inertiaWeight=inertiaWeight)
    output <- rbind(output, cv)
  }
}

#save the vector with kernel-acquisition-pair as name
name<-paste0("modelResults/", batchName, "/", kernelnames[model[[1]]], acqnames[model[[2]]], subjectId)
write.csv(output, paste0(name, '.csv'))

##############################################################################################################
#THE END
##############################################################################################################
