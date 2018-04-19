#Functions for Kalman filter, Gaussian Process Kernels, Matrix Inversion, and Acquisition Functions
#Charley Wu and Eric Schulz 2018
require("MASS")
##############################################################################################################
#KERNELS
##############################################################################################################

#Radial Basis Kernel
rbf <- function(X1,X2,theta){
  #transfer to matrices
  X1 <- as.matrix(X1)
  X2 <- as.matrix(X2)
  #check dimensions
  if(ncol(X1) != ncol(X2)){
    stop("X1 and X2 must contain input values of the same dimension.")
  } else if(!all(theta>=0)){
    stop("All parameters must be >= 0.")
  }
  #get dimensions
  N1 <- nrow(X1)
  N2 <- nrow(X2)
  d <- ncol(X1)
  #initialize sigma
  sigma <-  matrix(rep(0, N1*N2),nrow=N1)
  #observational variance
  sf <- theta[d+1]
  #noise variance
  sn <- theta[d+2]
  #loop through
  for(i in 1:d){
    #length scale
    l <- theta[i] #Note: assumes a unique length scale for each dimension
    #x-diff
    xdiff <- (outer(X1[,i],X2[,i],function(x,y) x - y)/l)^2
    sigma <- sigma + xdiff
  }
  #RBF function
  if(identical(X1,X2)){
    id <- diag(rep(1,N1))
    sigma.final <- sf*exp(-0.5*sigma) + sn*id
  } else {
    sigma.final <- sf*exp(-0.5*sigma)
  }
  #return final covariance matrix
  return(sigma.final)
}
class(rbf)<- c(class(rbf), "GP") #identify the rbf kernel as a gp model

#Ornstein-Uhlenbeck
oru <- function(X1,X2,theta){
  #matrices
  X1 <- as.matrix(X1)
  X2 <- as.matrix(X2)
  #check dimensions
  if(ncol(X1) != ncol(X2)){
    stop("X1 and X2 must contain input values of the same dimension.")
  } else if(!all(theta>=0)){
    stop("All parameters must be >= 0.")
  }
  #get dimensions
  N1 <- nrow(X1)
  N2 <- nrow(X2)
  d <- ncol(X1) 
  #intialize matrix
  sigma <-  matrix(rep(0, N1*N2),nrow=N1) 
  #observation variance
  sf <- theta[d+1]
  #noise variance
  sn <- theta[d+2]
  #loop through
  for(i in 1:d){
    #length scale
    l <- theta[i] #Note: assumes a unique length scale for each dimension
    #x dash
    xdiff <- abs(outer(X1[,i],X2[,i],function(x,y) x - y)/l)
    sigma <- sigma + xdiff
  }
  #apply Ornstein-Uhlenbeck
  if(identical(X1,X2)){
    id <- diag(rep(1,N1))
    sigma.final <- sf*exp(-0.5*sigma) + sn*id
  } else {
    sigma.final <- sf*exp(-0.5*sigma)
  }
  #return covariance matrix
  return(sigma.final)
}

class(oru)<- c(class(oru), "GP") #identify the ornstein-uhlenbeck kernel as a gp model

#Matern 3/2
matern <- function(X1,X2,theta){
  #matrices
  X1 <- as.matrix(X1)
  X2 <- as.matrix(X2)
  #check dimensions
  if(ncol(X1) != ncol(X2)){
    stop("X1 and X2 must contain input values of the same dimension.")
  } else if(!all(theta>=0)){
    stop("All parameters must be >= 0.")
  }
  #get dimensions
  N1 <- nrow(X1) 
  N2 <- nrow(X2) 
  d <- ncol(X1)
  #intialize matrix
  sigma <-  matrix(rep(0, N1*N2),nrow=N1)
  #observation variance
  sf <- theta[d+1]
  #noise variance
  sn <- theta[d+2]
  #loop through
  for(i in 1:d){
    #length scale
    l <- theta[i] #Note: assumes a unique length scale for each dimension
    #x dash
    xdiff <- abs(outer(X1[,i],X2[,i],function(x,y) x - y)/l)
    sigma <- sigma + xdiff
  }
  #apply Mater 3/2
  if(identical(X1,X2)){
    id <- diag(rep(1,N1))
    sigma.final <- sf*(1+(sqrt(3)*sigma))*exp(-sqrt(3)*sigma) + sn*id
  } else {
    sigma.final <- sf*(1+(sqrt(3)*sigma))*exp(-sqrt(3)*sigma)
  }
  #return covariance matrix
  return(sigma.final)
}

class(matern)<- c(class(matern), "GP") #identify the rbf kernel as a gp model


##############################################################################################################
#INERTIA MODEL
##############################################################################################################
#The inertial model only uses inverse manhattan distance from the previous click in order to model behavior
#Here we use it like a posterior, where the means are the inverse manhattan distance, and the variance is always zero
#Combined with the UCB model, these posteriors will yield a behavior akin to inertia
#since the variance is always zero, it will be a greedy mean strategy
#model posteriors are read in from CSV files that were pre-computed using blockdistance/blockdistance.R

readBlockDistance <- function(participantid, dataFolder = 'blockdistance'){#read precomputed .csv files of manhattan block distance for a specific participant
  fileNumbers <- seq((((participantid-1) * 8) + 1),(((participantid-1) * 8) + 8))
  files <- paste0(dataFolder,'/block',fileNumbers, '.csv')
  blocks <- lapply(files, FUN=function(x) read.csv(x, header=FALSE))
  return(blocks)
}

inertia<-function(blocks, trialnumber, roundnumber){
  distance <- blocks[[roundnumber]][,trialnumber]
  predictions <- data.frame(mu=distance, sig = rep(0,121))
  colnames(predictions) <- c("mu", "sig")
  return(predictions)
}
class(inertia) <- c(class(inertia), "inertia")


##############################################################################################################
#WIN STAY LOSE SHIFT MODEL
##############################################################################################################
#win := reward_t >= reward* ; else lose
#stay := repeat or cardinal neighbors
#go := any unrevealed tile with equal probability
#WSLG takes a vector of previously observed reweards Yt at inputs Xt (matrix), where the last entry in both is the previous action
WSLS <- function(Yt, Xt){
  mu <- rep(0,121) #initialize matrix storing the value of each tile for t+1
  allopts<-expand.grid(0:10, 0:10) #expand matrix of bivariate function to a vector of 121 unique locations
  y_t <- Yt[length(Yt)] #most recent observation
  x_t <- Xt[length(Yt),]
  y_star <- max(Yt[1:(length(Yt)-1)]) #previous best reward
  if (y_t >= y_star){ #WIN
    #STAY
    #loop through all neighboring tiles as well as current tile, and set value to 1
    for(i in 1:-1){
      for(j in 1:-1){
        mu[which(x_t[1]+i==allopts$Var1 & x_t[2]+j==allopts$Var2)] <- 1 
      }
    }
  }else{ #LOSE
    #GO
    revealed <- apply(Xt,1, FUN=function(x) which(x[1]==allopts$Var1 & x[2]==allopts$Var2))
    mu[-unlist(revealed)] <- 1 #assign all unrevealed tiles a value of 1
  }
  predictions <- data.frame(mu=mu, sig = rep(0,121))
  colnames(predictions) <- c("mu", "sig") #WSLG returns the same kind of data frame as GP regression, but the variance is always set to 0
  return(predictions)
}

class(WSLS) <- c(class(WSLS), "WSLS")
##############################################################################################################
#MATRIX INVERSION
##############################################################################################################

#calculate inverse of the cov-function using sigular value decomposition
cov.inverse.svd <- function(X, tol = sqrt(.Machine$double.eps)){
  # Generalized Inverse of a Matrix
  dnx <- dimnames(X)
  if(is.null(dnx)) dnx <- vector("list", 2)
  #singular value decomposition
  s <- svd(X)
  nz <- s$d > tol * s$d[1]
  #inverse
  K.inv <- structure(
    if(any(nz)) s$v[, nz] %*% (t(s$u[, nz])/s$d[nz]) else X,
    dimnames = dnx[2:1])
  #logarithm of determinant.
  log.K.det <- sum(log(s$d))
  #return inverse plus log determinant
  return(list(Inv = K.inv,lDet = log.K.det))
}

#calculate inverse of the cov-function using Cholesky
cov.inverse.chol <- function(X){
  #cholseky decomposition
  R <- chol(X)
  #complex conjugate
  Rt <- Conj(t(R))
  #invert
  R.inv <- solve(R)
  #invert
  Rt.inv <- solve(Rt)
  #multiply matrices
  X.inv <- R.inv %*% Rt.inv
  #log determinant
  log.X.det <- 2*sum(log(diag(R))) 
  #return both
  return(list(Inv = X.inv, lDet = log.X.det))
}

##############################################################################################################
#GAUSSIAN PROCESS
##############################################################################################################

#Gaussian Process function
#X.test: matrix for predcitions
#theta: vector of hyper-parameters (l1, l2, Sf, Sn)
#X; matrix of observations
#y: vector of observed outcomes
#kernel: used kernel function, can be "rbf", "oru", or "mat"
gpr<- function(X.test, theta, X, Y, k){
  #make it a matrix
  Xstar <- as.matrix(X.test)
  #dimensions
  d <- ncol(as.matrix(X))
  #calculate capital K
  K <- k(X,X,theta) 
  #Check if matrix is positive semi-definite
  if (is.positive.definite(K)){
    #KK <- cov.inverse.chol(K) #use Eric's Cholesky function
    KK.inv <- chol2inv(chol(K)) #MASS implementation of Cholesky
  } else {
    KK.inv <- cov.inverse.svd(K) #use SVD
  }
  #times y
  Ky <- KK.inv %*% Y
  #apply the kernel
  result <- apply(Xstar, 1, function(x){
    XX <- matrix(x,nrow=1) 
    Kstar <- k(X, XX, theta)
    Kstarstar <- k(XX,XX,theta)
    #get mean vector
    mu <- t(Kstar) %*% Ky
    #get covariance
    cv <- Kstarstar - (t(Kstar) %*% KK.inv %*% Kstar) #BUG: sometimes cv<0, leading to NaN when we return sqrt(cv)
    #DEBUG
    if (cv<0){ cv <- abs(cv)} #TEMPORARY SOLUTION: MANUALLY SET CV TO BE POSITIVE IF NEGATIVE
    #return means and variance
    return(c(mu, cv))
  })
  #as a data frame with names mu and sig
  prediction <- as.data.frame(t(result))
  colnames(prediction) <- c("mu", "sig")
  return(prediction)
}


##############################################################################################################
## KALMAN FILTER
##############################################################################################################

#Kalman Filter
#X.test: matrix for predcitions
#theta: vector of hyper-parameters
#X; matrix of observations
#y: vector of observed outcomes
KalmanFilter<- function(X.test, theta, X, Y, k=NULL){
  library('FKF')
  #set of choices to make predictions about
  Xstar <- as.matrix(X.test)
  # parameters
  mu0 <- 0 #par[1]
  var0 <- 5 #exp(par[1])
  vart <- theta[1] #innovation variance
  vare <- theta[2] #error varriance
  #Which of the 121 options were chosen at each time?
  allopts<-expand.grid(0:10, 0:10)
  chosen <- apply(X, 1, FUN=function(x) which(allopts$Var1==x[1] & allopts$Var2==x[2]))
  #beta <- exp(theta[2]) 
  #alpha <- exp(par[3])
  #lambda <- 1
  # setup for fkf
  a0 <- rep(mu0,121) #initial estimates of state variable
  P0 <- var0*diag(121) #covariance matrix of a0
  dt <- matrix(rep(0,121),ncol=1) #intercept of the transition equation
  ct <- matrix(rep(0,121),ncol=1) #intercept of the measurement equation
  Tt <- array(diag(121),dim=c(121,121,1)) # Factor of the transition equation
  Zt <- array(diag(121),dim=c(121,121,1)) # factor of  the measurement equation
  HHt <- array(vart*diag(121),dim=c(121,121,1)) #variance of the innovation and transition equatinos
  GGt <- array(vare*diag(121),dim=c(121,121,1)) #disturbances of the measurement equation
  yt <- matrix(NA,ncol=121,nrow=length(Y)) #observations
  #Loop over observations and add payoffs to yt
  for(i in 1:length(Y)) {
    #for each round, fill in the yt that was chosen with the payoff in Y
    yt[i,chosen[i]] <- Y[i]
  } 
  filt <- fkf(a0,P0,dt,ct,Tt,Zt,HHt,GGt,t(yt))
  E <- t(filt$at) #expectation
  S <- sapply(1:nrow(E), FUN=function(x) sqrt(diag(filt$Pt[,,x])))
  result <- rbind(E[length(Y) +1,], t(S)[length(Y)+1,]) #return expectation and variance for each of the 121 options
  #as a data frame with names mu and sig
  prediction <- as.data.frame(t(result))
  colnames(prediction) <- c("mu", "sig")
  return(prediction)
}
class(KalmanFilter)<- c(class(KalmanFilter), "KalmanFilter")

#Faster version of the kalman filter that computes a single update
fastKalmanFilter <- function(x, y, theta=c(1,1), prevPost=NULL){ #Updates the previous posterior based on a single observation
  #parameters
  mu0 <- 0 #prior mean
  var0 <- 5 #prior variance
  vart <- theta[1] #innovation variance
  vare <- theta[2] #error varriance
  if (is.null(prevPost)){#if no posterior prior, assume it is the first observation
    predictions <- data.frame(mu=rep(mu0,121), sig=rep(var0,121))
  }else{#if previous posterior is provided, update
    predictions <- prevPost
    }
  #Which of the 121 options were chosen at time?
  allopts<-expand.grid(0:10, 0:10)
  chosen <- which(allopts$Var1==x[1] & allopts$Var2==x[2])
  #Kalman gain
  kGain <- (predictions$sig[chosen] + vart^2) / (predictions$sig[chosen] + vart^2 + vare^2)
  #update mean
  predictions$mu[chosen] <- predictions$mu[chosen] + (kGain * (y-predictions$mu[chosen]))
  #update variance
  predictions$sig <- predictions$sig + vart^2
  predictions$sig[chosen] <- predictions$sig[chosen] * (1 - kGain)
  #return output
  return(predictions)
}
class(fastKalmanFilter)<- c(class(fastKalmanFilter), "KalmanFilter")

#Like a kalman filter, but with stable-mean and no innovation variance
bayesianMeanTracker <- function(x, y, theta=c(1), prevPost=NULL){ #Updates the previous posterior based on a single observation
  #parameters
  mu0 <- 0 #prior mean
  var0 <- 5 #prior variance
  vare <- theta[1] #error varriance
  if (is.null(prevPost)){#if no posterior prior, assume it is the first observation
    predictions <- data.frame(mu=rep(mu0,121), sig=rep(var0,121))
  }else{#if previous posterior is provided, update
    predictions <- prevPost
  }
  #Which of the 121 options were chosen at time?
  allopts<-expand.grid(0:10, 0:10)
  chosen <- which(allopts$Var1==x[1] & allopts$Var2==x[2])
  #Kalman gain
  kGain <- predictions$sig[chosen] / (predictions$sig[chosen] + vare^2)
  #update mean
  predictions$mu[chosen] <- predictions$mu[chosen] + (kGain * (y-predictions$mu[chosen]))
  #update variance for observed arm
  predictions$sig[chosen] <- predictions$sig[chosen] * (1 - kGain)
  #return output
  return(predictions)
}
class(bayesianMeanTracker)<- c(class(bayesianMeanTracker), "KalmanFilter")

##############################################################################################################
#ACQUISITION FUNCTIONS
##############################################################################################################

#Upper Confidence Bound Sampling
ucb<-function(out, pars, refactor=F){
  if (refactor==TRUE){
    gamma <- pars[1]
    beta_star<-pars[2]
    #calulate all the upper confidence bounds
    outtotal<-(gamma*out$mu)+(beta_star*sqrt(out$sig)) #refactored parameters in combination with softmax tau, where gamma = 1/tau and beta_star = beta/tau
    #avoid borderline cases
    outtotal[outtotal<0]<-0.0001
    outtotal[outtotal>100]<-100
    outtotal<-matrix(outtotal, nrow=nrow(out)/121, byrow=TRUE)
  }else{
    beta <- pars[1]
    #calulate all the upper confidence bounds
    outtotal<-out$mu+(beta*sqrt(out$sig)) #refactored parameters in combination with softmax tau, where gamma = 1/tau and beta_star = beta/tau
    #avoid borderline cases
    outtotal[outtotal<0]<-0.0001
    outtotal[outtotal>100]<-100
    outtotal<-matrix(outtotal, nrow=nrow(out)/121, byrow=TRUE)
  }
  #return them
  return(outtotal)
}
#add "UCB" to the class of the ucb function, so that modelFit can recognize that it has a longer parameter array
class(ucb)<- c(class(ucb), "UCB")

#Probability of Maximum Utility (the model formerly known as Thompson Sampling)
pmu<-function(out){
  #initialize
  mus<-matrix(out$mu, nrow=nrow(out)/121, byrow=TRUE)
  sigs<-matrix(sqrt(out$sig), nrow=nrow(out)/121, byrow=TRUE) #convert variance to standard deviation
  outtotal<-mus
  #loop through
  for (i in 1:nrow(mus)){
    #simulate from normals given the means and standard deviation
    simout<-apply(cbind(mus[i,], sigs[i,]),1, function(x){rnorm(n=1000,mean=x[1],sd=x[2])})
    #collect mean simulated wins per arm
    outtotal[i,]<-apply(simout==apply(simout,1,max),2,mean) 
  }
  #return them
  return(outtotal+0.0001)
}

#Probability of Improvement
poi<-function(out, y.star){
  #calulate the probabilities of improvement
  #pi<-pnorm((out$mu-y.star)/out$sig)
  laplace <- sqrt(.Machine$double.eps) #Laplacian smoothing to avoid NaN values
  pi<-pnorm((out$mu-y.star+laplace)/(sqrt(out$sig)+laplace)) 
  outtotal<-matrix(pi, nrow(out)/121, byrow=TRUE)
  #return them
  return(outtotal)
}
#add "IMP" to the class of the probofimp function, so that modelFit can recognize that it takes y.star as an argument
class(poi)<- c(class(poi), "Imp")

#Expected improvement
exofimp<-function(out, y.star){
  #calulate z-scores first, then expected improvements
  laplace <- sqrt(.Machine$double.eps) #Laplacian smoothing to avoid NaN values
  z<-(out$mu-y.star+laplace)/(sqrt(out$sig)+laplace)
  ei <-(out$mu-y.star)*pnorm(z)+sqrt(out$sig)*dnorm(z)
  outtotal<-matrix(ei, nrow(out)/121, byrow=TRUE)
  #return them
  return(outtotal)
}
#add "IMP" to the class of the exofimp function, so that modelFit can recognize that it takes y.star as an argument
class(exofimp)<- c(class(exofimp), "Imp")

#Greedy Mean
greedyMean <- function(out){
  outtotal<-out$mu #the value of each choice is solely based on the expectation of reward
  #avoid borderline cases
  outtotal[outtotal<0]<-0.0001
  outtotal[outtotal>100]<-100
  outtotal<-matrix(outtotal, nrow=nrow(out)/121, byrow=TRUE)
}

#Greedy Variance
greedyVar <- function(out){
  outtotal<- sqrt(out$sig) #the value of each choice is solely based on the expected uncertainty
  #avoid borderline cases
  outtotal[outtotal<0]<-0.0001
  outtotal[outtotal>100]<-100
  outtotal<-matrix(outtotal, nrow=nrow(out)/121, byrow=TRUE)
}

