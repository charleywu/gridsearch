#Functions for Data recovery of GP models
#Charley Wu 2018

packages <- c('plyr', 'jsonlite', 'reshape', 'DEoptim', "matrixcalc", "fields", 'gplots', 'RColorBrewer', 'ggplot2', 'gridExtra')
lapply(packages, require, character.only = TRUE)
source("Models.R")

##############################################################################################################
#Data recovery 
##############################################################################################################

###Multiplot function #######

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
#Data recovery of subject data
#Choose a subject id, a round, and the specific trial to recover
#Specify, kernel type, acquisition function, and parameters
#yields the revealed rewards, the posterior beliefs about the space, and the softmax probability surface of the acquisition function

GPPosterior <- function(subject=1, r=1, trial = 20, kernel=rbf,kernelName="RBF", acq = ucb, acqName="UCB", parVec= c(1,1,1,1), beta=1, tau = 1){ #pargp <-(lambda1, lambda2, signal variance, noise variance)
  #Subset a single subject
  subjD <- subset(d, id==subject) #NOTE: d is loaded from the .rmd file
  #subset data for round
  roundD <- subset(subjD, round==r)
  #short or long horizon round? Selects the integer value of the horizon length
  horizon = roundD$horizon[1]
  #Compile observations until the given trial
  y  <- roundD$z[1:trial]
  x1 <- roundD$x[1:trial]
  x2 <- roundD$y[1:trial]
  #create observation matrix
  X<-as.matrix(cbind(x1,x2))
  #initialize Xtest
  Xnew<-as.matrix(expand.grid(0:10,0:10))
  #make sure X is a matrix
  X<-as.matrix(X)
  #Compute kernel posterior predictions
  #rescale y just like in modelComparison.R
  yscaled <- (y - 50)/100
  if (inherits(kernel, "KalmanFilter")){ #kalman filter
    out <- KalmanFilter(X.test=Xnew, theta=parVec, X=X, Y=yscaled)
  }else if (inherits(kernel, "linGP")){ #linear GP
    out <- lingp(subject, trial, r) #linear kernel GP posterior
  }else{#GP
    out <- gpr(X.test=Xnew, theta=parVec, X=X, Y=yscaled, k=kernel)
  }
  #Compute softmax of acquition functions choice
  #Slightly different function calls for each acquisition function
  if (inherits(acq, "UCB")){ #UCB takes a beta parameter
    utilityVec <- acq(out, beta)
  }else if(inherits(acq, "Imp")){ #ProbImp or ExpectedImp
    y.star <- max(yscaled) #Best revealed solution
    utilityVec <- acq(out, y.star)
  }else{ #PMU or any other
    utilityVec <- acq(out)
  }
  p <-utilityVec/max(utilityVec) #rescale to max value of 1 to prevent inf or -inf
  p <- exp(p/tau)
  p <- p/rowSums(p)
  #avoid underflow by setting a floor and a ceiling
  p <- (pmax(p, 0.00001))
  p <- (pmin(p, 0.99999))
  #Next choice at t = trial + 1
  print(sum(p))
  nextChoice <- c(roundD$x[trial+1], roundD$y[trial+1])
  #Plot 1. Subject's search space, which has been revealed until the specific trial number defined by "trial"
  #construct dataframe
  d1 <- data.frame(cbind(Xnew,rep(NA, 121))) # fill all 121 spots with NA
  colnames(d1) <- c("x1", "x2", "y")
  #replace NA with observations
  for (row in seq(1:length(y))){ 
    obs <- X[row,]
    d1$y[d1$x1==obs[1] & d1$x2==obs[2]] <- y[row]
  }
  p1 <- ggplot(d1, aes(x = x1+1, y = x2+1, fill = y)) + 
    geom_tile(color='black', width=1, height=1) +
    theme_bw() +
    coord_equal() +
    xlim(0.5,11.5) +
    ylim(0.5,11.5) + 
    ggtitle(paste0('Revealed (t=',trial,")")) +
    #scale_fill_gradientn(name='Payoff', colours = hm.palette(100), values = seq(0, 100, length=9),  rescaler = function(x, ...) x, oob = identity) +
    scale_fill_distiller(palette = "Spectral",limits=c(-5,105), na.value = 'white', breaks=c(0,25,50,75,100))+
    labs(fill="Payoff     ")+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border=element_blank(), 
          plot.margin=unit(c(0,0,0,0),"mm"))
    
  #Plot 2. Posterior prediction
  d2 <- melt(matrix((out[,1]*100) + 50, nrow=11, ncol=11))
  names(d2) <- c('X1', 'X2', 'value')
  if (inherits(kernel, "KalmanFilter")){ #kalman filter
    plotTitle <- bquote(paste(.(kernelName), " Posterior (", sigma["i"],"=", .(round(parVec[1],2)), "; ", sigma["n"], "=", .(round(parVec[1],2)),")"))
  }else if(inherits(kernel, "GP")){ #GP
    plotTitle <- bquote(paste(.(kernelName), " Posterior (", lambda, "=", .(round(parVec[1],2)), ")"))
  }else{
    plotTitle <- bquote(paste(.(kernelName), " Posterior "))
  }
  p2<-ggplot(d2, aes(x = X1, y = X2, fill = value)) + 
    geom_tile(color='black', width=1, height=1) +
    theme_bw() +
    xlim(0.5,11.5) +
    ylim(0.5,11.5) + 
    coord_equal() +
    ggtitle(plotTitle) +
    #scale_fill_gradientn(name = "Exp. Payoff", colours = hm.palette(100),values = seq(0, 100, length=9)) +
    scale_fill_distiller(palette = "Spectral",limits=c(-5,105), na.value = 'white',  breaks=c(0,25,50,75,100))+
    labs(fill="E(Payoff)")+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border=element_blank(), 
          plot.margin=unit(c(0,0,0,0),"mm")) 
  
  #Plot 3. Softmax surface of acquisition function
  d3<- melt(matrix(p, nrow=11, ncol=11))
  names(d3) <- c('X1', 'X2', 'value')
  d3<- round(d3, 3)
  #plot title
  if (!acqName=="UCB"){
    roundedTau <- round(tau,2)
    if (acqName == "ProbOfImp") {
      acqName <- "POI"
    }
    plottitle<- bquote(paste(.(acqName)," Prediction (", tau, "=", .(roundedTau), ")"))
  }else{
    roundedTau <- round(tau,2)
    roundedBeta <- round(beta, 2)
    plottitle<- bquote(paste(.(acqName), " Prediction (", beta, "=", .(roundedBeta), "; ", tau, "=", .(roundedTau), ")"))
  }
  p3<-ggplot(d3, aes(x = X1, y = X2, fill = value)) + 
    geom_tile(color='black', width=1, height=1) +
    theme_bw() +
    xlim(0.5,11.5) +
    ylim(0.5,11.5) + 
    coord_equal() +
    ggtitle(plottitle) +
    #scale_fill_gradientn(name='P(choice)',colours = hm.palette(100),values = seq(0, 1, length=9)) +
    scale_fill_distiller(palette = "Spectral", na.value = 'white' )+
    labs(fill="P(Choice)")+
    annotate("text", x = nextChoice[1] + 1, y = nextChoice[2] + 1, label = "X") +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border=element_blank(), 
          plot.margin=unit(c(0,0,0,0),"mm"))
  return(multiplot(p1, p2, p3, cols=1))
  #return(grid.arrange(p1,p2,p3,cols=3))
}


#Create animated gifs, all models for a particular subject in a particular round
animatedGifs <- function(subjId, round, kernel, acq, kernelName, acqName, plotname){
  subDF <- subset(paramEstimates, participant==subjId & acq==acqName & kernel==kernelName)
  #Median parameter estimates
  lambda <- median(subDF$lambda)
  tau <- median(subDF$tau)
  beta <- median(subDF$beta)
  beta[is.na(beta)] <- 0
  #loop through trials
  for (t in 1:horizon){
    graphics.off()
    #png(filename = paste0(getwd(),"/plots/gifs/",plotname,sprintf("%03d", t),'.png',sep=""))
    pdf(file= paste0(getwd(),"/plots/gifs/",plotname,sprintf("%03d", t),'.pdf',sep=""))
    p<-GPPosterior(subject=subjId, r=round, trial=t, kernel=kernel, kernelName = kernelName, acq=acq,acqName=acqName, parVec =c(lambda,lambda,1,.1), beta=beta, tau = tau)
  }
  graphics.off()
}


#Predictions for custom data
#X is the 2 * t observation matrix and y is the vector of payoffs
customPredictions <- function(X, y, kernel=rbf,kernelName="RBF", acq = ucb, acqName="UCB", parVec= c(1,1,1,1), beta=1, tau = 1){ #parVec <-(lambda1, lambda2, signal variance, noise variance)
  #initialize Xtest
  Xnew<-as.matrix(expand.grid(0:10,0:10))
  #make sure X is a matrix
  X<-as.matrix(X)
  #Compute kernel posterior predictions
  #rescale y just like in modelComparison.R
  yscaled <- (y - 50)/100
  if (inherits(kernel, "KalmanFilter")){ #kalman filter
    out <- KalmanFilter(X.test=Xnew, theta=parVec, X=X, Y=yscaled)
  }else if (inherits(kernel, "linGP")){ #linear GP
    out <- lingp(subject, trial, r) #linear kernel GP posterior
  }else{#GP
    out <- gpr(X.test=Xnew, theta=parVec, X=X, Y=yscaled, k=kernel)
  }
  #Compute softmax of acquition functions choice
  #Slightly different function calls for each acquisition function
  if (inherits(acq, "UCB")){ #UCB takes a beta parameter
    utilityVec <- acq(out, beta)
  }else if(inherits(acq, "Imp")){ #ProbImp or ExpectedImp
    y.star <- max(yscaled) #Best revealed solution
    utilityVec <- acq(out, y.star)
  }else{ #PMU or any other
    utilityVec <- acq(out)
  }
  p <-utilityVec - max(utilityVec) #rescale to max value of 1 to prevent inf or -inf
  p <- exp(p/tau)
  p <- p/rowSums(p)
  #avoid underflow by setting a floor and a ceiling
  p <- (pmax(p, 0.00001))
  p <- (pmin(p, 0.99999))
  #Next choice at t = trial + 1
  #Plot 1. Posterior prediction
  d1 <- melt(matrix((out[,1]*100) + 50, nrow=11, ncol=11))
  names(d1) <- c('X1', 'X2', 'value')
  if (inherits(kernel, "KalmanFilter")){ #kalman filter
    plotTitle <- bquote(paste(.(kernelName), " Posterior (", sigma["i"],"=", .(round(parVec[1],2)), "; ", sigma["n"], "=", .(round(parVec[1],2)),")"))
  }else if(inherits(kernel, "GP")){ #GP
    plotTitle <- bquote(paste(.(kernelName), " Posterior (", lambda, "=", .(round(parVec[1],2)), ")"))
  }else{
    plotTitle <- bquote(paste(.(kernelName), " Posterior "))
  }
  #DEBUG
  plotTitle <- "Inferred Value"
  p1<-ggplot(d1, aes(x = X1, y = X2, fill = value)) + 
    geom_tile(color='black', width=1, height=1) +
    theme_bw() +
    xlim(0.5,11.5) +
    ylim(0.5,11.5) + 
    coord_equal() +
    ggtitle(plotTitle) +
    #scale_fill_gradientn(name = "Exp. Payoff", colours = hm.palette(100),values = seq(0, 100, length=9)) +
    scale_fill_distiller(palette = "Spectral",limits=c(5,90), na.value = 'white',  breaks=c(25,50,75))+
    guides(ticks=element_line(color='black'))+
    labs(fill= bquote(paste("E","(x)")))+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border=element_blank(), 
          plot.margin=unit(c(0,0,0,0),"mm"))
  
  #Plot 2. Softmax surface of acquisition function
  d2<- melt(matrix((utilityVec*100)+50, nrow=11, ncol=11)) #IMPORTANT: check that rescaling is correct
  names(d2) <- c('X1', 'X2', 'value')
  #maxValue <- subset(d2, value==max(d2$value))
  #d2<- round(d2, 3)
  #plot title
  if (!acqName=="UCB"){
    roundedTau <- round(tau,2)
    if (acqName == "ProbOfImp") {
      acqName <- "POI"
    }
    plottitle <-  bquote(paste(.(acqName)," Prediction"))
  }else{ #other model
    roundedTau <- round(tau,2)
    roundedBeta <- round(beta, 2)
    plottitle<- bquote(paste(.(acqName), " Prediction (", beta, "=", .(roundedBeta),")"))
  }
  #DEBUG
  plottitle <- "Upper Confidence Bound"
  p2<-ggplot(d2, aes(x = X1, y = X2, fill = value)) + 
    geom_tile(color='black', width=1, height=1) +
    theme_bw() +
    xlim(0.5,11.5) +
    ylim(0.5,11.5) + 
    coord_equal() +
    ggtitle(plottitle) +
    #scale_fill_gradientn(name='P(choice)',colours = hm.palette(100),values = seq(0, 1, length=9)) +
    scale_fill_distiller(palette = "Spectral", na.value = 'white' )+
    labs(fill="UCB(x)")+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border=element_blank(), 
          plot.margin=unit(c(0,0,0,0),"mm"))
    #p2 <- p2 + annotate("text", x =maxValue$X1, y = maxValue$X2, label = "X") 
  d3 <- melt(matrix(p, nrow=11, ncol=11))
  names(d3) <- c('X1', 'X2', 'value')
  plotTitle <- bquote(paste(" Softmax Choice Probabilities(", tau, "=", .(roundedTau),")"))
  #DEBUG
  plotTitle <- "Choice Probability"
  p3 <- ggplot(d3, aes(x = X1, y = X2, fill = value)) + 
    geom_tile(color='black', width=1, height=1) +
    theme_bw() +
    xlim(0.5,11.5) +
    ylim(0.5,11.5) + 
    coord_equal() +
    ggtitle(plotTitle) +
    #scale_fill_gradientn(name='P(choice)',colours = hm.palette(100),values = seq(0, 1, length=9)) +
    scale_fill_distiller(palette = "Spectral", na.value = 'white' )+
    labs(fill="P(x)")+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border=element_blank(), 
          plot.margin=unit(c(0,0,0,0),"mm"))
  
  return(multiplot(p1, p2,p3, cols=3))
  #return(grid.arrange(p1,p2,p3,cols=3))
}

#corresponding to screenshots from experiment
params <- read.csv("rationalModels/parameters/gp.csv") #Load participant parameter estimates
lambdaMed <- median(params$lambda)
betaMed <- median(params$beta)
tauMed <- median(params$tau)


x <- c(2, 7, 5, 3, 7, 9, 10, 7, 2,1)
y <- c(6, 2, 3, 2, 6, 5, 3, 8, 8,10)
z <- c(43, 29, 24, 15, 18, 31, 26, 20, 75,53)
Figure2 <- customPredictions(as.matrix(cbind(x,y)),z, kernel=rbf, kernelName="RBF", acq=ucb,acqName="UCB", parVec=c(lambdaMed,lambdaMed, 1, .0001), beta=betaMed, tau=tauMed)

ggsave("plots/modelOverview.pdf", Figure2, width = 7.64, height = 4.72, units = 'in')

