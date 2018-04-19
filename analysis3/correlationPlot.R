#Correlation plot of Gaussian process priors
#Charley Wu, March 2018

#house keeping
rm(list=ls())
theme_set(theme_bw(base_size=16))# use the b&w theme

#load packages
packages <- c('plyr', 'ggplot2', 'MASS', 'psych')
lapply(packages, require, character.only = TRUE)

#source models
source('models.R')

################################################################################################################
#Begin Simulation
################################################################################################################

reps <- 10000 #number of environments to generate for each length scale
kernelFun <- rbf #choose kernel function

#pre-compute all pairwise distances between points in a 11 x 11 matrix
grid <- as.matrix(expand.grid(1:11,1:11))
pairwiseDistances <- apply(grid, 1, function(i) apply(grid, 1, function(j) dist(rbind(i, j))))
uniqueDistances <- unique(c(pairwiseDistances))

#create data frame
corData <- data.frame(lambda = numeric(), correlation = numeric(), distance = numeric())

#Loop over different length scales
for (lambda in c(0.5, 1, 2)){
  #generate samples from GP-rbf priors
  sampledEnvs <- mvrnorm(n = reps, mu = rep(0,121), kernelFun(as.matrix(expand.grid(1:11,1:11)),as.matrix(expand.grid(1:11,1:11)),c(lambda,lambda,1,.0001)))
  #calculate the strength of correlations for each distance
  correlations <- c() #placeholder
  for (d in uniqueDistances){
    pairs <- which(pairwiseDistances==d, arr.ind=T) # select all pairs of points where the pairwise distance is equal to d
    valuesA <- sampledEnvs[,c(pairs[,1])] #rewards
    valuesB <- sampledEnvs[,c(pairs[,2])]
    correlations <- rbind(correlations, cor(c(valuesA),c(valuesB)))
  }
  newDat <- data.frame(lambda = rep(lambda,length(uniqueDistances)), correlation= correlations, distance = uniqueDistances)
  corData <- rbind(corData, newDat)
}

corData$Lambda <- factor(corData$lambda)

################################################################################################################
# Plot
################################################################################################################

 p <- ggplot(corData, aes(x=distance, y = correlation, color = Lambda, shape = Lambda, linetype = Lambda)) + 
  geom_line(size = 1) +
  #geom_point(size = 2) +
  ylab("Pearson Correlation")+xlab("Distance")+
  theme_classic()+
  xlim(c(0,10))+
  scale_color_brewer(palette="Set1", direction=-1)+
  #scale_color_manual(values=c("#E69F00", "#56B4E9", "#009E73","#F0E442", "#0072B2", "#D55E00", "#CC79A7"))+
  theme(text = element_text(size=14,  family="sans"), legend.justification = c(1, 1), legend.position = c(1, 1))+
  theme(strip.background=element_blank(), legend.key=element_rect(color=NA)) +
  ggtitle("Decay of Correlation")
p

ggsave(filename = 'plots/kernelCorrelation.pdf', p, height = 4, width = 5)
