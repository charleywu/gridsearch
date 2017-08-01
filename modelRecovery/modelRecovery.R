#Analysis of Model Recovery Results
#Charley Wu, 2017

#house keeping
rm(list=ls())
theme_set(theme_classic(base_size=16))

#load packages
packages <- c('plyr', 'jsonlite', 'ggplot2', 'gridExtra', 'reshape2', 'stargazer', 'coefplot', "grid", 'parallel', 'corrplot', 'rprojroot', 'MASS')
lapply(packages, require, character.only = TRUE)

rootDir <- rprojroot::find_rstudio_root_file()
setwd(rootDir)
########################################################################################################################
#Load dataframes
########################################################################################################################
#1D
setwd("analysis1D/modelResults/")

#Function learning model
FLmodelFit <- read.csv('GPUCBrecovery.csv')
FLparams <- read.csv('GPUCBparamEstimates.csv')

#Associative Learning model
ALmodelFit <- read.csv('BMTUCBrecovery.csv')
ALparams <- read.csv('BMTUCBparamEstimates.csv')

setwd(".."); setwd("..")

generatingParams1D <- read.csv("analysis1D/rationalModels/parameters/gp.csv")
#2D
setwd("analysis2D/modelResults/")

#Function learning model
localFLmodelFit <- read.csv('localGPUCBrecovery.csv')
localFLparams <- read.csv('localGPUCBparamEstimates.csv')

#Associative Learning model
localALmodelFit <- read.csv('localBMTUCBrecovery.csv')
localALparams <- read.csv('localBMTUCBparamEstimates.csv')

setwd(".."); setwd("..")

generatingParams2D <- read.csv("analysis2D/rationalModels/parameters/localgp.csv")
setwd("paper")

########################################################################################################################
#Join together dataframes and prepare for plotting
########################################################################################################################

#Add generating model
FLmodelFit$generatingModel <- "Function Learning"
ALmodelFit$generatingModel <- "Associative Learning"
localFLmodelFit$generatingModel <- "Function Learning"
localALmodelFit$generatingModel <- "Associative Learning"

#Add experiment
FLmodelFit$experiment <- "Experiment 1"
ALmodelFit$experiment <- "Experiment 1"
localFLmodelFit$experiment <- "Experiment 2"
localALmodelFit$experiment <- "Experiment 2"

modelRecovery <- Reduce(function(x, y) merge(x, y, all=TRUE), list(FLmodelFit, ALmodelFit, localFLmodelFit,localALmodelFit ))

#Factor
modelRecovery$generatingModel <- factor(modelRecovery$generatingModel)
modelRecovery$experiment <- factor(modelRecovery$experiment)

########################################################################################################################
#Predictive Accuracy plot
########################################################################################################################

p1 <- ggplot(modelRecovery, aes(x = ModelName, y = R2, fill=ModelName)) +
  stat_summary(fun.y = mean, geom = "bar", position = "dodge", color='black', width = 0.8) + 
  stat_summary(fun.data = mean_se, geom = "errorbar", position = position_dodge(width = 0.90), width = 0.2 ) +
  facet_grid(~generatingModel, scales="free_x", space="free_x") +
  scale_fill_manual(values = c(  "#F0E442", "#009E73", "#E69F00",  "#56B4E9")) +
  xlab("") + ylab("Predictive Accuracy")+ 
  ylim(c(0,0.5))+
  facet_wrap(experiment~generatingModel, nrow=1 , scales = 'free_x')+
  theme(text = element_text(size=16,  family="serif"), strip.background=element_blank(), legend.key=element_rect(color=NA), legend.position="right")
p1

ggsave(filename = "recovery.pdf", p1,height = 6, width = 12, units ="in")

########################################################################################################################
#Parameter Correlation plot
########################################################################################################################
source('multiplot.R')
#1D-------------------
generatingParams <- ddply(generatingParams1D, ~participant, summarize, lambda1=median(lambda), beta1 = median(beta), tau1 = median(tau))
recoveredParams <- ddply(subset(FLparams, kernel=="RBF" &  acq == "UCB"), ~participant, summarize, lambda2 = median(lambda), beta2= median(beta), tau2=median(tau))
cor1D <- merge(generatingParams, recoveredParams, by="participant")

#Lambda
lambdaCor <- cor.test(cor1D$lambda1, cor1D$lambda2) #correlation test
cor1D$lambdadensity <- fields::interp.surface(MASS::kde2d(cor1D$lambda1, cor1D$lambda2), cor1D[,c("lambda1", "lambda2")]) #bivariate point density
p2a <- ggplot(cor1D, aes(lambda1, lambda2, fill=lambdadensity, alpha = 1/lambdadensity)) +
  geom_point(pch=21, size =3, color='black')+
  geom_smooth(method='lm', color='black', se=FALSE, linetype=2, fullrange=TRUE)+
  scale_fill_gradient(low = "#0091ff", high = "#f0650e") +
  scale_alpha(range = c(.4, 1)) + guides(alpha="none", fill="none") +
  expand_limits(x = 0, y = 0) +
  xlab(expression(lambda[1]))+
  ylab(expression(lambda[2])) +
  annotate("text", x = 1, y = 6, label = "r == 0.62 * ~~'p<.001'", parse=TRUE, size=6) +
  theme(title = element_text(family = 'serif'))
p2a

#beta
betaCor <- cor.test(cor1D$beta1, cor1D$beta2) #correlation test
cor1D$betadensity <- fields::interp.surface(MASS::kde2d(cor1D$beta1, cor1D$beta2), cor1D[,c("beta1", "beta2")]) #bivariate point density
p2b <- ggplot(cor1D, aes(beta1, beta2, fill=betadensity, alpha = 1/betadensity)) +
  geom_point(pch=21, size =3, color='black')+
  geom_smooth(method='lm', color='black', se=FALSE, linetype=2, fullrange=TRUE)+
  scale_fill_gradient(low = "#0091ff", high = "#f0650e") +
  scale_alpha(range = c(.4, 1)) + guides(alpha="none", fill="none") +
  expand_limits(x = 0, y = 0) +
  xlab(expression(beta[1]))+
  ylab(expression(beta[2])) +
  annotate("text", x = 1, y = 0.4, label = "r == 0.62 * ~~'p<.001'", parse=TRUE, size=6)+
  theme(title = element_text(family = 'serif'))
p2b
  
#tau
tauCor <- cor.test(cor1D$tau1, cor1D$tau2) #correlation test
cor1D$taudensity <- fields::interp.surface(MASS::kde2d(cor1D$tau1, cor1D$tau2), cor1D[,c("tau1", "tau2")]) #bivariate point density
p2c <- ggplot(cor1D, aes(tau1, tau2, fill=taudensity, alpha = 1/taudensity)) +
  geom_point(pch=21, size =3, color='black')+
  geom_smooth(method='lm', color='black', se=FALSE, linetype=2, fullrange=TRUE)+
  scale_fill_gradient(low = "#0091ff", high = "#f0650e") +
  scale_alpha(range = c(.4, 1)) + guides(alpha="none", fill="none") +
  expand_limits(x = 0, y = 0) +
  xlab(expression(tau[1]))+
  ylab(expression(tau[2]))+
  annotate("text", x = 0.1, y = 0.5, label = "r == 0.91 * ~~'p<.001'", parse=TRUE, size=6)+
  theme(title = element_text(family = 'serif'))
p2c

p2 <- multiplot(p2a,p2b, p2c, cols=3)
ggsave(filename='paramRecovery1D.pdf', p2, height = 6, width = 12, units='in')

#2D -----------------
generatingParams <- ddply(generatingParams2D, ~participant, summarize, lambda1=median(lambda), beta1 = median(beta), tau1 = median(tau))
recoveredParams <- ddply(subset(localFLparams, kernel=="LRBF" &  acq == "UCB"), ~participant, summarize, lambda2 = median(lambda), beta2= median(beta), tau2=median(tau))
cor2D <- merge(generatingParams, recoveredParams, by="participant")

#Lambda
lambdaCor <- cor.test(cor2D$lambda1, cor2D$lambda2) #correlation test
cor2D$lambdadensity <- fields::interp.surface(MASS::kde2d(cor2D$lambda1, cor2D$lambda2), cor2D[,c("lambda1", "lambda2")]) #bivariate point density
p3a <- ggplot(cor2D, aes(lambda1, lambda2, fill=lambdadensity, alpha = 1/lambdadensity)) +
  geom_point(pch=21, size =3, color='black')+
  geom_smooth(method='lm', color='black', se=FALSE, linetype=2, fullrange=TRUE)+
  scale_fill_gradient(low = "#0091ff", high = "#f0650e") +
  scale_alpha(range = c(.4, 1)) + guides(alpha="none", fill="none") +
  expand_limits(x = 0, y = 0) +
  xlab(expression(lambda[1]))+
  ylab(expression(lambda[2])) +
  annotate("text", x = 1.5, y = 3, label = "r == 0.91 * ~~'p<.001'", parse=TRUE, size=6)+
  theme(title = element_text(family = 'serif'))
p3a

#beta
betaCor <- cor.test(cor2D$beta1, cor2D$beta2) #correlation test
cor2D$betadensity <- fields::interp.surface(MASS::kde2d(cor2D$beta1, cor2D$beta2), cor2D[,c("beta1", "beta2")]) #bivariate point density
p3b <- ggplot(cor2D, aes(beta1, beta2, fill=betadensity, alpha = 1/betadensity)) +
  geom_point(pch=21, size =3, color='black')+
  geom_smooth(method='lm', color='black', se=FALSE, linetype=2, fullrange=TRUE)+
  scale_fill_gradient(low = "#0091ff", high = "#f0650e") +
  scale_alpha(range = c(.4, 1)) + guides(alpha="none", fill="none") +
  expand_limits(x = 0, y = 0) +
  xlab(expression(beta[1]))+
  ylab(expression(beta[2])) +
  annotate("text", x = 0.4, y = 1.2, label = "r == 0.77 * ~~'p<.001'", parse=TRUE, size=6)+
  theme(title = element_text(family = 'serif'))
p3b

#tau
tauCor <- cor.test(cor2D$tau1, cor2D$tau2) #correlation test
cor2D$taudensity <- fields::interp.surface(MASS::kde2d(cor2D$tau1, cor2D$tau2), cor2D[,c("tau1", "tau2")]) #bivariate point density
p3c <- ggplot(cor2D, aes(tau1, tau2, fill=taudensity, alpha = 1/taudensity)) +
  geom_point(pch=21, size =3, color='black')+
  geom_smooth(method='lm', color='black', se=FALSE, linetype=2, fullrange=TRUE)+
  scale_fill_gradient(low = "#0091ff", high = "#f0650e") +
  scale_alpha(range = c(.4, 1)) + guides(alpha="none", fill="none") +
  expand_limits(x = 0, y = 0) +
  xlab(expression(tau[1]))+
  ylab(expression(tau[2]))+
  annotate("text", x = 0.05, y = 0.25, label = "r == 0.76 * ~~'p<.001'", parse=TRUE, size=6)+
  theme(title = element_text(family = 'serif'))
p3c

p3 <- multiplot(p3a,p3b, p3c, cols=3)
ggsave(filename='paramRecovery2D.pdf', p3, height = 6, width = 12, units='in')


