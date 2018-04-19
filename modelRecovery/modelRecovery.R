#Analysis of Model Recovery Results
#Charley Wu, 2018

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
localFLmodelFit2 <- read.csv('localGPUCBrecovery.csv')
localFLparams2 <- read.csv('localGPUCBparamEstimates.csv')

#Associative Learning model
localALmodelFit2 <- read.csv('localBMTUCBrecovery.csv')
localALparams2 <- read.csv('localBMTUCBparamEstimates.csv')

setwd(".."); setwd("..")

generatingParams2D <- read.csv("analysis2D/rationalModels/parameters/localgp.csv")

#Experiment 3
setwd("analysis3/modelResults/")

#Function learning model
localFLmodelFit3 <- read.csv('localGPUCBrecovery.csv')
localFLparams3 <- read.csv('localGPUCBparamEstimates.csv')

#Associative Learning model
localALmodelFit3 <- read.csv('localBMTUCBrecovery.csv')
localALparams3 <- read.csv('localBMTUCBparamEstimates.csv')

setwd(".."); setwd("..")

generatingParams3 <- read.csv("analysis3/rationalModels/parameters/localgp.csv")
setwd("paper")

########################################################################################################################
#Join together dataframes and prepare for plotting
########################################################################################################################

#Add generating model
FLmodelFit$generatingModel <- "Function Learning"
ALmodelFit$generatingModel <- "Option Learning"
localFLmodelFit2$generatingModel <- "Function Learning"
localALmodelFit2$generatingModel <- "Option Learning"
localFLmodelFit3$generatingModel <- "Function Learning"
localALmodelFit3$generatingModel <- "Option Learning"

#Add experiment
FLmodelFit$experiment <- "Experiment 1"
ALmodelFit$experiment <- "Experiment 1"
localFLmodelFit2$experiment <- "Experiment 2"
localALmodelFit2$experiment <- "Experiment 2"
localFLmodelFit3$experiment <- "Experiment 3"
localALmodelFit3$experiment <- "Experiment 3"

modelRecovery <- Reduce(function(x, y) merge(x, y, all=TRUE), list(FLmodelFit, ALmodelFit, localFLmodelFit2,localALmodelFit2,localFLmodelFit3,localALmodelFit3 ))

#Factor
modelRecovery$generatingModel <- factor(modelRecovery$generatingModel)
modelRecovery$experiment <- factor(modelRecovery$experiment)

#report means

modelRecovery %>% group_by(experiment, generatingModel, kernel) %>% summarize(avgR2= mean(R2))

########################################################################################################################
#Predictive Accuracy plot
########################################################################################################################

se<-function(x){sd(x)/sqrt(length(x))}
mtmDF <- ddply(modelRecovery, ~kernel+generatingModel+experiment+participant, summarise, R2 =mean(R2))
boxplotDF <- ddply(modelRecovery,  ~kernel+generatingModel+experiment, summarise,d_ymin = max(min(R2), quantile(R2, 0.25) - 1.5 * IQR(R2)), d_ymax = min(max(R2), quantile(R2, 0.75) + 1.5 * IQR(R2)),
                   d_lower = quantile(R2, 0.25),  d_middle = median(R2), d_upper = quantile(R2, 0.75),
                   mu=mean(R2))

set.seed(2018)
jitterVal <- runif(81, max = 0.3)

p1 <- ggplot(boxplotDF) +
  geom_line(data=mtmDF, aes(x = as.numeric(kernel) + 0.1 + jitterVal, group=participant, y=R2), color = 'black', alpha = 0.1)+
  #stat_summary(fun.y = mean, geom = "bar", position = "dodge", color='black') + 
  geom_boxplot(aes(x =as.numeric(kernel)-0.2, ymin = d_lower, ymax = d_upper, lower = d_lower, middle = d_middle, upper = d_upper, width = 2 * 0.2, fill = kernel), stat = "identity", color='black') +
  #whiskers
  geom_segment(aes(x = as.numeric(kernel), y = d_ymin, xend = as.numeric(kernel), yend = d_ymax)) +
  geom_segment(aes(x = as.numeric(kernel) - 0.1,  y = d_ymax - 0.0025, xend = as.numeric(kernel),  yend = d_ymax -0.0025)) + #top
  geom_segment(aes(x = as.numeric(kernel) - 0.1, y = d_ymin +0.0025, xend = as.numeric(kernel), yend = d_ymin + 0.0025)) + #bottom
  geom_point(aes(x = as.numeric(kernel)-0.2, y = mu), size = 1.3, shape=23, fill='white') +
  #geom_jitter(data=mtmDF, aes(x = as.numeric(kernel)+.2,  y = R2,  color = kernel), width = 0.2 - 0.25 * 0.2, height = 0, size=1, alpha = 0.8)+
  geom_point(data=mtmDF, aes(x = as.numeric(kernel)+0.1+jitterVal,  y = R2,  color = kernel), size=1, alpha = 0.8)+
  facet_grid(experiment~generatingModel, scales="free", space="free_x")+
  #geom_jitter(data = mainTextModels, aes(x=acq, y = R2, fill= kernel),  color='grey', size = 1, shape=21, alpha=0.2, position=position_jitterdodge(dodge.width=0.9, jitter.width = 0.2))+
  #geom_errorbar(aes(ymin=newR2 - se, ymax=newR2 + se),color='black', width = .4, position=position_dodge((width=0.9))) +
  #scale_fill_manual(values = c("#7F0000", "#00BCE2",  "#37045F" ))+
  scale_fill_manual(values = c( "#F0E442",  "#009E73", "#E69F00", "#56B4E9"))+
  scale_color_manual(values = c( "#F0E442",  "#009E73", "#E69F00", "#56B4E9"))+
  ylab("Predictive Accuracy")+ 
  coord_cartesian(ylim=c(-0.1, 1))+
  theme_classic()+
  theme(text = element_text(size=12,  family="sans"),axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), 
        strip.background=element_blank(),
        #strip.text = element_blank(),
        legend.key=element_rect(color=NA),
        panel.spacing.x=unit(0.2, "lines"),
        panel.spacing.y=unit(1, "lines"),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        legend.position="right")
p1

ggsave(filename = "recovery.pdf", p1,height = 5, width = 7, units ="in")

########################################################################################################################
#Parameter Correlation plot
########################################################################################################################
source('multiplot.R')
axisScale <- function(vec1, vec2){ #Finds the larger of the two axis limits that omits outliers more than 1.5IQR from the 75th percentile
  max1 <- min(max(vec1), median(vec1, 0.75) + 1.5 * IQR(vec1))
  max2 <- min(max(vec2), quantile(vec2, 0.75) + 1.5 * IQR(vec2))
  return(max(max1, max2)+0.1)
}

#1D-------------------
generatingParams <- ddply(generatingParams1D, ~participant, summarize, lambda1=median(lambda), beta1 = median(beta), tau1 = median(tau))
recoveredParams <- ddply(subset(FLparams, kernel=="RBF" &  acq == "UCB"), ~participant, summarize, lambda2 = median(lambda), beta2= median(beta), tau2=median(tau))
cor1D <- merge(generatingParams, recoveredParams, by="participant")

#Lambda
lambdaCor <- cor.test(cor1D$lambda1, cor1D$lambda2, method='kendall') #correlation test
cor1D$lambdadensity <- fields::interp.surface(MASS::kde2d(cor1D$lambda1, cor1D$lambda2), cor1D[,c("lambda1", "lambda2")]) #bivariate point density
p2a <- ggplot(cor1D, aes(lambda1, lambda2, fill=lambdadensity, alpha = 1/lambdadensity)) +
  geom_smooth(data=subset(cor1D, lambda1 <= axisScale(cor1D$lambda1, cor1D$lambda2) & lambda2 <= axisScale(cor1D$lambda1, cor1D$lambda2)), method='lm', color='black', se=FALSE, linetype=2, fullrange=TRUE)+
  geom_point(pch=21, size =3, color='black')+
  scale_fill_gradient(low = "#0091ff", high = "#f0650e") +
  scale_alpha(range = c(.4, 1)) + guides(alpha="none", fill="none") +
  expand_limits(x = 0, y = 0) +
  xlab(expression(lambda[1]))+
  ylab(expression(lambda[2])) +
  coord_cartesian(xlim=c(0,axisScale(cor1D$lambda1, cor1D$lambda2)), ylim=c(0,axisScale(cor1D$lambda1, cor1D$lambda2)))+
  annotate("text", x = 0.5, y = 1.5, label = "r[tau] == '.66;' * ~~'p<.001'", parse=TRUE, size=6) +
  theme(title = element_text(family = 'serif'))
#p2a <- ggExtra::ggMarginal(p2a, type = "boxplot", fill='transparent', size = 10,  outlier.shape=NA)
p2a

#beta
betaCor <- cor.test(cor1D$beta1, cor1D$beta2, method='kendall') #correlation test
cor1D$betadensity <- fields::interp.surface(MASS::kde2d(cor1D$beta1, cor1D$beta2), cor1D[,c("beta1", "beta2")]) #bivariate point density
p2b <- ggplot(cor1D, aes(beta1, beta2, fill=betadensity, alpha = 1/betadensity)) +
  geom_point(pch=21, size =3, color='black')+
  geom_smooth(data=subset(cor1D, beta1 <= axisScale(cor1D$beta1, cor1D$beta2) & beta2 <= axisScale(cor1D$beta1, cor1D$beta2)), method='lm', color='black', se=FALSE, linetype=2, fullrange=TRUE)+
  scale_fill_gradient(low = "#0091ff", high = "#f0650e") +
  scale_alpha(range = c(.4, 1)) + guides(alpha="none", fill="none") +
  expand_limits(x = 0, y = 0) +
  xlab(expression(beta[1]))+
  ylab(expression(beta[2])) +
  coord_cartesian(xlim=c(0,axisScale(cor1D$beta1, cor1D$beta2)), ylim=c(0,axisScale(cor1D$beta1, cor1D$beta2)))+
  annotate("text", x = .25, y = 0.6, label = "r[tau] == '.30;' * ~~'p<.001'", parse=TRUE, size=6)+
  theme(title = element_text(family = 'serif'))
#p2b <- ggExtra::ggMarginal(p2b, type = "boxplot", fill='transparent', size = 10, outlier.shape=NA)
p2b
  
#tau
tauCor <- cor.test(cor1D$tau1, cor1D$tau2, method='kendall') #correlation test
cor1D$taudensity <- fields::interp.surface(MASS::kde2d(cor1D$tau1, cor1D$tau2), cor1D[,c("tau1", "tau2")]) #bivariate point density
p2c <- ggplot(cor1D, aes(tau1, tau2, fill=taudensity, alpha = 1/taudensity)) +
  geom_point(pch=21, size =3, color='black')+
  geom_smooth(data=subset(cor1D, tau1 <= axisScale(cor1D$tau1, cor1D$tau2) & tau2 <= axisScale(cor1D$tau1, cor1D$tau2)), method='lm', color='black', se=FALSE, linetype=2, fullrange=TRUE)+
  scale_fill_gradient(low = "#0091ff", high = "#f0650e") +
  scale_alpha(range = c(.4, 1)) + guides(alpha="none", fill="none") +
  expand_limits(x = 0, y = 0) +
  xlab(expression(tau[1]))+
  ylab(expression(tau[2]))+
  coord_cartesian(xlim=c(0,axisScale(cor1D$tau1, cor1D$tau2)), ylim=c(0,axisScale(cor1D$tau1, cor1D$tau2)))+
  annotate("text", x = 0.05, y = 0.17, label = "r[tau] == '.54;' * ~~'p<.001'", parse=TRUE, size=6)+
  theme(title = element_text(family = 'serif'))
#p2c <- ggExtra::ggMarginal(p2c, type = "boxplot", fill='transparent', size = 10, outlier.shape=NA)
p2c

p2 <- plot_grid(p2a,p2b, p2c, labels=NA, ncol=3)
ggsave(filename='paramRecovery1D.pdf', p2, height = 4, width = 12, units='in')

#2D -----------------
generatingParams <- ddply(generatingParams2D, ~participant, summarize, lambda1=median(lambda), beta1 = median(beta), tau1 = median(tau))
recoveredParams <- ddply(subset(localFLparams2, kernel=="LRBF" &  acq == "UCB"), ~participant, summarize, lambda2 = median(lambda), beta2= median(beta), tau2=median(tau))
cor2D <- merge(generatingParams, recoveredParams, by="participant")

#Lambda
lambdaCor <- cor.test(cor2D$lambda1, cor2D$lambda2, method='kendall') #correlation test
cor2D$lambdadensity <- fields::interp.surface(MASS::kde2d(cor2D$lambda1, cor2D$lambda2), cor2D[,c("lambda1", "lambda2")]) #bivariate point density
p3a <- ggplot(cor2D, aes(lambda1, lambda2, fill=lambdadensity, alpha = 1/lambdadensity)) +
  geom_point(pch=21, size =3, color='black')+
  geom_smooth(data=subset(cor2D, lambda1 <= axisScale(cor2D$lambda1, cor2D$lambda2) & lambda2 <= axisScale(cor2D$lambda1, cor2D$lambda2)), method='lm', color='black', se=FALSE, linetype=2, fullrange=TRUE)+
  scale_fill_gradient(low = "#0091ff", high = "#f0650e") +
  scale_alpha(range = c(.4, 1)) + guides(alpha="none", fill="none") +
  expand_limits(x = 0, y = 0) +
  xlab(expression(lambda[1]))+
  ylab(expression(lambda[2])) +
  coord_cartesian(xlim=c(0,axisScale(cor2D$lambda1, cor2D$lambda2)), ylim=c(0,axisScale(cor2D$lambda1, cor2D$lambda2)))+
  annotate("text", x = .5, y = 1.2, label = "r[tau] == '.77;' * ~~'p<.001'", parse=TRUE, size=6)+
  theme(title = element_text(family = 'serif'))
p3a

#beta
betaCor <- cor.test(cor2D$beta1, cor2D$beta2, method='kendall') #correlation test
cor2D$betadensity <- fields::interp.surface(MASS::kde2d(cor2D$beta1, cor2D$beta2), cor2D[,c("beta1", "beta2")]) #bivariate point density
p3b <- ggplot(cor2D, aes(beta1, beta2, fill=betadensity, alpha = 1/betadensity)) +
  geom_point(pch=21, size =3, color='black')+
  geom_smooth(data=subset(cor2D, beta1 <= axisScale(cor2D$beta1, cor2D$beta2) & beta2 <= axisScale(cor2D$beta1, cor2D$beta2)), method='lm', color='black', se=FALSE, linetype=2, fullrange=TRUE)+
  scale_fill_gradient(low = "#0091ff", high = "#f0650e") +
  scale_alpha(range = c(.4, 1)) + guides(alpha="none", fill="none") +
  scale_x_continuous(expand=c(0,0), limits=c(0,axisScale(cor2D$beta1, cor2D$beta2))) +
  scale_y_continuous(expand=c(0,0), limits=c(0,axisScale(cor2D$beta1, cor2D$beta2))) +
  xlab(expression(beta[1]))+
  ylab(expression(beta[2])) +
  coord_cartesian(xlim=c(0,axisScale(cor2D$beta1, cor2D$beta2)), ylim=c(0,axisScale(cor2D$beta1, cor2D$beta2)))+
  annotate("text", x = 0.3, y = .9, label = "r[tau] == '.59;' * ~~'p<.001'", parse=TRUE, size=6)+
  theme(title = element_text(family = 'serif'))
p3b

#tau
tauCor <- cor.test(cor2D$tau1, cor2D$tau2, method='kendall') #correlation test
cor2D$taudensity <- fields::interp.surface(MASS::kde2d(cor2D$tau1, cor2D$tau2), cor2D[,c("tau1", "tau2")]) #bivariate point density
p3c <- ggplot(cor2D, aes(tau1, tau2, fill=taudensity, alpha = 1/taudensity)) +
  geom_point(pch=21, size =3, color='black')+
  geom_smooth(data=subset(cor2D, tau1 <= axisScale(cor2D$tau1, cor2D$tau2) & tau2 <= axisScale(cor2D$tau1, cor2D$tau2)), method='lm', color='black', se=FALSE, linetype=2, fullrange=TRUE)+
  scale_fill_gradient(low = "#0091ff", high = "#f0650e") +
  scale_alpha(range = c(.4, 1)) + guides(alpha="none", fill="none") +
  scale_x_continuous(expand=c(0,0), limits=c(0,axisScale(cor2D$tau1, cor2D$tau2))) +
  scale_y_continuous(expand=c(0,0), limits=c(0,axisScale(cor2D$tau1, cor2D$tau2))) +
  coord_cartesian(xlim=c(0,axisScale(cor2D$tau1, cor2D$tau2)), ylim=c(0,axisScale(cor2D$tau1, cor2D$tau2)))+
  xlab(expression(tau[1]))+
  ylab(expression(tau[2]))+
  coord_cartesian(xlim=c(0,axisScale(cor2D$tau1, cor2D$tau2)), ylim=c(0,axisScale(cor2D$tau1, cor2D$tau2)))+
  annotate("text", x = 0.1, y = 0.2, label = "r[tau] == '.61;' * ~~'p<.001'", parse=TRUE, size=6)+
  theme(title = element_text(family = 'serif'))
p3c

p3 <- plot_grid(p3a,p3b, p3c, labels=NA, ncol=3)
ggsave(filename='paramRecovery2D.pdf', p3, height = 4, width = 12, units='in')

#Experiment3
generatingParams <- ddply(generatingParams3, ~participant, summarize, lambda1=median(lambda), beta1 = median(beta), tau1 = median(tau))
recoveredParams <- ddply(subset(localFLparams3, kernel=="LRBF" &  acq == "UCB"), ~participant, summarize, lambda2 = median(lambda), beta2= median(beta), tau2=median(tau))
cor3 <- merge(generatingParams, recoveredParams, by="participant")

#Lambda
lambdaCor <- cor.test(cor3$lambda1, cor3$lambda2, method='kendall') #correlation test
cor3$lambdadensity <- fields::interp.surface(MASS::kde2d(cor3$lambda1, cor3$lambda2), cor3[,c("lambda1", "lambda2")]) #bivariate point density
p4a <- ggplot(cor3, aes(lambda1, lambda2, fill=lambdadensity, alpha = 1/lambdadensity)) +
  geom_point(pch=21, size =3, color='black')+
  geom_smooth(data=subset(cor3, lambda1 <= axisScale(cor3$lambda1, cor3$lambda2) & lambda2 <= axisScale(cor3$lambda1, cor3$lambda2)), method='lm', color='black', se=FALSE, linetype=2, fullrange=TRUE)+
  scale_fill_gradient(low = "#0091ff", high = "#f0650e") +
  scale_alpha(range = c(.4, 1)) + guides(alpha="none", fill="none") +
  expand_limits(x = 0, y = 0) +
  xlab(expression(lambda[1]))+
  ylab(expression(lambda[2])) +
  coord_cartesian(xlim=c(0,axisScale(cor3$lambda1, cor3$lambda2)), ylim=c(0,axisScale(cor3$lambda1, cor3$lambda2)))+
  annotate("text", x = .3, y = 1, label = "r[tau] == '.70;' * ~~'p<.001'", parse=TRUE, size=6)+
  theme(title = element_text(family = 'serif'))
p4a

#beta
betaCor <- cor.test(cor3$beta1, cor3$beta2, method='kendall') #correlation test
cor3$betadensity <- fields::interp.surface(MASS::kde2d(cor3$beta1, cor3$beta2), cor3[,c("beta1", "beta2")]) #bivariate point density
p4b <- ggplot(cor3, aes(beta1, beta2, fill=betadensity, alpha = 1/betadensity)) +
  geom_point(pch=21, size =3, color='black')+
  geom_smooth(data=subset(cor3, beta1 <= axisScale(cor3$beta1, cor3$beta2) & beta2 <= axisScale(cor3$beta1, cor3$beta2)), method='lm', color='black', se=FALSE, linetype=2, fullrange=TRUE)+
  scale_fill_gradient(low = "#0091ff", high = "#f0650e") +
  scale_alpha(range = c(.4, 1)) + guides(alpha="none", fill="none") +
  expand_limits(x = 0, y = 0) +
  xlab(expression(beta[1]))+
  ylab(expression(beta[2])) +
  coord_cartesian(xlim=c(0,axisScale(cor3$beta1, cor3$beta2)), ylim=c(0,axisScale(cor3$beta1, cor3$beta2)))+
  annotate("text", x = 0.4, y = 1.2, label = "r[tau] == '.76;' * ~~'p<.001'", parse=TRUE, size=6)+
  theme(title = element_text(family = 'serif'))
p4b

#tau
tauCor <- cor.test(cor3$tau1, cor3$tau2, method='kendall') #correlation test
cor3$taudensity <- fields::interp.surface(MASS::kde2d(cor3$tau1, cor3$tau2), cor3[,c("tau1", "tau2")]) #bivariate point density
p4c <- ggplot(cor3, aes(tau1, tau2, fill=taudensity, alpha = 1/taudensity)) +
  geom_point(pch=21, size =3, color='black')+
  geom_smooth(data=subset(cor3, tau1 <= axisScale(cor3$tau1, cor3$tau2) & tau2 <= axisScale(cor3$tau1, cor3$tau2)), method='lm', color='black', se=FALSE, linetype=2, fullrange=TRUE)+
  scale_fill_gradient(low = "#0091ff", high = "#f0650e") +
  scale_alpha(range = c(.4, 1)) + guides(alpha="none", fill="none") +
  expand_limits(x = 0, y = 0) +
  xlab(expression(tau[1]))+
  ylab(expression(tau[2]))+
  theme(title = element_text(family = 'serif'))+
  coord_cartesian(xlim=c(0,axisScale(cor3$tau1, cor3$tau2)), ylim=c(0,axisScale(cor3$tau1, cor3$tau2)))+
  annotate("text", x = 0.1, y = 0.35, label = "r[tau] =='.79;' * ~~'p<.001'", parse=TRUE, size=6)
p4c

p4 <- plot_grid(p4a,p4b, p4c, ncol=3, labels=NA)
ggsave(filename='paramRecovery3.pdf', p4, height = 4, width = 12, units='in')
