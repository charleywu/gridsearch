#Analysis of learning curves
#Charley Wu 2018

#house keeping
rm(list=ls())

#load packages
packages <- c('plyr', 'jsonlite', 'ggplot2', 'reshape2', "grid", 'matrixcalc', 'tools', 'data.table')
lapply(packages, require, character.only = TRUE)


##############################################################################################################
#Load simulation data
#saved version of dataframe #
##############################################################################################################

setwd('rationalModels/simulatedData/')

# Get the files names
files = list.files(pattern="*.csv")
df = do.call(rbind, lapply(files, fread))

setwd('..')
setwd('..')

df <- ddply(df, ~trial+scenario+horizon+Model, summarise, meanReward=mean(z), meanSE= sd(z)/sqrt(length(z)),  maxReward=mean(zmax), maxSE= sd(zmax)/sqrt(length(zmax)))
colnames(df) <- c("trial", "PayoffCondition", "Horizon", "Model", "meanReward", "meanSE", "maxReward", "maxSE")
df$Environment <- "Natural"
df$PayoffCondition <- factor(df$PayoffCondition)
levels(df$PayoffCondition) <- c("Cumulative", "Best")
df$Horizon <- factor(df$Horizon)
df$Model <- factor(df$Model)
df <- df[c("trial", "PayoffCondition", "Horizon",  "meanReward", "meanSE", "maxReward", "maxSE", "Model", "Environment")]

#simulatedData <- df
##############################################################################################################
#Add Human Data
##############################################################################################################
#add human data
source('dataMunging.R') #source data import function
d <- dataImport(normalize=FALSE)
humanData <- ddply(d, ~trial+scenario+horizon, summarise, meanReward=mean(z), meanSE= sd(z)/sqrt(length(z)),  maxReward=mean(zmax), maxSE= sd(zmax)/sqrt(length(zmax)))
humanData$Model <- rep("Human", nrow(humanData))
colnames(humanData) <- c("trial", "PayoffCondition", "Horizon",  "meanReward", "meanSE", "maxReward", "maxSE", "Model")
humanData$PayoffCondition <- factor(humanData$PayoffCondition)
levels(humanData$PayoffCondition) <- c("Cumulative", "Best")
humanData$Environment <- "Natural"
df <- rbind(df, humanData)

##############################################################################################################
#Add Random data
##############################################################################################################
randomDF <- read.csv("rationalModels/random.csv")
randomDF <- randomDF[ , (names(randomDF) != "X")]
randomDF$Horizon <- 40
randomDF$PayoffCondition <- "Cumulative"
randomDF2 <- randomDF
randomDF2$PayoffCondition <- "Best"
randomDF <- rbind(randomDF, randomDF2)
randomDF$Environment <- "Natural"
randomDF<- randomDF[c("trial", "PayoffCondition", "Horizon",  "meanReward", "meanSE", "maxReward", "maxSE", "Model", "Environment")] #double check this is correct
df<- rbind(df, randomDF)


##############################################################################################################
#Plot Data
##############################################################################################################
levels(df$Model) <- c("Option Learning", "Function Learning",  "Option Learning*", "Function Learning*","Human", "Random")
df$Model <- factor(df$Model, level=c("Random", "Option Learning", "Option Learning*","Function Learning", "Function Learning*",  "Human" ))
levels(df$Horizon) <- c("Short", "Long")
levels(df$PayoffCondition) <- c("Accumulation", "Maximization") 
#5 trial average
df$trial5 <- round((df$trial+1)/5)*5
df$trial5 <- ifelse(df$trial5<5,0,df$trial5)
plotdf <- ddply(df, ~trial5+Model+Environment+Horizon+PayoffCondition, summarise, meanReward=mean(meanReward), meanSE=sd(meanReward)/sqrt(length(meanReward)), maxReward=mean(maxReward), maxSE=sd(maxSE)/sqrt(length(maxSE)))

p1 <- ggplot(subset(plotdf, Model!="Local Search"), aes(x = trial5, y = meanReward, color = Model,fill=Model, linetype = Horizon))+
  stat_summary(fun.y=mean, geom='line', size=.7)+
  #geom_ribbon(aes(ymin=meanReward-meanSE, ymax=meanReward+meanSE),alpha=0.1, color=NA) +
  theme_classic()+
  xlab('Trial')+
  #coord_cartesian(ylim=c(45,85))+
  ylab('Average Reward')+
  scale_linetype_manual(values=c("dashed", "solid"))+
  scale_color_manual(values=c("black",  "#F0E442", "#E69F00", "#009E73", "#56B4E9", "#fb9a99"))+
  scale_fill_manual(values=c("black",  "#F0E442", "#E69F00", "#009E73", "#56B4E9", "#fb9a99"))+
  facet_grid(~PayoffCondition)+
  theme(legend.position='none',strip.background=element_blank(), legend.key=element_rect(color=NA))
p1
ggsave(filename = 'plots/separatedLearningCurve3Avg.pdf', p1,  height = 2.5, width = 4, unit='in')


p2 <- ggplot(subset(plotdf, Model!="Local Search"), aes(x = trial5, y = maxReward, color = Model,fill=Model, linetype = Horizon))+
  stat_summary(fun.y=mean, geom='line', size=.7)+
  #geom_ribbon(aes(ymin=maxReward-meanSE, ymax=maxReward+meanSE),alpha=0.1, color=NA) +
  theme_classic()+
  xlab('Trial')+
  #coord_cartesian(ylim=c(45,100))+
  ylab('Maximum Reward')+
  scale_linetype_manual(values=c("dashed", "solid"))+
  scale_color_manual(values=c("black",  "#F0E442", "#E69F00", "#009E73", "#56B4E9", "#fb9a99"))+
  scale_fill_manual(values=c("black",  "#F0E442", "#E69F00", "#009E73", "#56B4E9", "#fb9a99"))+
  facet_grid( ~PayoffCondition)+
  theme(legend.position='none',strip.background=element_blank(), legend.key=element_rect(color=NA))
p2
ggsave(filename = 'plots/separatedLearningCurve3Max.pdf', p2,  height = 2.5, width = 4, unit='in')


p3 <- ggplot(subset(plotdf, Model!="Local Search"), aes(x = trial5, y = meanReward, color = Model,fill=Model, shape=Model))+
  stat_summary(fun.y=mean, geom='line', size=.7)+
  stat_summary(fun.y=mean, geom='point', size=2)+
  #geom_ribbon(aes(ymin=meanReward-meanSE, ymax=meanReward+meanSE),alpha=0.1, color=NA) +
  theme_classic()+
  coord_cartesian(ylim=c(45,85))+
  xlab('Trial')+
  ylab('Average Reward')+
  scale_linetype_manual(values=c("dashed", "solid"))+
  scale_color_manual(values=c("black",  "#F0E442", "#E69F00", "#009E73", "#56B4E9", "#fb9a99"))+
  scale_fill_manual(values=c("black",  "#F0E442", "#E69F00", "#009E73", "#56B4E9", "#fb9a99"))+
  facet_grid(~Environment)+
  theme(strip.background=element_blank(), legend.key=element_rect(color=NA))
p3
ggsave(filename = 'plots/LearningCurve3Avg.pdf', p3, height = 2.5, width = 4.2, unit='in', useDingbats=FALSE)


################################

source('dataMunging.R') #source data import function
d <- dataImport(normalize=FALSE)
humanData <- ddply(d, ~id+trial+scenario+horizon, summarise, meanReward=mean(z), meanSE= sd(z)/sqrt(length(z)),  maxReward=mean(zmax), maxSE= sd(zmax)/sqrt(length(zmax)))
humanData$Model <- rep("Human", nrow(humanData))
colnames(humanData) <- c("id", "trial", "PayoffCondition", "Horizon",  "meanReward", "meanSE", "maxReward", "maxSE","Model")
humanData$PayoffCondition <- factor(humanData$PayoffCondition)
levels(humanData$PayoffCondition) <- c("Accumulation", "Maximization")


#5 trial average
humanData$trial<- humanData$trial - 1 #starting at trial 0
humanData$trial5 <- round((humanData$trial+1)/5)*5
humanData$trial5 <- ifelse(df$trial5<5,0,df$trial5)
plotdf <- ddply(humanData, ~trial5+id+Horizon+PayoffCondition, summarise, meanReward=mean(meanReward), meanSE=sd(meanReward)/sqrt(length(meanReward)), maxReward=mean(maxReward), maxSE=sd(maxSE)/sqrt(length(maxSE)))
plotdf$Horizon<- factor(plotdf$Horizon)
levels(plotdf$Horizon) <- c("Short", "Long")

p4<- ggplot(plotdf, aes(x=trial5, y = meanReward,group=interaction(id, Horizon),fill =Horizon, color=Horizon))+
  geom_line(alpha = 0.5)+
  #geom_ribbon(aes(ymin=meanReward-meanSE, ymax=meanReward+meanSE),alpha=0.05, color = NA) +
  #stat_summary(aes(group=Horizon), fun.y=mean, geom='line', size=.7)+
  theme_classic()+
  xlab('Trial')+
  ylab('Average Reward')+
  coord_cartesian(ylim=c(30,100))+
  scale_linetype_manual(values=c("dashed", "solid"))+
  scale_color_manual(values=c("#E69F00", "#56B4E9"))+
  scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
  facet_grid(PayoffCondition~.)+
  theme(legend.position='right',strip.background=element_blank(), legend.key=element_rect(color=NA))
p4

ggsave(filename = 'plots/indCurves3avg.pdf', p4, height = 4.3, width = 2.8, unit='in', useDingbats=FALSE)

p5<- ggplot(plotdf, aes(x=trial5, y = maxReward,group=interaction(id, Horizon),fill =Horizon, color=Horizon))+
  geom_line(alpha=.5)+
  #geom_ribbon(aes(ymin=meanReward-meanSE, ymax=meanReward+meanSE),alpha=0.05, color = NA) +
  #stat_summary(aes(group=Horizon), fun.y=mean, geom='line', size=.7)+
  theme_classic()+
  xlab('Trial')+
  ylab('Maximum Reward')+
  coord_cartesian(ylim=c(30,100))+
  #coord_cartesian(ylim=c(45,90))+
  scale_linetype_manual(values=c("dashed", "solid"))+
  scale_color_manual(values=c("#E69F00", "#56B4E9"))+
  scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
  facet_grid(PayoffCondition~.)+
  theme(legend.position='right',strip.background=element_blank(), legend.key=element_rect(color=NA))
p5
ggsave(filename = 'plots/indCurves3max.pdf', p5, height =  4.3, width = 2.8, unit='in', useDingbats=FALSE)
