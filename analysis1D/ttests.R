#house keeping
rm(list=ls())

#load packages
packages <- c('plyr', 'jsonlite', 'DEoptim', "matrixcalc", "fields", "lsr")
lapply(packages, require, character.only = TRUE)


####################################################################################################################################################################################
# IMPORT DATA
####################################################################################################################################################################################
source('dataMunging.R') #source data import function

d <- dataImport(normalize=FALSE)


############################################################################################################
# T.TESTS USED IN PAPER
############################################################################################################

#Locality of sampling

#locality of sampling
sampleSize <- 100000
randomDF <- data.frame(x=sample(x = seq(0:29), size = sampleSize, replace=TRUE), kernel=c(rep("Rough",sampleSize/2), rep("Smooth", sampleSize/2)), scenario = rep(c('Accumulators', 'Maximizers'), sampleSize))
randomDF <- randomDF %>%
  mutate(delta_x = abs(x - lag(x, default = NA)))

d$scenario <- factor(d$scenario)
levels(d$scenario) <-c("Avg. Reward Condition", "Max. Reward Condition")

aggregatedByIndividuals <- ddply(na.omit(d), .(id, scenario), summarize, meanDelta_x = mean(delta_x))
t.test(aggregatedByIndividuals$meanDelta_x, mu=mean(na.omit(randomDF$delta_x)))
cohensD(aggregatedByIndividuals$meanDelta_x, mu=mean(na.omit(randomDF$delta_x)))
#comparison to random

t.test(subset(aggregatedByIndividuals, scenario=="Avg. Reward Condition")$meanDelta_x, subset(aggregatedByIndividuals, scenario=="Max. Reward Condition")$meanDelta_x, var.equal=T)
cohensD(subset(aggregatedByIndividuals, scenario=="Avg. Reward Condition")$meanDelta_x, subset(aggregatedByIndividuals, scenario=="Max. Reward Condition")$meanDelta_x)


#LEARNING EFFECTS
#Correlate trial vs. performance
dm<-ddply(d, ~trial, summarize, mu=mean(y)) #avg. reward
cor.test(dm$trial, dm$mu)

dm<-ddply(d, ~trial, summarize, mu=mean(ymax)) #max reward
cor.test(dm$trial, dm$mu)

#learning over trials
dm<-ddply(d, ~id + trial)
dcor <- ddply(dm, ~id, summarize, cor = cor(y,trial))
t.test(dcor$cor, mu=0)
cohensD(dcor$cor, 0)

#learning over rounds
dm<-ddply(d, ~id + round)
dcor <- ddply(dm, ~id, summarize, cor = cor(y,as.numeric(round)))
t.test(dcor$cor, mu=0)
cohensD(dcor$cor, 0)

#Round vs. performance
dm<-ddply(d, ~id + round, summarize, mu=mean(y)) #avg. reward
cor.test(dm$mu, as.numeric(dm$round))

dm<-ddply(d, ~round, summarize, mu=mean(ymax))
cor.test(dm$mu, as.numeric(dm$round))

#PAY OFF CONDITIONS
#variance of locations sampled across payoffs
dvar<-ddply(d, ~id+scenario+round, summarise, var=var(x))
dvar<- ddply(dvar, ~id+scenario, summarize, var=mean(var)) #mean variance over all 8 rounds
t.test(subset(dvar, scenario=="Avg. Reward Condition")$var, subset(dvar, scenario=="Max. Reward Condition")$var, var.equal=T)

#avg rewards across payoffs 
dm<-ddply(d, ~id+scenario, summarize, mu=mean(y))
t.test(subset(dm, scenario=="Avg. Reward Condition")$mu, subset(dm, scenario=="Max. Reward Condition")$mu, var.equal=T)

#max rewards across payoffs 
dm<-ddply(d, ~id+scenario, summarize, mu=mean(ymax))
t.test(subset(dm, scenario=="Avg. Reward Condition")$mu, subset(dm, scenario=="Max. Reward Condition")$mu, var.equal=T)
cohensD(subset(dm, scenario=="Avg. Reward Condition")$mu, subset(dm, scenario=="Max. Reward Condition")$mu)

#SMOOTH+LONG only: avg rewards across payoffs 
dm<-ddply(subset(d,kernel=="Smooth" & horizon==10), ~id+scenario, summarize, mu=mean(y))
t.test(subset(dm, scenario=="Avg. Reward Condition")$mu, subset(dm, scenario=="Max. Reward Condition")$mu, var.equal=T)

#SMOOTH+LONG only: max rewards across payoffs 
dm<-ddply(subset(d,kernel=="Smooth" & horizon==10), ~id+scenario, summarize, mu=max(y))
t.test(subset(dm, scenario=="Avg. Reward Condition")$mu, subset(dm, scenario=="Max. Reward Condition")$mu, var.equal=T)

#ENVIRONMENT
#comparison of environment
dm<-ddply(d, ~id+kernel, summarize, mu=mean(y))
t.test(subset(dm, kernel=="Smooth")$mu, subset(dm, kernel=="Rough")$mu, var.equal=T)

dm<-ddply(d, ~id+kernel, summarize, mu=mean(ymax))
t.test(subset(dm, kernel=="Smooth")$mu, subset(dm, kernel=="Rough")$mu, var.equal=T)

#SEARCH HORIZON

#Average reward

dm<-ddply(subset(d), ~id+horizon, summarize, mu=mean(y))
t.test(subset(dm, horizon==10)$mu, subset(dm, horizon==5)$mu, var.equal=T, paired=T)

dm<-ddply(subset(d, scenario=="Avg. Reward Condition" & kernel=="Rough"), ~id+horizon, summarize, mu=mean(z))
t.test(subset(dm, horizon==20)$mu, subset(dm, horizon==40)$mu, var.equal=T)


#MaxReward
dm<-ddply(subset(d, scenario=="MaxReward"), ~id+horizon, summarize, mu=mean(ymax))
t.test(subset(dm, horizon=="Long")$mu, subset(dm, horizon=="Short")$mu, var.equal=T)



dm<-ddply(subset(d, scenario=="MaxReward"), ~id+horizon, summarize, mu=mean(y))
t.test(subset(dm, horizon=="Short")$mu, subset(dm, horizon=="Long")$mu, var.equal=T)


dm<-ddply(d, ~id+scenario, summarize, mu=mean(y))
t.test(subset(dm, scenario=="Avg. Reward Condition")$mu, subset(dm, scenario=="Max. Reward Condition")$mu, var.equal=T)

dm<-ddply(d, ~id+scenario, summarize, mu=mean(ymax))
t.test(subset(dm, scenario=="Avg. Reward Condition")$mu, subset(dm, scenario=="Max. Reward Condition")$mu, var.equal=T)


dm<- ddply(dm, ~id+scenario, summarize, mu=mean(mu))
t.test(subset(dm, scenario=="Avg. Reward Condition")$mu, subset(dm, scenario=="Max. Reward Condition")$mu, var.equal=T)


