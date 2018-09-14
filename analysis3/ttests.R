#house keeping
rm(list=ls())

#load packages
packages <- c('psych', 'plyr', 'jsonlite', 'DEoptim', "matrixcalc", "fields", "lsr")
lapply(packages, require, character.only = TRUE)


####################################################################################################################################################################################
# IMPORT DATA
####################################################################################################################################################################################
source('dataMunging.R') #source data import function

d <- dataImport(normalize=FALSE)


############################################################################################################
# T.TESTS USED IN PAPER
############################################################################################################

#locality of sampling

sampleSize <- 100000
randomDF <- data.frame(x=sample(x = seq(0:10), size = sampleSize, replace=TRUE), y=sample(x = seq(0:10), size = sampleSize, replace=TRUE))
randomDF <- randomDF %>%
  mutate(delta_x = abs(x - lag(x, default = NA)) + abs(y - lag(y, default = NA)) ) 

d$scenario <- factor(d$scenario)
levels(d$scenario) <-c("Avg. Reward Condition", "Max. Reward Condition")

aggregatedByIndividuals <- ddply(na.omit(d), .(id, scenario), summarize, meanDelta_x = mean(delta_x))
t.test(aggregatedByIndividuals$meanDelta_x, mu=mean(na.omit(randomDF$delta_x)))
psych::cohen.d.ci(cohensD(aggregatedByIndividuals$meanDelta_x, mu=mean(na.omit(randomDF$delta_x))), n1 = 80)

#comparison to random

t.test(subset(aggregatedByIndividuals, scenario=="Avg. Reward Condition")$meanDelta_x, subset(aggregatedByIndividuals, scenario=="Max. Reward Condition")$meanDelta_x, var.equal=T)
psych::cohen.d.ci(cohensD(subset(aggregatedByIndividuals, scenario=="Avg. Reward Condition")$meanDelta_x, subset(aggregatedByIndividuals, scenario=="Max. Reward Condition")$meanDelta_x),n1 =40, n2 = 40) 

#LEARNING EFFECTS
#Correlate trial vs. performance
dm<-ddply(d, ~trial, summarize, mu=mean(z)) #avg. reward
cor.test(dm$trial, dm$mu)

dm<-ddply(d, ~trial, summarize, mu=mean(zmax)) #max reward
cor.test(dm$trial, dm$mu)


#learning over trials
dm<-ddply(d, ~id + trial)
dcor <- ddply(dm, ~id, summarize, cor = cor(z,trial))
t.test(dcor$cor, mu=0)
psych::cohen.d.ci(cohensD(dcor$cor, 0), n=80)

#learning over rounds
dm<-ddply(d, ~id + round)
dcor <- ddply(dm, ~id, summarize, cor = cor(z,as.numeric(round)))
t.test(dcor$cor, mu=0)
psych::cohen.d.ci(cohensD(dcor$cor, 0), n=80)

#Round vs. performance
dm<-ddply(d, ~round, summarize, mu=mean(z)) #avg. reward
cor.test(dm$mu, as.numeric(dm$round))

dm<-ddply(d, ~round, summarize, mu=mean(zmax))
cor.test(dm$mu, as.numeric(dm$round))

#PAY OFF CONDITIONS
#variance of locations sampled across payoffs
dvar<-ddply(d, ~id+scenario+round, summarise, var=var(x)+var(y))
dvar<- ddply(dvar, ~id+scenario, summarize, var=mean(var)) #mean variance over all 8 rounds
t.test(subset(dvar, scenario=="Avg. Reward Condition")$var, subset(dvar, scenario=="Max. Reward Condition")$var, var.equal=T)

#avg rewards across payoffs 
dm<-ddply(d, ~id+scenario, summarize, mu=mean(z))
t.test(subset(dm, scenario=="Avg. Reward Condition")$mu, subset(dm, scenario=="Max. Reward Condition")$mu, var.equal=T)
psych::cohen.d.ci(cohensD(subset(dm, scenario=="Avg. Reward Condition")$mu, subset(dm, scenario=="Max. Reward Condition")$mu), n1 = 40, n2 = 40)
#max rewards across payoffs 
dm<-ddply(d, ~id+scenario, summarize, mu=mean(zmax))
t.test(subset(dm, scenario=="Avg. Reward Condition")$mu, subset(dm, scenario=="Max. Reward Condition")$mu, var.equal=T)
psych::cohen.d.ci(cohensD(subset(dm, scenario=="Avg. Reward Condition")$mu, subset(dm, scenario=="Max. Reward Condition")$mu), n1 =40, n2 = 40)
#SMOOTH+LONG only: avg rewards across payoffs 
dm<-ddply(subset(d,kernel=="Smooth" & horizon==40), ~id+scenario, summarize, mu=mean(z))
t.test(subset(dm, scenario=="Avg. Reward Condition")$mu, subset(dm, scenario=="Max. Reward Condition")$mu, var.equal=T)

#ENVIRONMENT
#comparison of environment
dm<-ddply(d, ~id+kernel, summarize, mu=mean(z))
t.test(subset(dm, kernel=="Smooth")$mu, subset(dm, kernel=="Rough")$mu, var.equal=T)
cohensD(subset(dm, kernel=="Smooth")$mu, subset(dm, kernel=="Rough")$mu)


dm<-ddply(d, ~id+kernel, summarize, mu=mean(zmax))
t.test(subset(dm, kernel=="Smooth")$mu, subset(dm, kernel=="Rough")$mu, var.equal=T)

#SEARCH HORIZON

#Average reward
dm<-ddply(subset(d), ~id+horizon, summarize, mu=mean(z))
t.test(subset(dm, horizon==20)$mu, subset(dm, horizon==40)$mu, var.equal=T, paired=T)
psych::cohen.d.ci(cohensD(subset(dm, horizon==20)$mu, subset(dm, horizon==40)$mu), n = 80)

#MaxReward
dm<-ddply(subset(d, scenario=="MaxReward"), ~id+horizon, summarize, mu=mean(zmax))
t.test(subset(dm, horizon=="Long")$mu, subset(dm, horizon=="Short")$mu, var.equal=T)





dm<-ddply(subset(d, scenario=="MaxReward"), ~id+horizon, summarize, mu=mean(z))
t.test(subset(dm, horizon=="Short")$mu, subset(dm, horizon=="Long")$mu, var.equal=T)


dm<-ddply(d, ~id+scenario, summarize, mu=mean(z))
t.test(subset(dm, scenario=="AvgReward")$mu, subset(dm, scenario=="MaxReward")$mu, var.equal=T)

dm<-ddply(d, ~id+scenario, summarize, mu=mean(zmax))
t.test(subset(dm, scenario=="AvgReward")$mu, subset(dm, scenario=="MaxReward")$mu, var.equal=T)


dm<- ddply(dm, ~id+scenario, summarize, mu=mean(mu))
t.test(subset(dm, scenario=="AvgReward")$mu, subset(dm, scenario=="MaxReward")$mu, var.equal=T)


