#house keeping
rm(list=ls())

#load packages
packages <- c('plyr', 'jsonlite', 'DEoptim', "matrixcalc", "fields")
lapply(packages, require, character.only = TRUE)


####################################################################################################################################################################################
# IMPORT DATA
####################################################################################################################################################################################
source('dataMunging.R') #source data import function

d <- dataImport(normalize=FALSE)


############################################################################################################
# T.TESTS USED IN PAPER
############################################################################################################

#LEARNING EFFECTS
#Correlate trial vs. performance
dm<-ddply(d, ~trial, summarize, mu=mean(z)) #avg. reward
cor.test(dm$trial, dm$mu)

dm<-ddply(d, ~trial, summarize, mu=mean(zmax)) #max reward
cor.test(dm$trial, dm$mu)

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

#max rewards across payoffs 
dm<-ddply(d, ~id+scenario, summarize, mu=max(z))
t.test(subset(dm, scenario=="Avg. Reward Condition")$mu, subset(dm, scenario=="Max. Reward Condition")$mu, var.equal=T)

#SMOOTH+LONG only: avg rewards across payoffs 
dm<-ddply(subset(d,kernel=="Smooth" & horizon==40), ~id+scenario, summarize, mu=mean(z))
t.test(subset(dm, scenario=="Avg. Reward Condition")$mu, subset(dm, scenario=="Max. Reward Condition")$mu, var.equal=T)

#SMOOTH+LONG only: max rewards across payoffs 
dm<-ddply(subset(d,kernel=="Smooth" & horizon==40), ~id+scenario, summarize, mu=max(z))
t.test(subset(dm, scenario=="Avg. Reward Condition")$mu, subset(dm, scenario=="Max. Reward Condition")$mu, var.equal=T)

#ENVIRONMENT
#comparison of environment
dm<-ddply(d, ~id+kernel, summarize, mu=mean(z))
t.test(subset(dm, kernel=="Smooth")$mu, subset(dm, kernel=="Rough")$mu, var.equal=T)

dm<-ddply(d, ~id+kernel, summarize, mu=mean(zmax))
t.test(subset(dm, kernel=="Smooth")$mu, subset(dm, kernel=="Rough")$mu, var.equal=T)

#SEARCH HORIZON

#Average reward
dm<-ddply(subset(d, scenario=="Avg. Reward Condition" & kernel=="Rough"), ~id+horizon, summarize, mu=mean(z))
t.test(subset(dm, horizon==20)$mu, subset(dm, horizon==40)$mu, var.equal=T)


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


