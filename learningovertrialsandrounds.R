#Learning over trials and rounds
#Eric Schulz and Charley Wu, 2018

rm(list=ls())
library(dplyr)
library(ggplot2)

alldata <- data.frame()
############################################################################################################
# 2D
############################################################################################################
setwd('analysis2D')
source('dataMunging.R')
d <- dataImport(normalize=FALSE) #load experiment data


#Linear model of effect of trials
co<-rep(0, max(d$id))
cose<-rep(0, max(d$id))
for(i in 1:max(d$id)){
  dummy<-subset(d, id==i)
  dummy$z<-scale(dummy$z)
  dummy$trial<-scale(dummy$trial)
  m<-lm(z~trial, data=dummy)
  co[i]<-coef(summary(m))[2,1]
  cose[i]<-coef(summary(m))[2,2]
}


dplot1<-data.frame(variable="Trial", estimate=co, lb=co-cose, ub=co+cose)
dplot1<-dplot1[order(-dplot1$estimate),]
dplot1$level<-order(dplot1$estimate)
dplot1$experiment<-"Experiment 2"

#Linear model of learning over rounds
co<-rep(0, max(d$id))
cose<-rep(0, max(d$id))
for(i in 1:max(d$id)){
  dummy<-subset(d, id==i)
  dummy$z<-scale(dummy$z)
  dummy$round<-scale(as.numeric(dummy$round))
  m<-lm(z~round, data=dummy)
  co[i]<-coef(summary(m))[2,1]
  cose[i]<-coef(summary(m))[2,2]
}

dplot2<-data.frame(variable="Round", estimate=co, lb=co-cose, ub=co+cose)
dplot2<-dplot2[order(-dplot2$estimate),]
dplot2$level<-order(dplot2$estimate)

dplot2$experiment<-"Experiment 2"

d$experiment <- "Experiment 2"
alldata <- rbind(alldata, d)
############################################################################################################
# 1D
############################################################################################################
setwd('..')
setwd('analysis1D')
source('dataMunging.R')
d <- dataImport(normalize=FALSE)

#learning over trials
co<-rep(0, max(d$id))
cose<-rep(0, max(d$id))
for(i in 1:max(d$id)){
  dummy<-subset(d, id==i)
  dummy$y<-scale(dummy$y)
  dummy$trial<-scale(dummy$trial)
  m<-lm(y~trial, data=dummy)
  co[i]<-coef(summary(m))[2,1]
  cose[i]<-coef(summary(m))[2,2]
}



dplot3<-data.frame(variable="Trial", estimate=co, lb=co-cose, ub=co+cose)
dplot3<-dplot3[order(-dplot3$estimate),]
dplot3$level<-order(dplot3$estimate)
dplot3$experiment<-"Experiment 1"

#Learning over  rounds
co<-rep(0, max(d$id))
cose<-rep(0, max(d$id))
for(i in 1:max(d$id)){
  dummy<-subset(d, id==i)
  dummy$y<-scale(dummy$y)
  dummy$round<-scale(as.numeric(dummy$round))
  m<-lm(y~round, data=dummy)
  co[i]<-coef(summary(m))[2,1]
  cose[i]<-coef(summary(m))[2,2]
}

dplot4<-data.frame(variable="Round", estimate=co, lb=co-cose, ub=co+cose)
dplot4<-dplot4[order(-dplot4$estimate),]
dplot4$level<-order(dplot4$estimate)
dplot4$experiment<-"Experiment 1"

d$experiment <- "Experiment 1"
alldata <- rbind.fill(alldata, d)
############################################################################################################
# Experiment 3
############################################################################################################
setwd('..')
setwd('analysis3')
source('dataMunging.R')
d <- dataImport(normalize=FALSE)

#learning over trials
co<-rep(0, max(d$id))
cose<-rep(0, max(d$id))
for(i in 1:max(d$id)){
  dummy<-subset(d, id==i)
  dummy$z<-scale(dummy$z)
  dummy$trial<-scale(dummy$trial)
  m<-lm(z~trial, data=dummy)
  co[i]<-coef(summary(m))[2,1]
  cose[i]<-coef(summary(m))[2,2]
}



dplot5<-data.frame(variable="Trial", estimate=co, lb=co-cose, ub=co+cose)
dplot5<-dplot5[order(-dplot5$estimate),]
dplot5$level<-order(dplot5$estimate)
dplot5$experiment<-"Experiment 3"

#Learning over  rounds
co<-rep(0, max(d$id))
cose<-rep(0, max(d$id))
for(i in 1:max(d$id)){
  dummy<-subset(d, id==i)
  dummy$z<-scale(dummy$z)
  dummy$round<-scale(as.numeric(dummy$round))
  m<-lm(z~round, data=dummy)
  co[i]<-coef(summary(m))[2,1]
  cose[i]<-coef(summary(m))[2,2]
}

dplot6<-data.frame(variable="Round", estimate=co, lb=co-cose, ub=co+cose)
dplot6<-dplot6[order(-dplot6$estimate),]
dplot6$level<-order(dplot6$estimate)
dplot6$experiment<-"Experiment 3"

d$experiment <- "Experiment 3"
alldata <- rbind.fill(alldata, d)

############################################################################################################
# Put it all together
############################################################################################################

dplot<-rbind(dplot3, dplot4, dplot1, dplot2, dplot5, dplot6)

dline<-data.frame(experiment=rep(c("Experiment 1", "Experiment 2", "Experiment 3"), each=2),
                  variable=rep(c("Trial", "Round"), 3),
                  xint=c(mean(dplot3$estimate),mean(dplot4$estimate),mean(dplot1$estimate), mean(dplot2$estimate), mean(dplot5$estimate),mean(dplot5$estimate)))


g1 <- ggplot(dplot, aes(estimate,level,xmin=lb,xmax=ub))+
  geom_errorbarh(height=0)+
  geom_vline(xintercept=0,lty=2)+
  geom_point(size=0.75)+
  geom_vline(data=dline, aes(xintercept=xint), lty=1, col="red")+
  facet_grid(variable~experiment, scales="free")+theme_classic()+xlab("Effect Size")+
  ylab("Participants")+ theme(text = element_text(size=18,  family="sans")) +
  theme(legend.position = "top", strip.background=element_blank(), legend.key=element_rect(color=NA), 
        axis.line=element_line()) 
g1

setwd("..")
setwd("paper")
ggsave("effectsize.pdf")



############################################################################################################
# T-tests
############################################################################################################
library(lsr)

#Experiment 1
t.test(subset(dplot, experiment =='Experiment 1' & variable=='Trial')$estimate, mu=0)
cohensD(subset(dplot, experiment =='Experiment 1' & variable=='Trial')$estimate, mu=0)

t.test(subset(dplot, experiment =='Experiment 1' & variable=='Round')$estimate, mu=0)
cohensD(subset(dplot, experiment =='Experiment 1' & variable=='Round')$estimate, mu=0)


#Experiment 2
t.test(subset(dplot, experiment =='Experiment 2' & variable=='Trial')$estimate, mu=0)
cohensD(subset(dplot, experiment =='Experiment 2' & variable=='Trial')$estimate, mu=0)

t.test(subset(dplot, experiment =='Experiment 2' & variable=='Round')$estimate, mu=0)
cohensD(subset(dplot, experiment =='Experiment 2' & variable=='Round')$estimate, mu=0)



#Experiment 3
t.test(subset(dplot, experiment =='Experiment 3' & variable=='Trial')$estimate, mu=0)
cohensD(subset(dplot, experiment =='Experiment 3' & variable=='Trial')$estimate, mu=0)

t.test(subset(dplot, experiment =='Experiment 3' & variable=='Round')$estimate, mu=0)
cohensD(subset(dplot, experiment =='Experiment 3' & variable=='Round')$estimate, mu=0)
  