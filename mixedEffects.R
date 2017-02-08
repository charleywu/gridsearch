#Mixed Effects Modeling

#house keeping
rm(list=ls())
#load packages
#packages <- c('plyr', 'jsonlite', "matrixcalc", "fields", "lme4", "lmerTest", "stargazer", "arm")
packages <- c('plyr', 'jsonlite', "matrixcalc", "fields", "lme4", "lmerTest", "stargazer", "arm", "xtable")
lapply(packages, require, character.only = TRUE)

maxton<-function(x){
  maxn<-rep(0,length(x))
  maxn[1]<-x[1]
  for (i in 2:length(x)){
    if (x[i]>maxn[i-1]){maxn[i]<-x[i]}
    else{maxn[i]<-maxn[i-1]}
  }
  return(maxn)
}



#read in data
dat<-read.csv("ExperimentData/newdat.csv")

#dummy data frame
d<-data.frame(id=numeric(), trial=numeric(), x=numeric(), y=numeric(), 
              z=numeric(), zmax=numeric(), kernel=numeric(), scenario=numeric(), 
              horizon=numeric(), round=numeric())

#loop through data
for (i in 1:nrow(dat)){
  #write to json file
  write(paste(dat$searchHistory[i]), file = "data.JSON")
  #load from json file just written
  dfinal <- fromJSON(txt="data.JSON")
  #sampled x value
  x<-unlist(dfinal$xcollect)
  #sampled y value
  y<-unlist(dfinal$ycollect)
  #sampled z value
  z<-unlist(dfinal$zcollect)
  #length
  len<-sapply(dfinal[1]$xcollect, length)
  #trial 1:20 or 1:40
  trial<-c(1:len[1],1:len[2],1:len[3],1:len[4])
  #round, 1:4
  round<-rep(1:8, len)
  #horizon, 20 or 40
  horizon<-rep(len-1, len)
  #maximum value or cumulative
  scenario<-rep(dat$scenario[i], sum(len))
  #smoothness
  kernel<-rep(dat$kernel[i], sum(len))
  #id number
  id<-rep(i, sum(len))
  zmax<-unlist(sapply(dfinal$zcollect, maxton))
  #dummy frame
  dummy<-data.frame(id, trial, x, y, z, zmax, kernel, scenario, horizon, round)
  #bind them together
  d<-rbind(d, dummy)
}

d$kernel<-ifelse(d$kernel==0, "Rough", "Smooth")
#recode horizon
d$horizon<-ifelse(d$horizon==20, "Short", "Long")
#recode reward function
d$scenario<-ifelse(d$scenario==0, "AvgReward", "MaxReward")

head(d)

meannull <- lmer(z ~ (1|round), data=d, REML=FALSE)
summary(meannull)


meanm1 <- lmer(z ~ trial+(1|round), data=d, REML=FALSE)
summary(meanm1)
anova(meannull, meanm1)

meanm2a<-lmer(z ~ trial+kernel+(1|round), data=d, REML=FALSE) #Trial and kernel
summary(meanm2a)
anova(meanm1, meanm2a)

meanm2b<-lmer(z ~ trial+scenario+(1|round), data=d, REML=FALSE) #Trial and scenario
summary(meanm2b)
anova(meanm1, meanm2b)


meanm3<-lmer(z ~ trial+kernel+horizon+(1|round), data=d, REML=FALSE)
summary(meanm3)
anova(meanm2a, meanm3)

meanm4<-lmer(z ~ trial+kernel+horizon+scenario+(1|round), data=d, REML=FALSE)
summary(meanm4)
anova(meanm3, meanm4)

meanm5<-lmer(z ~ trial+kernel+horizon+scenario+(1|round)+(1|id), data=d, REML=FALSE)
summary(meanm5)
anova(meanm4, meanm5)
coef(meanm5)


maxnull <- lmer(zmax ~ (1|round), data=d, REML=FALSE)

maxm1 <- lmer(zmax ~ trial+(1|round), data=d, REML=FALSE)
summary(maxm1)
anova(maxnull, maxm1)

maxm2<-lmer(zmax ~ trial+kernel+(1|round), data=d, REML=FALSE)
summary(maxm2)
anova(maxm1, maxm2)

maxm3<-lmer(zmax ~ trial+kernel+horizon+(1|round), data=d, REML=FALSE)
summary(maxm3)
anova(maxm2, maxm3)

maxm4<-lmer(zmax ~ trial+kernel+horizon+scenario+(1|round), data=d, REML=FALSE)
summary(maxm4)
anova(maxm3, maxm4)

maxm5<-lmer(zmax ~ trial+kernel+horizon+scenario+(1|round)+(1|id), data=d, REML=FALSE)
summary(maxm5)
anova(maxm4, maxm5)

dvar<-ddply(d, ~id+kernel+scenario+horizon+round, summarise, var=var(x)+var(y))
mvar<-glm(var~scenario+kernel+horizon, data=dvar)
summary(mvar)



##Modes to include in paper
meanm5
arm::display(meanm5)
summary(meanm5, digits=2)

maxm5
arm::display(maxm5)
summary(maxm5, digits=2)

mvar
arm::display(mvar)
stargazer(mvar, type='text')