#script to produce plots for GP grid experiment
#Eric Schulz, July 2016

#house keeping
rm(list=ls())
theme_set(theme_bw(base_size=16))# use the b&w theme

#load packages
packages <- c('plyr', 'ggplot2', 'jsonlite', 'MASS', 'gridExtra', 'zoo')
lapply(packages, require, character.only = TRUE)

####################################################################################################################################################################################
# IMPORT DATA
####################################################################################################################################################################################
source('dataMunging.R') #source data import function

d <- dataImport(normalize=FALSE)

####################################################################################################################################################################################
# AGGREGATE DATA (compute mean, max, resp. standard errors; repeat for 5-trial average)
####################################################################################################################################################################################
#standard error function
se<-function(x){return(sd(x)/sqrt(length(x)))}

#summarize results for different average reward and max reward performance measures
dplota<-ddply(d,~trial+Horizon+scenario+kernel,summarise,mean=mean(z), se=se(z))
dplota$Measure<-"Average"
dplotb<-ddply(d,~trial+Horizon+scenario+trial+kernel,summarise,mean=mean(zmax), se=se(zmax))
dplotb$Measure<-"Maximum found"

#repeat, but calculate 5-trial average
d$trial5<-round(d$trial/5)*5
d$trial5<-ifelse(d$trial5<5,0,d$trial5)
dplot5a<-ddply(d,~trial5+scenario+Horizon+kernel,summarise,mean=mean(z), se=se(z))
dplot5a$Measure <- "Avg. Reward"
dplot5b<-ddply(d,~scenario+trial5+Horizon+kernel,summarise,mean=mean(zmax), se=se(zmax))
dplot5b$Measure <- "Max. Reward"


####################################################################################################################################################################################
# AGGREGATE DATA (compute mean, max, resp. standard errors; repeat for 5-trial average)
####################################################################################################################################################################################

#Set into parent folder to read kernel files
setwd("..")
setwd("experiment")

#A. Simulate rough environments
dfinal <- fromJSON("kernel1.json", flatten=TRUE)
mu<-mx<-matrix(0, 40, 10000)
for(i in 1:10000){
  m<-as.numeric(unlist(dfinal[[sample(1:20,1)]]))[seq(2,3*121,3)]*100
  s<-sample(1:121, 40, replace=TRUE)
  mu[,i]<-m[s]
  mx[,i]<-maxton(m[s])
}

#save results
dplotRanda<-data.frame(trial=rep(1:40, 4),
                    Horizon=rep("Long", 160), 
                    scenario=rep("Random Model",160),
                    kernel=rep("Rough", 160), 
                    mean=rep(c(apply(mu,1,mean), apply(mx,1,mean)),2),
                    se=rep(c(apply(mu,1,se), apply(mx,1,se)),2),
                    Measure=c(rep("Avg. Reward", 40), rep("Max. Reward", 40),rep("Avg. Reward", 40), rep("Max. Reward", 40)))

#B. simulate smooth environments
dfinal <- fromJSON("kernel2.json", flatten=TRUE)
mu<-mx<-matrix(0, 40, 10000)
for(i in 1:10000){
  m<-as.numeric(unlist(dfinal[[sample(1:20,1)]]))[seq(2,3*121,3)]*100
  s<-sample(1:121, 40, replace=TRUE)
  mu[,i]<-m[s]
  mx[,i]<-maxton(m[s])
}


dplotRandb<-data.frame(trial=rep(1:40, 4),
                    Horizon=rep("Long", 160), 
                    scenario=rep("Random Model",160),
                    kernel=rep("Smooth", 160),
                    mean=rep(c(apply(mu,1,mean), apply(mx,1,mean)),2),
                    se=rep(c(apply(mu,1,se), apply(mx,1,se)),2),
                    Measure=c(rep("Avg. Reward", 40), rep("Max. Reward", 40),rep("Avg. Reward", 40), rep("Max. Reward", 40)))


#Step back into analysis folder
setwd("..")
setwd("analysis")

#Join data frames together
dplotRand<-rbind(dplotRanda, dplotRandb)
dplot<-rbind(dplota, dplotb, dplotRand)

dplot$Horizon <- factor(dplot$Horizon, levels=c("Short", "Long", "Random"))
dplot$kernel <- factor(dplot$kernel, levels=c("Smooth", "Rough"))
dplot$Measure <- factor(dplot$Measure, levels=c("Average", "Maximum found"))

#Repeat aggregation for 5-trial average
dplotRanda$trial5<-round(dplotRanda$trial/5)*5
dplotRandb$trial5<-round(dplotRandb$trial/5)*5
dplotRanda$trial5<-ifelse(dplotRanda$trial5<5,0,dplotRanda$trial5)
dplotRandb$trial5<-ifelse(dplotRandb$trial5<5,0,dplotRandb$trial5)
#Join to existing dataframes
dplot5RandSuma<- ddply(dplotRanda,~scenario+trial5+Horizon+kernel+Measure,summarise,mean=mean(mean), se=se(mean))
dplot5RandSumb<- ddply(dplotRandb,~scenario+trial5+Horizon+kernel+Measure,summarise,mean=mean(mean), se=se(mean))
dplot5 <- rbind(dplot5a, dplot5b, dplot5RandSuma, dplot5RandSumb)


################################################################################################################
# MAIN PLOT
#################################################################################################################position
pd <- position_dodge(.1)
#Average Reward as Measure
p1<-ggplot(subset(dplot5,Measure=="Avg. Reward"), aes(x=trial5, y=mean, colour=scenario)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.1, size=1, position=pd) +
  #geom_point(size=1.5)+
  geom_line(aes(linetype=Horizon),position=pd, size=1) +
  ylab("Avg. Reward (±SE)")+xlab("Trial")+
  theme_bw()+
  ylim(40,80)+
  scale_color_manual(values = c("#A20000", "#0082C1", "black"))+
  theme(text = element_text(size=18,  family="serif"), legend.position="top")+
  facet_wrap(~kernel)+
  theme(legend.position="bottom", strip.background=element_blank(), legend.key=element_rect(color=NA))
p1

ggsave(filename = "plots/avgReward.pdf", plot = p1, height =3.66, width = 4.25, units = "in")


p2<-ggplot(subset(dplot5,Measure=="Max. Reward"), aes(x=trial5, y=mean, colour=scenario)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.1, size=1, position=pd) +
  #geom_point(size=1.5)+
  geom_line(aes(linetype=Horizon),position=pd, size=1) +
  ylab("Max. Reward (±SE)")+xlab("Trial")+
  theme_bw()+
  coord_cartesian(ylim=c(60,100))+
  scale_color_manual(values = c("#A20000", "#0082C1", "black"))+
  theme(text = element_text(size=18,  family="serif"), legend.position="top")+
  facet_wrap(~kernel)+
  theme(legend.position="bottom", strip.background=element_blank(), legend.key=element_rect(color=NA))
p2
ggsave(filename = "plots/maxReward.pdf", plot = p2, height =3.66, width = 4.25, units = "in")

################################################################################################################
# OTHER SUMMARY PLOTS
################################################################################################################
#summarize over all subjects, the number of unique tiles clicked and the number of repeats
uniqueTiles <- function(x,y){ return(nrow(unique(cbind(x,y))))}
summarydf <- ddply(d, .(id, scenario, kernel, horizon, round), summarize, uniqueTiles = uniqueTiles(x,y)) 
summarydf$repeats <- (summarydf$horizon +1) - summarydf$uniqueTiles #total minus unique

p3 <- ggplot(summarydf, aes(x=factor(horizon), y = uniqueTiles, fill=scenario))+
  stat_summary(fun.y = mean, geom = "bar", position = "dodge") + 
  stat_summary(fun.data = mean_se, geom = "errorbar", position = position_dodge(width = 0.90), width = 0.2 ) +
  facet_wrap(~kernel) +
  scale_fill_manual(values = c("#A20000", "#0082C1")) +
  ylab("Number of Unique Tiles Revealed") +
  xlab("Search Horizon") +
  ggtitle("Unique tiles clicked (per round)") +
  theme(text = element_text(size=18,  family="serif"))+
  theme(legend.position="bottom", strip.background=element_blank(), legend.key=element_rect(color=NA))
p3
ggsave(filename = "plots/unique.pdf", plot = p3,height =3.66, width = 4.25, units = "in")

p4 <- ggplot(summarydf, aes(x=factor(horizon), y = repeats, fill=scenario))+
  stat_summary(fun.y = mean, geom = "bar", position = "dodge") + 
  stat_summary(fun.data = mean_se, geom = "errorbar", position = position_dodge(width = 0.90), width = 0.2 ) +
  facet_wrap(~kernel) +
  scale_fill_manual(values = c("#A20000", "#0082C1")) +
  ylab("Number of Repeat Clicks") +
  xlab("Search Horizon") +
  ggtitle("Repeat clicks (per round)") +
  theme(text = element_text(size=18,  family="serif"))+
  theme(legend.position="bottom", strip.background=element_blank(), legend.key=element_rect(color=NA))
p4
ggsave(filename = "plots/repeat.pdf", plot = p4, height =3.66, width = 4.25, units = "in")

g <- ggplot(d, aes(x=z, color=scenario, fill=scenario)) +
  geom_density(alpha=.1, adjust=1) +
  facet_wrap(kernel~horizon, nrow=1) +
  ggtitle("Distribution of Rewards") +
  scale_fill_manual(values = c("#A20000", "#0082C1")) +
  scale_color_manual(values = c("#A20000", "#0082C1")) +
  ylab("Density") +
  xlab("Payoff value") +
  theme(text = element_text(size=18,  family="serif"))+
  theme(legend.position="bottom", strip.background=element_blank(), legend.key=element_rect(color=NA))

ggsave(filename = "plots/density.pdf", plot = p5, height =3.66, width = 8.5, units = "in")



################################################################################################################
# T-Tests
################################################################################################################
#unique tiles
uniqueDF <- ddply(summarydf, .(id, scenario, kernel), summarize, uniqueTiles = mean(uniqueTiles))

################################################################################################################
# OLD PLOTS
################################################################################################################

#standard error function
se<-function(x){return(sd(x)/sqrt(length(x)))}
#aggregate over 5 trial rounds
d$trial2<-round(d$trial/5)*5
d$trial2<-ifelse(d$trial2<5,0,d$trial2)
#recode smoothness
d$kernel<-ifelse(d$kernel==0, "Rough", "Smooth")
#recode horizon
d$horizon<-ifelse(d$horizon==20, "Short", "Long")
#recode reward function
d$scenario<-ifelse(d$scenario==0, "AvgReward", "MaxReward")

#aggregate for plot
dplot<-ddply(d,~scenario+trial2+horizon+kernel,summarise,AvgReward=mean(z), se=se(z))

#position
pd <- position_dodge(.1)



#plot mean over trials per time horizon and reward function
ggplot(dplot, aes(x=trial2, y=AvgReward, colour=horizon, linetype = kernel)) +
  geom_errorbar(aes(ymin=AvgReward-se, ymax=AvgReward+se), width=.7, position=pd) +
  #geom_point(size=2.5)+
  #scale_shape_discrete(solid=F) +
  geom_line(position=pd, lwd=.7) +
  scale_linetype_manual(values = c("dashed","solid"))+
  ylab("Avg. Reward (±SE)")+ggtitle("Average Reward")+xlab("Trial")+
  facet_wrap(Measure~ scenario, ncol = 2) +
  theme_bw()+
  ylim(0,80)+
  scale_colour_manual(values = c("#A20000", "#0082C1"))+
  theme_bw() +
  guides(linetype = guide_legend("Environment"), colour = guide_legend("Horizon")) +
  #theme(text = element_text(size=20,  family="serif"), legend.position="bottom")+
  theme(legend.position="bottom", strip.background=element_blank(), legend.key=element_rect(color=NA))
p

#write to pdf

ggsave("plots/AvgReward.pdf", plot = p, height =4.33, width = 7, units = "in")


#create the maximum value seen
zmax<-numeric()
#loop through ids
for (i in 1:max(d$id)){
  #loop through rounds
  for (j in 1:max(round)){
    #select round per participant
    dummy<-subset(d, id==i & round==j)
    #get first value
    start<-rep(dummy$z[1], nrow(dummy))
    #loop through values
    for (k in 2:nrow(dummy)){
      #if bigger than prior value, accept, else keep prior value
      start[k]<-ifelse(dummy$z[k]>start[k-1], dummy$z[k], start[k-1])
    }
    #concatenate
    zmax<-c(zmax,start)
  }
}

#add to frame
d$zmax<-zmax

#aggregate over trials, reward, horizon, and kernel
dplot<-ddply(d,~scenario+trial2+horizon+kernel,summarise, MaxReward=mean(zmax), se=se(zmax))

#plot maximum seen value
pd <- position_dodge(.1)

#plot mean over trials per time horizon and reward function
p<-ggplot(dplot, aes(x=trial2, y=MaxReward, colour=horizon, linetype = kernel)) +
  geom_errorbar(aes(ymin=MaxReward-se, ymax=MaxReward+se), width=.7, position=pd) +
  #geom_point(size=2.5)+
  #scale_shape_discrete(solid=F) +
  geom_line(position=pd, lwd=.7) +
  scale_linetype_manual(values = c("dashed","solid"))+
  ylab("Max Reward (±SE)")+ggtitle("Max Reward")+xlab("Trial")+
  facet_wrap(~ scenario, ncol = 2) +
  theme_bw()+
  ylim(0,100)+
  scale_colour_manual(values = c("#A20000", "#0082C1"))+
  theme_bw() +
  guides(linetype = guide_legend("Environment"), colour = guide_legend("Horizon")) +
  #theme(text = element_text(size=20,  family="serif"), legend.position="bottom")+
  theme(legend.position="bottom", strip.background=element_blank(), legend.key=element_rect(color=NA))

#write to pdf
ggsave("plots/MaxReward.pdf", plot = p, height =4.33, width = 7, units = "in")


#select all trial up to 20
d20<-subset(d, trial<=20)

#mean reward for different scenarios
d1<-ddply(d20,~scenario, summarise,mean=mean(z), se=se(z))
#mean reward for different smoothness
d2<-ddply(d20,~kernel, summarise,mean=mean(z), se=se(z))
#mean reward for different horizons
d3<-ddply(d20,~horizon, summarise,mean=mean(z), se=se(z))

#dummy function for approximate 95% confidence interval
limits <- aes(ymax = mean + 2*se, ymin=mean - 2*se)

#plot for scenarios
p1 <- ggplot(d1, aes(y=mean, x=scenario)) + 
  #bars
  geom_bar(stat="identity")+
  #0 to 1
  ylim(c(0,100)) + 
  #golden ratio error bars
  geom_errorbar(limits, position="dodge", width=0.31)+
  #point size
  geom_point(size=3)+
  #title
  ggtitle("Reward")+theme_classic() +xlab("\nScenario")+ylab("Mean Reward\n")+
  #change fonts
  theme(text = element_text(size=20, family="serif"))

p1

#plot for smoothness
p2 <- ggplot(d2, aes(y=mean, x=kernel)) + 
  #bars
  geom_bar(stat="identity")+
  #0 to 1
  ylim(c(0,100)) + 
  #golden ratio error bars
  geom_errorbar(limits, position="dodge", width=0.31)+
  #point size
  geom_point(size=3)+
  #title
  ggtitle("Smoothness")+theme_classic() +xlab("\nEnvironment")+ylab("Mean Reward\n")+
  #change font
  theme(text = element_text(size=20, family="serif"))

p2

#plot for horizon
p3 <- ggplot(d3, aes(y=mean, x=horizon)) + 
  #bars
  geom_bar(stat="identity")+
  #0 to 1
  ylim(c(0,100)) + 
  #golden ratio error bars
  geom_errorbar(limits, position="dodge", width=0.31)+
  #point size
  geom_point(size=3)+
  #title
  ggtitle("Horizon")+theme_classic() +xlab("\nHorizon")+ylab("Mean Reward\n")+
  #change font
  theme(text = element_text(size=20, family="serif"))

p3

#plot all three in one pdf in one row
pdf("plots/comparisonreward.pdf", width=12, height=6)
grid.arrange(p1,p2,p3, nrow=1)
dev.off()



#aggregate maximum seen over reward function
d1<-ddply(d20,~scenario, summarise,mean=mean(zmax), se=se(zmax))
#aggregate maximum seen over kernel
d2<-ddply(d20,~kernel, summarise,mean=mean(zmax), se=se(zmax))
#aggregate maximum seen over horizon
d3<-ddply(d20,~horizon, summarise,mean=mean(zmax), se=se(zmax))

#95%CI
limits <- aes(ymax = mean + 2*se, ymin=mean - 2*se)

#plot for scenarios
p1 <- ggplot(d1, aes(y=mean, x=scenario)) + 
  #bars
  geom_bar(stat="identity")+
  #0 to 1
  ylim(c(0,100)) + 
  #golden ratio error bars
  geom_errorbar(limits, position="dodge", width=0.31)+
  #point size
  geom_point(size=3)+
  #title
  ggtitle("Reward")+theme_classic() +xlab("\nScenario")+ylab("Mean Best\n")+
  #change font
  theme(text = element_text(size=20, family="serif"))
p1

#plot for smoothness
p2 <- ggplot(d2, aes(y=mean, x=kernel)) + 
  #bars
  geom_bar(stat="identity")+
  #0 to 1
  ylim(c(0,100)) + 
  #golden ratio error bars
  geom_errorbar(limits, position="dodge", width=0.31)+
  #point size
  geom_point(size=3)+
  #title
  ggtitle("Smoothness")+theme_classic() +xlab("\nEnvironment")+ylab("Mean Best\n")+
  #change font
  theme(text = element_text(size=20, family="serif"))

p2

#plot for horizon
p3 <- ggplot(d3, aes(y=mean, x=horizon)) + 
  #bars
  geom_bar(stat="identity")+
  #0 to 1
  ylim(c(0,100)) + 
  #golden ratio error bars
  geom_errorbar(limits, position="dodge", width=0.31)+
  #point size
  geom_point(size=3)+
  #title
  ggtitle("Horizon")+theme_classic() +xlab("\nHorizon")+ylab("Mean Best\n")+
  #change font
  theme(text = element_text(size=20, family="serif"))

p3

#write all together in one pdf in one row
pdf("plots/comparisonbestfound.pdf", width=12, height=6)
grid.arrange(p1,p2,p3, nrow=1)
dev.off()
#THE END