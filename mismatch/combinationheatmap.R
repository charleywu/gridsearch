#Mismatch lambda simulations
#Charley Wu and Eric Schulz, 2018

library(scales)
library(plyr)
library(ggplot2)
se <- function(x) sqrt(var(x)/length(x))
rootDir <- rprojroot::find_rstudio_root_file()
setwd(rootDir)
############################################################################################################
## Generalized case                                                                                       ##
############################################################################################################
setwd('mismatch')
#Load data
dr<-read.csv("regret.csv")
head(dr)
dr<-subset(dr, setting=="gp_ucb") #select only gp-ucb simulations
dr$X
dr<-subset(dr, X %in% c(0,1,2,4))


dp1<-data.frame(trial=dr$X, lteach=dr$l0, llearn=dr$l1, Score=dr$m) #construct new dataframe for plotting
dp1$Score<- rescale(-log(dp1$Score)) #recale score to log
dp1$trial
dp1$trial<-mapvalues(dp1$trial, c(0,1,2,4), c("t=1","t=5","t=10","t=20")) #remap x to trial numbers
dp1$trial<-factor(dp1$trial,levels= c("t=1","t=5","t=10","t=20"))
dp1$sim<-"Generalized"

dl1<-ddply(dp1, ~lteach+trial, summarize, y=median(Score[lteach>llearn])-median(Score[lteach<=llearn])) #compute optimal score
dl1$llearn<-dl1$lteach-dl1$y
dl1$llearn<-ifelse(dl1$llearn<0.1,0.1,dl1$llearn)
dl1$llearn<-ifelse(dl1$llearn>0.9,0.9,dl1$llearn)

nrow(dp2)

nrow(dp1)*11


dd<-subset(dp1, lteach==llearn)

p1<-ggplot(dp1, aes(x = lteach, y =llearn, fill = Score)) +
  geom_tile() +
  coord_equal()+
  #scale_fill_gradient(trans = 'log')+
  scale_fill_distiller(palette = "Spectral", direction = -1)+
  #geom_point(data=dl1, col="black", fill="black", size=1.2, shape=4) +
  facet_grid(~trial)+
  scale_x_continuous(expression("Teacher"~lambda[0]), expand = c(0, 0), breaks=seq(0.2,0.8, 0.2))+ 
  scale_y_continuous(expression("Student"~lambda[1]), expand = c(0, 0), breaks=seq(0.2,0.8, 0.2))+
  theme_bw()+
  ggtitle(expression("Generalized effect of mismatch"))+
  geom_abline(slope=1, intercept=0, linetype = 'dotted', size=1.5)+
  theme(legend.key.width = unit(0.2, "cm"),
        legend.key.height = unit(0.4, "cm"))+
  theme(text = element_text(size=14, family="sans"))+
  theme(strip.background = element_blank())
p1

ggsave(plot = p1, filename='generalizedPlot.pdf', height =4.33, width = 7, units = "in", useDingbats=F)
############################################################################################################
## Experiment 1                                                                                           ##
############################################################################################################
setwd('mismatch')
L<-expand.grid(seq(0.1,3,0.1), seq(0.1,3,0.1))

d<-read.csv("mismatch1d.csv")
d$lteach<-rep(L[,1], each=11)
d$llearn<-rep(L[,2], each=11)
d <- d%>% group_by(lteach, llearn) %>% mutate(avgCumReward = cumsum(m)/seq_along(m))
subset(d, lteach==1 & llearn==1)
dp2<-subset(d, trial %in% c(1,3,5,10))
#dp2$m<- rescale(dp2$m, from=c(-1.8, 1.8), to=c(0,100))
dp2$avgCumReward<- rescale(dp2$avgCumReward, from=c(-1.8, 1.8), to=c(0,100))
dp2$trial<-mapvalues(dp2$trial, c(1,3,5,10), c("t=1","t=3","t=5","t=10"))
dp2$trial<-factor(dp2$trial,levels= c("t=1","t=3","t=5","t=10"))
dp2$Score<-dp2$avgCumReward
dp2$X<-NULL
dp2$m<-NULL


#Human estimates
setwd('..')
setwd('analysis1D')
gpParams <- read.csv('rationalModels/parameters/gp.csv')
setwd('..')
gpParams$lteach <- ifelse(gpParams$environment=='Rough', 1, 2)
ldf <- gpParams %>% group_by(environment) %>%
  dplyr::summarise(medLambda = median(lambda),  meanLambda = mean(lambda[!lambda %in% boxplot.stats(lambda)$out]), lower=quantile(lambda, .25), upper=quantile(lambda, .75)) 
  
ldf<- ldf %>%  mutate(lteach = ifelse(environment=='Smooth', 2, 1))

p2<-ggplot(dp2, aes(x = lteach, y =llearn)) +
  geom_tile(aes(fill = Score)) +
  coord_fixed(ratio = 1, xlim = c(0,3), ylim = c(0,3), expand = FALSE)+
  geom_boxplot(data=gpParams, aes(group=environment, x=lteach, y=lambda), width = 0.2, fill=NA, outlier.shape = NA, color='black')+
  #geom_errorbar(data=ldf, aes(x=lteach, ymin=lower, ymax=upper, y= medLambda), fill=NA, width=0.1)+
  geom_point(data=ldf , aes(group=environment,x=lteach, y=meanLambda, shape=environment), fill='black')+
  #geom_point(data=dl2, aes(shape=factor(lteach)), col="black", fill=NA, size=2) +
  #scale_fill_gradient(trans = 'log')+
  #scale_fill_distiller(palette = "Spectral", direction = -1, limits = c(40,105), breaks=seq(40,100,20), labels=seq(40,100,20))+
  scale_fill_distiller(palette = "Spectral", direction = -1)+
  facet_grid(~trial)+
  scale_x_continuous(expression("Teacher"~lambda[0]), expand = c(0, 0), breaks=seq(0.5,2.5,1))+ 
  scale_y_continuous(expression("Student"~lambda[1]), expand = c(0, 0),  breaks=seq(0.5,2.5,1))+
  theme_bw()+
  scale_shape_manual(labels = c(expression(hat(lambda)["Rough"]),expression(hat(lambda)["Smooth"])), values=c(2,1))+
  labs(fill="Reward", shape="Human \nEstimates")+
  ggtitle(expression("Experiment 1"))+
  geom_abline(slope=1, intercept=0, linetype = 'dotted', size=1)+
  theme(legend.key.width = unit(0.2, "cm"),
        legend.key.height = unit(0.4, "cm"))+
  theme(text = element_text(size=14, family="sans"))+
  theme(strip.background = element_blank())
p2

setwd('mismatch')
ggsave(plot = p2, filename='Exp1Plot.pdf', height =4.33, width = 7, units = "in", useDingbats=F)
setwd('..')

############################################################################################################
## Experiment 2                                                                                        ##
############################################################################################################
setwd('mismatch')
d<-read.csv("mismatch2d.csv")
d$lteach<-rep(L[,1], each=40)
d$llearn<-rep(L[,2], each=40)
d <- d%>% group_by(lteach, llearn) %>% mutate(avgCumReward = cumsum(m)/seq_along(m))
dp3<-subset(d, trial %in% c(5,10,20,40))
#dp3$m<-rescale(dp3$m, from=c(-1.8, 1.8), to=c(0,100))
dp3$avgCumReward<-rescale(dp3$avgCumReward, from=c(-1.8, 1.8), to=c(0,100))
dp3$trial<-mapvalues(dp3$trial, c(5,10,20,40), c("t=5","t=10","t=20","t=40"))
#dp3$trial<-mapvalues(dp3$trial, c(1,5,20,40), c("t=1","t=5","t=20","t=40"))
dp3$trial<-factor(dp3$trial,levels= c("t=5","t=10","t=20","t=40"))
#dp3$trial<-factor(dp3$trial,levels= c("t=1","t=5","t=20","t=40"))
dp3$Score<-dp3$avgCumReward
dp3$X<-NULL
dp3$m<-NULL
dp3$sim<-"Empirical-2D"


setwd('..')
setwd('analysis2D')
gpParams <- read.csv('rationalModels/parameters/gp.csv')
setwd('..')
gpParams$lteach <- ifelse(gpParams$environment=='Rough', 1, 2)
ldf <- gpParams %>% group_by(environment) %>%
  dplyr::summarise(medLambda = median(lambda), meanLambda = mean(lambda[!lambda %in% boxplot.stats(lambda)$out]), lower=quantile(lambda, .25), upper=quantile(lambda, .75)) 

ldf<- ldf %>%  mutate(lteach = ifelse(environment=='Smooth', 2, 1))


p3<-ggplot(dp3, aes(x = lteach, y =llearn)) +
  geom_tile(aes(fill = Score)) +
  coord_fixed(ratio = 1, xlim = c(0,3), ylim = c(0,3), expand = FALSE)+
  geom_boxplot(data=gpParams, aes(group=environment, x=lteach, y=lambda), width = 0.2, fill=NA, outlier.shape = NA, color='black')+
  #geom_errorbar(data=ldf, aes(x=lteach, ymin=lower, ymax=upper, y= medLambda), fill=NA, width=0.1)+
  geom_point(data=ldf , aes(group=environment,x=lteach, y=meanLambda, shape=environment), fill='black')+
  #scale_fill_gradient(trans = 'log')+
  #scale_fill_distiller(palette = "Spectral", direction = -1,  limits = c(40,105), breaks=seq(40,100,20), labels=seq(40,100,20)) +  
  scale_fill_distiller(palette = "Spectral", direction = -1) +  
  facet_grid(~trial)+
  scale_x_continuous(expression("Teacher"~lambda[0]), expand = c(0, 0), breaks=seq(0.5,2.5,1))+ 
  scale_y_continuous(expression("Student"~lambda[1]), expand = c(0, 0), breaks=seq(0.5,2.5,1))+
  theme_bw()+
  scale_shape_manual(labels = c(expression(hat(lambda)["Rough"]),expression(hat(lambda)["Smooth"])), values=c(2,1))+
  labs(fill="Reward", shape="Human \nEstimates")+
  ggtitle(expression("Experiment 2"))+
  geom_abline(slope=1, intercept=0, linetype = 'dotted', size=1)+
  theme(legend.key.width = unit(0.2, "cm"),
        legend.key.height = unit(0.4, "cm"))+
  theme(text = element_text(size=14, family="sans"))+
  theme(strip.background = element_blank())
p3

setwd('mismatch')
ggsave(plot = p3, filename='Exp2Plot.pdf', height =4.33, width = 7, units = "in", useDingbats=F)
setwd('..')

############################################################################################################
## Experiment 3                                                                                       ##
############################################################################################################
setwd('analysis3')
targetFolder <- 'modelResults/mismatch/'

df3 <- data.frame()
for (i in 1:31){
  filename <- paste0(targetFolder, i, ".csv") #read output file
  if (file.exists(filename)){
    dp<-read.csv(filename)
    dp$avgCumReward <- cumsum(dp$meanReward)/seq_along(dp$meanReward)
    df3 <- rbind(df3, dp)
    
  }
}
#calculate average cumulative reward

gpPars <- read.csv('rationalModels/parameters/gp.csv')

timepoints <- c(5, 10, 20,40)
median <- median(gpPars$lambda)
lowerBound <- quantile(gpPars$lambda, 0.25)
upperBound <- quantile(gpPars$lambda, 0.75)

df3 <- df3 %>% filter(trial %in% timepoints) %>% group_by(lambda, trial) %>% dplyr::select(lambda, trial, avgCumReward, maxReward) %>% melt(id.vars=c("lambda", "trial"))

p4 <- ggplot(df3, aes(x=lambda, y = value, color=variable))+
  geom_vline(xintercept=lowerBound, linetype='dashed')+
  geom_vline(xintercept=upperBound, linetype='dashed')+
  geom_line()+
  facet_grid(~trial) +
  theme_classic() +
  ylab('Reward')+
  scale_fill_manual(values =  c("#7F0000","#00BCE2")) +
  scale_color_manual(values =  c("#7F0000","#00BCE2")) +
  theme(text = element_text(size=14, family="sans"))+
  theme(strip.background = element_blank())
p4
setwd('..')
setwd('mismatch')
ggsave(plot = p4, filename='Exp3Plot.pdf', height =2, width = 7, units = "in", useDingbats=F)

lambdaMax <- ddply(df3, ~trial, summarize, lambda=lambda[which.max(avgCumReward)])

ggplot(df3, aes(x=trial, y = lambda, fill=avgCumReward))+
  geom_tile() +
  geom_point(data=lambdaMax, aes(x=trial, y = lambda), fill='black', shape=4)+
  #scale_fill_gradient(trans = 'log')+
  geom_hline(yintercept=upperBound, linetype='dashed')+
  geom_hline(yintercept=lowerBound, linetype='dashed')+
  scale_fill_distiller(palette = "Spectral", direction = -1) +  
  theme_classic()+
  theme(text = element_text(size=14, family="sans"))+
  theme(strip.background = element_blank())



#Local GP
targetFolder <- 'modelResults/mismatch/local'

df4 <- data.frame()
for (i in 1:31){
  filename <- paste0(targetFolder, i, ".csv") #read output file
  if (file.exists(filename)){
    dp<-read.csv(filename)
    dp$avgCumReward <- cumsum(dp$meanReward)/seq_along(dp$meanReward)
    df4 <- rbind(df4, dp)
    
  }
}


localgpPars <- read.csv('rationalModels/parameters/localgp.csv')
lowerBound <- quantile(localgpPars$lambda, 0.25)
upperBound <- quantile(localgpPars$lambda, 0.75)
timepoints <- c(1, 10, 20,40)

ggplot(subset(df4, trial %in% timepoints), aes(x=lambda, y = avgCumReward))+
  geom_line()+
  facet_grid(~trial) +
  theme_classic()+
  geom_vline(xintercept=lowerBound, linetype='dashed')+
  geom_vline(xintercept=upperBound, linetype='dashed')+
  theme(text = element_text(size=14, family="sans"))+
  theme(strip.background = element_blank())


lambdaMax <- ddply(df4, ~trial, summarize, lambda=lambda[which.max(avgCumReward)])

ggplot(df4, aes(x=trial, y = lambda, fill=avgCumReward))+
  geom_tile() +
  #scale_fill_gradient(trans = 'log')+
  geom_hline(yintercept=upperBound, linetype='dashed')+
  geom_hline(yintercept=lowerBound, linetype='dashed')+
  geom_point(data=lambdaMax, aes(x=trial, y = lambda), fill='black', shape=4)+
  scale_fill_distiller(palette = "Spectral", direction = -1) +  
  theme_classic()+
  theme(text = element_text(size=14, family="sans"))+
  theme(strip.background = element_blank())



