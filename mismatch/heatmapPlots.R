#Mismatch lambda simulations
#Erich Schulz and Charley Wu, Dec 2017

library(scales)
library(plyr)
library(ggplot2)


############################################################################################################
## Generalized case                                                                                       ##
############################################################################################################
#Load data computed using 'generalizedMismatch.py'
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

d<-read.csv("mismatch1d.csv") #load computed using '1dmismatchempirical.R'
L<-expand.grid(seq(0.1,3,0.1), seq(0.1,3,0.1))
d$lteach<-rep(L[,1], each=11)
d$llearn<-rep(L[,2], each=11)
dp2<-subset(d, trial %in% c(1,3,5,10))
dp2$m<- rescale(dp2$m, from=c(-1.8, 1.8), to=c(0,100))
dp2$trial<-mapvalues(dp2$trial, c(1,3,5,10), c("t=1","t=3","t=5","t=10"))
dp2$trial<-factor(dp2$trial,levels= c("t=1","t=3","t=5","t=10"))
dp2$Score<-dp2$m
dp2$X<-NULL
dp2$m<-NULL


#Human estimates
dl2<-subset(dp2, lteach %in% c(1,2))
dl2$llearn<-ifelse(dl2$lteach==1,0.78, 0.82)

p2<-ggplot(dp2, aes(x = lteach, y =llearn, fill = Score)) +
  geom_tile() +
  coord_equal()+
  geom_point(data=dl2, aes(shape=factor(lteach)), col="black", fill=NA, size=2) +
  #annotate("text", label = "hat(lambda)[rough]", x = 1, y = 0.5, color = "black", parse = TRUE)+
  #scale_fill_gradient(trans = 'log')+
  scale_fill_distiller(palette = "Spectral", direction = -1, limits = c(40,105), breaks=seq(40,100,20), labels=seq(40,100,20))+
  facet_grid(~trial)+
  scale_x_continuous(expression("Teacher"~lambda[0]), expand = c(0, 0), breaks=seq(0.5,2.5,1))+ 
  scale_y_continuous(expression("Student"~lambda[1]), expand = c(0, 0),  breaks=seq(0.5,2.5,1))+
  theme_bw()+
  scale_shape_manual(labels = c(expression(hat(lambda)["Rough"]),expression(hat(lambda)["Smooth"])), values=c(2,1))+
  labs(fill="Reward", shape="Human \nEstimates")+
  ggtitle(expression("Experiment 1"))+
  geom_abline(slope=1, intercept=0, linetype = 'dotted', size=1.5)+
  theme(legend.key.width = unit(0.2, "cm"),
        legend.key.height = unit(0.4, "cm"))+
  theme(text = element_text(size=14, family="sans"))+
  theme(strip.background = element_blank())
p2

ggsave(plot = p2, filename='Exp1Plot.pdf', height =4.33, width = 7, units = "in", useDingbats=F)


############################################################################################################
## Experiment 2                                                                                        ##
############################################################################################################

d<-read.csv("mismatch2d.csv") #load computed using '1dmismatchempirical.R'
d$lteach<-rep(L[,1], each=40)
d$llearn<-rep(L[,2], each=40)
dp3<-subset(d, trial %in% c(1,3,5,40))
#dp3<-subset(d, trial %in% c(1,5,20,40))
dp3$m<-rescale(dp3$m, from=c(-1.8, 1.8), to=c(0,100))
dp3$trial<-mapvalues(dp3$trial, c(1,3,5,40), c("t=1","t=3","t=5","t=40"))
#dp3$trial<-mapvalues(dp3$trial, c(1,5,20,40), c("t=1","t=5","t=20","t=40"))
dp3$trial<-factor(dp3$trial,levels= c("t=1","t=3","t=5","t=40"))
#dp3$trial<-factor(dp3$trial,levels= c("t=1","t=5","t=20","t=40"))
dp3$Score<-dp3$m
dp3$X<-NULL
dp3$m<-NULL
dp3$sim<-"Empirical-2D"


dl3<-subset(dp3, lteach %in% c(1,2))
dl3$llearn<-ifelse(dl3$lteach==1,0.78, 0.92)


p3<-ggplot(dp3, aes(x = lteach, y =llearn, fill = Score)) +
  geom_tile() +
  coord_equal()+
  geom_point(data=dl3, aes(shape=factor(lteach)), col="black", fill=NA, size=2) +
  #scale_fill_gradient(trans = 'log')+
  scale_fill_distiller(palette = "Spectral", direction = -1,  limits = c(40,105), breaks=seq(40,100,20), labels=seq(40,100,20)) +  
  facet_grid(~trial)+
  scale_x_continuous(expression("Teacher"~lambda[0]), expand = c(0, 0), breaks=seq(0.5,2.5,1))+ 
  scale_y_continuous(expression("Student"~lambda[1]), expand = c(0, 0), breaks=seq(0.5,2.5,1))+
  theme_bw()+
  scale_shape_manual(labels = c(expression(hat(lambda)["Rough"]),expression(hat(lambda)["Smooth"])), values=c(2,1))+
  labs(fill="Reward", shape="Human \nEstimates")+
  ggtitle(expression("Experiment 2"))+
  geom_abline(slope=1, intercept=0, linetype = 'dotted', size=1.5)+
  theme(legend.key.width = unit(0.2, "cm"),
        legend.key.height = unit(0.4, "cm"))+
  theme(text = element_text(size=14, family="sans"))+
  theme(strip.background = element_blank())
p3

ggsave(plot = p3, filename='Exp2Plot.pdf', height =4.33, width = 7, units = "in", useDingbats=F)



library(gridExtra)
pdf("mismatchgrid.pdf", useDingbats=F)
grid.arrange(p2,p3,p1, nrow=3)
dev.off()
