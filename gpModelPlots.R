#Script to plot model fitting results from GPs.R (negative log likelihoods)
rm(list=ls())

#read in data
dat<-read.csv("newdat.csv")
#Model parameter data
dataFolder <- "output4/"
#load packages
packages <- c('plyr', 'ggplot2', 'jsonlite', 'MASS', 'gridExtra', 'parallel', "reshape2")
lapply(packages, require, character.only = TRUE)

#dummy data frame
data<-data.frame(id=numeric(), trial=numeric(), x=numeric(), y=numeric(), 
                 z=numeric(), kernel=numeric(), scenario=numeric(), 
                 horizon=numeric(), round=numeric())

#loop through data
for (i in 1:nrow(dat)){
  #write to json file
  write(paste(dat$searchHistory[i]), file = "json/data.JSON")
  #load from json file just written
  dfinal <- fromJSON(txt="json/data.JSON")
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
  #round, 1:8
  round<-rep(1:8, len)
  #horizon, 20 or 40
  horizon<-rep(len-1, len)
  #maximum value or cumulative
  scenario<-rep(dat$scenario[i], sum(len))
  #smoothness
  kernel<-rep(dat$kernel[i], sum(len))
  #id number
  id<-rep(i, sum(len))
  #dummy frame
  dummy<-data.frame(id, trial, x, y, z, kernel, scenario, horizon, round)
  #bind them together
  data<-rbind(data, dummy)
}

allopts<-expand.grid(0:10, 0:10)
data$chosen<--99
for (i in 1:nrow(data)){
  data$chosen[i]<-which(data$x[i]==allopts$Var1 & data$y[i]==allopts$Var2)
}
data$overall<-1:nrow(data)-with(data, ave(id, paste(id, round), FUN = seq_along))
data<-transform(data,overall=as.numeric(factor(overall)))



#kernels <-  c("RBF", "Oru") #GP kernels
kernels <-  c("RBF") #GP kernels
# acqFuncs <- c("UCB", "Thompson", "ProbOfImp","ExpectedImp") #Acquisition functions
acqFuncs <- c("UCB", "ProbOfImp") #Acquisition functions

modelFit <- data.frame(participant=numeric(), reward=numeric(), environment=numeric(), nLL=numeric(),kernel=numeric(), acq=numeric(), lambda=numeric(), noise=numeric(), tau=numeric(), beta_s=numeric(), beta_l=numeric()) 

#loop through data
for (k in kernels){
	for (a in acqFuncs){
		for (i in 1:80){ #subjects
			filename <- paste0(dataFolder, k, a, i, ".csv") #read output file
			if (file.exists(filename)){
			  dp<-read.csv(filename)
			  dummy<-subset(data, id==i) #subset participant in the dataframe
			  if (dummy$kernel[1] == 0 ){ #Environment 
			    environment <- "Rough"
			  } else { environment <- "Smooth"}
			  if (dummy$scenario[1] == 0){
			    reward <- "TotalScore"
			  }else{ reward <- "BestReward"}
			  nLL<-dp$nLL #retrieve AIC from the .csv file data
			  participant<-dummy$id[i]  #replicate ID
			  kernel<-k
			  acq<-a
			  lambda <- exp(dp[4])
			  noise <- exp(dp[5])
			  tau <- exp(dp[length(dp)])
			  if (length(dp)>6){
			    beta_s <- exp(dp[6])
			    beta_l <- exp(dp[7])
			  }else {
			    beta_s <- 0
			    beta_l <- 0
			  }
			  dadd<-data.frame(participant=participant,reward=reward, environment=environment, nLL=nLL, kernel=kernel, acq = acq, lambda=lambda, noise=noise, tau=tau, beta_s=beta_s, beta_l=beta_l)
			  names(dadd) <- names(modelFit)
			  modelFit<-rbind(modelFit, dadd)
			}
		}
	}
}

#BIC,AIC, AICCc
ucb <- subset(modelFit, acq=="UCB")
ucb$AIC <- (ucb$nLL * 2) + (2 * 5)
ucb$AICc <- (ucb$nLL * 2) + (((2 * 5) * (5+1))/(240-5-1))
ucb$BIC <- (ucb$nLL * 2) + (5 * log(240))

probimp <- subset(modelFit, acq=="ProbOfImp")
probimp$AIC <- (probimp$nLL * 2) +(2*3)
probimp$AICc <- (((2 * 3) * (3+1))/(240-3-1))
probimp$BIC <- (probimp$nLL * 2) + (3 * log(240))
modelFit <- rbind(ucb, probimp)

#RANDOM MODEL
randomModel <- -2 * sum(rep(log(1/121),240))

#McFadden R^2
UCB_r2 <- 1 - (mean(ucb$BIC)/randomModel)
probimp_r2 <- 1 - (mean(probimp$BIC)/randomModel)
##Model comparison plot

p1 <- ggplot(modelFit, aes(x=environment, y=BIC, fill = acq)) + 
  geom_bar(position="dodge", stat="identity") + 
  facet_wrap( ~ reward) + geom_hline(yintercept=randomModel)
    
print(p1)
#Count of which model best explains which subject

#Subsets
#ucbSub <- ucb
#probimpSub <- probimp
ucbSub <- subset(ucb, environment=="Rough")
probimpSub <- subset(probimp, environment=="Rough")

ucbCount <- 0
probImpCount <- 0
for (id in probimpSub$participant){
  if (match(id, ucbSub$participant, nomatch=FALSE) && match(id, probimpSub$participant, nomatch=FALSE)){
    if (ucbSub$BIC[ucbSub$participant == id] > probimpSub$BIC[probimpSub$participant==id]){
      probImpCount <- probImpCount + 1
    } else { ucbCount <- ucbCount + 1}
  }
}
barplot(height = c(ucbCount/(ucbCount+probImpCount) ,probImpCount/(ucbCount+probImpCount) ),  names.arg=c("UCB", "ProbImp"), ylab=("Prop. Subjects Best Explained"))


##Parameter estimates
#LAMBDA length scale
p2<- ggplot(modelFit, aes(y=lambda, x=environment,  fill=environment)) +
  geom_violin() +
  geom_boxplot(width=0.1) +
  xlab("Environment") +
  ylab("Smoothness (lambda)")+ 
  ggtitle("Lambda (How Smooth?)") +
  theme_bw()+
  theme(legend.position = "none")

print(p2)

ggsave(p2, filename = "plots/lambda.pdf", height=6, width=5, units="in" )

#NOISE
p3<- ggplot(modelFit, aes(y=noise, x=acq, color=acq, shape = acq)) +
  geom_jitter() +
  facet_wrap(~environment)+
  xlab("Model") +
  ylab("Noise Variance")+ 
  ggtitle("NoiseSigma (How noisey?)") +
  theme_bw()+
  theme(legend.position = "none")
print(p3)
ggsave(p3, filename = "plots/noise.pdf", height=6, width=5, units="in" )

#TAU
p4<- ggplot(modelFit, aes(y=tau, x=acq, fill=acq, color=acq, shape=acq)) +
  geom_jitter() +
  facet_wrap(~environment)+
  xlab("Model") +
  ylab("Tau")+ 
  ggtitle("Tau (How Predictable?)") +
  theme_bw()+
  theme(legend.position = "none")
print(p4)
ggsave(p4, filename = "plots/tau.pdf", height=6, width=5, units="in" )
#BETA
betaDF<- melt(data = subset(modelFit,acq=="UCB"), id.vars = c("participant", "reward", "environment", "kernel", "acq"), measure.vars = c("beta_s", "beta_l"), variable.name = "Parameter", value.name = "Beta")
p5<- ggplot(betaDF, aes(y=Beta, x=Parameter, color=Parameter, shape = Parameter)) +
  geom_jitter() +
  facet_wrap( ~ reward) +
  xlab("Horizon") +
  scale_x_discrete(labels=c("Short", "Long"))+ 
  ylab("Beta")+ 
  ggtitle("Beta (How Much To Explore?)") +
  theme_bw()+
  theme(legend.position = "none")
print(p5)
ggsave(p5, filename = "plots/beta.pdf", height=6, width=5, units="in" )

modelFit$deltaBeta <- modelFit$beta_l - modelFit$beta_s
p6<- ggplot(modelFit, aes(y=deltaBeta, x=environment, color=environment, shape = environment)) +
  geom_jitter() +
  facet_wrap( ~ reward)
print(p6)

#Model fit
#Factors
p7<- ggplot(modelFit, aes(y=BIC, x=acq, color=acq, shape=acq)) +
  geom_jitter() + 
  geom_hline(yintercept = randomModel) + facet_wrap(environment ~ reward)
print(p7)


#GENERAL
p8<- ggplot(modelFit, aes(y=BIC, x=acq, color=acq, shape=acq)) +
  geom_jitter() + 
  geom_hline(yintercept = randomModel, linetype = "dotdash" )+
  xlab("Model") +
  ylab("BIC")+ 
  theme_bw()+
  theme(legend.position = "none") +
  ylim(c(500, 2500))
print(p8)
ggsave(p8, filename = "plots/BIC.pdf", height=4, width = 6, units = "in")

#Bayes Factor
 
ggsave(p9, filename = "plots/BayesFactor.pdf", height=4, width = 6, units = "in")






#Subset of dataframe with short horizon
daic20<-subset(daic, horizon==20)
aic_random = -log(1/121) * 20 * 4 #log loss of a random model
daic20$cox<-1-daic20$aic/aic_random
dplot20<-ddply(daic20,~horizon+acq+kernel+reward,summarise,mu=mean(cox))

#subset of dataframe with long horizon
daic40<-subset(daic, horizon==40)
aic_random = -log(1/121) * 40 * 4 #log loss of a random model
daic40$cox<-1-daic40$aic/aic_random
dplot40<-ddply(daic40,~horizon+acq+kernel+reward,summarise,mu=mean(cox))

dplot<-rbind(dplot20, dplot40)
dplot$Kernel<-dplot$kernel

dplot$model<-paste0(dplot$kernel, "\n", dplot$acq)
dplot$horizon<-ifelse(dplot$horizon==20, "Short", "Long")
dplot$reward<-ifelse(dplot$reward==0, "Average Reward", "Max Reward")

p <- ggplot(dplot, aes(y=mu, x=acq, fill=Kernel)) + 
  #bars
  geom_bar(position="dodge", stat="identity")+
  #golden ratio error bars
  ggtitle("Model Comparison")+theme_classic() +xlab("\nModel")+ylab("McFadden r2\n")+
  scale_fill_manual(values=c("#A20000", "#0082C1"))+
  #scale_y_continuous(breaks = round(seq(0, 1, 0.5), 1), limits=c(-0.1,1.2))+
  #geom_text(aes(label=mu), position=position_dodge(width=0.9), vjust=-0.25)+
  #fill with Wes Anderson colors  scale_fill_manual(values = c("#A20000", "#0082C1"))

  facet_grid(reward ~ horizon)+
  theme(legend.position="top",
        text = element_text(size=16, family="serif"),
        plot.title = element_text(size=26))
pdf("comparison.pdf", width=12, height=8)  
p
dev.off()
