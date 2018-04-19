#Charley Wu, March 2018
#Interpret Experiment 3 modeling results

#house keeping
rm(list=ls())
theme_set(theme_bw(base_size=16))# use the b&w theme

#load packages
packages <- c('plyr', 'jsonlite', 'lsr', 'ggplot2', 'gridExtra', 'reshape2', 'stargazer', 'coefplot', "grid", 'corrplot')
lapply(packages, require, character.only = TRUE)

#############################################################################################################################
# IMPORT PARTICIPANT DATA 
#############################################################################################################################

#Participant data
data<-read.csv("ExperimentData/exp3.csv",  sep=";")

#add factors
data$scenario <- as.factor(data$scenario)
levels(data$scenario) <- c("Cumulative", "Best")
data$horizon <- as.factor(data$horizon)

#Scale variables
data$meanScale<- sapply(as.character(data$scale), FUN = function(x) mean(fromJSON(x)))
data$varScale<- sapply(as.character(data$scale), FUN = function(x) var(fromJSON(x)))

#Sync ids to trial data
data$id <- seq(1,80)

#############################################################################################################################
# IMPORT MODELING RESULTS
#############################################################################################################################

importModelResults <- function(dataFolder, kernels, acqFuncs){
  #initialize data frames
  modelFit <- data.frame(participant=numeric(), reward=numeric(), nLL=numeric(), kernel=numeric(), acq=numeric(), R2=numeric(), bonus=numeric(), kError=numeric(), lambda=numeric(), beta=numeric(), tau=numeric()) 
  paramEstimates <- data.frame(participant=numeric(), reward=numeric(),horizon=numeric(),leaveoutindex=numeric(), nLL=numeric(),R2=numeric(), kernel=numeric(), acq=numeric(), kError=numeric(), lambda=numeric(), beta = numeric(), tau=numeric(), roundnLL=numeric(), bonus=numeric())
  #loop through data
  for (k in kernels){
    for (a in acqFuncs){
      for (i in 1:80){ #subjects
        filename <- paste0(dataFolder, k, a, i, ".csv") #read output file
        if (file.exists(filename)){
          dp<-read.csv(filename)
          #interpret parameters based on mode combination
          if (k=="IMD"| k=="LocalSearch"| k=="WSLS" | k=="LWSLS"){#inverse manhattan distance aka Inertia
            colnames(dp) <- c("", "Horizon", "leaveoutindex", "nLL", "tau")
          }else if (k=="BMT" | k=="LBMT"){#Bayesian mean tracker
            ifelse(a=='UCB', colnames(dp)<- c("", "Horizon", "leaveoutindex", "nLL", "kError", "beta","tau"), colnames(dp)<- c("", "Horizon", "leaveoutindex", "nLL","kError", "tau"))
          }else if (k=="LIN" | k=="LLIN"){ #linear kernel
            ifelse(a=='UCB', colnames(dp)<- c("", "Horizon", "leaveoutindex", "nLL", "beta","tau"), colnames(dp)<- c("", "Horizon", "leaveoutindex", "nLL", "tau"))
          }else { #normal GP kernel
            ifelse(a=='UCB', colnames(dp)<- c("", "Horizon", "leaveoutindex", "nLL", "lambda", "beta","tau"), colnames(dp)<- c("", "Horizon", "leaveoutindex", "nLL", "lambda", "tau"))
          }
          #demographics
          dummy<-subset(data, id==i) #subset participant in the dataframe
          reward <- dummy$scenario
          participant<-dummy$id  #replicate ID
          kernel<-k
          acq<-a
          bonus<- dummy$reward
          #Total out of sample nLL
          nLL <- sum(dp$nLL)
          randomModel240 <- -log(1/121)*240
          R2 <- 1 - (nLL/randomModel240)
          #blank mean parameter estimates
          kErrorMean <- NA
          lambdaMean <- NA
          betaMean <- NA
          tauMean <- median(exp(dp$tau)) #all models should have tau
          #mean parameter estimates for UCB RBF
          if (a=="UCB" & k != "IMD"){ #UCB when not used for IMD should have beta
            betaMean <- median(exp(dp$beta))
          }
          if (k=="RBF" | k=="LRBF"){
            lambdaMean <- median(exp(dp$lambda))
          }
          if (k=="BMT" | k=="LBMT"){ #BMT
            kErrorMean <- median(exp(dp$kError))
          }
          #save modelFit
          dadd<-data.frame(participant=participant,reward=reward, nLL=nLL, kernel=kernel, acq = acq, R2 =R2, bonus = bonus, kError = kErrorMean, lambda=lambdaMean, beta = betaMean, tau=tauMean)
          names(dadd) <- names(modelFit)
          modelFit<-rbind(modelFit, dadd)
          #Loop through short and long horizons to save parameter estimates
          for (h in c(20,40)){
            #loop through leave out index to save all 8 parameter estimates for each subject
            for (j in c(1,2,3,4)){
              subDP <- subset(dp, Horizon == h & leaveoutindex == j)
              roundnLL <- subDP$nLL
              #exponentiation of all parameters
              kError <- ifelse("kError" %in% colnames(subDP), exp(subDP$kError), NA)
              lambda <- ifelse("lambda" %in% colnames(subDP), exp(subDP$lambda), NA)
              beta <- ifelse("beta" %in% colnames(subDP), exp(subDP$beta),  NA)
              tau <- exp(subDP$tau)
              dadd<-data.frame(participant=participant,reward=reward, horizion=h, leaveoutindex =j, nLL=nLL, R2 =R2,  kernel=kernel, acq = acq, kError = kError, lambda=lambda, beta=beta,  tau=tau, roundnLL=roundnLL, bonus=bonus)
              names(dadd) <- names(paramEstimates)
              paramEstimates<-rbind(paramEstimates, dadd)
            }}}}}}
  return(list(modelFit, paramEstimates))
}

#############################################################################################################################
# COMPILE MODELING RESULTS FROM POTENTIALLY MULTIPLE SOURCES
# Load all models so far:
#modelFit <- read.csv('modelResults/modelFit.csv')
#paramEstimates <- read.csv('modelResults/paramEstimates.csv')
#############################################################################################################################
modelResults <- importModelResults('modelResults/first/', c("LocalSearch",  "WSLS", "LWSLS",  "BMT",  "LBMT", "RBF", "LRBF"), c("", "GV", "GM", "UCB", "EXI", "POI", "PMU"))
#modelResults <- importModelResults('modelResults/recoveryLocalGP/', c("LocalSearch",  "WSLS", "LWSLS",  "BMT",  "LBMT", "RBF", "LRBF"), c("", "GV", "GM", "UCB", "EXI", "POI", "PMU"))
#modelResults <- importModelResults('modelResults/recoveryLocalBMT/', c("LocalSearch",  "WSLS", "LWSLS",  "BMT",  "LBMT", "RBF", "LRBF"), c("", "GV", "GM", "UCB", "EXI", "POI", "PMU"))

#separate into overall per participant and individual cross validation blocks
modelFit <- modelResults[[1]]
paramEstimates <- modelResults[[2]]

#reorder acqisition function levels and add "ModelName" to identify unique models
modelFit$acq <-factor(modelFit$acq, levels = c("", "GV", "GM", "UCB", "EXI", "POI", "PMU"))
modelFit$ModelName <- paste(modelFit$kernel,modelFit$acq, sep = "-")
modelFit$ModelName <- factor(modelFit$ModelName)


#SAVE MODEL RESULTS TO DISK
#write.csv(modelFit,'modelResults/modelFit.csv')
#write.csv(paramEstimates,'modelResults/paramEstimates.csv')

#write.csv(modelFit,'modelResults/localGPUCBrecovery.csv')
#write.csv(paramEstimates,'modelResults/localGPUCBparamEstimates.csv')

#write.csv(modelFit,'modelResults/localBMTUCBrecovery.csv')
#write.csv(paramEstimates,'modelResults/localBMTUCBparamEstimates.csv')
#############################################################################################################################
# SAVE PARAMETER ESTIMATES FOR RATIONAL MODELS
#############################################################################################################################

#localPars <- subset(paramEstimates, kernel=="LocalSearch")
#write.csv(localPars[,c("participant", "leaveoutindex", "reward","horizon", "tau")], file = "rationalModels/parameters/local.csv", row.names = FALSE)

#gpUCBpars <- subset(paramEstimates, kernel=="RBF" & acq=="UCB")
#write.csv(gpUCBpars[,c("participant", "leaveoutindex","reward","horizon","lambda", "beta", "tau")], file = "rationalModels/parameters/gp.csv", row.names = FALSE)

#localgpUCBpars <- subset(paramEstimates, kernel=="LRBF" & acq=="UCB")
#write.csv(localgpUCBpars[,c("participant", "leaveoutindex", "reward","horizon","lambda", "beta", "tau")], file = "rationalModels/parameters/localgp.csv", row.names = FALSE)

#BMTpars <- subset(paramEstimates, kernel=="BMT" & acq=="UCB")
#write.csv(BMTpars[,c("participant", "leaveoutindex", "reward","horizon","kError", "beta", "tau")], file = "rationalModels/parameters/BMT.csv", row.names = FALSE)

#localBMTpars <- subset(paramEstimates, kernel=="LBMT" & acq=="UCB")
#write.csv(localBMTpars[,c("participant", "leaveoutindex", "reward","horizon","kError", "beta", "tau")], file = "rationalModels/parameters/localBMT.csv", row.names = FALSE)

#############################################################################################################################
# Summary of Modeling Results
#############################################################################################################################
#ALL MODELS

#Number of participants best described
models <-  rep(0, length(levels(modelFit$ModelName)))
names(models) <- levels(modelFit$ModelName)
for (pid in 1:80){
  subDF <- subset(modelFit, participant==pid)
  best <- subDF$ModelName[which(subDF$R2==max(subDF$R2))]
  models[best] <- models[best] + 1
}
#add best described to modelFit 
modelFit$bestDescribed <- NA
for (mod in names(models)){
  modelFit[modelFit$ModelName == mod,]$bestDescribed <- models[mod]
}


#Summary of results
cat("All conditions")
for (k in levels(modelFit$kernel)){
  for (a in levels(modelFit$acq)){
    if (paste(k,a,sep="-") %in% names(models)){
      sub <- subset(modelFit, acq == a & kernel == k) #sub of model fit
      psub <- subset(paramEstimates, acq ==a & kernel==k)
      randomModel <- -log(1/121) * 240
      r_2 <- 1 - (mean(sub$nLL)/randomModel)
      bestDescribed <- models[paste(k,a,sep="-")]
      cat( k,"-",a, " (n=",nrow(sub),")", "; R2 = ", round(r_2, digits=2),"; Best Described = ", bestDescribed , "; Lambda = ", round(median(psub$lambda), digits=2),"; Beta = ", round(median(psub$beta), digits=2), "; kError = ", round(median(psub$kError), digits=2), "; Tau = ", round(median(psub$tau), digits=2), "\n", sep="")
    }
  }
}


#MAIN TEXT COMPARISON
modelList <- c("BMT-GV", "BMT-GM", "BMT-UCB", "LBMT-GV", "LBMT-GM", "LBMT-UCB","RBF-GV", "RBF-GM", "RBF-UCB", "LRBF-GV", "LRBF-GM", "LRBF-UCB") #subset models
mainTextModels <- modelFit[modelFit$ModelName %in% modelList,]
mainTextModels$ModelName <- factor(mainTextModels$ModelName) #refactor levels of subset
#create a vector to count number best described
modelSub <-  rep(0, length(levels(mainTextModels$ModelName)))
names(modelSub) <- levels(mainTextModels$ModelName)
for (pid in 1:80){
  subDF <- subset(mainTextModels, participant==pid)
  best <- subDF$ModelName[which(subDF$R2==max(subDF$R2))]
  modelSub[best] <- modelSub[best] + 1
}
modelSub
#add best described to modelFit 
modelFit$mainTextBest <- NA
for (mod in names(modelSub)){
  modelFit[modelFit$ModelName == mod,]$mainTextBest <- modelSub[mod] 
}

#SAVE MAINTEXT RESULTS TO DISK
#write.csv(mainTextModels,'modelResults/maintext.csv')

#############################################################################################################################
# Model Comparison Plot
#############################################################################################################################
#MAIN TEXT
mainTextModels <- modelFit[modelFit$ModelName %in% modelList,]
mainTextModels$ModelName <- factor(mainTextModels$ModelName) #refactor levels of subset'
mainTextModels$kernel <- factor(mainTextModels$kernel, levels = c("BMT", "LBMT", "RBF", "LRBF"))
mainTextModels$acq <- factor(mainTextModels$acq, levels=c("UCB", "GM", "GV"))
se<-function(x){sd(x)/sqrt(length(x))}
mtmDF <- ddply(mainTextModels, ~kernel+acq, summarise, newR2 =mean(R2), se=se(R2), best=mean(mainTextBest))
boxplotDF <- ddply(mainTextModels, ~kernel+acq, summarise,d_ymin = max(min(R2), quantile(R2, 0.25) - 1.5 * IQR(R2)), d_ymax = min(max(R2), quantile(R2, 0.75) + 1.5 * IQR(R2)),
                   d_lower = quantile(R2, 0.25),  d_middle = median(R2), d_upper = quantile(R2, 0.75),
                   mu=mean(R2))

main <- ggplot(boxplotDF) +
  #stat_summary(fun.y = mean, geom = "bar", position = "dodge", color='black') + 
  geom_boxplot(aes(x =as.numeric(kernel)-0.2, ymin = d_lower, ymax = d_upper, lower = d_lower, middle = d_middle, upper = d_upper, width = 2 * 0.2, fill = kernel), stat = "identity", color='black') +
  #whiskers
  geom_segment(aes(x = as.numeric(kernel), y = d_ymin, xend = as.numeric(kernel), yend = d_ymax)) +
  geom_segment(aes(x = as.numeric(kernel) - 0.1,  y = d_ymax - 0.0025, xend = as.numeric(kernel),  yend = d_ymax -0.0025)) + #top
  geom_segment(aes(x = as.numeric(kernel) - 0.1, y = d_ymin +0.0025, xend = as.numeric(kernel), yend = d_ymin + 0.0025)) + #bottom
  geom_point(aes(x = as.numeric(kernel)-0.2, y = mu), size = 1.3, shape=23, fill='white') +
  geom_jitter(data=mainTextModels, aes(x = as.numeric(kernel)+.2,  y = R2,  color = kernel), 
              width = 0.2 - 0.25 * 0.2, height = 0, size=0.4, alpha = 0.4)+
  facet_grid(~acq)+
  #geom_jitter(data = mainTextModels, aes(x=acq, y = R2, fill= kernel),  color='grey', size = 1, shape=21, alpha=0.2, position=position_jitterdodge(dodge.width=0.9, jitter.width = 0.2))+
  #geom_errorbar(aes(ymin=newR2 - se, ymax=newR2 + se),color='black', width = .4, position=position_dodge((width=0.9))) +
  #scale_fill_manual(values = c("#7F0000", "#00BCE2",  "#37045F" ))+
  scale_fill_manual(values = c( "#F0E442", "#E69F00", "#009E73", "#56B4E9"))+
  scale_color_manual(values = c( "#F0E442", "#E69F00", "#009E73", "#56B4E9"))+
  ylab("Predictive Accuracy")+ 
  coord_cartesian(ylim=c(-0.22, 0.8))+
  theme_classic()+
  theme(text = element_text(size=12,  family="sans"),axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), 
        strip.background=element_blank(),
        strip.text = element_blank(),
        legend.key=element_rect(color=NA),
        panel.spacing.x=unit(0.2, "lines"),
        panel.spacing.y=unit(1, "lines"),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        legend.position="none")+
  guides(color=FALSE, shape=FALSE)+
  ggtitle("Experiment 3")
main
ggsave(filename = "plots/mainTextR2.pdf", main,height =2.5, width = 4, units ="in")
#FULL MODEL COMPARISON IN SI
#Model order
modelFit$ModelName <- factor(modelFit$ModelName, levels = c("BMT-UCB", "BMT-GM", "BMT-GV", "BMT-EXI", "BMT-POI", "BMT-PMU", "LBMT-UCB", "LBMT-GM", "LBMT-GV", "LBMT-EXI", "LBMT-POI", "LBMT-PMU", "RBF-UCB", "RBF-GM", "RBF-GV", "RBF-EXI", "RBF-POI", "RBF-PMU", "LRBF-UCB", "LRBF-GM", "LRBF-GV", "LRBF-EXI", "LRBF-POI", "LRBF-PMU", "WSLS-", "LWSLS-", "LocalSearch-"))
modelFit$reward <- factor(modelFit$reward, levels=c("Cumulative", "Best"))
 
#Recode simple models
modelFit$acq <- as.character(modelFit$acq) #unfactor for now
modelFit$kernel <- as.character(modelFit$kernel)
modelFit[modelFit$kernel=="WSLS",]$acq <- rep("WSLS", 80) 
modelFit[modelFit$kernel=="WSLS",]$kernel <- rep("Simple Strategies", 80)
modelFit[modelFit$kernel=="LWSLS",]$acq <- rep("WSLS*", 80) 
modelFit[modelFit$kernel=="LWSLS",]$kernel <- rep("Simple Strategies", 80)
modelFit[modelFit$kernel=="LocalSearch",]$acq <- rep("Local Search", 80) 
modelFit[modelFit$kernel=="LocalSearch",]$kernel <- rep("Simple Strategies", 80)

#Refactor
modelFit$acq <- factor(modelFit$acq, levels = c("UCB", "GM", "GV", "EXI", "POI", "PMU", "WSLS", "WSLS*", "Local Search"))
levels(modelFit$acq) <- c("UCB", "Exploit", "Explore",  "EXI", "POI", "PMU", "WSLS", "WSLS*", "Local Search") 
modelFit$kernel <- factor(modelFit$kernel, levels = c("BMT", "LBMT", "RBF", "LRBF", "Simple Strategies"))
levels(modelFit$kernel) <- c("Mean Tracker", "Mean Tracker*", "Function Learning", "Function Learning*", "Simple Strategies")

full <- ggplot(modelFit , aes(y=R2, x=acq, fill = reward)) + #fill=interaction(environment, reward)
  stat_summary(fun.y = mean, geom = "bar", position = "dodge") + 
  geom_jitter(color='grey', shape=21, alpha=0.4, size = .8, position=position_jitterdodge(dodge.width=0.9, jitter.width=0.2))+
  stat_summary(fun.data = mean_se, geom = "errorbar", position = position_dodge(width = 0.90), width = 0.2 ) +
  #geom_boxplot(position="dodge")+
  #geom_text(aes(label=bestDescribed, y = 0.5), family="serif") +
  xlab("") +
  #scale_fill_manual(values = c("#7F0000", "#DA4233" , "#005F8D", "#00BCE2"))+
  scale_fill_manual(values = c( "#DA4233" , "#00BCE2"))+
  ylab("Predictive Accuracy")+ 
  scale_y_continuous(breaks = c(0, .2, .4, .6, .8))+
  theme_classic()+
  theme(text = element_text(size=16,  family="sans"), strip.background=element_blank(), legend.key=element_rect(color=NA), legend.position="none")+
  coord_cartesian(ylim=c(-.25, 0.8)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.title=element_blank())+
  ggtitle("Experiment 3") +
  guides(shape=FALSE) +
  facet_grid(~kernel,  scales="free_x", space="free_x")
full

ggsave(filename = "plots/fullComparison2.pdf", plot = full, height =5, width = 12, units = "in") 

#RANDOM MODEL BIC
randomModel240 <- -log(1/121)*240 #for full task

##Model comparison plots
ggplot(modelFit, aes(y=nLL, x=acq, color=acq, shape=environment)) +
  geom_jitter() + 
  geom_hline(yintercept = randomModel240, linetype = "dotdash" )+
  xlab("Model") +
  ylab(" Out of Sample nLL")+ 
  theme_bw()+
  theme(legend.position = "right") +
  #ylim(c(500, 2500)) +
  facet_wrap(~kernel) +
  ggtitle("Model Comparison")


######tests########
#model comparison
t.test(subset(modelFit, ModelName=="BMT-UCB")$R2, subset(modelFit, ModelName=="RBF-UCB")$R2, paired=TRUE)
cohensD(subset(modelFit, ModelName=="BMT-UCB")$R2, subset(modelFit, ModelName=="RBF-UCB")$R2)

t.test(subset(modelFit, ModelName=="LBMT-UCB")$R2, subset(modelFit, ModelName=="LRBF-UCB")$R2, paired = TRUE)
wilcox.test(subset(modelFit, ModelName=="LBMT-UCB")$R2, subset(modelFit, ModelName=="LRBF-UCB")$R2, paired=TRUE)
cohensD(subset(modelFit, ModelName=="LBMT-UCB")$R2, subset(modelFit, ModelName=="LRBF-UCB")$R2)

t.test(subset(modelFit, ModelName=="LBMT-UCB")$R2, subset(modelFit, ModelName=="BMT-UCB")$R2, paired=TRUE)
cohensD(subset(modelFit, ModelName=="LBMT-UCB")$R2, subset(modelFit, ModelName=="BMT-UCB")$R2)

t.test(subset(modelFit, ModelName=="LRBF-UCB")$R2, subset(modelFit, ModelName=="RBF-UCB")$R2, paired=TRUE)
cohensD(subset(modelFit, ModelName=="LRBF-UCB")$R2, subset(modelFit, ModelName=="RBF-UCB")$R2)

t.test(subset(modelFit, ModelName=="LRBF-UCB")$R2, subset(modelFit, ModelName=="LRBF-GM")$R2, paired=TRUE)
cohensD(subset(modelFit, ModelName=="LRBF-UCB")$R2, subset(modelFit, ModelName=="LRBF-GM")$R2)

t.test(subset(modelFit, ModelName=="LRBF-UCB")$R2, subset(modelFit, ModelName=="LRBF-GV")$R2, paired=TRUE)
cohensD(subset(modelFit, ModelName=="LRBF-UCB")$R2, subset(modelFit, ModelName=="LRBF-GV")$R2)
#separated by 
t.test(subset(modelFit, acq=="UCB" & kernel=="Local Search")$R2, subset(modelFit, acq=="UCB" & kernel=="RBF")$R2, paired=TRUE)
t.test(subset(modelFit, acq=="UCB" & kernel=="LRBF")$R2, subset(modelFit, acq=="UCB" & kernel=="RBF")$R2, var.equal=TRUE)
t.test(subset(modelFit, acq=="UCB" & kernel=="Local GP")$R2, subset(modelFit, acq=="UCB" & kernel=="Local Search")$R2, var.equal=TRUE)

#Payoff conditions
t.test(subset(modelFit, reward=="Cumulative" &  ModelName=="LBMT-UCB")$R2, subset(modelFit, reward=="Cumulative" &  ModelName=="LRBF-UCB")$R2, paired=T)
t.test(subset(modelFit, reward=="Best" &  ModelName=="LBMT-UCB")$R2, subset(modelFit, reward=="Best" &  ModelName=="LRBF-UCB")$R2, paired=T)
#############################################################################################################################
# Parameter Estimates
#############################################################################################################################

#local GP-UCB
mDF <- melt(subset(modelFit, kernel=='LRBF' & acq == 'UCB'), id.vars = c("acq", "kernel", "reward", "R2", "bonus"), measure.vars = c("lambda", "beta", "tau"))
boxplotDF <- ddply(mDF, ~variable+acq+kernel, summarise,
                   d_ymin = max(min(value), quantile(value, 0.25) - 1.5 * IQR(value)), 
                   d_ymax = min(max(value), quantile(value, 0.75) + 1.5 * IQR(value)),
                   d_lower = quantile(value, 0.25),  d_middle = median(value), d_upper = quantile(value, 0.75),
                   mu=mean(value))

p<- ggplot(boxplotDF) +
  #stat_summary(fun.y = mean, geom = "bar", position = "dodge", color='black') + 
  geom_boxplot(aes(x =as.numeric(variable)-0.2, ymin = d_lower, ymax = d_upper, lower = d_lower, middle = d_middle, upper = d_upper, width = 2 * 0.2), fill="#56B4E9", stat = "identity", color='black', alpha  =0.8) +
  #whiskers
  geom_segment(aes(x = as.numeric(variable), y = d_ymin, xend = as.numeric(variable), yend = d_ymax)) +
  geom_segment(aes(x = as.numeric(variable) - 0.1,  y = d_ymax, xend = as.numeric(variable),  yend = d_ymax)) + #top
  geom_segment(aes(x = as.numeric(variable) - 0.1, y = d_ymin, xend = as.numeric(variable), yend = d_ymin )) + #bottom
  geom_point(aes(x = as.numeric(variable)-0.2, y = mu), size = 1.3, shape=23, fill='white') +
  geom_jitter(data=mDF, aes(x = as.numeric(variable)+.2,  y = value), color="#56B4E9", size=1, alpha = 0.6, width = 0.2 - 0.25 * 0.2)+
  theme_classic()+
  ylab("Estimate")+
  theme(text = element_text(size=12,  family="sans"),axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), 
        strip.background=element_blank(),
        strip.text = element_blank(),
        legend.key=element_rect(color=NA),
        panel.spacing.x=unit(0.2, "lines"),
        panel.spacing.y=unit(1, "lines"),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        legend.position="none")+
  guides(color=FALSE, shape=FALSE)+
  scale_y_continuous(trans="log10", breaks = c(0.01, 0.1, 1, 10)) +
  annotation_logticks(base=10, sides='l')+
  coord_cartesian(ylim=c(0.005,50))
#ggtitle("Parameter Estimates: Function Learning") +
print(p)

min(mDF$value)
max(mDF$value)
ggsave(filename = "plots/paramEstimates.pdf", plot = p, height =2, width = 4, units = "in") 

exp3Lambda <- subset(mDF, variable=='lambda')$value
exp3beta <- subset(mDF, variable=='beta')$value
exp3tau <- subset(mDF, variable=='tau')$value

#BMT-UCB

#Melt to separate out parameters
mDF <- melt(modelFit, id.vars = c("acq", "kernel", "reward", "environment", "R2", "bonus"), measure.vars = c("kError", "beta", "tau"))
#Select which acq and which kernel to display
k<- "LBMT"
a<- "UCB"

p<- ggplot(subset(mDF, kernel==k &acq==a), aes(y=value, x=variable, color=variable)) +
  geom_jitter(size=2, alpha=.9) +
  geom_boxplot(width=0.2, color="black", outlier.shape=NA) +
  xlab("Parameter") +
  ylab("Estimate (log scale)")+ 
  scale_y_continuous(trans="log10") +
  annotation_logticks(base=10, sides='l')+
  ggtitle("Parameter Estimates: Local BMT-UCB") +
  theme_bw()+
  theme(text = element_text(size=18,  family="serif")) +
  scale_x_discrete("",labels=c("Error Variance", bquote(paste(beta, " (Exploration Bonus)")),bquote(paste(tau, " (Temperature)"))))
#facet_wrap(~environment)
print(p)

bestDF<- subset(modelFit, kernel=="LRBF" & acq=="UCB")

summarizedBetaDF <- ddply(bestDF, ~participant+reward+environment, summarize, meanBeta=mean(beta))
#t-tests
t.test(subset(bestDF, environment=="Rough")$lambda, mu=1) #lambda
t.test(subset(bestDF, environment=="Smooth")$lambda, mu=2) #lambda
t.test(subset(paramEstimates, kernel=="LRBF" & acq=="UCB" & environment=="Smooth")$lambda, mu=2)
t.test(subset(paramEstimates, kernel=="LRBF" & acq=="UCB" & environment=="Rough")$lambda, mu=1) 

#BETA - each individual estimate
t.test(subset(bestDF, environment=="Smooth")$beta, subset(bestDF, environment=="Rough")$beta, var.equal=TRUE) #beta
t.test(subset(bestDF, reward=="Cumulative")$beta, subset(bestDF, reward=="Best")$beta, var.equal=TRUE) #beta
t.test(subset(paramEstimates, kernel=="LRBF" & acq=="UCB" & reward=="Cumulative")$beta, subset(paramEstimates, kernel=="LRBF" & acq=="UCB" & reward=="Best")$beta, var.equal=TRUE) #beta
t.test(subset(paramEstimates, kernel=="LRBF" & acq=="UCB" &  environment=="Rough")$beta, subset(paramEstimates, kernel=="LRBF" & acq=="UCB" &  environment=="Smooth")$beta, var.equal=TRUE) #beta
t.test(subset(paramEstimates, kernel=="LRBF" & acq=="UCB" &  horizon==20)$beta, subset(paramEstimates, kernel=="LRBF" & acq=="UCB" &  horizon==40)$beta, var.equal=TRUE) #beta

#BETA - Averaged over subject
t.test(subset(summarizedBetaDF, environment=="Smooth")$meanBeta, subset(summarizedBetaDF, environment=="Rough")$meanBeta, var.equal=TRUE) #beta
t.test(subset(summarizedBetaDF, reward=="Cumulative")$meanBeta, subset(summarizedBetaDF, reward=="Best")$meanBeta, var.equal=TRUE) #beta
#average beta for each participant for each horizon length
horizonBeta <- ddply(subset(paramEstimates, kernel=="LRBF" & acq=="UCB"), ~participant+reward+environment+horizon, summarize, meanBeta=mean(beta)) 
t.test(subset(horizonBeta, horizon==20)$meanBeta, subset(horizonBeta, horizon==40)$meanBeta, paired=TRUE)

t.test(subset(bestDF, reward=="Cumulative")$tau, subset(bestDF, reward=="Best")$tau, var.equal=TRUE) #beta
t.test(subset(paramEstimates, kernel=="LRBF" & acq=="UCB" &  reward=="Best")$tau, subset(paramEstimates, kernel=="LRBF" & acq=="UCB" &  reward=="Cumulative")$tau, var.equal=TRUE) #tau
t.test(subset(paramEstimates, kernel=="LRBF" & acq=="UCB" &  horizon==20)$tau, subset(paramEstimates, kernel=="LRBF" & acq=="UCB" &  horizon==40)$tau, var.equal=TRUE) 
t.test(subset(paramEstimates, kernel=="LRBF" & acq=="UCB" &  environment=="Smooth")$tau, subset(paramEstimates, kernel=="LRBF" & acq=="UCB" &  environment=="Rough")$tau, var.equal=TRUE) 



###Wilcoxon tests

bestDF<- subset(paramEstimates, kernel=="LRBF" & acq=="UCB")
bestDF <- ddply(bestDF, ~participant+reward, summarize, lambda = median(lambda), beta=median(beta), tau = median(tau))

median(bestDF$lambda)

#Compare with results from exp2
setwd("..")
exp2 <- paramEstimates <- read.csv('analysis2D/modelResults/paramEstimates.csv')
exp2 <- subset(exp2, kernel=="LRBF" & acq=="UCB")
exp2 <- ddply(exp2, ~participant+reward+environment, summarize, lambda = median(lambda), beta=median(beta), tau = median(tau))


median(bestDF$lambda)
median(exp2$lambda)
wtest <- wilcox_test(exp2$lambda~bestDF$lambda)
wtest
qnorm(wtest$p.value)
qnorm(wtest$p.value)/sqrt(length(bestDF$lambda) + length(exp2$lambda))


median(bestDF$beta)
median(exp2$beta)
wtest <- wilcox.test(exp2$beta, bestDF$beta, conf.int=TRUE, conf.level=0.95)
wtest
qnorm(wtest$p.value)
abs(qnorm(wtest$p.value))/sqrt(length(bestDF$lambda) + length(exp2$lambda))


median(bestDF$tau)
median(exp2$tau)
wtest <- wilcox.test(exp2$tau, bestDF$tau, conf.int=TRUE, conf.level=0.95)
wtest
qnorm(wtest$p.value)
abs(qnorm(wtest$p.value))/sqrt(length(bestDF$tau) + length(exp2$tau))


wtest <- wilcox.test(bestDF$lambda, subset(exp2, environment=='Smooth')$lambda)
wtest
qnorm(wtest$p.value)
qnorm(wtest$p.value)/sqrt(length(bestDF$lambda) + length( subset(exp2, environment=='Smooth')$lambda))

wtest <- wilcox.test(bestDF$lambda, subset(exp2, environment=='Rough')$lambda)
wtest
qnorm(wtest$p.value)
qnorm(wtest$p.value)/sqrt(length(bestDF$lambda) + length( subset(exp2, environment=='Rough')$lambda))


median(subset(bestDF, environment=='Rough')$lambda)
wtest <- wilcox.test(subset(bestDF, environment=='Rough')$lambda, mu=1)
wtest
qnorm(wtest$p.value)
qnorm(wtest$p.value)/sqrt(length(subset(bestDF, environment=='Rough')$lambda))


median(bestDF$beta)
wtest <- wilcox_test(bestDF$beta, mu=exp(-5))
wtest
qnorm(wtest$p.value)
qnorm(wtest$p.value)/sqrt(80)

            