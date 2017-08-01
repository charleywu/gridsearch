#Charley Wu, July 2017
#Interpret 2D modeling results

#house keeping
rm(list=ls())
theme_set(theme_bw(base_size=16))# use the b&w theme

#load packages
packages <- c('plyr', 'jsonlite', 'ggplot2', 'gridExtra', 'reshape2', 'stargazer', 'coefplot', "grid", 'corrplot')
lapply(packages, require, character.only = TRUE)

#############################################################################################################################
# IMPORT PARTICIPANT DATA 
#############################################################################################################################

#Participant data
data<-read.csv("ExperimentData/experimentData2D.csv",  sep="\t")
#remove null rows
data <- subset(data, MTurkID != "NULL")

#add factors
data$scenario <- as.factor(data$scenario)
levels(data$scenario) <- c("Cumulative", "Best")
data$kernel <- as.factor(data$kernel)
levels(data$kernel) <- c("Rough", "Smooth")
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
  modelFit <- data.frame(participant=numeric(), reward=numeric(), environment=numeric(), nLL=numeric(), kernel=numeric(), acq=numeric(), R2=numeric(), bonus=numeric(), kError=numeric(), lambda=numeric(), beta=numeric(), tau=numeric()) 
  paramEstimates <- data.frame(participant=numeric(), reward=numeric(), environment=numeric(), horizon=numeric(),leaveoutindex=numeric(), nLL=numeric(),R2=numeric(), kernel=numeric(), acq=numeric(), kError=numeric(), lambda=numeric(), beta = numeric(), tau=numeric(), roundnLL=numeric(), bonus=numeric())
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
          environment <- dummy$kernel
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
          dadd<-data.frame(participant=participant,reward=reward, environment=environment, nLL=nLL, kernel=kernel, acq = acq, R2 =R2, bonus = bonus, kError = kErrorMean, lambda=lambdaMean, beta = betaMean, tau=tauMean)
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
              dadd<-data.frame(participant=participant,reward=reward, environment=environment, horizion=h, leaveoutindex =j, nLL=nLL, R2 =R2,  kernel=kernel, acq = acq, kError = kError, lambda=lambda, beta=beta,  tau=tau, roundnLL=roundnLL, bonus=bonus)
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
modelResults <- importModelResults('modelResults/final/', c("LocalSearch",  "WSLS", "LWSLS",  "BMT",  "LBMT", "RBF", "LRBF"), c("", "GV", "GM", "UCB", "EXI", "POI", "PMU"))
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
#write.csv(localPars[,c("participant", "leaveoutindex", "reward","environment","horizon", "tau")], file = "rationalModels/parameters/local.csv", row.names = FALSE)

#gpUCBpars <- subset(paramEstimates, kernel=="RBF" & acq=="UCB")
#write.csv(gpUCBpars[,c("participant", "leaveoutindex","reward","environment","horizon","lambda", "beta", "tau")], file = "rationalModels/parameters/gp.csv", row.names = FALSE)

#localgpUCBpars <- subset(paramEstimates, kernel=="LRBF" & acq=="UCB")
#write.csv(localgpUCBpars[,c("participant", "leaveoutindex", "reward","environment","horizon","lambda", "beta", "tau")], file = "rationalModels/parameters/localgp.csv", row.names = FALSE)

#BMTpars <- subset(paramEstimates, kernel=="BMT" & acq=="UCB")
#write.csv(BMTpars[,c("participant", "leaveoutindex", "reward","environment","horizon","kError", "beta", "tau")], file = "rationalModels/parameters/BMT.csv", row.names = FALSE)

#localBMTpars <- subset(paramEstimates, kernel=="LBMT" & acq=="UCB")
#write.csv(localBMTpars[,c("participant", "leaveoutindex", "reward","environment","horizon","kError", "beta", "tau")], file = "rationalModels/parameters/localBMT.csv", row.names = FALSE)

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

main <- ggplot(mtmDF, aes(y=newR2, x=acq, fill=kernel)) +
  stat_summary(fun.y = mean, geom = "bar", position = "dodge", color='black') + 
  geom_errorbar(aes(ymin=newR2 - se, ymax=newR2 + se), width = .2, position=position_dodge((width=0.9))) +
  #geom_text(aes(label=best, y = newR2 + se + 0.02),position = position_dodge(width=0.9)) +
  xlab("Sampling Strategy") +
  #scale_fill_manual(values = c("#7F0000", "#00BCE2",  "#37045F" ))+
  scale_fill_manual(values = c(  "#F0E442", "#E69F00", "#009E73", "#56B4E9"))+
  ylab("Predictive Accuracy")+ 
  theme_classic()+
  ylim(c(0,.45)) +
  theme(text = element_text(size=16,  family="serif"), strip.background=element_blank(), legend.key=element_rect(color=NA), legend.position="none")+
  guides(color=FALSE, shape=FALSE)+
  ggtitle("Experiment 2")
main
ggsave(filename = "plots/mainTextComparison.pdf", plot = main, height =2.82, width = 5.38, units = "in") 

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

full <- ggplot(modelFit , aes(y=R2, x=acq, fill=interaction(environment, reward))) +
  stat_summary(fun.y = mean, geom = "bar", position = "dodge") + 
  stat_summary(fun.data = mean_se, geom = "errorbar", position = position_dodge(width = 0.90), width = 0.2 ) +
  #geom_text(aes(label=bestDescribed, y = 0.5), family="serif") +
  xlab("") +
  scale_fill_manual(values = c("#7F0000", "#DA4233" , "#005F8D", "#00BCE2"))+
  ylab("Predictive Accuracy")+ 
  theme_classic()+
  theme(text = element_text(size=16,  family="sans"), strip.background=element_blank(), legend.key=element_rect(color=NA), legend.position="none")+
  coord_cartesian(ylim=c(0, 0.5)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.title=element_blank())+
  ggtitle("Experiment 2") +
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


#############################################################################################################################
# Parameter Estimates
#############################################################################################################################

#local GP-UCB

#Melt to separate out parameters
mDF <- melt(modelFit, id.vars = c("acq", "kernel", "reward", "environment", "R2", "bonus"), measure.vars = c("lambda", "beta", "tau"))
#Select which acq and which kernel to display
k<- "LRBF"
a<- "UCB"

p<- ggplot(subset(mDF, kernel==k &acq==a), aes(y=value, x=NA, color=variable)) +
  geom_jitter(aes(color=variable), size=1, alpha=.6, width = 0.2) +
  geom_boxplot(aes(fill=variable),width=0.2, color="black", outlier.shape=NA, lwd=0.5, alpha=0) +
  stat_summary(color='black',fun.y=mean, geom="point", shape=5, size=2, stroke = 1) + #mean
  xlab("Parameter") +
  ylab("Estimate")+ 
  scale_fill_brewer(palette="Set1") +
  scale_color_brewer(palette="Set1") +
  #scale_y_continuous(trans="log10") +
  #annotation_logticks(base=10, sides='l')+
  coord_cartesian(ylim=c(0.005,3.5))+
  #ggtitle("Parameter Estimates: Function Approximatoin") +
  theme_classic()+
  theme(text = element_text(size=14,  family="serif"),strip.background=element_blank(), legend.key=element_rect(color=NA), legend.position="None") +
  scale_x_discrete("",labels=c("")) + 
  facet_grid(~variable)
print(p)

ggsave(filename = "plots/paramEstimates.pdf", plot = p, height =2.7, width = 3, units = "in") 

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
            