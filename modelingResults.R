#Charley Wu, December 2016
#Interpret modeling results

#house keeping
rm(list=ls())
theme_set(theme_bw(base_size=16))# use the b&w theme

#load packages
packages <- c('plyr', 'jsonlite', 'ggplot2', 'gridExtra', 'reshape2', 'stargazer', 'coefplot', "grid")
lapply(packages, require, character.only = TRUE)

#############################################################################################################################
# IMPORT PARTICIPANT DATA 
#############################################################################################################################

#Participant data
data<-read.csv("ExperimentData/fullExperiment.csv",  sep="\t")
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
  modelFit <- data.frame(participant=numeric(), reward=numeric(), environment=numeric(), nLL=numeric(), kernel=numeric(), acq=numeric(), R2=numeric(), bonus=numeric(), lambda=numeric(), beta=numeric(), tau=numeric()) 
  paramEstimates <- data.frame(participant=numeric(), reward=numeric(), environment=numeric(), horizon=numeric(),leaveoutindex=numeric(), nLL=numeric(),R2=numeric(), kernel=numeric(), acq=numeric(), kError=numeric(), lambda=numeric(), beta = numeric(), tau=numeric(), roundnLL=numeric(), bonus=numeric())
  #loop through data
  for (k in kernels){
    for (a in acqFuncs){
      for (i in 1:80){ #subjects
        filename <- paste0(dataFolder, k, a, i, ".csv") #read output file
        if (file.exists(filename)){
          dp<-read.csv(filename)
          #interpret parameters based on mode combination
          if (k=="IMD"| k=="WSLG"){#inverse manhattan distance aka Inertia
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
          #save modelFit
          dadd<-data.frame(participant=participant,reward=reward, environment=environment, nLL=nLL, kernel=kernel, acq = acq, R2 =R2, bonus = bonus, lambda=lambdaMean, beta = betaMean, tau=tauMean)
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
#############################################################################################################################
modelResults <- importModelResults('modelResults/paper/', c("IMD","WSLG", "RBF", "LRBF"), c("GV", "GM", "UCB"))
#modelResults <- importModelResults('modelResults/NestedUCB/', c("IMD", "RBF", "LRBF"), c("UCB", "GM", "GV"))

modelFit <- modelResults[[1]]
paramEstimates <- modelResults[[2]]

#reorder acqisition function levels
modelFit$acq <-factor(modelFit$acq, levels = c("GV", "GM", "UCB"))

#############################################################################################################################
# SAVE PARAMETER ESTIMATES FOR RATIONAL MODELS
#############################################################################################################################
#localPars <- subset(paramEstimates, kernel=="IMD")
#gpUCBpars <- subset(paramEstimates, kernel=="RBF" & acq=="UCB")
#localgpUCBpars <- subset(paramEstimates, kernel=="LRBF" & acq=="UCB")

#write.csv(localPars$tau, file = "rationalModels/parameters/local.csv", row.names = FALSE)
#write.csv(cbind(gpUCBpars$lambda, gpUCBpars$beta, gpUCBpars$tau), file = "rationalModels/parameters/gp.csv", row.names = FALSE)
#write.csv(cbind(localgpUCBpars$lambda, localgpUCBpars$beta, localgpUCBpars$tau), file = "rationalModels/parameters/localgp.csv", row.names = FALSE)

#############################################################################################################################
# PLOTS
#############################################################################################################################
#RANDOM MODEL BIC
#randomModel20 <- -2 * sum(rep(log(1/121),20)) #equal probability for each 121 options for short horizon
#randomModel40<- -2 * sum(rep(log(1/121),20)) #equal probability for each 121 options for long horizon
randomModel240 <- -log(1/121)*240 #for full task

##Model comparison plots
ggplot(modelFit, aes(y=nLL, x=acq, color=acq, shape=environment)) +
  geom_jitter(aes(alpha=.5)) + 
  geom_hline(yintercept = randomModel240, linetype = "dotdash" )+
  xlab("Model") +
  ylab(" Out of Sample nLL")+ 
  theme_bw()+
  theme(legend.position = "none") +
  #ylim(c(500, 2500)) +
  facet_wrap(~kernel) +
  ggtitle("Model Comparison")


#Some rather unsavoury manipulations to put simple models and computational models on the same plot
#1. Rename the kernel 

levels(modelFit$kernel) <-c("Local Search", "WSLS","GP", "Local GP")
p1 <- ggplot(subset(modelFit, kernel=="Local Search"| kernel=="WSLS"), aes(x=kernel, y = R2, fill=reward))+
  stat_summary(fun.y = mean, geom = "bar", position = "dodge") + 
  stat_summary(fun.data = mean_se, geom = "errorbar", position = position_dodge(width = 0.90), width = 0.2 ) +
  xlab("Model") +
  coord_cartesian(ylim=c(0,.45))+
  #geom_hline(yintercept = 0.068, linetype = "dotdash" )+ #Inertia model
  ylab(expression(paste("McFadden's R"^"2"," (" %+-% "SE)") ))+ 
  theme_bw()+
  scale_fill_manual(values = c("#A20000", "#0082C1"),name="Reward Condition", labels=c("Average Reward","Maximum Reward")) +
  ggtitle("Model Comparison") +
  #scale_x_discrete("",labels=c("UCB", "PMU", "POI"))+
  theme(text = element_text(size=24,  family="serif"), strip.background=element_blank(), legend.key=element_rect(color=NA), legend.position="bottom")
p1
ggsave(filename = "plots/R2plota.pdf", plot = p1, height =4.33, width = 3, units = "in") 

p2  <- ggplot(subset(modelFit, kernel=="WSLS" | kernel=="GP" |kernel=="Local GP") , aes(x=acq, y = R2, fill=reward))+
  stat_summary(fun.y = mean, geom = "bar", position = "dodge") + 
  stat_summary(fun.data = mean_se, geom = "errorbar", position = position_dodge(width = 0.90), width = 0.2 ) +
  facet_wrap(~kernel) +
  xlab("Model") +
  coord_cartesian(ylim=c(0,.45))+
  #geom_hline(yintercept = 0.068, linetype = "dotdash" )+ #Inertia model
  ylab(expression(paste("McFadden's R"^"2"," (" %+-% "SE)") ))+ 
  theme_bw()+
  scale_fill_manual(values = c("#A20000", "#0082C1"),name="Reward Condition", labels=c("Average Reward","Maximum Reward")) +
  ggtitle("Model Comparison") +
  #scale_x_discrete("",labels=c("UCB", "PMU", "POI"))+
  theme(text = element_text(size=24,  family="serif"), strip.background=element_blank(), legend.key=element_rect(color=NA), legend.position="bottom")
p2
ggsave(filename = "plots/R2plotb.pdf", plot = p2, height =4.33, width = 9, units = "in") 

#McFadden R^2
cat("All conditions")
for (k in levels(modelFit$kernel)){
  for (a in levels(modelFit$acq)){
    sub <- subset(modelFit, acq == a & kernel == k)
    randomModel <- -log(1/121) * 240
    r_2 <- 1 - (mean(sub$nLL)/randomModel)
    cat( k,"-",a, " (n=",nrow(sub),")", " = ", round(r_2, digits=2), "\n", sep="")
  }
}
  

######T-tests########
#model comparison
t.test(subset(modelFit, acq=="GM" & kernel=="WSLS")$R2, mu=0)
t.test(subset(modelFit, acq=="UCB" & kernel=="Local Search")$R2, subset(modelFit, acq=="UCB" & kernel=="GP")$R2, var.equal=TRUE)
t.test(subset(modelFit, acq=="UCB" & kernel=="Local GP")$R2, subset(modelFit, acq=="UCB" & kernel=="GP")$R2, var.equal=TRUE)
t.test(subset(modelFit, acq=="UCB" & kernel=="Local GP")$R2, subset(modelFit, acq=="UCB" & kernel=="Local Search")$R2, var.equal=TRUE)
  
#Payoff conditions and environment type
t.test(subset(modelFit, reward=="Cumulative")$R2, subset(modelFit, reward=="Best")$R2, var.equal=TRUE)
t.test(subset(modelFit, environment=="Smooth")$R2, subset(modelFit, environment=="Rough")$R2, var.equal=TRUE)

#Melt to separate out parameters
mDF <- melt(modelFit, id.vars = c("acq", "kernel", "reward", "environment", "R2", "bonus"), measure.vars = c("lambda", "beta", "tau"))
#Select which acq and which kernel to display
k<- "Local GP"
a<- "UCB"

p<- ggplot(subset(mDF, kernel==k &acq==a), aes(y=value, x=variable)) +
  geom_jitter(aes(color=variable), size=2, alpha=.9) +
  geom_boxplot(width=0.2, color="black", outlier.shape=NA) +
  xlab("Parameter") +
  ylab("Estimate (log scale)")+ 
  scale_y_continuous(trans="log10") +
  annotation_logticks(base=10, sides='l')+
  ggtitle("Parameter Estimates: Local GP-UCB") +
  theme_bw()+
  theme(text = element_text(size=24,  family="serif"),strip.background=element_blank(), legend.key=element_rect(color=NA), legend.position="None") +
  scale_x_discrete("",labels=c(bquote(paste(lambda, " (Length-Scale)")), bquote(paste(beta, " (Exploration Bonus)")),bquote(paste(tau, " (Temperature)"))))
#facet_wrap(~environment)
print(p)

ggsave(filename = "plots/paramEstimates.pdf", plot = p, height =4.33, width = 9, units = "in") 

bestDF<- subset(modelFit, kernel=="Local GP" & acq=="UCB")

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
            