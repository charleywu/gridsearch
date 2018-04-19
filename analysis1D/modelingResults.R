#Charley Wu, December 2016
#Interpret 1D modeling results
#1D case

#house keeping
rm(list=ls())
theme_set(theme_bw(base_size=16))# use the b&w theme

#load packages
packages <- c('plyr', 'jsonlite', 'ggplot2', 'gridExtra', 'reshape2', 'stargazer', 'coefplot', "grid", 'parallel', 'corrplot')
lapply(packages, require, character.only = TRUE)


#############################################################################################################################
# IMPORT PARTICIPANT DATA 
#############################################################################################################################

#read in data
data<-read.csv("ExperimentData/experimentData1D.csv")

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
data$id <- seq(1,81)


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
      for (i in 1:81){ #subjects
        filename <- paste0(dataFolder, k, a, i, ".csv") #read output file
        if (file.exists(filename)){
          dp<-read.csv(filename)
          #interpret parameters based on mode combination
          if (k=="LocalSearch"| k=="WSLS" | k=="LWSLS"){#inverse manhattan distance aka Inertia
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
          randomModel120 <- -log(1/30)*120
          R2 <- 1 - (nLL/randomModel120)
          #blank mean parameter estimates
          lambdaMean <- NA
          kErrorMean <- NA
          betaMean <- NA
          tauMean <- median(exp(dp$tau)) #all models should have tau
          #mean parameter estimates for UCB RBF
          if (a=="UCB"){ #UCB should have beta
            betaMean <- median(exp(dp$beta))
          }
          if (k=="RBF" | k=="LRBF"){
            lambdaMean <- median(exp(dp$lambda))
          }
          if (k=="BMT" | k=="LBMT"){
            kErrorMean <- median(exp(dp$kError))
          }
          #save modelFit
          dadd<-data.frame(participant=participant,reward=reward, environment=environment, nLL=nLL, kernel=kernel, acq = acq, R2 =R2, bonus = bonus, kError = kErrorMean, lambda=lambdaMean, beta = betaMean, tau=tauMean)
          names(dadd) <- names(modelFit)
          modelFit<-rbind(modelFit, dadd)
          #Loop through short and long horizons to save parameter estimates
          for (h in c(5,10)){
            #loop through leave out index to save all 8 parameter estimates for each subject
            for (j in 1:8){
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
# modelResults <- importModelResults('modelResults/recovery2/', c("LocalSearch",  "WSLS", "LWSLS",  "BMT",  "LBMT", "RBF", "LRBF"), c("", "GV", "GM", "UCB", "EXI", "POI", "PMU"))
# modelResults <- importModelResults('modelResults/recoveryBMT2/', c("LocalSearch",  "WSLS", "LWSLS",  "BMT",  "LBMT", "RBF", "LRBF"), c("", "GV", "GM", "UCB", "EXI", "POI", "PMU"))

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

#write.csv(modelFit,'modelResults/GPUCBrecovery.csv')
#write.csv(paramEstimates,'modelResults/GPUCBparamEstimates.csv')

#write.csv(modelFit,'modelResults/BMTUCBrecovery.csv')
#write.csv(paramEstimates,'modelResults/BMTUCBparamEstimates.csv')

#############################################################################################################################
# SAVE PARAMETER ESTIMATES FOR RATIONAL MODELS
#############################################################################################################################

#localPars <- subset(paramEstimates, kernel=="LocalSearch")
#write.csv(localPars[,c("participant", "leaveoutindex","reward","environment","horizon", "tau")], file = "rationalModels/parameters/local.csv", row.names = FALSE)

#gpUCBpars <- subset(paramEstimates, kernel=="RBF" & acq=="UCB")
#write.csv(gpUCBpars[,c("participant", "leaveoutindex","reward","environment","horizon","lambda", "beta", "tau")], file = "rationalModels/parameters/gp.csv", row.names = FALSE)

#localgpUCBpars <- subset(paramEstimates, kernel=="LRBF" & acq=="UCB")
#write.csv(localgpUCBpars[,c("participant", "leaveoutindex","reward","environment","horizon","lambda", "beta", "tau")], file = "rationalModels/parameters/localgp.csv", row.names = FALSE)

#BMTpars <- subset(paramEstimates, kernel=="BMT" & acq=="UCB")
#write.csv(BMTpars[,c("participant", "leaveoutindex","reward","environment","horizon","kError", "beta", "tau")], file = "rationalModels/parameters/BMT.csv", row.names = FALSE)

#localBMTpars <- subset(paramEstimates, kernel=="LBMT" & acq=="UCB")
#write.csv(localBMTpars[,c("participant", "leaveoutindex","reward","environment","horizon","kError", "beta", "tau")], file = "rationalModels/parameters/localBMT.csv", row.names = FALSE)


#############################################################################################################################
# NUMBER OF PARTICIPANTS BEST DESCRIBED
#############################################################################################################################
#ALL MODELS
models <-  rep(0, length(levels(modelFit$ModelName)))
names(models) <- levels(modelFit$ModelName)
for (pid in 1:81){
  subDF <- subset(modelFit, participant==pid)
  best <- subDF$ModelName[which(subDF$R2==max(subDF$R2))]
  models[best] <- models[best] + 1
}
#add best described to modelFit 
modelFit$bestDescribed <- NA
for (mod in names(models)){
  modelFit[modelFit$ModelName == mod,]$bestDescribed <- models[mod] / 81
}

#McFadden R^2
cat("All conditions")
for (k in levels(modelFit$kernel)){
  for (a in levels(modelFit$acq)){
    if (paste(k,a,sep="-") %in% names(models)){
      sub <- subset(modelFit, acq == a & kernel == k)
      randomModel <- -log(1/30) * 120
      r_2 <- 1 - (mean(sub$nLL)/randomModel)
      bestDescribed <- models[paste(k,a,sep="-")]
      #subset(modelFit, acq == a & kernel == k)$BestDescribed <- bestDescribed
      cat( k,"-",a, " (n=",nrow(sub),")", "; R2 = ", round(r_2, digits=2),"; Best Described = ", bestDescribed , "; Lambda = ", round(median(sub$lambda), digits=2),"; Beta = ", round(median(sub$beta), digits=2), "; kError = ", round(median(sub$kError), digits=2), "; Tau = ", round(median(sub$tau), digits=2), "\n", sep="")
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
for (pid in 1:81){
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

##############################################################################################################################################
# Correlation of R2
###############################################################################################################################################
#correlation test for matrix
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

#create matrix with models as row names and 
dummyDF <- data.frame()
for (m in levels(mainTextModels$ModelName)){
  dummyDF <- rbind(dummyDF,  mainTextModels[mainTextModels$ModelName==m,]$R2)
}
rownames(dummyDF) <- levels(mainTextModels$ModelName)
rownames(dummyDF)<- c("Assoc.Learn Exploit", "Assoc.Learn Explore", "Assoc.Learn UCB", "Assoc.Learn* Exploit", "Assoc.Learn* Explore", "Assoc.Learn* UCB", "Func.Learn* Exploit", "Func.Learn* Explore", "Func.Learn* UCB", "Func.Learn Exploit", "Func.Learn Explore", "Func.Learn UCB")
colnames <- NA

# matrix of correlations and of p-value of the correlation
M <- cor(t(dummyDF))
p.mat <- cor.mtest(dummyDF)

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
#open output pdf
pdf("plots/mainCorrelations1.pdf",height =8, width = 8, useDingbats=FALSE) 
corrplot(M, method="circle", col=col(200),  
         order = "hclust", type='upper',
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45,
         p.mat = p.mat, sig.level = 0.05, insig = "blank")

dev.off()
#############################################################################################################################
# PLOTS
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
  ggtitle("Experiment 1")
main
ggsave(filename = "plots/mainTextR2.pdf", main,height =2.5, width = 4, units ="in")

#FULL MODEL COMPARISON IN SI
#Model order
modelFit$ModelName <- factor(modelFit$ModelName, levels = c("BMT-UCB", "BMT-GM", "BMT-GV", "BMT-EXI", "BMT-POI", "BMT-PMU", "LBMT-UCB", "LBMT-GM", "LBMT-GV", "LBMT-EXI", "LBMT-POI", "LBMT-PMU", "RBF-UCB", "RBF-GM", "RBF-GV", "RBF-EXI", "RBF-POI", "RBF-PMU", "LRBF-UCB", "LRBF-GM", "LRBF-GV", "LRBF-EXI", "LRBF-POI", "LRBF-PMU", "WSLS-", "LWSLS-", "LocalSearch-"))
modelFit$reward <- factor(modelFit$reward, levels=c("Cumulative", "Best"))

#Recode simple models
modelFit$acq <- as.character(modelFit$acq) #unfactor for now
modelFit$kernel <- as.character(modelFit$kernel)
modelFit[modelFit$kernel=="WSLS",]$acq <- rep("WSLS", 81) 
modelFit[modelFit$kernel=="WSLS",]$kernel <- rep("Simple Strategies", 81)
modelFit[modelFit$kernel=="LWSLS",]$acq <- rep("WSLS*", 81) 
modelFit[modelFit$kernel=="LWSLS",]$kernel <- rep("Simple Strategies", 81)
modelFit[modelFit$kernel=="LocalSearch",]$acq <- rep("Local Search", 81) 
modelFit[modelFit$kernel=="LocalSearch",]$kernel <- rep("Simple Strategies", 81)

#Refactor
modelFit$acq <- factor(modelFit$acq, levels = c("UCB", "GM", "GV", "EXI", "POI", "PMU", "WSLS", "WSLS*", "Local Search"))
levels(modelFit$acq) <- c("UCB", "Exploit", "Explore", "EXI", "POI", "PMU", "WSLS", "WSLS*", "Local Search") 
modelFit$kernel <- factor(modelFit$kernel, levels = c("BMT", "LBMT", "RBF", "LRBF", "Simple Strategies"))
levels(modelFit$kernel) <- c("Mean Tracker", "Mean Tracker*", "Function Learning", "Function Learning*", "Simple Strategies")

full <- ggplot(modelFit , aes(y=R2, x=acq, fill=interaction(environment, reward))) +
  stat_summary(fun.y = mean, geom = "bar", position = "dodge") + 
  geom_jitter(color='grey', shape=21, alpha=0.4, size = .8, position=position_jitterdodge(dodge.width=0.9, jitter.width=0.2))+
  stat_summary(fun.data = mean_se, geom = "errorbar", position = position_dodge(width = 0.90), width = 0.4 ) +
  #geom_boxplot(position="dodge")+
  #geom_text(aes(label=bestDescribed, y = 0.5), family="serif") +
  xlab("") +
  scale_fill_manual(values = c("#7F0000", "#DA4233" , "#005F8D", "#00BCE2"))+
  ylab("Predictive Accuracy")+ 
  theme_classic()+
  scale_y_continuous(breaks = c(0, .2, .4, .6, .8))+
  theme(text = element_text(size=16,  family="sans"), strip.background=element_blank(), legend.key=element_rect(color=NA), legend.position="none")+
  coord_cartesian(ylim=c(-0.25, 0.8)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.title=element_blank())+
  ggtitle("Experiment 1") +
  guides(shape=FALSE) +
  facet_grid(~kernel,  scales="free_x", space="free_x")
full

ggsave(filename = "plots/fullComparison1.pdf", plot = full, height =5, width = 12, units = "in") 

#RANDOM MODEL BIC
#randomModel20 <- -2 * sum(rep(log(1/121),20)) #equal probability for each 121 options for short horizon
#randomModel40<- -2 * sum(rep(log(1/121),20)) #equal probability for each 121 options for long horizon
randomModel120 <- -log(1/30)*120 #for full task

##Model comparison plots
ggplot(modelFit, aes(y=nLL, x=ModelName, color=reward, shape=environment)) +
  geom_jitter(alpha=1) + 
  geom_hline(yintercept = randomModel120, linetype = "dotdash" )+
  xlab("Model") +
  ylab(" Out of Sample nLL")+ 
  theme_bw()+
  #theme(legend.position = "none") +
  #ylim(c(500, 2500)) +
  #facet_wrap(~kernel) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggtitle("Model Comparison")

#COGSCI PLOT
#Some rather unsavoury manipulations to put simple models and computational models on the same plot
#1. Rename the kernel 

levels(modelFit$kernel) <-c("BMT", "Local BMT","GP", "Local GP")
p1 <- ggplot(subset(modelFit, kernel=="Local Search"| kernel=="WSLS"), aes(x=kernel, y = R2, fill=reward))+
  stat_summary(fun.y = mean, geom = "bar", position = "dodge") + 
  stat_summary(fun.data = mean_se, geom = "errorbar", position = position_dodge(width = 0.90), width = 0.2 ) +
  xlab("Model") +
  #geom_hline(yintercept = 0.068, linetype = "dotdash" )+ #Inertia model
  ylab(expression(paste("McFadden's R"^"2"," (" %+-% "SE)") ))+ 
  theme_bw()+
  scale_fill_manual(values = c("#A20000", "#0082C1"),name="Reward Condition", labels=c("Average Reward","Maximum Reward")) +
  ggtitle("Model Comparison") +
  #scale_x_discrete("",labels=c("UCB", "PMU", "POI"))+
  theme(text = element_text(size=24,  family="serif"), strip.background=element_blank(), legend.key=element_rect(color=NA), legend.position="bottom")
p1
ggsave(filename = "plots/R2plota.pdf", plot = p1, height =4.33, width = 3, units = "in") 

p2  <- ggplot(subset(modelFit, acq=="GV" | acq=="GM" |acq=="UCB") , aes(x=acq, y = R2))+
  stat_summary(fun.y = mean, geom = "bar", position = "dodge") + 
  stat_summary(fun.data = mean_se, geom = "errorbar", position = position_dodge(width = 0.90), width = 0.2 ) +
  facet_wrap(~kernel, nrow = 1) +
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
ggsave(filename = "plots/R2pnas.pdf", plot = p2, height =4.33, width = 12, units = "in") 



######T-tests########
#model comparison
t.test(subset(modelFit, acq=="GM" & kernel=="WSLS")$R2, mu=0)
t.test(subset(modelFit, acq=="UCB" & kernel=="Local Search")$R2, subset(modelFit, acq=="UCB" & kernel=="RBF")$R2, var.equal=TRUE)
t.test(subset(modelFit, acq=="UCB" & kernel=="LRBF")$R2, subset(modelFit, acq=="UCB" & kernel=="RBF")$R2, var.equal=TRUE)
t.test(subset(modelFit, acq=="UCB" & kernel=="Local GP")$R2, subset(modelFit, acq=="UCB" & kernel=="Local Search")$R2, var.equal=TRUE)
  
#Payoff conditions and environment type
t.test(subset(modelFit, reward=="Cumulative")$R2, subset(modelFit, reward=="Best")$R2, var.equal=TRUE)
t.test(subset(modelFit, environment=="Smooth")$R2, subset(modelFit, environment=="Rough")$R2, var.equal=TRUE)


#############################################################################################################################
# Parameter Estimates
#############################################################################################################################

#GP-UCB
#Melt to separate out parameters
mDF <- melt(subset(modelFit, kernel=='RBF' & acq == 'UCB'), id.vars = c("acq", "kernel", "reward", "R2", "bonus"), measure.vars = c("lambda", "beta", "tau"))
boxplotDF <- ddply(mDF, ~variable+acq+kernel, summarise,
                   d_ymin = max(min(value), quantile(value, 0.25) - 1.5 * IQR(value)), 
                   d_ymax = min(max(value), quantile(value, 0.75) + 1.5 * IQR(value)),
                   d_lower = quantile(value, 0.25),  d_middle = median(value), d_upper = quantile(value, 0.75),
                   mu=mean(value))

p<- ggplot(boxplotDF) +
  #stat_summary(fun.y = mean, geom = "bar", position = "dodge", color='black') + 
  geom_boxplot(aes(x =as.numeric(variable)-0.2, ymin = d_lower, ymax = d_upper, lower = d_lower, middle = d_middle, upper = d_upper, width = 2 * 0.2, fill = variable), stat = "identity", color='black', fill = "#009E73") +
  #whiskers
  geom_segment(aes(x = as.numeric(variable), y = d_ymin, xend = as.numeric(variable), yend = d_ymax)) +
  geom_segment(aes(x = as.numeric(variable) - 0.1,  y = d_ymax, xend = as.numeric(variable),  yend = d_ymax )) + #top
  geom_segment(aes(x = as.numeric(variable) - 0.1, y = d_ymin, xend = as.numeric(variable), yend = d_ymin)) + #bottom
  geom_point(aes(x = as.numeric(variable)-0.2, y = mu), size = 1.3, shape=23, fill='white') +
  geom_jitter(data=mDF, aes(x = as.numeric(variable)+.2,  y = value), color="#009E73", size=1, alpha = 0.6, width = 0.2 - 0.25 * 0.2)+
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

min(boxplotDF$d_ymin)
max(boxplotDF$d_ymax)
ggsave(filename = "plots/paramEstimates.pdf", plot = p, height =2, width = 4, units = "in") 

#individual variance
#Melt to separate out parameters
individualDF <- melt(paramEstimates, id.vars = c("participant", "acq", "kernel", "reward", "environment", "R2", "bonus"), measure.vars = c("lambda", "beta", "tau"))
#Select which acq and which kernel to display
k<- "RBF"
a<- "UCB"
p<- ggplot(subset(individualDF, kernel==k &acq==a), aes(y=value, x=participant, color=variable)) +
  geom_boxplot(fill=NA,width=0.5, color="black", outlier.shape=NA, lwd=0.5) +
  stat_summary(fun.data = "mean_cl_boot", alpha=.2, size = 0.4, lwd = 1) + #mean
  #geom_jitter(aes(fill=variable), shape=21, size=1,color='black', alpha=.4, width = 0.2) +
  #geom_boxplot(aes(fill=variable), alpha = 0.3, outlier.shape=NA) +
  #stat_summary(fun.y=mean, geom="point", 21, size=1, alpha=0.2) + #mean
  #xlab("Parameter") +
  ylab("Estimate")+ 
  scale_fill_brewer(palette="Set1") +
  scale_color_brewer(palette="Set1") +
  #scale_y_continuous(trans="log10") +
  #annotation_logticks(base=10, sides='l')+
  coord_cartesian(ylim=c(0.005,3.5))+
  theme_classic()+
  theme(text = element_text(size=16,  family="sans"),strip.background=element_blank(), legend.key=element_rect(color=NA), legend.position="None") +
  #scale_x_discrete("",labels=c("")) + 
  facet_grid(~variable)
print(p)

ggsave(filename = "plots/individualparamEstimates.pdf", plot = p, height =2.82, width = 5.38, units = "in") 

#BMT-UCB

#Melt to separate out parameters
mDF <- melt(modelFit, id.vars = c("acq", "kernel", "reward", "environment", "R2", "bonus"), measure.vars = c("kError", "beta", "tau"))
#Select which acq and which kernel to display
k<- "BMT"
a<- "UCB"

p<- ggplot(subset(mDF, kernel==k &acq==a), aes(y=value, x=variable, fill=environment, color = environment)) +
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


#Mean-based estimates
bestDF<- subset(paramEstimates, kernel=="RBF" & acq=="UCB")
bestDF <- ddply(bestDF, ~participant+reward+environment, summarize, beta=mean(beta), lambda = mean(lambda), tau=mean(tau))
#t-tests
t.test(subset(bestDF, environment=="Rough")$lambda, mu=1) #lambda
cohensD(subset(bestDF, environment=="Rough")$lambda, mu=1)
t.test(subset(bestDF, environment=="Smooth")$lambda, mu=2) #lambda
cohensD(subset(bestDF, environment=="Smooth")$lambda, mu=2)

#BETA - each individual estimate
t.test(bestDF$beta, mu=exp(-5))
cohensD(bestDF$beta, mu=exp(-5))
t.test(subset(bestDF, environment=="Smooth")$beta, subset(bestDF, environment=="Rough")$beta, var.equal=TRUE) #beta
t.test(subset(bestDF, reward=="Cumulative")$beta, subset(bestDF, reward=="Best")$beta, var.equal=TRUE) #beta
t.test(subset(paramEstimates, kernel=="RBF" & acq=="UCB" & reward=="Cumulative")$beta, subset(paramEstimates, kernel=="RBF" & acq=="UCB" & reward=="Best")$beta, var.equal=TRUE) #beta
t.test(subset(paramEstimates, kernel=="RBF" & acq=="UCB" &  environment=="Rough")$beta, subset(paramEstimates, kernel=="RBF" & acq=="UCB" &  environment=="Smooth")$beta, var.equal=TRUE) #beta
t.test(subset(paramEstimates, kernel=="RBF" & acq=="UCB" &  horizon==5)$beta, subset(paramEstimates, kernel=="RBF" & acq=="UCB" &  horizon==10)$beta, var.equal=TRUE) #beta

#BETA - Averaged over subject
t.test(subset(summarizedBetaDF, environment=="Smooth")$meanBeta, subset(summarizedBetaDF, environment=="Rough")$meanBeta, var.equal=TRUE) #beta
t.test(subset(summarizedBetaDF, reward=="Cumulative")$meanBeta, subset(summarizedBetaDF, reward=="Best")$meanBeta, var.equal=TRUE) #beta
#average beta for each participant for each horizon length
horizonBeta <- ddply(subset(paramEstimates, kernel=="RBF" & acq=="UCB"), ~participant+reward+environment+horizon, summarize, meanBeta=mean(beta)) 
t.test(subset(horizonBeta, horizon==5)$meanBeta, subset(horizonBeta, horizon==10)$meanBeta, paired=TRUE)

#tau
mean(bestDF$tau)
t.test(subset(bestDF, reward=="Cumulative")$tau, subset(bestDF, reward=="Best")$tau, var.equal=TRUE) #beta
t.test(subset(paramEstimates, kernel=="RBF" & acq=="UCB" &  reward=="Best")$tau, subset(paramEstimates, kernel=="RBF" & acq=="UCB" &  reward=="Cumulative")$tau, var.equal=TRUE) #tau
t.test(subset(paramEstimates, kernel=="RBF" & acq=="UCB" &  horizon==5)$tau, subset(paramEstimates, kernel=="RBF" & acq=="UCB" &  horizon==10)$tau, var.equal=TRUE) 
t.test(subset(paramEstimates, kernel=="RBF" & acq=="UCB" &  environment=="Smooth")$tau, subset(paramEstimates, kernel=="RBF" & acq=="UCB" &  environment=="Rough")$tau, var.equal=TRUE) 

###Median based tests

bestDF<- subset(paramEstimates, kernel=="RBF" & acq=="UCB")
bestDF <- ddply(bestDF, ~participant+reward+environment, summarize, beta=median(beta), lambda = median(lambda), tau=median(tau))

#lambda
wtest <- wilcox.test(subset(bestDF, environment=="Smooth")$lambda,subset(bestDF, environment=="Rough")$lambda)
wtest
qnorm(wtest$p.value) #z-statistic
qnorm(wtest$p.value)/sqrt(81) #effect size; r value
median(subset(bestDF, environment=="Smooth")$lambda)
median(subset(bestDF, environment=="Rough")$lambda)

wtest<- wilcox.test(subset(bestDF, environment=="Smooth")$lambda, mu = 2)
wtest
qnorm(wtest$p.value) #z-statistic
qnorm(wtest$p.value)/sqrt(length(subset(bestDF, environment=="Smooth")$lambda)) #effect size; r value

wtest<- wilcox.test(subset(bestDF, environment=="Rough")$lambda, mu = 1)
wtest
qnorm(wtest$p.value) #z-statistic
qnorm(wtest$p.value)/sqrt(length(subset(bestDF, environment=="Rough")$lambda)) #effect size; r value
            
#beta
median(bestDF$beta)
wtest<- wilcox.test(bestDF$beta, mu = exp(-5))
wtest
qnorm(wtest$p.value) #z-statistic
qnorm(wtest$p.value)/sqrt(length(bestDF$beta)) #effect size; r value

wtest<- wilcox.test(subset(bestDF, reward=="Cumulative")$beta, subset(bestDF, reward="Best")$beta)
wtest

wtest<- wilcox.test(subset(bestDF,environment=="Rough")$beta, subset(bestDF, environment=="Smooth")$beta)
wtest
qnorm(wtest$p.value) #z-statistic
qnorm(wtest$p.value)/sqrt(81) #effect size; r value



