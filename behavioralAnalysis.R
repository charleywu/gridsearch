#Script to import experiment data
#Charley Wu August 2016

#house keeping
rm(list=ls())

#load packages
packages <- c('plyr', 'jsonlite', 'ggplot2')
lapply(packages, require, character.only = TRUE)

##############################################################################################################
##  1. Subject Level Data                                                                                   ##
##############################################################################################################

#read in data
data<-read.csv("ExperimentData/fullExperiment.csv",  sep=",")
#remove null rows
data <- subset(data, searchHistory != "NULL")

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


##############################################################################################################
##  2. Trial Level Data                                                                                     ##
##############################################################################################################
tdata <- read.csv('ExperimentData/newdat.csv')
#dummy data frame
d<-data.frame()

#function to test if a search decision is performing a local exploitation or a global exploration
exploreOrExploit<- function(x, y){
  output <- rep(NA, length(x_round)) # dummy output vector
  for (i in 2:length(x)){
    current <- c(x[i], y[i])
    hist <- cbind(x_round[1:(i-1)],y_round[1:(i-1)] )
    #if x and y in history --> repeat
    if (length(which(apply(hist,1,function(z) identical(z,current))))>0){
      output[i] <- "Repeat"
    }else{ #if x and y are new:
      #if x and y are a neighbor of what is in the history --> exploit
      nearestNeighbor <- min(sqrt((current[1] - hist[,1])^2 + (current[2] - hist[,2])^2))
      if (nearestNeighbor <= sqrt(2)){
        output[i] <- "Exploit"
      }else{#if x and y are not a neighbor -- > explore
        output[i] <- "Explore"
      }
    } 
  }
  return(output)
}

#loop through each subject
for (i in 1:nrow(tdata)){
  searchHistory <- fromJSON(as.character(tdata$searchHistory[i]))
  #sampled x value
  x<-unlist(searchHistory$xcollect)
  #sampled y value
  y<-unlist(searchHistory$ycollect)
  #sampled z value
  z<-unlist(searchHistory$zcollect)
  #length
  len<-sapply(searchHistory[1]$xcollect, length)
  #trial 1:20 or 1:40
  trial<-c(1:len[1],1:len[2],1:len[3],1:len[4], 1:len[5],1:len[6],1:len[7],1:len[8])
  #round, 1:8
  round<-rep(1:8, len)
  #horizon, 20 or 40
  horizon<-rep(len-1, len)
  #maximum value or cumulative
  scenario<-rep(tdata$scenario[i], sum(len))
  #smoothness
  kernel<-rep(tdata$kernel[i], sum(len))
  #id number
  id<-rep(i, sum(len))
  #Environment Vector
  environments <- fromJSON(as.character(tdata$envOrder[i]))
  environment <- sapply(round, FUN=function(x) environments[x])
  #Payoff Scaling
  scales <- fromJSON(as.character(tdata$scale[i]))
  payoffScale <- sapply(round, FUN=function(x) scales[x])
  #Distance from previous location, change in reward, and change in reward for previous trial
  distance <- c()
  deltaReward <- c()
  deltaRewardPrev <- c()
  #Explore or Exploit
  exploreExploit <- c()
  #improvement on best
  improvement <- c()
  #loop through rounds
  for (j in 1:8){
    end <- sum(len[1:j]) #trial indicator for end of round
    start <- end  - horizon[sum(len[1:j])]
    x_round <- x[end:start]
    y_round <- y[end:start]
    z_round <- z[end:start]
    # sqrt((xnew - xold)^2 + (ynew - yold)^2)
    dist <- sqrt((x_round[-1] - x_round[1:horizon[sum(len[1:j])]])^2 + (y_round[-1] - y_round[1:horizon[sum(len[1:j])]])^2)
    #add a NA to the first value, since it was not a choice
    dist <- append(NA,dist)
    distance <-append(distance,dist)
    #new reward minus old
    dReward <- z_round[-1] - z_round[1:horizon[sum(len[1:j])]]
    #add a NA to the first value, since it was not a choice
    dReward <- append(NA,dReward)
    deltaReward <- append(deltaReward,dReward)
    #calculate previous trial's delta reward
    prevdReward <- append(c(NA,NA), dReward[2:(length(dReward)-1)])
    deltaRewardPrev <- append(deltaRewardPrev, prevdReward)
    #round explore or exploit
    roundExpExp <- exploreOrExploit(x_round, y_round)
    exploreExploit <- append(exploreExploit, roundExpExp)
    #improvement on best
    improv <- sapply(2:length(z_round), FUN = function(z_i) ifelse(z_round[z_i]>max(z_round[1:(z_i - 1)]), "Improvement", "No Improvement"))
    improv <- append(NA,improv) #last round had no choice
    improvement <- append(improvement, improv)
  }
  #Stay or go
  stay <- ifelse(distance==0, "Stay", "Go")
  #proportion of repeat/explore/exploit clicks
  propRepeat <- length(na.omit(exploreExploit[exploreExploit=="Repeat"]))/length(na.omit(exploreExploit))
  propExplore <- length(na.omit(exploreExploit[exploreExploit=="Explore"]))/length(na.omit(exploreExploit))
  propExploit <- length(na.omit(exploreExploit[exploreExploit=="Exploit"]))/length(na.omit(exploreExploit))
  #given improvement or no improvement, proportion of repeat/explore/exploit clicks
  improvements <- exploreExploit[improvement=="Improvement"]
  impPropRepeat <- length(na.omit(improvements[improvements=="Repeat"]))/length(na.omit(improvements))
  impPropExplore <- length(na.omit(improvements[improvements=="Explore"]))/length(na.omit(improvements))
  impPropExploit <- length(na.exclude(improvements[improvements=="Exploit"]))/length(na.omit(improvements))
  noImprovements <- exploreExploit[improvement=="No Improvement"]
  noimpPropRepeat <- length(na.omit(noImprovements[noImprovements=="Repeat"]))/length(na.omit(noImprovements))
  noimpPropExplore <- length(na.omit(noImprovements[noImprovements=="Explore"]))/length(na.omit(noImprovements))
  noimpPropExploit <- length(na.exclude(noImprovements[noImprovements=="Exploit"]))/length(na.omit(noImprovements))
  #dummy frame
  dummy<-data.frame(id, trial, x, y, z, kernel, scenario, horizon, round, environment, payoffScale, distance,propRepeat, propExplore, propExploit, deltaReward,deltaRewardPrev, stay, exploreExploit, improvement, impPropRepeat, impPropExplore, impPropExploit, noimpPropRepeat, noimpPropExplore, noimpPropExploit)
  #bind them together
  d<-rbind(d, dummy)
}

#create a unique listing for each space in the grid
allopts<-expand.grid(0:10, 0:10)
d$chosen<--99
for (i in 1:nrow(d)){ #loop through data and assign the unique listing from allopts to data$chosen based on the chosen x and y values
  d$chosen[i]<-which(d$x[i]==allopts$Var1 & d$y[i]==allopts$Var2)
}
#trial level data
tdata <- d


##############################################################################################################
##  3. Join together                                                                                        ##
##############################################################################################################
#Drop overlapping columns
#trial level data
drops <- c("kernel","scenario")
tdata <- tdata[ , !(names(tdata) %in% drops)]
#subject data
drops <- c("searchHistory", "horizon")
data <- data[ , !(names(data) %in% drops)]
#Merge
df <- merge(tdata, data, by = 'id')
#factor horizon
df$horizon<- as.factor(df$horizon)
#factor environment
df$environment <- as.factor(df$environment)

#add interaction of stay and kernel
df$stay.kernel <- interaction(df$stay, df$kernel)
#add interaction of kernel and scenario
df$kernel.scenario <- interaction(df$kernel, df$scenario)
#add interaction of kernel and environment
df$kernel.environment <- interaction(df$environment, df$kernel)
#add interaction of improvement and exploreExploit
df$improvExploreExploit <- interaction(df$exploreExploit, df$improvement)


##############################################################################################################
#Test Graphs
##############################################################################################################


#Explore Exploit behavior
ggplot(na.omit(df), aes(x = exploreExploit, fill = exploreExploit)) +  geom_bar(position='dodge') + facet_wrap(kernel.scenario~horizon)

# Stay - Go Behavior
#proportion of go and stay moves
bw = 10 #binwidth
ggplot(na.omit(df), aes(x = stay.kernel, fill = stay.kernel)) +  geom_bar(position='dodge') + facet_wrap(scenario~horizon)

#previous change in reward and stay/go behavior
ggplot(na.omit(df),aes(x=stay.kernel, y = deltaRewardPrev, fill = stay.kernel)) + geom_violin() + geom_boxplot(width = 0.1)   + facet_wrap(scenario~horizon) + theme_bw()

#distance moved
#scatter
ggplot(df, aes(x=trial, y = distance, col = kernel)) + geom_jitter(alpha = 0.1) + geom_smooth() + facet_wrap(scenario~horizon) + theme_bw()
#Density
ggplot(df, aes(x = distance, col = kernel)) + geom_density() + facet_wrap(scenario~horizon) + theme_bw() + ggtitle('Density Plot of Euclidian Distance Moved between trials')

#change in reward
ggplot(df, aes(x=deltaRewardPrev, y = distance, col = kernel)) + geom_jitter(alpha = 0.1) + geom_smooth() + facet_wrap(horizon~scenario) + theme_bw()

#change in reward x distance moved
ggplot(df, aes(x=distance, y = deltaReward, col = kernel)) + geom_jitter(alpha = 0.1) + geom_smooth() + facet_wrap(horizon~scenario) + theme_bw()

#X-y plot
ggplot(df, aes(x = x, y = y, col = scenario)) + geom_jitter(alpha = 0.1) + facet_wrap(~kernel.environment, nrow = 5) + theme_bw()

#reward over time
ggplot(df, aes(x=trial, y = z, col = kernel)) + geom_jitter(alpha = 0.1) + geom_smooth() + facet_wrap(scenario~horizon) + theme_bw()

#reward by environment
ggplot(df, aes(x = environment, y = reward, fill = environment)) + geom_violin() + geom_boxplot(width = 0.2) + facet_wrap(~kernel, nrow = 2)

#rewards histogram
ggplot(df, aes(x=reward, fill=horizon)) + geom_histogram(bins=10) + facet_wrap(kernel~scenario)

#Age
ggplot(data,aes(x=age, y=reward,fill = age)) + geom_boxplot() + facet_wrap(kernel~scenario)

#gender
ggplot(data,aes(x=gender, y=reward,fill = gender)) + geom_boxplot() + facet_wrap(kernel~scenario)

#scale
ggplot(df, aes(x=as.numeric(scale), y = z, col = scenario)) + geom_jitter(alpha = 0.1) + geom_smooth() + facet_wrap(kernel~horizon)

ggplot(data,aes(x=meanScale, y=reward)) + geom_point() + facet_wrap(kernel~scenario)
ggplot(data,aes(x=varScale, y=reward, col=horizon)) + geom_point() + facet_wrap(kernel~scenario)

