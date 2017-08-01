#Read and process participant data
#Charley Wu 2017

#load packages
packages <- c('plyr', 'jsonlite')
lapply(packages, require, character.only = TRUE)

#Main data import function
dataImport <- function(normalize=TRUE){
  #read in data
  dat<-read.csv("ExperimentData/experimentData1D.csv")
  
  #dummy data frame
  data<-data.frame(id=numeric(), trial=numeric(), x=numeric(), y=numeric(), 
                  ymax=numeric(), kernel=numeric(), scenario=numeric(), 
                   horizon=numeric(), round=numeric(), reward = numeric(), experimentDuration = numeric(), instructionDuration=numeric(), env=numeric())
  
  #Compile experiment  data
  for (i in 1:nrow(dat)){
    #parse JSON
    searchHistory <- fromJSON(as.character(dat$searchHistory[i]))
    #sampled x value
    x<-unlist(searchHistory$xcollect)
    #sampled z value (UNSCALED!)
    y<-unlist(searchHistory$ycollect)
    #max reward found at each time point
    ymax<-unlist(sapply(searchHistory$ycollect, maxton))
    if (normalize==TRUE){
      #normalize z and zmax
      y <- (y-50)/100
      ymax <- (ymax-50)/100
    }
    #length of trial
    len<-sapply(searchHistory[1]$xcollect, length)
    #trial lengths (i.e., horizon)
    trial<-c(c(1:len[1],1:len[2],1:len[3],1:len[4],1:len[5],1:len[6],1:len[7],1:len[8]),c(1:len[1],1:len[2],1:len[3],1:len[4],1:len[5],1:len[6],1:len[7],1:len[8]))
    #round, 1:8
    round<-rep(1:16, len)
    #Env number
    env <- rep(fromJSON(as.character(dat$envOrder[i])),len)
    #horizon, 20 or 40
    horizon<-rep(len-1, len)
    #payoff structure: 0 = cumulative reward; 1= maximum value
    scenario<-rep(dat$scenario[i], sum(len))
    #smoothness: 0 = rough (length scale = 1); 1 = smooth (length scale = 2)
    kernel<-rep(dat$kernel[i], sum(len))
    #id number for each subject
    id<-rep(i, sum(len))
    #performance bonus
    reward <- rep(dat$reward[i], sum(len))
    #timing
    experimentDuration <- rep(dat$ExperimentDuration[i], sum(len))
    instructionsDuration <- rep(dat$instructionDuration[i], sum(len))
    #dummy frame
    dummy<-data.frame(id, trial, x, y, ymax, kernel, scenario, horizon, round, reward, experimentDuration, instructionsDuration, env)
    #bind them together
    data<-rbind(data, dummy)
  }
  
  #create an overall id, unique for each subject for each round
  data$overall<-1:nrow(data)-with(data, ave(id, paste(id, round), FUN = seq_along))
  data<-transform(data,overall=as.numeric(factor(overall)))
  data$round <- factor(data$round) #convert round into factor
  #Rename variables
  data$kernel<-ifelse(data$kernel==0, "Rough", "Smooth")
  #recode horizon
  data$Horizon<-ifelse(data$horizon==5, "Short", "Long")
  #recode reward function
  data$scenario<-ifelse(data$scenario==0, "Avg. Reward Condition", "Max. Reward Condition")
  
  #return it
  return(data)
}

#Max reward up until a certain trial number
maxton<-function(x){
  maxn<-rep(0,length(x))
  maxn[1]<-x[1]
  for (i in 2:length(x)){
    if (x[i]>maxn[i-1]){maxn[i]<-x[i]}
    else{maxn[i]<-maxn[i-1]}
  }
  return(maxn)
}