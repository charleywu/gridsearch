#Charley Wu 2018
#Read and process participant data

#load packages
packages <- c('tidyr', 'dplyr', 'plyr', 'jsonlite')
lapply(packages, require, character.only = TRUE)


#Main function
dataImport <- function(normalize=TRUE){
  
  #read in data
  dat<-read.csv("ExperimentData/experimentData2D.csv")
  
  #dummy data frame
  data<-data.frame(id=numeric(), trial=numeric(), x=numeric(), y=numeric(), 
                   z=numeric(), zmax=numeric(), kernel=numeric(), scenario=numeric(), 
                   horizon=numeric(), round=numeric(), env=numeric(), delta_x=numeric())
  
  #Compile experiment  data
  for (i in 1:nrow(dat)){
    #parse JSON
    searchHistory <- fromJSON(as.character(dat$searchHistory[i]))
    #sampled x value
    x<-unlist(searchHistory$xcollect)
    #sampled y value
    y<-unlist(searchHistory$ycollect)
    #sampled z value (UNSCALED!)
    z<-unlist(searchHistory$zcollect)
    #max reward found at each time point
    zmax<-unlist(sapply(searchHistory$zcollect, maxton))
    if (normalize==TRUE){
      #normalize z and zmax
      z <- (z-50)/100
      zmax <- (zmax-50)/100
    }
    #length of trial
    len<-sapply(searchHistory[1]$xcollect, length)
    #trial 1:20 or 1:40
    trial<-c(1:len[1],1:len[2],1:len[3],1:len[4])
    #round, 1:8
    round<-rep(1:8, len)
    #env number
    env <- rep(fromJSON(as.character(dat$envOrder[i])),len)
    #horizon, 20 or 40
    horizon<-rep(len-1, len)
    #payoff structure: 0 = cumulative reward; 1= maximum value
    scenario<-rep(dat$scenario[i], sum(len))
    #smoothness: 0 = rough (length scale = 1); 1 = smooth (length scale = 2)
    kernel<-rep(dat$kernel[i], sum(len))
    #id number for each subject
    id<-rep(i, sum(len))
    #dummy frame
    dummy<-data.frame(id, trial, x, y, z, zmax, kernel, scenario, horizon, round, env)
    #calculate manhattan distance between clicks
    #calculate distance between clicks
    dummy <- dummy %>%
      group_by(round) %>%
      mutate(delta_x = abs(x - lag(x, default = NA)) + abs(y - lag(y, default = NA)) ) 
    dummy <- as.data.frame(dummy)
    dummy$delta_x[dummy$trial==1]<-NA #set as NA for all first clicks, since it was randomly selected
    #bind them together
    data<-rbind(data, dummy)
  }
  
  #create a unique listing for each space in the grid
  allopts<-expand.grid(0:10, 0:10)
  data$chosen<--99
  for (i in 1:nrow(data)){ #loop through data and assign the unique listing from allopts to data$chosen based on the chosen x and y values
    data$chosen[i]<-which(data$x[i]==allopts$Var1 & data$y[i]==allopts$Var2)
  }
  #create an overall id, unique for each subject for each round
  data$overall<-1:nrow(data)-with(data, ave(id, paste(id, round), FUN = seq_along))
  data<-transform(data,overall=as.numeric(factor(overall)))
  data$round <- factor(data$round) #convert round into factor
  #Rename variables
  data$kernel<-ifelse(data$kernel==0, "Rough", "Smooth")
  #recode horizon
  data$Horizon<-ifelse(data$horizon==20, "Short", "Long")
  #recode reward function
  data$scenario<-ifelse(data$scenario==0, "Avg. Reward Condition", "Max. Reward Condition")
  
  #return it
  return(data)
}

#Maximum reward up until trial x
maxton<-function(x){
  maxn<-rep(0,length(x))
  maxn[1]<-x[1]
  for (i in 2:length(x)){
    if (x[i]>maxn[i-1]){maxn[i]<-x[i]}
    else{maxn[i]<-maxn[i-1]}
  }
  return(maxn)
}