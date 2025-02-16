##### libraries ####
{
library(tidyverse)
library(readr)
library(sf)
library(dplyr)
library(raster)
library(RColorBrewer)
}
#####

#### import ####
{
main_path <- ".\\LSLS"
setwd(main_path)

source(".\\scenario_fnct.r") #functions

env_table <- read_csv(".\\tables\\env100_raster_table_.csv") #read pixel table with env variables
max_h <-  read_csv(".\\tables\\max_height_PA.csv") #pixel table with development zones
protected <- read_csv(".\\tables\\protected_.csv")
TAD <-  read_csv(".\\tables\\inTAD.csv")
stats <- read_csv(".\\tables\\stat_areas.csv") #statistical areas


#2050 projection (5 year interval)
trgts <- read_csv(".\\tables\\pop_projection_2050.csv")

#determine number of scenario iterations
k <- 10

#which scenarios to run
sce_names <- c("LShH","LSpH", "LSpF/30", "LSpL/30")

#define year time steps
years <-  c("y2025","y2030","y2035","y2040","y2045","y2050")

}
######

#### data prep ####
{
trgts <- column_to_rownames(trgts, var = "year")

env_table <- left_join(env_table, max_h[,c("FID", "max_h")], by="FID")  #assign max h to pixels 
env_table <- left_join(env_table, protected[,c("OBJECTID", "protected")],  by="OBJECTID")   
env_table <- left_join(env_table, TAD, by="OBJECTID")  
env_table <- left_join(env_table, stats[c("X","Y","stat","res_stat")], by=c("X","Y"))

#bldg index
env_table["bldg_vol_res"][env_table["bldg_vol_res"] == 0] <- NA # replace 0 with NA for calculation
env_table$bldg_grp <- ntile(env_table$bldg_vol_res,10) #divide built cells into 10 quantiles
env_table[is.na(env_table)] <- 0 # replace all NA with 0


}
#####


##### scenario run

#### LShH - Sharing
{
  scendata <- data.frame(scen = rep("LShH", k*length(years)), 
                              k = rep(1:k,each=length(years)),
                              year = rep(years,k))

for (i in 1:k){
  print(paste("starting LSh",i))
  
  sim_data <- data.frame()
  change <- numeric()
  #create path
  scen_path <- paste0(main_path,"/scenarios/LShH/iter_",i)
  dir.create(scen_path, showWarnings = FALSE)
  
   #simulate scenarios 
 env_LShH <- list(y2020=env_table)

 for (y in 1:length(years)){
   print(paste("starting",years[y]))
   
   res_list <- LShH (env_LShH[[y]] , trgts[y+1,"bldg_vol"]) #call the simulation function
   
   #store results
   env_LShH[[y+1]] <- res_list[["env_new"]]
   names(env_LShH)[y+1] <- years[y]
   sim_data <- rbind(sim_data,data.frame(year=years[y],res_list$sim.data))
   change <- append(change, sum(env_LShH[[y+1]]$change)-sum(change))
   }
 
 #metadata
 scendata[which(scendata$k == i),names(sim_data)] <- as.list(sim_data)
 scendata[which(scendata$k == i),"changed"] <- change
 
 #summarize
 sum_df <- scen_sum(env_LShH)
 write_csv(sum_df, paste0(scen_path,"/summary.csv"))

 #export
 env_LShH <- env_LShH[-1]
 for(l in 1:length(env_LShH)){
     scen_export(env_LShH[l], scen_path)
 }
 
}
#export metadata
write_csv(scendata, paste0(main_path,"/scenarios/LShH/scen_metadata.csv"))
}

##### LSpH - Sparing
{
#increase max height by 10% to allow meeting housing targets
env_table2 <- mutate(env_table, max_h=max_h*1.1)

scendata <- data.frame(scen = rep("LSpH", k*length(years)), 
                            k = rep(1:k,each=length(years)))

for (i in 1:k){
  print(paste("starting LSP",i))
  sim_data <- data.frame()
  change <- numeric()
  #create path
  scen_path <- paste0(main_path,"/scenarios/LSpH/iter_",i)
  dir.create(scen_path, showWarnings = FALSE)
  
  #simulate scenarios 
  env_LSpH <- list(y2020=env_table2)
  
  for (y in 1:length(years)){
    print(paste("starting",years[y]))
    
    res_list <- LSpH (env_LSpH[[y]] , trgts[y+1,"bldg_vol"],0) #call the simulation function
    
    #store results
    env_LSpH[[y+1]] <- res_list$env_new 
    names(env_LSpH)[y+1] <- years[y]
    sim_data <- rbind(sim_data,data.frame(year=years[y],res_list$sim.data))
    change_sum <- sum(subset(env_LSpH[[y+1]], spare == 0)$change)
    change <- append(change, change_sum-sum(change))
  }

  #metadata
  scendata[which(scendata$k == i),names(sim_data)] <- as.list(sim_data)
  scendata[which(scendata$k == i),"changed"] <- change
  
  #summarize
sum_df <- scen_sum(env_LSpH)
write_csv(sum_df, paste0(scen_path,"/summary.csv"))

#export
env_LSpH <- env_LSpH[-1]
for(l in 1:length(env_LSpH)){
  scen_export(env_LSpH[l], scen_path)
}
}
#export metadata
write_csv(scendata, paste0(main_path,"/scenarios/LSpH/scen_metadata.csv"))
}

##### LSpF - Full Sparing (Regional)

## 30% sparing by 2050
{

  scendata <- data.frame(scen = rep("LSpF/30", k*length(years)), k = rep(1:k,each=length(years)))

for (i in 1:k){
  print(paste("starting LSF/30",i))
  sim_data <- data.frame()
  change <- numeric()
  #create path
  scen_path <- paste0(main_path,"/scenarios/LSpF/30/iter_",i)
  dir.create(scen_path, recursive = TRUE,showWarnings = TRUE)
  
  #simulate scenarios 
  env_LSpF <- list(y2020=env_table)
  
  for (y in 1:length(years)){
    print(paste("starting",years[y]))
    
    res_list <- LSpF (env_LSpF[[y]] , trgts[y+1,"bldg_vol"],0.05) #call the simulation function
    
    #store results
    env_LSpF[[y+1]] <- res_list$env_new 
    names(env_LSpF)[y+1] <- years[y]
    sim_data <- rbind(sim_data,data.frame(year=years[y],res_list$sim.data))
    change_sum <- sum(subset(env_LSpF[[y+1]], spare == 0)$change)
    change <- append(change, change_sum-sum(change))
  }
  
  #metadata
  scendata[which(scendata$k == i),names(sim_data)] <- as.list(sim_data)
  scendata[which(scendata$k == i),"changed"] <- change
  
  #summarize
  sum_df <- scen_sum(env_LSpF)
  write_csv(sum_df, paste0(scen_path,"/summary.csv"))
  
  #export
  env_LSpF <- env_LSpF[-1]
  for(l in 1:length(env_LSpF)){
    scen_export(env_LSpF[l], scen_path)
  }
}

#export metadata
write_csv(scendata, paste0(main_path,"/scenarios/LSpF/30/scen_metadata.csv"))
}


##### LSpL - Local Sparing (Neighbourhood)
## 30% sparing by 2050
{
  scendata <- data.frame(scen = rep("LSpL/30", k*length(years)), k = rep(1:k,each=length(years)))

for (i in 1:k){
  print(paste("starting LSpL/30",i))
  sim_data <- data.frame()
  change <- numeric()
  #create path
  scen_path <- paste0(main_path,"/scenarios/LSpL/30/iter_",i)
  dir.create(scen_path, recursive = TRUE,showWarnings = TRUE)
  
  #simulate scenarios 
  env_LSpL <- list(y2020=env_table)
  
  for (y in 1:length(years)){
    print(paste("starting",years[y]))
    
    res_list <- LSpL (env_LSpL[[y]] , trgts[y+1,"bldg_vol"],0.05) #call the simulation function
    
    #store results
    env_LSpL[[y+1]] <- res_list$env_new 
    names(env_LSpL)[y+1] <- years[y]
    sim_data <- rbind(sim_data,data.frame(year=years[y],res_list$sim.data))
    change_sum <- sum(subset(env_LSpL[[y+1]], spare == 0)$change)
    change <- append(change, change_sum-sum(change))
  }
  
  #metadata
  scendata[which(scendata$k == i),names(sim_data)] <- as.list(sim_data)
  scendata[which(scendata$k == i),"changed"] <- change
  
  #summarize
  sum_df <- scen_sum(env_LSpL)
  write_csv(sum_df, paste0(scen_path,"/summary.csv"))
  
  #export
  env_LSpL <- env_LSpL[-1]
  for(l in 1:length(env_LSpL)){
    scen_export(env_LSpL[l], scen_path)
  }
}

#export metadata
write_csv(scendata, paste0(main_path,"/scenarios/LSpL/30/scen_metadata.csv"))
}

