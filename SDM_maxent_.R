
##### Libraries #####
{
library(dismo)
library(readr)
library(sf)
library(tidyverse)
library(ENMeval)
library(ecospat)
library(blockCV)
library(terra)
library(raster)
library(trend)
  
}

######

##### IMPORT #####
{
setwd(".\\LSLS")

#functions
source(".\\SDM_fnct.r") #source for get_focals function

# load species list
sp_data <- read_csv(".\\BirdGroups_Final.csv")

# access presence data and locations (and identify absences too using the global pc database)
pc_data <- read_csv(".\\birdsurvey_alldata_clean.csv")
pc_xy <- st_read(".\\survey_points_new.shp")

# predictor layers
predictors <- stack(list.files(".\\100_grid_variables", pattern='tif$', full.names=TRUE ))
r_mask <-  raster(".\\100_grid_variables\\bldg_cover.tif")
  
# number of model folds
k_fold <- 5

#cross validation iterations
k_cross <- 5
}
#####

##### Preprocess ####
{
#add 500m predictors (focal)
focal_list <- c( "bldg_cover", "low_veg", "tree_cover")
predictors <- get_focals(predictors,focal_list)

#import xy coordinates to pc_data
pc_xy <- cbind(pc_xy, st_coordinates(pc_xy))
pc_data <- subset(pc_data, select=-c(x,y))%>% 
  left_join(x=., y=pc_xy[,c("pointid","X","Y")], by=c("PointID"="pointid"))%>%
  filter(!is.na(X))

#exclude points outside boundaries
temp <- cbind(pc_data,data.frame(raster::extract(predictors$bldg_cover_500,pc_data[,c("X","Y")]))) 
temp <- temp[!is.na(temp[,ncol(temp)]), ] 

pc_data <- filter(pc_data, PointID %in% temp$PointID)

# exclude rare species (<15 points)

xspecies <- pc_data %>% group_by(Sp)%>% 
  summarise(count = n_distinct(PointID))%>%
  filter(count>=15) 

}
#######

##### SDM run ####

### cross-validation loop ###
for(i in 1:length(xspecies$Sp)){
  print(paste("Starting sp",i, xspecies$Sp[i],"..."))
  dir.create(paste0("./Maxent/y2020/",xspecies$Sp[i]), showWarnings = FALSE)
  
    # get occurrence locations:
  # filter for the species
  pc_spec <- subset(pc_data, Sp == xspecies$Sp[i])  
  #only unique cases (to avoid pseudoreps)
  sp_ID <- pc_spec[,c("PointID","X","Y")] %>% distinct(PointID, .keep_all = TRUE)
  
  #cross-validation repetitions
  boyce_list <- list() 
  for (c in 1:k_cross){
    sp_boyce <- c()
    #create SpatialFolds
    #create pa file with coordinates of presences and random background points
  occ <- data.frame(x=sp_ID$X,y=sp_ID$Y)
  temp <- length(data.frame(rasterToPoints(predictors[[1]], na.rm=TRUE))$x) #how many pixels are there in the study area
  acc <- as.data.frame(randomPoints(predictors[[1]],10000,occ)) #select 10000 random background pixels
  
  occ$presence <- rep(1,nrow(occ))
  acc$presence <- rep(0,nrow(acc))
  
  pa <-rbind(occ,acc)
  
  #add environmental information to pa file
  pa <- cbind(pa,data.frame(raster::extract(predictors,pa[,c("x","y")])))
  
  #create the spatial folds
  pa_sv <- vect(pa, geom = c("x", "y"),crs=raster::crs(predictors)) #change to SpatVector
  block_size <- 5000

  #blocks <- spatialBlock(speciesData = as(pa_sv, "Spatial"), species = "presence", theRange = block_size, k = k_fold, selection = "random", showBlocks = FALSE)
  blocks <-  cv_spatial(x = as(pa_sv, "Spatial"), column = "presence", size = block_size, k = k_fold, selection = "random", plot = FALSE)
  
    # add spatial block ID:
  pa$foldID <- blocks$folds_ids
  
  
  # split data into training and test subset
  for(k in 1:k_fold){   #5-fold
    print(paste("model fold",k))
    
    dir.create(paste0("./Maxent/y2020/",xspecies$Sp[i],"/rep_",c,"/iter_",k), recursive = TRUE, showWarnings = FALSE)

    pa_train <- pa[pa$foldID != k, ]
    pa_test <- pa[pa$foldID == k, ]
    
    #go to the next fold if there are no presences in the test set
    if (sum(pa_test$presence[which(pa_test$presence == 1 )]) == 0){
      sp_boyce <- append(sp_boyce,NA)
      next
    } else {
    
    ### do the SDM
    ## with maxent
    mxnt_env <- pa_train[, names(pa_train) %in% names(predictors)] #only predictor columns
    mxnt_p <-pa_train$presence
    mxnt_sp <- maxent(mxnt_env, p = mxnt_p,
                      path = paste0(getwd(),"/Maxent/y2020/",xspecies$Sp[i],"/rep_",c,"/iter_",k))
    p_clog <- predict(object = mxnt_sp, x= predictors, args=c("outputformat=cloglog"))


    ## evaluation with p/a
    pa_test$hs <- raster::extract(p_clog,pa_test[c(1:2)])
    mxnt_eval_p <- pa_test[pa_test$presence ==1,]$hs
    mxnt_eval_a <- pa_test[pa_test$presence ==0,]$hs
    e2 <- evaluate(mxnt_sp,p=mxnt_eval_p, a=mxnt_eval_a)
    TSS <- max(e2@TPR-e2@FPR)

    ## Boyce index
    boyce_k <- ecospat.boyce(fit=na.omit(getValues(p_clog)), obs=na.omit(mxnt_eval_p),
                              rm.duplicate = T, method='kendall')
    boyce_df <- data.frame( HS = boyce_k$HS, F.ratio = boyce_k$F.ratio)
    res <- rbind(data.frame(mxnt_sp@results), TSS = TSS, boyce = boyce_k$cor)
    sp_boyce <- append(sp_boyce,boyce_k$cor)

    ## Export
    #results table
    write.csv(res, paste0("./Maxent/y2020/",xspecies$Sp[i],"/rep_",c,"/iter_",k,"/maxent_results.csv"))
    #Boyce values
    write.csv(boyce_df, paste0("./Maxent/y2020/",xspecies$Sp[i],"/rep_",c,"/iter_",k,"/boyce.csv"))
    #evaluation object
    saveRDS(e2, paste0("./Maxent/y2020/",xspecies$Sp[i],"/rep_",c,"/iter_",k,"/evaluate.rda"))
    #model object
    saveRDS(mxnt_sp, paste0("./Maxent/y2020/",xspecies$Sp[i],"/rep_",c,"/iter_",k,"/maxent.rda"))
    #prediction raster
    p_clog <- raster::crop(p_clog, extent(r_mask)) 
    p_clog <- raster::mask(p_clog, r_mask) #crop and mask to match research area boundaries
    raster::writeRaster(p_clog, paste0("./Maxent/y2020/",xspecies$Sp[i],"/rep_",c,"/iter_",k,"/prediction.grd"),overwrite=TRUE)
      }
  }
  boyce_c <-  data.frame(sp = xspecies$Sp[i],
                         rep=c,
                           "X1" = sp_boyce[1],
                           "X2" = sp_boyce[2],
                           "X3" = sp_boyce[3],
                           "X4" = sp_boyce[4],
                           "X5" = sp_boyce[5],
                           mean=mean(na.omit(sp_boyce)),
                           sd=sd(na.omit(sp_boyce)))
  write.csv(boyce_c,paste0("./Maxent/y2020/",xspecies$Sp[i],"/rep_",c,"/boyce_sum.csv") )
  
  boyce_list[[c]]  <- boyce_c
  }
  
  boyce_sum <-  do.call("rbind",boyce_list)
  write.csv(boyce_sum, paste0("./Maxent/y2020/",xspecies$Sp[i],"/boyce_sum.csv"))
}

# dataframe for list of species
{columns <-  c("Sp",c(paste0("M",as.character(1:25))),"exclude") 
sp_df <-  data.frame(matrix(nrow = 0, ncol = length(columns))) 
colnames(sp_df) <-  columns}

#### evaluation loop ####
for(i in 1:length(xspecies$Sp)){

   #get summary of cross validation Boyce index
  boyce_sum <-  read.csv(paste0("./Maxent/y2020/",xspecies$Sp[i],"/boyce_sum.csv"))
   
  # Calculate means and sd as the runs go
  x_values <- boyce_sum[,-1] %>% select(starts_with("X")) 
  x_values <-  c(t(x_values))
  
  cumulative_means <- sapply(1:length(x_values), function(x) mean(x_values[1:x],na.rm=TRUE))
  cumulative_sd <- sapply(1:length(x_values), function(x) sd(x_values[1:x],na.rm=TRUE))
  
  # Vector indicating whether mean is larger than sd
  new_vector <- sapply(1:length(cumulative_means), function(x) {
    if (is.na(cumulative_sd[x])) {
      return(NA)
    } else if (cumulative_means[x] >  cumulative_sd[x]) {
      return(1)
    } else {
      return(0)
    }})
  
  cumulative_target <- sapply(1:length(x_values), function(x) mean(new_vector[1:x], na.rm = TRUE))
  
  # Graphical representation
  {
    # Create a data frame for plotting
    cumulative_data <- data.frame(
      Runs = 1:length(cumulative_means),
      CumulativeMeans = cumulative_means,
      CumulativeSD = cumulative_sd,
      NewVector = new_vector
    )
    
    # Plot cumulative means and SD with y-axis limits set between -1 and 1
    p1 <- ggplot(cumulative_data, aes(x = Runs)) +
      geom_line(aes(y = CumulativeMeans, color = "Mean")) +
      geom_line(aes(y = CumulativeSD, color = "SD")) +
      geom_point(aes(y = CumulativeMeans, color = ifelse(NewVector == 1, "MeanPoint", "MeanPoint"))) +
      geom_point(aes(y = CumulativeSD, color = ifelse(NewVector == 1, "SDPoint", "SDPoint"))) +
      scale_color_manual(name = "Legend",
                         values = c("Mean" = "red", "SD" = "black", "MeanPoint" = "red", "SDPoint" = "black"),
                         labels = c("Mean", "SD")) +
      labs(title = paste("Cumulative Means and SD vs Number of Runs for", xspecies$Sp[i]),
           x = "Number of Runs",
           y = "Value") +
      ylim(-1, 1) +
      theme_minimal() +
      theme(legend.position = "bottom")
    
    # Save the first plot
    ggsave(filename = file.path(paste0("./Maxent/y2020/",xspecies$Sp[i],"/Cumulative_Means_SD.jpg")), plot = p1, device = "jpeg")
    
    # Create a data frame for plotting
    cumulative_data <- data.frame(
      Runs = 1:length(cumulative_target),
      CumulativeTarget = cumulative_target
    )
    
    # Plot cumulative target with y-axis limits set between 0 and 1
    p2 <- ggplot(cumulative_data, aes(x = Runs, y = CumulativeTarget)) +
      geom_line(color = "blue") +
      geom_point(color = "red") +
      labs(title = paste("Cumulative Target vs Number of Runs for", xspecies$Sp[i]),
           x = "Number of Runs",
           y = "Cumulative Target") +
      ylim(0, 1) +
      theme_minimal()
    
    # Save the second plot
    ggsave(filename = file.path(paste0("./Maxent/y2020/",xspecies$Sp[i],"/Cumulative_Target.jpg")), plot = p2, device = "jpeg")
  }
  
  # Remove NaN/NA values
  cumulative_target_nonNA <- cumulative_target[!is.na(cumulative_target)]
  
  # Check for variation and perform Mann-Kendall Trend Test
  if (all(cumulative_target_nonNA == 1)) {
    exclude <- 0
  } else if (all(cumulative_target_nonNA == 0)) {
    exclude <- 1
  } else {
    mk_test <- mk.test(cumulative_target_nonNA)
    p_value <- mk_test$p.value
    tau <- mk_test$estimate[3]
    
    if (p_value < 0.05 && tau > 0) {
      exclude <- 0
    } else {
      exclude <- 1
    }
  }
  
    sp_df <- rbind(sp_df,c(xspecies$Sp[i], x_values, exclude))
    names(sp_df) <- columns
}

write.csv(sp_df, paste0("./Maxent/results/boyce_summary_",
                                                format(Sys.time(), "%Y-%m-%d"),".csv"  ))

#### Final model loop ####

sp_df2 <- read.csv(paste0("./Maxent/results/boyce_summary_",
                                    "2024-08-14",".csv"  ))#change if needed
sp_list <- filter(sp_df2, exclude == 0 )
length(sp_list$Sp)

for(i in 1:length(sp_list$Sp)){
    print(paste("Sp",i, sp_list$Sp[i],"final model..."))
  
  #1. fit the model
   # get sp occurrence locations:
   pc_spec <- subset(pc_data, Sp == sp_list$Sp[i])  
   sp_ID <- pc_spec[,c("PointID","X","Y")] %>% 
     distinct(PointID, .keep_all = TRUE) #only unique cases (to avoid pseudoreps)
   occ <- data.frame(x=sp_ID$X,y=sp_ID$Y)
   temp <- length(data.frame(rasterToPoints(predictors[[1]], na.rm=TRUE))$x) #how many pixels are there in the study area
   acc <- as.data.frame(randomPoints(predictors[[1]],10000,occ)) #select 10000 random background pixels
   occ$presence <- rep(1,nrow(occ))
   acc$presence <- rep(0,nrow(acc))
   pa <- rbind(occ,acc) 
   
   #add environmental information to pa file
   pa <- cbind(pa,data.frame(raster::extract(predictors,pa[,c("x","y")])))
   
    mxnt_env <- pa[, names(pa) %in% names(predictors)] #only predictor columns
    mxnt_p <-pa$presence
    mxnt_sp <- maxent(mxnt_env, p = mxnt_p,
                      path = paste0(getwd(),"/Maxent/y2020/",sp_list$Sp[i]))
    p_clog <- predict(object = mxnt_sp, x= predictors, args=c("outputformat=cloglog"))
    
    
    #3. evaluation metrics
    pa$hs <- raster::extract(p_clog,pa[c(1:2)])
    mxnt_eval_p <- pa[pa$presence ==1,]$hs
    mxnt_eval_a <- pa[pa$presence ==0,]$hs
    e2 <- evaluate(mxnt_sp,p=mxnt_eval_p, a=mxnt_eval_a)
    TSS <- max(e2@TPR-e2@FPR)
    
    boyce_k <- ecospat.boyce(fit=na.omit(getValues(p_clog)), obs=na.omit(mxnt_eval_p),
                             rm.duplicate = T, method='kendall')
    boyce_df <- data.frame( HS = boyce_k$HS, F.ratio = boyce_k$F.ratio)
    res <- rbind(data.frame(mxnt_sp@results), TSS = TSS, boyce = boyce_k$cor)
    sp_df2[sp_df2$Sp == sp_list$Sp[i],"final"] <- boyce_k$cor

    #4. Export
     #results table
    write.csv(res, paste0("./Maxent/y2020/",sp_list$Sp[i],"/maxent_results.csv"))
     #Boyce values
    write.csv(boyce_df, paste0("./Maxent/y2020/",sp_list$Sp[i],"/boyce.csv"))
     #evaluation object
    saveRDS(e2, paste0("./Maxent/y2020/",sp_list$Sp[i],"/evaluate.rda"))
     #model object
    saveRDS(mxnt_sp, paste0("./Maxent/y2020/",sp_list$Sp[i],"/maxent.rda"))
     #prediction raster
    p_clog <- raster::crop(p_clog, extent(r_mask)) 
    p_clog <- raster::mask(p_clog, r_mask) #crop and mask to match research area boundaries
    raster::writeRaster(p_clog, paste0("./Maxent/y2020/",sp_list$Sp[i],"/prediction.grd"),overwrite=TRUE)
    
   #5. response curve data
    dir.create(paste0("./Maxent/y2020/",sp_list$Sp[i],"/responses"), showWarnings = FALSE)
    for(o in 1:length(names(predictors))){
      r <- response(mxnt_sp, var = names(predictors)[o],  expand = 0)
      write.csv(r, paste0("./Maxent/y2020/",sp_list$Sp[i],"/responses","/response_",names(predictors)[o],".csv"))
    }

    print(paste("species",i, sp_list$Sp[i] ,"complete, Boyce index=", boyce_k$cor))

  }

write.csv(sp_df2, paste0("./Maxent/results/final_boyce_summary_",
                        format(Sys.time(), "%Y-%m-%d"),".csv"  ))

#####


############# Scenario SDMs ##################

##### produce predictions ######
sp_list <- filter(sp_df2, exclude == 0 )

sce_names <- c("LSpH", "LShH","LSpF/30","LSpL/30")
years <- c("y2025","y2030","y2035","y2040","y2045","y2050")

k_scen <- 10
#focal_list <- c( "bldg_cover", "UGS","tree_cover","low_veg" )
focal_list <- c( "bldg_cover")


# loop through scenarios -> iterations -> years -> species to produce new predictions

for (i in 1:length(sce_names)){ #scenario types
  dir.create(paste0("./Maxent/",sce_names[i]), recursive = TRUE, showWarnings = F)
  
  for (j in 1:length(years)){ #years
    dir.create(paste0("./Maxent/",sce_names[i],"/",years[j]), showWarnings = FALSE)
    
    for (k in 1:k_scen){ #scenario & SDM iterations
      print(paste("Starting", sce_names[i],years[j],"iterarion",k ))
      
      path_k <- paste0("./Maxent/",sce_names[i],"/",years[j],"/iter_",k)
      dir.create(path_k, showWarnings = FALSE)
      
      ### import and process predictors for scenario k-th iteration  
      #create predictor list
      pred_list <-  list.files(paste0("./scenarios/",sce_names[i],"/iter_",k,"/",years[j]), pattern='tif$', full.names=FALSE )
      pred_list <-  tools::file_path_sans_ext(pred_list)
      pred_list <- pred_list[pred_list!="bldg_vol"]
      pred_list <- pred_list[pred_list!="UGS"] #exclude UGS from predictors
      
      
      #read and add buffer to allow 500m calculation
      predictors_s <- stack()
      for (r in 1:length(pred_list)){
        r_name <- pred_list[r]
        rast <- raster(paste0("./scenarios/",sce_names[i],"/iter_",k,"/",years[j],"/",r_name,".tif"))
        raster::crs(rast) <- raster::crs(predictors[[r_name]]) #match crs
        rast <- raster::extend(rast,extent(predictors[[r_name]])) #extend extent
        rast <- raster::cover(rast, predictors[[r_name]]) #get buffer values
        predictors_s <- stack(predictors_s, rast) #stack raster
        names(predictors_s)[r] <- r_name #set variable name
      }
      
      #add coast
      predictors_s <- stack(predictors_s, predictors$coast)
      #add focals
      predictors_s <- get_focals(predictors_s,focal_list)
      #crop and mask to match research area boundaries
      predictors_s <- raster::crop(predictors_s, extent(r_mask)) 
      predictors_s <- raster::mask(predictors_s, r_mask)
      
      if (any(!names(predictors) %in% names(predictors_s))) { 
        print ("predictors incomplete! missing:")
        print (setdiff(names(predictors), names(predictors_s)))
      }
     
       for (s in 1:length(sp_list$Sp)){ #species
        #for (m in 1:k_fold){ #maxent iterations
          
         #import model 
         mxnt_sp <-  readRDS(paste0("./Maxent/y2020/",sp_list$Sp[s],"/maxent.rda"))
         
        #predict for new predictors
         p_clog_s <- predict(object = mxnt_sp, x= predictors_s, args=c("outputformat=cloglog"))
         raster::crs(p_clog_s) <- raster::crs(predictors_s)
        
         #export prediction map
         raster::writeRaster(p_clog_s, paste0(path_k,"/",sp_list$Sp[s],".grd"),overwrite=TRUE)
         
         print(paste("species",s,sp_list$Sp[s],"exported"))
       }
      }
    print(paste(sce_names[i],years[j], "completed"))
  }
}





