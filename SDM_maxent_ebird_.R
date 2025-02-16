##### clear ####
rm(list = ls())

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
library(GeoThinneR)  
}

######

##### IMPORT #####
{
setwd(".\\ebird")

#functions
source(".\\SDM_fnct.r") #source for get_focals function

# load species list
sp_data <- read_csv(".\\BirdGroups_Final.csv")


ebird_pc <- read_csv(".\\ebird_2024.csv")
ebird_xy <- st_read(".\\eBird_points_filtered_ITM.shp")
  
# predictor layers
predictors <- stack(list.files(".\\100_grid_variables", pattern='tif$', full.names=TRUE ))
r_mask <-  raster(".\\100_grid_variables\\New\\bldg_cover.tif")
bias_extended <- raster(".\\all_bias_log10.tif")
  
# number of model folds
k_fold <- 5

#cross validation iterations
k_cross <- 5
}
#####

ebird_final_pc <- read_csv(".\\Bird datasets\\ebird\\2024\\ebird_2024_thinned_final.csv")

##### Preprocess ####
{
  
## adjust bias layer

bias_proj <- project(rast(bias_extended), rast(predictors[[1]]))
bias_proj <- raster(bias_proj)
bias_crop <- raster::crop(bias_proj,predictors[[1]] )
bias_TAD <- raster::mask(bias_crop,predictors[[1]] )
  
#add 500m predictors (focal)
focal_list <- c( "bldg_cover", "low_veg", "tree_cover")
predictors <- get_focals(predictors,focal_list)

}


#### ebird data prep ###
{
#import xy coordinates to pc_data
ebird_xy <-  subset(ebird_xy, select = -c(X,Y)) %>%
  cbind(st_coordinates(ebird_xy))%>%
  mutate(IDSourc = as.numeric(IDSourc))
ebird_pc <- subset(ebird_pc, select=-c(X,Y))%>%
  left_join(x=., y=ebird_xy[,c("IDSourc","X","Y")], by="IDSourc")%>%
  filter(!is.na(X)) %>%
  rename("ID" = "IDSourc")

#exclude points outside boundaries
temp <- cbind(ebird_pc,data.frame(raster::extract(predictors$bldg_cover_500,ebird_pc[,c("X","Y")])))
temp <- temp[!is.na(temp[,ncol(temp)]), ]
ebird_pc <- filter(ebird_pc, ID %in% temp$ID)

## filter observations:

#get spp list
ebird_sp <- unique(ebird_pc$Binomil)
ebird_pc_filtered <- data.frame()

# Reduce observations to one per grid cell
for(i in 1:length(ebird_sp)){
  pc_spec <- filter(ebird_pc, Binomil == ebird_sp[i])

  # Extract the cell number for each point using the raster
  coords <- cbind(pc_spec$X, pc_spec$Y)
  pc_spec$cell <- raster::cellFromXY(predictors[[1]], coords)
  # Filter to keep only one point per unique cell
  sp_ID <- pc_spec %>%
    group_by(cell) %>%
    slice(1) %>%
    ungroup() %>%
    select(-cell)
  # Accumulate the results
  ebird_pc_filtered  <- rbind(ebird_pc_filtered, sp_ID)
}

# Thin occurrences in non-neighboring grid cells

thinned_data_list <- list()

# Loop over each species
for (i in seq_along(ebird_sp)) {
  species <- ebird_sp[i]
  cat("Processing species:", species, "(", i, "of", length(ebird_sp), ")\n")
  # Subset the data for the current species
  species_data <- ebird_pc_filtered[ebird_pc_filtered$Binomil == species, ]
  # Prepare the data for thinning
  species_data <- species_data[, c("X", "Y", "Binomil")]

# Perform thinning using the thin_points function
  thinned_species_data <- thin_points(
    data = species_data,
    long_col = "X",
    lat_col = "Y",
    method = "grid",
    raster_obj = rast(predictors[[1]]),
    trials = 50,
    all_trials = FALSE,
    seed = 123
  )[[1]]

  # Store the thinned data in the list
  thinned_data_list[[species]] <- thinned_species_data
}

# Combine all thinned data into a single data frame
all_thinned_data <- do.call(rbind, thinned_data_list)
rownames(all_thinned_data) <- NULL

# Filter out species with fewer than 15 occurrences
ebird_final_pc <- all_thinned_data %>%
  group_by(Binomil) %>%
  filter(n() >= 15) %>%
  ungroup()

write_csv(ebird_final_pc,".\\ebird_2024_thinned_final.csv")
}

#sp list
xspecies_ebird <- unique(ebird_final_pc$Binomil)


#######

##### SDM run ####

### cross-validation loop ###

for(i in 1:length(xspecies_ebird)){
  sp = xspecies_ebird[i]
  print(paste("Starting sp", i, sp, "..."))
  dir.create(paste0("./y2020/", sp), showWarnings = FALSE)
  
  # Get occurrence locations:
  # Filter for the species
  sp_ID <- filter(ebird_final_pc, Binomil == sp) 
  
  # Cross-validation repetitions
  boyce_list <- list() 
  for (c in 1:k_cross){
    sp_boyce <- c()
    # Create Spatial Folds
    # Create pa file with coordinates of presences and random background points
    occ <- data.frame(x = sp_ID$X, y = sp_ID$Y)
    temp <- length(data.frame(rasterToPoints(predictors[[1]], na.rm = TRUE))$x) # Number of pixels in the study area
    acc <- as.data.frame(randomPoints(predictors[[1]], 10000, occ)) # Select 10000 random background pixels
    
    occ$presence <- rep(1, nrow(occ))
    acc$presence <- rep(0, nrow(acc))
    
    pa <- rbind(occ, acc)
    
    # Add environmental information to pa file
    pa <- cbind(pa, data.frame(raster::extract(predictors, pa[, c("x", "y")])))
    
    # Create the spatial folds
    pa_sv <- vect(pa, geom = c("x", "y"), crs = raster::crs(predictors)) # Convert to SpatVector
    block_size <- 5000
    
    blocks <- cv_spatial(x = as(pa_sv, "Spatial"), column = "presence", size = block_size, k = k_fold, selection = "random", plot = FALSE)
    
    # Add spatial block ID:
    pa$foldID <- blocks$folds_ids
    
    # Split data into training and test subsets
    for(k in 1:k_fold){   # 5-fold
      print(paste("Model fold", k))
      
      dir.create(paste0("./y2020/", sp, "/rep_", c, "/iter_", k), recursive = TRUE, showWarnings = FALSE)
      
      pa_train <- pa[pa$foldID != k, ]
      pa_test <- pa[pa$foldID == k, ]
      
      # Go to the next fold if there are no presences in the test set
      if (sum(pa_test$presence == 1) == 0){
        sp_boyce <- append(sp_boyce, NA)
        next
      } else {
        
        ### Perform the Species Distribution Modeling
        ## With MaxEnt
        
        # Extract occurrence and background points for training
        occ_train <- pa_train[pa_train$presence == 1, c("x", "y")]
        acc_train <- pa_train[pa_train$presence == 0, c("x", "y")]
        
        # Convert to SpatialPoints
        occ_train_points <- SpatialPoints(occ_train, proj4string = raster::crs(predictors))
        acc_train_points <- SpatialPoints(acc_train, proj4string = raster::crs(predictors))
        
        # Run MaxEnt with bias file
        mxnt_sp <- maxent(
          x = predictors,                # RasterStack of predictors
          p = occ_train_points,          # SpatialPoints of occurrence training data
          a = acc_train_points,          # SpatialPoints of background training data
          biasfile = bias_TAD,           # Specify bias file
          path = paste0(getwd(), "/y2020/", sp, "/rep_", c, "/iter_", k)
        )
        
        # Predict using the fitted model
        p_clog <- predict(
          object = mxnt_sp,
          x = predictors,
          args = c("outputformat=cloglog")
        )
        
        ## Evaluation with presence/absence data
        # Extract predictions at test points
        pa_test$hs <- raster::extract(p_clog, pa_test[, c("x", "y")])
        
        mxnt_eval_p <- pa_test[pa_test$presence == 1, ]$hs
        mxnt_eval_a <- pa_test[pa_test$presence == 0, ]$hs
        
        # Remove NA values if any
        mxnt_eval_p <- mxnt_eval_p[!is.na(mxnt_eval_p)]
        mxnt_eval_a <- mxnt_eval_a[!is.na(mxnt_eval_a)]
        
        # Evaluate the model
        e2 <- evaluate(p = mxnt_eval_p, a = mxnt_eval_a)
        TSS <- max(e2@TPR - e2@FPR)
        
        ## Boyce index
        boyce_k <- ecospat.boyce(fit=na.omit(getValues(p_clog)), obs=na.omit(mxnt_eval_p),
                              rm.duplicate = T, method='kendall')
        boyce_df <- data.frame( HS = boyce_k$HS, F.ratio = boyce_k$F.ratio)
        res <- rbind(data.frame(mxnt_sp@results), TSS = TSS, boyce = boyce_k$cor)
        sp_boyce <- append(sp_boyce,boyce_k$cor)

        ## Export
        #results table
        write.csv(res, paste0("./y2020/",sp,"/rep_",c,"/iter_",k,"/maxent_results.csv"))
        #Boyce values
        write.csv(boyce_df, paste0("./y2020/",sp,"/rep_",c,"/iter_",k,"/boyce.csv"))
        #evaluation object
        saveRDS(e2, paste0("./y2020/",sp,"/rep_",c,"/iter_",k,"/evaluate.rda"))
        #model object
        saveRDS(mxnt_sp, paste0("./y2020/",sp,"/rep_",c,"/iter_",k,"/maxent.rda"))
        #prediction raster
        p_clog <- raster::crop(p_clog, extent(r_mask)) 
        p_clog <- raster::mask(p_clog, r_mask) #crop and mask to match research area boundaries
        raster::writeRaster(p_clog, paste0("./y2020/",sp,"/rep_",c,"/iter_",k,"/prediction.grd"),overwrite=TRUE)
          }
  }
  boyce_c <-  data.frame(sp = sp,
                         rep=c,
                           "X1" = sp_boyce[1],
                           "X2" = sp_boyce[2],
                           "X3" = sp_boyce[3],
                           "X4" = sp_boyce[4],
                           "X5" = sp_boyce[5],
                           mean=mean(na.omit(sp_boyce)),
                           sd=sd(na.omit(sp_boyce)))
  write.csv(boyce_c,paste0("./y2020/",sp,"/rep_",c,"/boyce_sum.csv") )
  
  boyce_list[[c]]  <- boyce_c
  }
  
  boyce_sum <-  do.call("rbind",boyce_list)
  write.csv(boyce_sum, paste0("./y2020/",sp,"/boyce_sum.csv"))
}

# dataframe for list of species
{columns <-  c("Sp",c(paste0("M",as.character(1:25))),"exclude") 
sp_df <-  data.frame(matrix(nrow = 0, ncol = length(columns))) 
colnames(sp_df) <-  columns}

#### evaluation loop ####
for(i in 1:length(xspecies_ebird)){
sp=xspecies_ebird[i]
   #get summary of cross validation Boyce index
  boyce_sum <-  read.csv(paste0("./y2020/",sp,"/boyce_sum.csv"))
   
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
      labs(title = paste("Cumulative Means and SD vs Number of Runs for", sp),
           x = "Number of Runs",
           y = "Value") +
      ylim(-1, 1) +
      theme_minimal() +
      theme(legend.position = "bottom")
    
    # Save the first plot
    ggsave(filename = file.path(paste0("./y2020/",sp,"/Cumulative_Means_SD.jpg")), plot = p1, device = "jpeg")
    
    # Create a data frame for plotting
    cumulative_data <- data.frame(
      Runs = 1:length(cumulative_target),
      CumulativeTarget = cumulative_target
    )
    
    # Plot cumulative target with y-axis limits set between 0 and 1
    p2 <- ggplot(cumulative_data, aes(x = Runs, y = CumulativeTarget)) +
      geom_line(color = "blue") +
      geom_point(color = "red") +
      labs(title = paste("Cumulative Target vs Number of Runs for", sp),
           x = "Number of Runs",
           y = "Cumulative Target") +
      ylim(0, 1) +
      theme_minimal()
    
    # Save the second plot
    ggsave(filename = file.path(paste0("./y2020/",sp,"/Cumulative_Target.jpg")), plot = p2, device = "jpeg")
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
  
    sp_df <- rbind(sp_df,c(sp, x_values, exclude))
    names(sp_df) <- columns
}

write.csv(sp_df, paste0("./results/boyce_summary_",
                                                format(Sys.time(), "%Y-%m-%d"),".csv"  ))

#### Final model loop ####
sp_df2 <- sp_df

sp_list <- filter(sp_df2, exclude == 0 )
length(sp_list$Sp)

for(i in 1:length(sp_list$Sp)){
    print(paste("Sp",i, sp_list$Sp[i],"final model..."))
  
  #1. fit the model
   # get sp occurrence locations:
  sp_ID <- subset(ebird_final_pc, Binomil == sp_list$Sp[i])
  occ <- data.frame(x = sp_ID$X, y = sp_ID$Y)
  
  # 2. Generate random background points
  acc <- as.data.frame(randomPoints(predictors[[1]], 10000, occ)) # 10000 random background points
  
  # 3. Convert occurrences to SpatialPoints
  occ_points <- SpatialPoints(occ, proj4string = raster::crs(predictors))
  acc_points <- SpatialPoints(acc, proj4string = raster::crs(predictors))
  
  # 4. Fit the MaxEnt model
  
  mxnt_sp <- maxent(
    x = predictors,               # RasterStack of predictors
    p = occ_points,               # SpatialPoints of occurrences
    a = acc_points,               # SpatialPoints of background points
    biasfile=bias_TAD, # Specify bias file
    path = paste0(getwd(), "/y2020/", sp_list$Sp[i]) # Output directory
  )
  
     p_clog<- predict(object = mxnt_sp, x= predictors, args=c("outputformat=cloglog"))
  
    #3. evaluation metrics
    # Extract predicted values at presence (occurrence) points
    mxnt_eval_p <- raster::extract(p_clog, occ_points)
    # Extract predicted values at evaluation background points
    mxnt_eval_a <- raster::extract(p_clog, acc_points)
    
    e2 <- evaluate(mxnt_sp,p=mxnt_eval_p, a=mxnt_eval_a)
    TSS <- max(e2@TPR-e2@FPR)
    
    boyce_k <- ecospat.boyce(fit=na.omit(getValues(p_clog)), obs=na.omit(mxnt_eval_p),
                             rm.duplicate = T, method='kendall')
    boyce_df <- data.frame( HS = boyce_k$HS, F.ratio = boyce_k$F.ratio)
    res <- rbind(data.frame(mxnt_sp@results), TSS = TSS, boyce = boyce_k$cor)
    sp_df2[sp_df2$Sp == sp_list$Sp[i],"final"] <- boyce_k$cor

    #4. Export
     #results table
    write.csv(res, paste0("./y2020/",sp_list$Sp[i],"/maxent_results.csv"))
     #Boyce values
    write.csv(boyce_df, paste0("./y2020/",sp_list$Sp[i],"/boyce.csv"))
     #evaluation object
    saveRDS(e2, paste0("./y2020/",sp_list$Sp[i],"/evaluate.rda"))
     #model object
    saveRDS(mxnt_sp, paste0("./y2020/",sp_list$Sp[i],"/maxent.rda"))
     #prediction raster
    p_clog <- raster::crop(p_clog, extent(r_mask)) 
    p_clog <- raster::mask(p_clog, r_mask) #crop and mask to match research area boundaries
    raster::writeRaster(p_clog, paste0("./y2020/",sp_list$Sp[i],"/prediction.grd"),overwrite=TRUE)
    
   #5. response curve data
    dir.create(paste0("./y2020/",sp_list$Sp[i],"/responses"), showWarnings = FALSE)
    for(o in 1:length(names(predictors))){
      r <- response(mxnt_sp, var = names(predictors)[o],  expand = 0)
      write.csv(r, paste0("./y2020/",sp_list$Sp[i],"/responses","/response_",names(predictors)[o],".csv"))
    }

    print(paste("species",i, sp_list$Sp[i] ,"complete, Boyce index=", boyce_k$cor))

  }

write.csv(sp_df2, paste0("./results/final_boyce_summary_",
                        format(Sys.time(), "%Y-%m-%d"),".csv"  ))



#################

for(i in 1:length(sp_list$Sp)){
  print(paste("Sp",i, sp_list$Sp[i],"2020 prediction..."))

  mxnt_sp <-  readRDS(paste0("./y2020/",sp_list$Sp[i],"/maxent.rda"))

  #predict
  p_clog<- predict(object = mxnt_sp, x= predictors, args=c("outputformat=cloglog"))
  raster::crs(p_clog) <- raster::crs(predictors)
  p_clog <- raster::crop(p_clog, extent(r_mask))
  p_clog <- raster::mask(p_clog, r_mask) #crop and mask to match research area boundaries
  raster::writeRaster(p_clog, paste0("./y2020/",sp_list$Sp[i],"/prediction.grd"),overwrite=TRUE)
  sp_ID <- subset(ebird_final_pc, Binomil == sp_list$Sp[i])
  occ <- data.frame(x = sp_ID$X, y = sp_ID$Y)
  
  # 2. Generate random background points
  acc <- as.data.frame(randomPoints(predictors[[1]], 10000, occ)) # 10000 random background points
  
  # 3. Convert occurrences to SpatialPoints
  occ_points <- SpatialPoints(occ, proj4string = raster::crs(predictors))
  acc_points <- SpatialPoints(acc, proj4string = raster::crs(predictors))
  
  #3. evaluation metrics
  # Extract predicted values at presence (occurrence) points
  mxnt_eval_p <- raster::extract(p_clog, occ_points)
  # Extract predicted values at evaluation background points
  mxnt_eval_a <- raster::extract(p_clog, acc_points)
  
  e2 <- evaluate(mxnt_sp,p=mxnt_eval_p, a=mxnt_eval_a)
  TSS <- max(e2@TPR-e2@FPR)
  
  boyce_k <- ecospat.boyce(fit=na.omit(getValues(p_clog)), obs=na.omit(mxnt_eval_p),
                           rm.duplicate = T, method='kendall')
  boyce_df <- data.frame( HS = boyce_k$HS, F.ratio = boyce_k$F.ratio)
  res <- rbind(data.frame(mxnt_sp@results), TSS = TSS, boyce = boyce_k$cor)
  sp_df2[sp_df2$Sp == sp_list$Sp[i],"final"] <- boyce_k$cor
  
  #4. Export
  #results table
  write.csv(res, paste0("./y2020/",sp_list$Sp[i],"/maxent_results.csv"))
  #Boyce values
  write.csv(boyce_df, paste0("./y2020/",sp_list$Sp[i],"/boyce.csv"))
  #evaluation object
  saveRDS(e2, paste0("./y2020/",sp_list$Sp[i],"/evaluate.rda"))
}

############# Scenario SDMs ##################



##### produce predictions ######

sp_list <- filter(sp_df2, exclude == 0 )

sce_names <- c("LSpH", "LShH","LSpF/30", "LSpL/30")
years <- c("y2025","y2030","y2035","y2040","y2045","y2050")

k_scen <- 10
#focal_list <- c( "bldg_cover", "UGS","tree_cover","low_veg" )
focal_list <- c( "bldg_cover")

# average bias layer for predictions
#bias_TAD_avg <-  bias_TAD * 0 + mean(raster::values(bias_TAD), na.rm=TRUE)


# loop through scenarios -> iterations -> years -> species to produce new predictions

for (i in 1:length(sce_names)){ #scenario types
  dir.create(paste0("./",sce_names[i]), recursive = TRUE, showWarnings = F)
  
  for (j in 1:length(years)){ #years
    dir.create(paste0("./",sce_names[i],"/",years[j]), showWarnings = FALSE)
    
    for (k in 1:k_scen){ #scenario & SDM iterations
      print(paste("Starting", sce_names[i],years[j],"iterarion",k ))
      
      path_k <- paste0("./",sce_names[i],"/",years[j],"/iter_",k)
      dir.create(path_k, showWarnings = FALSE)
      
      ### import and process predictors for scenario k-th iteration  
      #create predictor list
      pred_list <-  list.files(paste0("./R/chapter 2/scenarios/",sce_names[i],"/iter_",k,"/",years[j]), pattern='tif$', full.names=FALSE )
      pred_list <-  tools::file_path_sans_ext(pred_list)
      pred_list <- pred_list[pred_list!="bldg_vol"]
      pred_list <- pred_list[pred_list!="UGS"] #exclude UGS from predictors
      pred_list <- pred_list[pred_list!="richness_500"]
      #read and add buffer to allow 500m calculation
      predictors_s <- stack()
      for (r in 1:length(pred_list)){
        r_name <- pred_list[r]
        rast <- raster(paste0("./R/chapter 2/scenarios/",sce_names[i],"/iter_",k,"/",years[j],"/",r_name,".tif"))
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
      # add bias layer
      #predictors_s <- stack(predictors_s,bias_TAD_avg)
      
      #crop and mask to match research area boundaries
      predictors_s <- raster::crop(predictors_s, extent(r_mask)) 
      predictors_s <- raster::mask(predictors_s, r_mask)
      
      if (any(!names(predictors) %in% names(predictors_s))) { 
        print ("predictors incomplete! missing:")
        print (setdiff(names(predictors), names(predictors_s)))
      }
     
       for (s in 135:length(sp_list$Sp)){ #species
        #for (m in 1:k_fold){ #maxent iterations
          
         #import model 
         mxnt_sp <-  readRDS(paste0("./y2020/",sp_list$Sp[s],"/maxent.rda"))
         
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


