
#### libraries ####

library(tidyverse)
library(sf)
library(dplyr)
library(raster)
#####

### Distance functions etc.
{

#focal variables
focal_500 <- function(rast, pad_rast){
  rast_e <- raster::extend(rast, pad_rast)
  rast_e <- raster::overlay(rast_e, pad_rast, fun = function(r, pad) {
    # Use value from raster1 if not NA, otherwise use green_pad
    return(ifelse(!is.na(r), r, pad))
  })
  rast_500 <- raster::focal(rast_e, w=matrix(1/25,nrow=5,ncol=5))
  rast_500 <- raster::crop(rast_500, rast)
  rast_500 <- raster::mask(rast_500,rast)
  names(rast) <- paste0(rast,"_500")
  return(rast_500)
}



#### get_pixeldist ####

# calculate distance to changed pixels                        
get_changedist <- function (env_df){  
  rast <- rasterFromXYZ(env_df[, c("X","Y", "change")]) #create res bldg vol raster
  rast[rast[]<= 0] <- NA #ignore cells with zeros

  if(all(is.na(rast[]))){ #if all cells are NA - return max value
    env_df$change_dist <- 100
    env_df$change_dist_w <- 0
     } else {
       
  rast <- raster::distance(rast) #calculate distance to non NA pixels
  df <-raster::extract(rast,env_df[, c("X","Y")],df=T) #extract back to df
  env_df$change_dist <- df[,"layer"]
  env_df$change_dist_w <-  0.99^env_df$change_dist #transform dist values - exponential decay
     }
  
  return(env_df)
}

# distance to built pixels
get_builtdist <- function (env_df){  
  rast <- rasterFromXYZ(env_df[, c("X","Y", "bldg_vol_res")]) #create res bldg vol raster
  bldg_1 <- max(env_df$bldg_vol_res[env_df$bldg_grp<2])
  rast[rast[]<= bldg_1] <- NA #ignore cells with zero or very low values
  rast <- raster::distance(rast) #calculate distance to non NA pixels
  df <-raster::extract(rast,env_df[, c("X","Y")],df=T) #extract back to df
  env_df$bldg_dist <- df[,"layer"] 
  return(env_df)
}

# distance to district's periphery
get_opentdist <- function (env_df){  
  
  stat_bldg <- env_df %>% 
    dplyr::select(stat,res_stat,bldg_cover) %>%
    filter(stat > 0) %>%
    group_by(stat,res_stat) %>%
    summarise(area = n(), 
              bldg_stat = sum(bldg_cover)/area)
  
  env_df_temp <- env_df %>%
    left_join(stat_bldg[,c("stat","bldg_stat")], by = "stat") 
  env_df_temp$bldg_stat[is.na(env_df_temp$bldg_stat)] <- 0
  
  rast <- rasterFromXYZ(env_df_temp[, c("X","Y", "bldg_stat")]) #create stat area raster
  rast[rast[]> 0.01] <- NA #convert built stat areas to NA
  rast <- raster::distance(rast) #calculate distance to non NA pixels = periphery pixels
  df <-raster::extract(rast,env_df[, c("X","Y")],df=T) #extract back to original df
  env_df$open_dist <- df[,"layer"] 
  
  return(env_df)
}

}


#### SCENARIOS ####

#### LShH - high land sharing ####

LShH <- function(env_table, vol_trgt){
 
  start_time <- Sys.time()
  vol_added <- 0
  env_new <- env_table
  i <- 0 #counter for all iterations
  Si <- 0 #counter for successful iterations
  dist_500 <- 0 #check if used built_dist = 200
 
  samp_size <- numeric() #vectors for documenting simulation data 
  rep_size <- numeric()

  if ((!"change" %in% names(env_new))){  #only if it's the first scenario run
    env_new$change <- 0
  }
  
  env_new$CFID <- env_new$FID
  
  while ( vol_added < vol_trgt ){ #loop until reaching bldg volume target

    if (i %% 20 == 0 | dist_500 != 0){ #recalculate distances every 20 iterations
    env_new <- get_builtdist(env_new) #get distances
    env_new <- get_changedist(env_new)
    
    #calculate w
    env_new <- env_new %>%  mutate( w = (inTAD+change_dist_w)/2) 
    }
    
    dist_500 <- 0
    i <- i+1
    
    #break after too many unsuccessful iterations
    if (i-Si>1000 ) {
      break
    } 
    
    samp_df <- filter(env_new, (protected == 0) & (bldg_cover<0.3) & (UGS<0.3) & 
                     (paved<0.2) & (bldg_grp<2) & (bldg_dist == 100))
    samp <- slice_sample(samp_df, n=1, weight_by = w)   #sample pixel to be replaced
    
    if (dim(samp)[1] == 0) { #allow larger distance if cannot find pixel
      samp_df <- filter(env_new, (protected == 0) & (bldg_cover<0.3) & (UGS<0.3) & 
                       (paved<0.2) & (bldg_grp<2) & (bldg_dist <= 500)) 
      samp <- slice_sample(samp_df, n=1, weight_by = w)  
      
      dist_500 <- 1
      
      if (dim(samp)[1] == 0) {
        samp_df <- filter(env_new, (protected == 0) & (bldg_cover<0.3) & (UGS<0.3) & 
                         (paved<0.2) & (bldg_grp<2) )  
        samp <- slice_sample(samp_df, n=1, weight_by = w)  
        dist_500 <- 2
      }
      if (dim(samp)[1] == 0) {
        break()
      }
    }
    
    samp_size <- append(samp_size, dim(samp_df)[1])

    #subset replacing pixels
     rep_df <- subset(env_new, 
                     (bldg_h_res >0) & (bldg_h_res <13) & 
                     (paved >= (samp$paved-0.05))& (bldg_cover >= (samp$bldg_cover-0.05)))  
 
    rep_size <- append(rep_size, dim(rep_df)[1])
 
     #next iteration if cannot find rep
     if (dim(rep_df)[1] == 0) { 
    
      next
    }
    
    rep <- slice_sample(rep_df, n=1) #sample one of possible reps
   
    rep <- rep %>% #keep some of the original cell values
      mutate(FID = samp$FID, X = samp$X, Y = samp$Y, max_h = samp$max_h, change = 1) 
    env_new [env_new$FID == samp$FID, ] <-  rep #replace in dataframe
    
    vol_added <-  vol_added + (rep$bldg_vol_res -  samp$bldg_vol_res) #recalculate bldg volume
    Si <- Si+1
    
}
 
  sim.data <- data.frame(i=length(samp_size),
                         samp_mean=mean(samp_size),
                         samp_sd = sd(samp_size),
                         rep_mean=mean(rep_size),
                         rep_sd = sd(rep_size))
  
  res.list <- list(env_new=env_new, sim.data=sim.data)
    
  time_diff <- Sys.time()- start_time
  print(paste("LShH run completed, ", round((vol_added/vol_trgt*100),0) ,
              "% of target volume added. run time:"))
  print(time_diff)
  
  return(res.list)
}
######

#### LSpH - high land sparing ####
LSpH <- function(env_table, vol_trgt, non_res_trgt){
  start_time <- Sys.time()
  env_new <- env_table
  
  vol_added <- 0
  non_res_added <- 0
  i <- 0 #counter for all iterations
  Si <- 0 #counter for successful iterations

  samp_size <- numeric() #vectors for documenting simulation data 
  rep_size <- numeric()
    
  if ((!"change" %in% names(env_new))){  #only if it's the first scenario run
    env_new$change <- 0
  }
  if ((!"spare" %in% names(env_new))){  #create an empty spare column if it doesn't exist
    env_new$spare <- 0
  }
  
 env_new$CFID <- env_new$FID

 if( sum(env_new$spare)>0){  # remove max_h restrictions in full sparing scenario
   max_h <- dplyr::select(env_new, FID, max_h)
   env_new$max_h <- max(env_new$max_h)
 }
 
 while ( vol_added < vol_trgt ) #loop until reaching bldg volume target
  {
    if (i-Si>1000 ) {
      break
    }
    i <- i+1
    
    samp_df <- filter(env_new, bldg_grp>0 & bldg_h_res<max_h & res_parcel>=0.1 & 
                     inTAD == 1 & spare == 0 )
      
    samp_size <- append(samp_size, dim(samp_df)[1])
    
    samp <- slice_sample(samp_df, n=1, weight_by = (max_h-bldg_h_res))  #sample pixel to be replaced
    
    rep_df <- filter(env_new,  bldg_h_res<samp$max_h &  
                    bldg_vol_res>samp$bldg_vol_res & 
                    res_parcel>(samp$res_parcel-0.05) & 
                    res_parcel<(samp$res_parcel+0.05 ))  
    
   if (dim(rep_df)[1] == 0) { #next iteration if cannot find rep
      next
    }
    
    rep_size <- append(rep_size, dim(rep_df)[1]) 
    
    rep <- slice_sample(rep_df, n=1)  %>% #sample replacing pixel
      #keep some of the original cell values
           mutate(FID = samp$FID, X = samp$X, Y = samp$Y, max_h = samp$max_h, 
             stat = samp$stat, res_stat = samp$res_stat, spare = 0, inTAD = 1,
             change = 1 )
    env_new [env_new$FID == samp$FID, ] <-  rep #replace in dataframe
    
    vol_added <- vol_added + (rep$bldg_vol_res -  samp$bldg_vol_res) #recalculate bldg volume
    non_res_added <- non_res_added + (rep$bldg_vol - rep$bldg_vol_res) - 
      (samp$bldg_vol - samp$bldg_vol_res)# non-residential building volume - only relevant in full sparing
    
    Si <- Si+1
    
 }
 
 # loop to achieve non_res_trgt - only in full/local sparing scenarios
 i <- 0
 Si <- 0
 if( sum(env_new$spare)>0){
 
   while (non_res_added < non_res_trgt){
   if (i-Si>10000 ) {
     break
   }
   i <- i+1
   
   samp <- filter(env_new, bldg_grp>0 & 
                    inTAD == 1 & spare == 0 ) %>% 
     slice_sample( n=1)  #sample pixel to be replaced
   rep <- filter(env_new,    
                 bldg_vol>samp$bldg_vol & 
                   bldg_vol_res>(samp$bldg_vol_res-0.05*samp$bldg_vol_res) & 
                   bldg_vol_res<(samp$bldg_vol_res+0.05*samp$bldg_vol_res ) &
                   res_parcel>(samp$res_parcel-0.05) & 
                   res_parcel<(samp$res_parcel+0.05 )&
                   paved>(samp$paved-0.05)) %>% 
     slice_sample( n=1) #sample replacing pixel
   
   if (dim(rep)[1] == 0) { #next iteration if cannot find rep
     next
   }
   
   rep <- rep %>% #keep some of the original cell values
     mutate(FID = samp$FID, X = samp$X, Y = samp$Y, max_h = samp$max_h, stat = samp$stat, 
            res_stat = samp$res_stat, spare = 0, inTAD = 1,
            change = 1 ) 
   env_new [env_new$FID == samp$FID, ] <-  rep #replace in dataframe
   
   vol_added <- vol_added + (rep$bldg_vol_res -  samp$bldg_vol_res) #recalculate bldg volume
   non_res_added <- non_res_added + (rep$bldg_vol - rep$bldg_vol_res) - 
     (samp$bldg_vol - samp$bldg_vol_res)# non-residential building volume - only relevant in full sparing
   
   Si <- Si+1
   }
   }
 
 #get max_h info back
 env_new <- env_new %>%
 left_join(dplyr::select(max_h, FID, new_max_h = max_h), by = "FID") %>% 
   mutate(max_h = coalesce(new_max_h, max_h)) %>%
   dplyr::select(-new_max_h)
 
 sim.data <- data.frame(res_i=length(samp_size),
                        samp_mean=mean(samp_size),
                        samp_sd = sd(samp_size),
                        rep_mean=mean(rep_size),
                        rep_sd = sd(rep_size),
                        non_res_i <- i)
 
 res.list <- list(env_new=env_new, sim.data=sim.data)
 
 if( sum(env_new$spare)==0){ #print only in simple sparing scenarios
 time_diff <- Sys.time()- start_time
  print(paste("LSpH scenario run complete with", round((vol_added/vol_trgt*100),0) ,"% of target volume added. run time:"))
  print(time_diff)
 }
 else{print(paste("LSp densification:",round((vol_added/vol_trgt*100),0) ,"% of residential target  volume added and",
                  round((non_res_added/non_res_trgt*100),0) ,"% of non-res target volume"))}
 
  return(res.list)
}



##### LSpF - Full Sparing (Regional) #######

LSpF <- function(env_table, vol_trgt, Spare_trgt){ #Spare_trgt is the % area that should be spared
  start_time <- Sys.time()
  env_new <- env_table
  env_new <- get_opentdist(env_new)
  
   # create table for statistical areas
{
  options(dplyr.summarise.inform = FALSE)
  
    stat_table <- env_new %>% 
    dplyr::select(stat,res_stat,open_dist,bldg_cover) %>%
    filter(res_stat > 0) %>%
    group_by(res_stat,stat ) %>%
    summarise(area = n(), 
              min_dist = min(open_dist),
              stat_bldg = sum(bldg_cover)/area) %>% 
     ungroup()%>%
    mutate( spare = 0)
  
  
 #check for spared areas from previous scenario runs
  if ("spare" %in% names(env_new))
  {
    stat_table <- stat_table %>%
      mutate(spare = env_new$spare[match(stat, env_new$stat)])
  }
   } 
 
   # select stat areas to spare
{
  stat_temp <- stat_table %>%
    filter(area>=5)  #exclude areas with less than 5 pixels for target calculations - to match requirements of LSPL
  
  Spare_trgt_area <- Spare_trgt*sum(stat_temp$area)
  spared_area <- 0
  
  while (spared_area < Spare_trgt_area){
    samp <- filter(stat_table, min_dist==0, res_stat>0, spare == 0  ) %>% 
      slice_sample( n=1)  #sample stat area to be spared
    if (dim(samp)[1] == 0) {   
      samp <- filter(stat_table, min_dist<150, res_stat>0, spare == 0  ) %>% 
        slice_sample( n=1 ) 
    }
    
    stat_table[stat_table$stat == samp$stat, "spare" ] <- 1
    area_value <- stat_table[stat_table$stat == samp$stat, "area"][[1]]
    spared_area <- spared_area+area_value
}

 # calculate added buildings from spared areas
  spared_bldg_res <- sum(env_new %>%
                       filter(stat %in% stat_table[stat_table$spare == 1, ]$stat) %>%
                       pull(bldg_vol_res), na.rm = TRUE)
  spared_bldg_nonres <- sum(env_new %>%
                           filter(stat %in% stat_table[stat_table$spare == 1, ]$stat) %>%
                           pull(bldg_vol), na.rm = TRUE) - spared_bldg_res
  
  #modify env table
  env_new <- env_new %>%
    mutate(spare = case_when(
      stat %in% stat_table[stat_table$spare == 1, ]$stat ~ 1,  # If stat matches, set spare to 1
      TRUE ~ 0 ))
} 
  
   # densify with LSpH
 { 
   lsp_res <- LSpH (env_new , vol_trgt+spared_bldg_res, spared_bldg_nonres)
   env_new1 <- lsp_res$env_new
   sim.data <- lsp_res$sim.data
 }
 
  # clear spared areas
{
  spare_list <- env_new1 %>% # list of all pixels to clear
    filter(spare == 1) %>%
    dplyr::select(FID, stat)
  
  rep_size <- numeric()
  
  for (i in spare_list$FID)
   {
    samp <- env_new1[env_new1$FID == i, ]
    
    if (samp$bldg_cover < 0.01 & samp$paved<0.2){ #skip replacement if pixel is open
      next
    }
    
    rep_df <- filter(env_new1,bldg_cover == 0, agri_field ==0, agri_orchard ==0, paved<0.2)
    
    rep_size <- append(rep_size, dim(rep_df)[1])
    
    rep <- slice_sample(rep_df, n=1) %>% #keep some of the original cell values
    mutate(FID = i, X = samp$X, Y = samp$Y, stat = samp$stat, res_stat = samp$res_stat, 
           change = 1, spare = 1) 
    env_new1 [env_new1$FID == i, ] <-  rep #replace in data frame
  }
}  
  
  
  
 # summarize scenario
{
  vol_added <-  sum(env_new1$bldg_vol_res)-sum(env_table$bldg_vol_res)
  
  sim.data$spared_res <- sum(stat_table$spare)
  sim.data$spared_i <- sum(env_new1$spare)
  sim.data$spare_rep_mean <- mean(rep_size)
  sim.data$spare_rep_sd <- sd(rep_size)
  
  res.list <- list(env_new=env_new1, sim.data=sim.data)
    
  time_diff <- Sys.time()- start_time
  print(paste("full sparing complete with", round((vol_added/vol_trgt)*100,0) ,"% of building volume target,",
              round((spared_area/Spare_trgt_area)*100,0),  "% of sparing target achieved. run time:"))
  print(time_diff)
}
    
  return(res.list)
}


##### LSpL - Local Sparing (Neighbourhood) #######

LSpL <- function(env_table, vol_trgt, spare_trgt){ #Spare_trgt is the % area that should be spared
  start_time <- Sys.time()
  env_new <- env_table

    if ((!"spare" %in% names(env_new))){  #only if it's the first scenario run
    env_new$spare <- 0 
    }

    # table for statistical areas
  {
    options(dplyr.summarise.inform = FALSE)

    stat_table <- env_new %>% 
     dplyr::select(res_stat,bldg_vol_res) %>%
     filter(res_stat > 0) %>%
     group_by(res_stat ) %>%
     summarise(area = n(), 
              bldg_vol_res = sum(bldg_vol_res))%>%
      filter(area>=5)  #exclude areas with less than 5 pixels - too small for some raster operations
    
    
    Spare_trgt_area <- spare_trgt*sum(stat_table$area)
    
    }
  
  # loop over statistical areas
  for (s in stat_table$res_stat){
   env_stat <- filter(env_new, res_stat == s)
  
   # define local targets
   {
   area <-  stat_table %>% filter(res_stat == s) %>% pull(area)
   stat_spare_trgt <- round(spare_trgt * area,0)
   spared_before <- sum(env_stat$spare)
   
   }
   
   # Only in first scenario run - identify largest UGS
   if (spared_before == 0)
   {
   UGS <- rasterFromXYZ(env_stat[,c("X","Y","UGS")]) #UGS raster
   UGS[UGS>=0.6] <- 1 #threshold for UGS pixel
   UGS[UGS<0.6] <- NA
   
   if( sum(raster::values(!is.na(UGS)))>0 )
     {
     UGS_clumps <- clump(UGS, directions=8) #identify clumps of UGS pixels
     UGS_polys <- rasterToPolygons(UGS_clumps, dissolve=TRUE) #convert to to polygons
     UGS_polys$area <- sapply(slot(UGS_polys, "polygons"), function(x) slot(x, "area"))#calculate area
     largest_UGS <- which.max(UGS_polys$area)#find largest UGS
     env_stat$clump <- raster::extract(UGS_clumps, env_stat[,c("X","Y")])#get clump ID to table
   
     env_stat <- mutate(env_stat,spare = ifelse(is.na(clump), spare, 
                                       ifelse(clump == largest_UGS, 1, spare)))
    }   else  { #start sparing in a pixel with lowest bldg_cover + paved
       spare_start <- env_stat %>%
       filter(paved<0.5 & bldg_cover<0.5)%>%
       filter((bldg_cover + paved) == min(bldg_cover + paved))%>%
       slice_sample(n=1)
     env_stat <-  mutate(env_stat,spare = ifelse(FID == spare_start$FID, 1, spare))
   }
    }
   
   # Sparing loop
   spare_added <- sum(env_stat$spare)- spared_before
   samp_size <- numeric()
   
   while(spare_added < stat_spare_trgt){
     
     #identify sparing neighbors
     spare_raster <- rasterFromXYZ(env_stat[, c("X", "Y", "spare")])
     raster::values(spare_raster)[raster::values(spare_raster) != 1] <- NA 
     adjacencyMatrix <- matrix(1, nrow = 3, ncol = 3) # Define the matrix to specify 8-directional adjacency
     count_spare_neighbors <- raster::focal(spare_raster, w = adjacencyMatrix, pad = TRUE, padValue = NA, fun = function(x) sum(x == 1, na.rm = TRUE))
     env_stat$spare_neighbor_count <- raster::extract(count_spare_neighbors, env_stat[, c("X", "Y")])
   
     #sample pixel to spare
     samp_pool <- env_stat %>%
       filter(spare == 0, spare_neighbor_count > 0, UGS > 0.1) 
     
     if(dim(samp_pool)[1]==0){ # if no such pixels exist, omit the UGS condition
           samp_pool <- env_stat %>%
         filter(spare == 0, spare_neighbor_count > 0) 
     }
     
    max_n <- max(samp_pool$spare_neighbor_count)
     samp_df <- samp_pool %>%
       filter(spare_neighbor_count == max_n) 
     
     samp_size <- append(samp_size, dim(samp_df)[1])

     samp <- slice_sample(samp_df, n = 1)
     
       if(dim(samp)[1]==0)
     {
       print(paste("Samp = 0 in stat",s))
       break
     }
     env_stat$spare[env_stat$FID == samp$FID] <- 1
     spare_added <- spare_added + 1
   }
   
   # validate sparing loop
   if(spare_added != sum(env_stat$spare) - spared_before)
   {
     print(paste("Error in sparing loop s=",s))
     break
   }
   
  #update spared pixels in env_new
   env_new$spare[match(env_stat$FID, env_new$FID)] <- env_stat$spare
   
  }
 
  # sum bldg vol in spared pixels
  spared_bldg_res <- sum(env_new %>%
                           filter(spare == 1) %>%
                           pull(bldg_vol_res), na.rm = TRUE)
  spared_bldg_nonres <- sum(env_new %>%
                              filter(spare == 1) %>%
                              pull(bldg_vol), na.rm = TRUE) - spared_bldg_res
  
  # densify with LSpH
  { 
    lsp_res <- LSpH (env_new , vol_trgt+spared_bldg_res, spared_bldg_nonres)
    env_new1 <- lsp_res$env_new
    sim.data <- lsp_res$sim.data
  }

  # clear spared pixels
  {
    spare_list <- env_new1 %>% # list of all pixels to clear
      filter(spare == 1) %>%
      dplyr::select(FID, stat)
    
    rep_size <- numeric()
    
    for (i in spare_list$FID)
    {
      samp <- env_new1[env_new1$FID == i, ]
      
      if (samp$bldg_cover < 0.01 & samp$paved<0.2){ #skip replacement if pixel is open
        next
      }
      
      rep_df <- filter(env_new1,bldg_cover < 0.01, agri_field ==0, agri_orchard ==0, paved<0.2, UGS>0.1)
        
      rep_size <- append(rep_size, dim(rep_df)[1])
      
      rep <- slice_sample(rep_df, n=1) %>% #keep some of the original cell values
        mutate(FID = samp$FID, X = samp$X, Y = samp$Y, stat = samp$stat, res_stat = samp$res_stat,
               change = 1, spare = 1) 
      env_new1 [env_new1$FID == i, ] <-  rep #replace in data frame
    }
  }   
  

  # summarize scenario
  {
    vol_added <-  sum(env_new1$bldg_vol_res)-sum(env_table$bldg_vol_res)
    spared_area <- sum(env_new$spare)- sum(env_table$spare)
    
    sim.data$spared_i <- sum(env_new1$spare)
    sim.data$spare_samp_mean <- mean(rep_size)
    sim.data$spare_samp_sd <- sd(rep_size)
    sim.data$spare_rep_mean <- mean(rep_size)
    sim.data$spare_rep_sd <- sd(rep_size)
    
    res.list <- list(env_new=env_new1, sim.data=sim.data)
    
    time_diff <- Sys.time()- start_time
    print(paste("Local sparing complete with", round((vol_added/vol_trgt)*100,0) ,"% of building volume target,",
                round((spared_area/Spare_trgt_area)*100,0),  "% of sparing target achieved. run time:"))    
    print(time_diff)
  }
  return(res.list)
}
###### summary and export ####

scen_sum <- function (env_list){
  sumlist <-  lapply(env_list, function (x){summarise(x,bldg_vol_sum=sum(bldg_vol),
                                                            bldg_vol_mean=mean(bldg_vol[which(bldg_vol!=0)]), #don't include zeros in mean
                                                            bldg_h_mean=mean(bldg_h[which(bldg_h!=0)]), #don't include zeros in mean
                                                            bldg_cover_mean=mean(bldg_cover),
                                                            bldg_cover_sum=sum(bldg_cover),
                                                            low_veg_mean=mean(low_veg),
                                                            low_veg_sum=sum(low_veg),
                                                            tree_cover_mean=mean(tree_cover),
                                                            tree_cover_sum=sum(tree_cover),
                                                            agri_field_sum=sum(agri_field),
                                                            agri_orchard_sum=sum(agri_orchard))})
  
  sum_df <- as.data.frame(t(do.call(rbind, sumlist))) %>%  #store summary in scenario list
    mutate(total_change = (.[[length(sumlist)]] - y2020)/y2020)
  sum_df <- rownames_to_column(sum_df, var="variable")
  
  return(sum_df)
}

scen_export <- function (env_list_i,path){ #get item from a single scenario list 
  env_df <- env_list_i[[1]]
  
  #create predictor raster stack
  stacklist <- c("agri_field", "agri_orchard", "bldg_cover" , "bldg_h","bldg_vol", "water", "paved",
                 "UGS","veg_h_sd")
  predictors <- stack()
  for(r in 1:length(stacklist)){
    raster1 <- rasterFromXYZ(env_df[ ,c("X","Y",stacklist[r])])  #Convert first two columns as lon-lat and third as value   
    predictors <- raster::stack(predictors , raster1 )
  }
  
  greenlist <- c( "tree_cover","low_veg")
  for(r in 1:length(greenlist)){
    green_raster <- rasterFromXYZ(env_df[, c("X","Y",greenlist[r])])  #Convert first two columns as lon-lat and third as value   
    predictors <- raster::stack(predictors , green_raster )
      #extend to allow 500m calculation
      green_pad <-  raster(paste0(".\\100_grid_variables\\wBuffer\\",greenlist[r],".tif"))
      green_500 <- focal_500(green_raster,green_pad)
      names(green_500) <- paste0(greenlist[r],"_500")
      predictors <- raster::stack(predictors , green_500 )
        }
  
  dir.create(paste0(path,"/",names(env_list_i)))
  raster::writeRaster(predictors, filename=paste0(path,"/",names(env_list_i),"/",names(predictors)), bylayer=TRUE,overwrite=TRUE,format="GTiff")
  write_csv(env_df, paste0(path,"/",names(env_list_i),"/","env_table.csv")) 
  print(paste(names(env_list_i),"exported"))
  
}


