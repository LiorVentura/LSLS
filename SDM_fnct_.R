##### Libraries ####
library(tidyverse)
########

##### get_focals ####
#get predictor raster stack and variable vector, return stack with added focal rasters 

get_focals <- function(predictors,focal_list ){  
  for(r in 1:length(focal_list)){
    raster1 <- raster::focal(subset(predictors,focal_list[r]), w=matrix(1/25,nrow=5,ncol=5))
    names(raster1) <- paste0(focal_list[r],"_500")
    predictors <- raster::stack( predictors , raster1 )
  }
  return(predictors)
}

