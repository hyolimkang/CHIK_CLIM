
## does this function checks whether the bg points overlap with current occurrence or FOI data locations?
## how to find out true absence points?

thin_bg_chik <- function(occ){
  
  chik_occ <- occ # create the same dataframe before transforming it into SpatialPoints
    
  
  trv <- 1:length(as.vector(template)) # unique pixel values
  trv = trv * !is.na(as.vector(template))
  pixelID = template
  values(pixelID) = trv
  

  
  bg <- data.frame(Longitude = runif(nrow(occ) * 100, min = -180, max = 180),
                   Latitude = runif(nrow(occ) * 100, min = -59, max = 75))
  
  
  # converting occ into SpatialPoints: now the occ is spatial form. 
  coordinates(occ) <- ~ long + lat
  crs(occ) <- crs(template)
  occRaster <- rasterize(occ, template, fun=function(...) 1, background=NA)
  
  
  # trim to those on land and trim to size of occurrence data
  onland <- !is.na(raster::extract(template, bg[, 1:2]))
  bg = bg[onland, ]
  
  # Exclude points in occurrence areas
  notInOcc <- is.na(raster::extract(occRaster, bg[, 1:2]))
  bg = bg[notInOcc, ]
  

  # thinning
  bg$uID <- raster::extract(pixelID, bg[, 1:2])
  # remove duplicates
  bg = bg[!duplicated(bg$uID), ]
  
  # trim to size of 10 times of FOI dataset
  
  if (nrow(bg) >= nrow(chik_occ)) {
    bg = bg[1:nrow(chik_foi), ]
  }
} 


