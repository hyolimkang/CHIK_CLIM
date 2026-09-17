library(randomForest)
library(dplyr)
library(tidyr)
library(ggplot2)
library(sf)
library(exactextractr)
library(raster)
library(snowfall)
library(viridisLite)
library(gbm)
library(spatialRF)
#library(MODIStsp)
library(sf)
library(geodata)
library(ecospat)
library(fields)
library(ncdf4)
library(png)
library(FNN)

load("MainData/chik_occ.RData")
source("Raster/fixNAs_adj.R")
source("Raster/thinning_chik.R")
source("D:/OneDrive - London School of Hygiene and Tropical Medicine/chik_mapping/plotRaster.R")
options(scipen = 999)
chik_foi<- read.csv("MainData/chik_foi.csv")


# distribution of raw data (foi) -----------------------------------------------
# https://rpubs.com/Alema/1000582
hist(chik_foi$mfoi)
ggplot(chik_foi, aes(x=mfoi)) + geom_histogram(color="black", fill="white")+
  geom_vline(aes(xintercept=mean(mfoi)),
             color="blue", linetype="dashed", size=1)+
  geom_density(alpha=.2, fill="#FF6666")+
  theme_bw()
par(mfrow = c(1,1))
boxplot(chik_foi$mfoi,
        ylab = "FOI"
)

out <- boxplot.stats(chik_foi$mfoi)$out
out_ind <- which(chik_foi$mfoi %in% c(out))
out_ind
outliers <- chik_foi[out_ind, ]

# ------------------------------------------------------------------------------
raster.scale <- function(ras, logT = TRUE){
  vec = as.vector(ras) # converting raster into a vector
  if(logT){
    zinf <- min(vec[vec > 0], na.rm = T) # finding the smallest possible values, w/o NA
    print(paste("zinf:", zinf)) # prints out smallest possible positive value
    vec = log(vec + 0.5 * zinf) # adds half of the smallest value to all raster values and takes log
    print(head(vec))
  }
  vec = scale(vec) # scale raster values to make mean = 0 and sd = 1
  values(ras) = vec # assigns scaled values back to the raster
  return(ras)
}
raster.scale.neg <- function(ras, logT = TRUE){
  vec = as.vector(ras)
  if(logT){
    min_val <- min(vec, na.rm = TRUE) # This is the minimum value in the data
    offset <- abs(min_val) + 1        # Calculate the offset to make all values positive
    print(paste("offset:", offset))
    vec <- log(vec + offset)          # Apply the offset and then the log transformation
    print(head(vec))
  }
  vec <- scale(vec)                   # Scale the transformed data
  values(ras) <- vec                  # Assign the scaled values back to the raster
  return(ras)
}
# 1. covariate1: temp suit------------------------------------------------------
## 2010-2020 
tsuit <- "01_Data/final/Dengue_temperature_suitaiblity_masked_.tif"
# inspect data distribution to determine if logT should be used 
tsuit.hist <- raster(tsuit)
hist(values(tsuit.hist)) # left skewed --> when highly skewed, logarithmic trasformation is needed
tsuit <- raster.scale(raster(tsuit))
#temp <- raster(temp)
plot(tsuit)

points(p_covs$Longitude, p_covs$Latitude, pch = 20, col = 'red')
temp_colors <- colorRampPalette(c("blue", "green", "yellow", "red"))(length(unique(p_covs$Temp)))
point_colors <- temp_colors[as.factor(p_covs$Tsuit)]
points(p_covs$Longitude, p_covs$Latitude, pch = 20, col = point_colors)
image.plot(legend.only = TRUE, zlim = c(min(p_covs$Temp), max(p_covs$Temp)), col = temp_colors, 
           horizontal = FALSE, legend.lab = "Temperature")

# 2. covariate1: temp min ------------------------------------------------------
## 2010-2020 
tmin <- "D:/OneDrive - London School of Hygiene and Tropical Medicine/chik_mapping/final/Tmin_TerraClim_2010_2020_005dg_masked_.tif"
tmin <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/chik_mapping/final/Tmin_TerraClim_2010_2020_005dg_masked_.tif"
# inspect data distribution to determine if logT should be used 
tmin.hist <- raster(tmin)
hist(values(tmin.hist)) # pretty normally distributed (no need for logT)
tmin <- raster.scale(raster(tmin), logT = FALSE)
#temp <- raster(temp)
plot(tmin)

points(p_covs$Longitude, p_covs$Latitude, pch = 20, col = 'red')
temp_colors <- colorRampPalette(c("blue", "green", "yellow", "red"))(length(unique(p_covs$Temp)))
point_colors <- temp_colors[as.factor(p_covs$Temp)]
points(p_covs$Longitude, p_covs$Latitude, pch = 20, col = point_colors)
image.plot(legend.only = TRUE, zlim = c(min(p_covs$Temp), max(p_covs$Temp)), col = temp_colors, 
           horizontal = FALSE, legend.lab = "Temperature")

# 3. covariate3: mean precipitation --------------------------------------------
## 2010-2020 
precip <- "D:/OneDrive - London School of Hygiene and Tropical Medicine/chik_mapping/final/PRCP_TerraClim_2010_2020_005dg_masked_.tif"
precip <- "D:/OneDrive - London School of Hygiene and Tropical Medicine/chik_mapping/covlist_upd/PRCP_TerraClim_2010_2020_005dg_masked_.tif"
precip <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/chik_mapping/covlist_upd/PRCP_TerraClim_2010_2020_005dg_masked_.tif"
# inspect data distribution to determine if logT should be used 
precip.hist <- raster(precip)
hist(values(precip.hist)) # right skewed -> need scaling
precip <- raster.scale(raster(precip)) # scaling 
#precip <- raster(precip) # not scaling
plot(precip)

prcp_colors <- colorRampPalette(c("blue", "green", "yellow", "red"))(length(unique(p_covs$PRCP)))
point_colors <- temp_colors[as.factor(p_covs$PRCP)]
points(p_covs$Longitude, p_covs$Latitude, pch = 20, col = point_colors)
image.plot(legend.only = TRUE, zlim = c(min(p_covs$PRCP), max(p_covs$PRCP)), col = temp_colors, 
           horizontal = FALSE, legend.lab = "PRCP")

# 4. covariate4: elevation -----------------------------------------------------
## 1970-2000 average WorldClim
elev <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/chik_mapping/final/wc2.1_5m_elev.tif"
elev.hist <- raster(elev) # skewed 
hist(values(elev.hist)) 
elev <- raster.scale.neg(raster(elev))
plot(elev)

elev_colors <- colorRampPalette(c("blue", "green", "yellow", "red"))(length(unique(p_covs$Elev)))
point_colors <- temp_colors[as.factor(p_covs$Elev)]
points(p_covs$Longitude, p_covs$Latitude, pch = 20, col = point_colors)
image.plot(legend.only = TRUE, zlim = c(min(p_covs$Elev), max(p_covs$Elev)), col = temp_colors, 
           horizontal = FALSE, legend.lab = "Elev")

# 5. covariate4: population density --------------------------------------------
## 2000, 2005, 2010, 2015, 2020 (all years combined average; SEDAC)
pop_dens <- "D:/OneDrive - London School of Hygiene and Tropical Medicine/chik_mapping/final/landscan_global_2022_masked.tif"
pop_dens <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/chik_mapping/final/landscan_global_2022_masked.tif"
pop_dens <- raster.scale(raster(pop_dens))
#pop_dens <- raster(pop_dens)
plot(pop_dens)

pop_colors <- colorRampPalette(c("blue", "green", "yellow", "red"))(length(unique(p_covs$Pop_dens)))
point_colors <- temp_colors[as.factor(p_covs$Pop_dens)]
points(p_covs$Longitude, p_covs$Latitude, pch = 20, col = point_colors)
image.plot(legend.only = TRUE, zlim = c(min(p_covs$Pop_dens), max(p_covs$Pop_dens)), col = temp_colors, 
           horizontal = FALSE, legend.lab = "pop")

# 6. covariate5: GDP national --------------------------------------------------
## 2009-2019
gdp_nat <- "D:/OneDrive - London School of Hygiene and Tropical Medicine/chik_mapping/final/GDP_2009_2019_National_masked.tif"
gdp_nat <- "D:/OneDrive - London School of Hygiene and Tropical Medicine/chik_mapping/covlist_upd/GDP_2009_2019_National_masked_.tif"
gdp_nat <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/chik_mapping/covlist_upd/GDP_2009_2019_National_masked_.tif"
gdp_nat.hist <- raster(gdp_nat) # skewed 
hist(values(gdp_nat.hist))
gdp_nat <- raster.scale(raster(gdp_nat))
#gdp_nat <- raster(gdp_nat)
plot(gdp_nat)


# 7. covariate7: Albopictus ----------------------------------------------------
## 2020
albo <- "D:/OneDrive - London School of Hygiene and Tropical Medicine/chik_mapping/final/Albopictus_mean_2020_rcp60_spreadXsuit_masked_.tif"
albo <- "D:/OneDrive - London School of Hygiene and Tropical Medicine/chik_mapping/covlist_upd/Albopictus_mean_2020_rcp60_spreadXsuit_masked_.tif"
albo <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/chik_mapping/covlist_upd/Albopictus_mean_2020_rcp60_spreadXsuit_masked_.tif"
albo.hist <- raster(albo) # skewed 
hist(values(albo.hist))
albo <- raster.scale(raster(albo), logT= FALSE)
#albo <- raster(albo)
plot(albo)

# 8. covariate8: Aegypti -------------------------------------------------------
## 2020
aegyp <- "D:/OneDrive - London School of Hygiene and Tropical Medicine/chik_mapping/final/Aegypti_mean_2020_rcp60_spreadXsuit_masked_.tif"
aegyp <- "D:/OneDrive - London School of Hygiene and Tropical Medicine/chik_mapping/covlist_upd/Aegypti_mean_2020_rcp60_spreadXsuit_masked_.tif"
aegyp <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/chik_mapping/covlist_upd/Aegypti_mean_2020_rcp60_spreadXsuit_masked_.tif"
aegyp.hist <- raster(aegyp) # skewed 
hist(values(aegyp.hist))
aegyp <- raster.scale(raster(aegyp), logT= FALSE)
#aegyp <- raster(aegyp)
plot(aegyp)

# 9. covariate9: NDVI ----------------------------------------------------------
## 2010-2020
ndvi <- "D:/OneDrive - London School of Hygiene and Tropical Medicine/chik_mapping/final/NDVI_2010_2020_005dg_masked_.tif"
ndvi <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/chik_mapping/final/NDVI_2010_2020_005dg_masked_.tif"
ndvi.hist <- raster(ndvi) 
hist(values(ndvi.hist))
ndvi <- raster.scale(raster(ndvi), logT = FALSE)
#ndvi <- raster(ndvi)
plot(ndvi)

# 10. covariate10: chikungunya risk map ----------------------------------------------------------
## 2010-2020
chikrisk <- "D:/OneDrive - London School of Hygiene and Tropical Medicine/chik_mapping/chik_data_ahyoung/CHIK_riskmap_wmean_masked.tif"
chikrisk <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/chik_mapping/chik_data_ahyoung/CHIK_riskmap_wmean_masked.tif"
chikrisk.hist <- raster(chikrisk) 
plot(chikrisk.hist)
hist(values(chikrisk.hist))
chikrisk <- raster.scale(raster(chikrisk))
plot(chikrisk)

#11. covariate11: urbanisation -------------------------------------------------
urban <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/chik_mapping/final/urban_2010_2020_1km_masked_.tif"
urban.hist <- raster(urban) 
hist(values(urban.hist))
urban <- raster.scale(raster(urban))
plot(urban)

#12. covariate11: GDP per capita -----------------------------------------------
gdp_per <- "D:/OneDrive - London School of Hygiene and Tropical Medicine/chik_mapping/final/GDP_2009_2019_1km_masked_.tif"
gdp_per <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/chik_mapping/final/GDP_2009_2019_1km_masked_.tif"
gdp_per.hist <- raster(gdp_per) 
hist(values(gdp_per.hist))
gdp_per <- raster.scale(raster(gdp_per))
plot(gdp_per)

#13. covariate12: DHI ----------------------------------------------------------
dhi <- "D:/OneDrive - London School of Hygiene and Tropical Medicine/chik_mapping/final/DHI_global_clusters_14c_sv20_masked_.tif"
dhi <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/chik_mapping/final/DHI_global_clusters_14c_sv20_masked_.tif"
dhi.hist <- raster(dhi) 
hist(values(dhi.hist))
dhi <- raster(dhi)
plot(dhi)


#13. combining histograms
par(mfrow = c(3,4))
hist(values(tsuit.hist), main = "Tsuit")
hist(values(tmin.hist), main = "Tmin") # pretty normally distributed (no need for logT)
hist(values(precip.hist), main = "Precip") # right skewed -> need scaling
hist(values(elev.hist), main = "Elevation")
hist(values(gdp_nat.hist), main = "National GDP")
hist(values(gdp_per.hist), main = "GDP per capita")
hist(values(albo.hist), main = "Albopictus")
hist(values(aegyp.hist), main = "Aegypti")
hist(values(ndvi.hist), main = "NDVI")
hist(values(chikrisk.hist), main = "CHIK risk")
hist(values(urban.hist), main = "Urban")


#14. aligning covariates -------------------------------------------------------
# extent, ncells, resolutions
tmin    <- projectRaster(tmin, tsuit, method="bilinear")
precip <- projectRaster(precip, tsuit, method="bilinear")
gdp_nat   <- projectRaster(gdp_nat, tsuit, method="bilinear")
albo   <- projectRaster(albo, tsuit, method="bilinear")
aegyp   <- projectRaster(aegyp, tsuit, method="bilinear")
ndvi   <- projectRaster(ndvi, tsuit, method="bilinear")
pop_dens <- projectRaster(pop_dens, tsuit, method="bilinear")
chikrisk <- resample(chikrisk, tsuit, method="bilinear")
urban <- resample(urban, tsuit, method="bilinear")
gdp_per <- resample(gdp_per, tsuit, method="bilinear")
dhi <- resample(dhi, tsuit, method="bilinear")
elev <- resample(elev, tsuit, method="bilinear")

#15. stack covariates ----------------------------------------------------------
covlist <- stack(tsuit, tmin, precip, pop_dens, gdp_nat, albo, aegyp, ndvi, elev, chikrisk)
names(covlist) <- c("Tsuit", "Tmin", "PRCP", "Pop_dens", "GDP", "Albo", "Aegyp", "NDVI", "Elev", "CHIKRisk")
plot(covlist)

covlist1 <- stack(tsuit, tmin, precip, gdp_nat, albo, aegyp, ndvi, chikrisk, urban, gdp_per, pop_dens)
names(covlist1) <- c("Tsuit", "Tmin", "PRCP", "GDP", "Albo", "Aegyp", "NDVI", "CHIKRisk", "Urban", "GDP_capita", "Pop_dens")
plot(covlist1)

covlist_upd <- stack(tsuit, tmin, precip, gdp_nat, albo, aegyp, chikrisk, gdp_per, pop_dens, ndvi, dhi)
names(covlist_upd) <- c("Tsuit", "Tmin", "PRCP", "GDP", "Albo", "Aegyp","CHIKRisk", "GDP_cap", "pop_dens", "NDVI", "DHI")
plot(covlist_upd)


#16. pcov data -----------------------------------------------------------------
# extract covariate values
# create only occurrence + foi
chik_foi <- chik_foi[,c(5:9, 15)]
new_order <- c("long", "lat", "mfoi")
chik_foi <- chik_foi[,new_order]
chik_occ <- chik_occ[,c(1:2)]
names(chik_occ) <- c("long", "lat")
chik_foi_merge <- chik_foi[,c(1:2)]
chik_occ <- rbind(chik_foi_merge, chik_occ)

template = tsuit

#chik_bg <- thin_bg_chik(chik_occ)
chik_bg <- thin_bg_buffer_chik_interpolated(chik_occ) # create buffer around foi and occ points by 5x5 km
#chik_bg <- thin_bg_buffer_chik(chik_occ) # create buffer around foi and occ points by 5x5 km

chik_bg <- chik_bg[,c(3:4)]
names(chik_bg) <- c("long", "lat")
chik_bg$mfoi <- 0

chik_foi_all <- rbind(chik_foi, chik_bg)

p_covs <- extract(covlist_upd, data.frame(chik_foi$long, chik_foi$lat), df= T, na.rm=TRUE)
p_covs$FOI <- chik_foi$mfoi
p_covs$FOI_lo <- chik_foi$mfoi_lo
p_covs$FOI_hi <- chik_foi$mfoi_hi
p_covs$Longitude <- chik_foi$long
p_covs$Latitude <- chik_foi$lat
p_covs$NumTested <- chik_foi$NumTested
p_covs$region <- chik_foi$region
p_covs$study_no <- chik_foi$study_no
p_covs$country <- chik_foi$country
p_covs <- fixNAs(p_covs, covlist_upd)
(point_idx <- which(apply(p_covs[, names(covlist_upd)], 1, function(row) any(is.na(row)))))


# fill out missing value in row 50 
p_covs_sf <- st_as_sf(p_covs, coords = c("Longitude", "Latitude"))

ggplot(data = p_covs_sf) +
  geom_sf() +  # Plot the spatial data
  geom_sf_text(aes(label = ID), check_overlap = TRUE, nudge_y = 0.01) +
  theme_minimal()

coords_matrix <- st_coordinates(p_covs_sf)

# Find the row with NA (for example, row 50)
na_row <- 50

# Find nearest neighbors (for example, the nearest 3 neighbors)
neighbors <- get.knnx(data = coords_matrix, query = coords_matrix[na_row, , drop = FALSE], k = 3)
mean_value_tsuit <- mean(p_covs_sf$Tsuit[neighbors$nn.index], na.rm = TRUE)
mean_value_tmin <- mean(p_covs_sf$Tmin[neighbors$nn.index], na.rm = TRUE)
mean_value_prcp <- mean(p_covs_sf$PRCP[neighbors$nn.index], na.rm = TRUE)
mean_value_gdp <- mean(p_covs_sf$GDP[neighbors$nn.index], na.rm = TRUE)
mean_value_albo <- mean(p_covs_sf$Albo[neighbors$nn.index], na.rm = TRUE)
mean_value_aegyp <- mean(p_covs_sf$Aegyp[neighbors$nn.index], na.rm = TRUE)
mean_value_risk <- mean(p_covs_sf$CHIKRisk[neighbors$nn.index], na.rm = TRUE)
mean_value_gdp_per <- mean(p_covs_sf$GDP_cap[neighbors$nn.index], na.rm = TRUE)
mean_value_pop_dens <- mean(p_covs_sf$pop_dens[neighbors$nn.index], na.rm = TRUE)
mean_value_ndvi <- mean(p_covs_sf$NDVI[neighbors$nn.index], na.rm = TRUE)
mean_value_dhi <- mean(p_covs_sf$DHI[neighbors$nn.index], na.rm = TRUE)

p_covs$Tsuit[50] <- mean_value_tsuit
p_covs$Tmin[50] <- mean_value_tmin
p_covs$PRCP[50] <- mean_value_prcp
p_covs$GDP[50] <- mean_value_gdp
p_covs$Albo[50] <- mean_value_albo
p_covs$Aegyp[50] <- mean_value_aegyp
p_covs$CHIKRisk[50] <- mean_value_risk
p_covs$GDP_cap[50] <- mean_value_gdp_per
p_covs$pop_dens[50] <- mean_value_pop_dens
p_covs$NDVI[50] <- mean_value_ndvi
p_covs$DHI[50] <- mean_value_dhi

# mapping foi and bg points 
chik_sf <- chik_foi_all
coordinates(chik_sf) <- ~long + lat
crs(chik_sf) <- crs(template)

# Convert the raster to a dataframe
template_df <- as.data.frame(template, xy=TRUE)
names(template_df) <- c("x", "y", "layer")
# Convert chik_sf to an sf object
chik_sf_obj <- st_as_sf(chik_sf)
chik_sf_df <- as.data.frame(chik_sf_obj)
chik_sf_obj_points <- st_coordinates(chik_sf_obj)
chik_sf_df$long <- chik_sf_obj_points[, "X"]
chik_sf_df$lat <- chik_sf_obj_points[, "Y"]

chik_sf_df <- chik_sf_df %>%
  mutate(mfoi_category = ifelse(mfoi == 0, "0", ">0"))

ggplot() +
  geom_raster(data = template_df, aes(x = x, y = y, fill = layer), alpha = 1) +
  geom_point(data = chik_sf_df, aes(x = long, y = lat, color = mfoi_category), size = 1) +
  scale_color_manual(values = c("0" = "red", ">0" = "blue")) +
  coord_fixed() +
  theme_bw()+
  guides(fill = "none")


#17. interpolating FOI in presence data using FOI data -------------------------
# number of bg points = nrow(chik_occ)

load("chik_occ.RData")
chik_occ <- chik_occ[,c(1:2)]
names(chik_occ) <- c("long", "lat")
chik_occ$mfoi <- NA # original occ data is a binary presence data (no FOI)
chik_foi_all_inter <- rbind(chik_foi, chik_occ) # first merge occ + foi data and then interpolate

pr_covs <- extract(covlist, data.frame(chik_foi_all_inter$long, chik_foi_all_inter$lat), df= T, na.rm=TRUE)
pr_covs$FOI <- chik_foi_all_inter$mfoi
pr_covs$Longitude <- chik_foi_all_inter$long
pr_covs$Latitude  <- chik_foi_all_inter$lat
pr_covs <- fixNAs(pr_covs, covlist)
pr_point_idx <- which(apply(pr_covs[, names(covlist)], 1, function(row) any(is.na(row))))

# train FOI model
foi_model <- randomForest(FOI ~ Temp + PRCP + Pop_dens + GDP + Albo + Aegyp + NDVI +Elev,
                          data = pr_covs,
                          na.action = na.omit)

pr_covs$FOI <- predict(foi_model, pr_covs) # interpolated FOI values in the occ data

# background point generation after merging foi + occ and interpolation
chik_occ <- pr_covs[, c("Longitude", "Latitude", "FOI")]
names(chik_occ) <- c("long", "lat", "mfoi") # extract occ + FOI data (total 4570 FOI values)
chik_bg_inter <- thin_bg_buffer_chik_interpolated(chik_occ) # then create 1:1 pseudo-absence points based on 4570 FOI points
chik_bg_inter <- chik_bg_inter[, c(3:4)]
names(chik_bg_inter) <- c("long", "lat")
chik_bg_inter$mfoi <- 0 
chik_foi_all_inter <- rbind(chik_occ, chik_bg_inter)

pr_covs <- extract(covlist, data.frame(chik_foi_all_inter$long, chik_foi_all_inter$lat), df= T, na.rm=TRUE)
pr_covs$FOI <- chik_foi_all_inter$mfoi
pr_covs$Longitude <- chik_foi_all_inter$long
pr_covs$Latitude  <- chik_foi_all_inter$lat
pr_covs <- fixNAs(pr_covs, covlist)
pr_point_idx <- which(apply(pr_covs[, names(covlist)], 1, function(row) any(is.na(row))))
pr_covs <- pr_covs[-pr_point_idx,]

# mapping
chik_foi_all_inter <- pr_covs[,c(10:12)]
chik_sf_all_inter <- chik_foi_all_inter
coordinates(chik_sf_all_inter) <- ~Longitude + Latitude
crs(chik_sf_all_inter) <- crs(template)

chik_sf_obj_inter <- st_as_sf(chik_sf_all_inter)
chik_sf_df_inter <- as.data.frame(chik_sf_obj_inter)
chik_sf_obj_points_inter <- st_coordinates(chik_sf_obj_inter)
chik_sf_df_inter$long <- chik_sf_obj_points_inter[, "X"]
chik_sf_df_inter$lat <- chik_sf_obj_points_inter[, "Y"]

chik_sf_df_inter <- chik_sf_df_inter %>%
  mutate(FOI_category = ifelse(FOI == 0, "0", ">0"))

ggplot() +
  geom_raster(data = template_df, aes(x = x, y = y, fill = layer), alpha = 0.8) +
  geom_point(data = chik_sf_df_inter, aes(x = long, y = lat, color = FOI_category), size = 0.8) +
  scale_color_manual(values = c("0" = "red", ">0" = "blue")) +
  coord_fixed() +
  theme_bw()+
  guides(fill = "none")

#17. save data------------------------------------------------------------------
save("tsuit", "tmin", "precip", "elev", "pop_dens", "gdp_nat", "albo", "aegyp", "ndvi", "chikrisk", "covlist", 
     "p_covs",
     file = "covariates.RData")

save("tsuit", "tmin", "precip", "elev", "pop_dens", "gdp_nat", "albo", "aegyp", "ndvi", "chikrisk", "covlist1", 
     "p_covs", "urban", "chikrisk", "gdp_per",
     file = "covariates1.RData")

save("tsuit", "tmin", "precip", "gdp_nat", "albo", "aegyp","chikrisk", "pop_dens", "gdp_per", "ndvi", "dhi", "covlist_upd", 
     "p_covs",
     file = "covariates_upd.RData")

#18. correlation analysis ------------------------------------------------------
cov_df <- as.data.frame(stack(tsuit, precip, gdp_nat, albo, aegyp), xy = TRUE)
colnames(cov_df) <- c("x", "y", "tsuit", "precip", "gdp", "albo", "aegyp")
coordinates(chik_foi) <- ~long + lat

chik_foi$tsuit <- extract(tsuit, chik_foi)
chik_foi$precip <- extract(precip, chik_foi)
chik_foi$gdp <- extract(gdp_nat, chik_foi)
chik_foi$albo <- extract(albo, chik_foi)
chik_foi$aegyp <- extract(aegyp, chik_foi)

model <- lm(mfoi ~ tsuit + precip + gdp + albo + aegyp, data = chik_foi)
summary(model)

chik_foi_df <- as.data.frame(chik_foi)

model1 <- gam(FOI ~ s(Tsuit), data = p_covs, weights = NumTested)
model2 <- gam(FOI ~ s(PRCP), data = p_covs, weights = NumTested)
model3 <- gam(FOI ~ s(GDP), data = p_covs, weights = NumTested)
model4 <- gam(FOI ~ s(Albo), data = p_covs, weights = NumTested)
model5 <- gam(FOI ~ s(Aegyp), data = p_covs, weights = NumTested)

(r2_model1 <- summary(model1)$r.sq)
(r2_model2 <- summary(model2)$r.sq)
(r2_model3 <- summary(model3)$r.sq)
(r2_model4 <- summary(model4)$r.sq)
(r2_model5 <- summary(model5)$r.sq)

pdf("gamcorr.pdf", height = 12, width = 20)

par(mfrow = c(3, 2))  # Set up the plotting area for multiple plots
plot(model1, main = "GAM for Temperature suitability", shade = TRUE, ylab = "Estimated Effect on FOI")
plot(model2, main = "GAM for Precipitation", shade = TRUE, ylab = "Estimated Effect on FOI")
plot(model3, main = "GAM for National GDP", shade = TRUE, ylab = "Estimated Effect on FOI")
plot(model4, main = "GAM for Albopictus distribution", shade = TRUE, ylab = "Estimated Effect on FOI")
plot(model5, main = "GAM for Aegypti distribution", shade = TRUE, ylab = "Estimated Effect on FOI")

dev.off()


GDP_nat <- values(gdp_nat)
GDP_capita <- values(gdp_per)
Tsuit <- values(tsuit)
Tmin <- values(tmin)
Aegyp <- values(aegyp)
Albo <- values(albo)
Pop_dens <- values(pop_dens)
Precip <- values(precip)
NDVI <- values(ndvi)
CHIK_risk <- values(chikrisk)
urban <- values(urban)
elev <- values(elev)
# Combine values into a data frame
data <- data.frame(GDP_nat, GDP_capita, Tsuit, Tmin, Aegyp, 
                   Albo, urban, Precip, NDVI, CHIK_risk)
data <- na.omit(data)
corr_matrix <- cor(data, use = "complete.obs")

ggcorrplot(corr_matrix, lab = TRUE)