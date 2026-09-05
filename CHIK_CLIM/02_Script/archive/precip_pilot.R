# ------------------------------------------------------------------------------
# Pilot prediction setup: train 100 RF models on current data
# (precipitation pilot — others identical to current)
# ------------------------------------------------------------------------------

# ---- Libraries ---------------------------------------------------------------
library(randomForest)
library(data.table)
library(dplyr)
library(foreach)
library(doParallel)
library(parallel)

# ---- Load tuning results -----------------------------------------------------
load("01_Data/hyperparam_results_raw_0625.RData")  # provides `model_results`, `train_list`
load("01_Data/mcmc_tot.RData")
# ---- Register parallel backend ----------------------------------------------
ncpus <- max(1, parallel::detectCores() - 1)
cl    <- parallel::makeCluster(ncpus)
doParallel::registerDoParallel(cl)
on.exit(parallel::stopCluster(cl), add = TRUE)


# 1. create 100 random FOI dataset and convert to sf ----------------------------------------
for(i in 1:100){
  
  current_dataset <- data.frame(foi = numeric(76), lat = numeric(76), long = numeric(76))
  
  for(k in 1:nrow(chik_foi)){
    # randomly sample logit scale FOI from each distribution
    #non_zero_values <- mcmc_tot[[k]][, 2][mcmc_tot[[k]][, 2] != 0]
    current_dataset$foi[k] <- sample(mcmc_tot[[k]][, 2], 1) # for raw foi model fitting
    #current_dataset$nonzero[k] <- sample(non_zero_values, 1)
    current_dataset$logfoi[k] <- sample(mcmc_tot[[k]][, "logFOI"], 1) # for logfoi model fitting
    current_dataset$logitfoi[k] <- sample(mcmc_tot[[k]][, "logitFOI"], 1) # for logitfoi model fitting
    # Assign corresponding latitude and longitude
    current_dataset$lat[k] <- mcmc_tot[[k]]$lat[1]
    current_dataset$long[k] <- mcmc_tot[[k]]$long[1]
    current_dataset$study_no[k] <- mcmc_tot[[k]]$study_no[1]
    current_dataset$ID[k]    <- mcmc_tot[[k]]$ID[1]
  }
  
  FOI[[i]] <- current_dataset
  
  # create blocks 
  FOI_sf <- st_as_sf(FOI[[i]], coords = c("long", "lat"))
  data_extent <- st_bbox(FOI_sf)
  FOI_sf_list[[i]] <- FOI_sf
}

FOI_df <- do.call(rbind, FOI_sf_list)

#2. Merge each of 100 FOI dataset with covariate dataset for a complete data ---------------

for(k in 1:length(FOI_sf_list)) {
  # Convert sf objects to regular data frames
  df1 <- as.data.frame(FOI_sf_list[[k]])
  
  # Merge the data frames based on "ID"
  merged_df <- left_join(df1, p_covs, by = "ID")
  
  merge_sf <- st_as_sf(merged_df, coords = c("Longitude", "Latitude"), crs = st_crs(FOI_sf_list[[k]]))
  
  merge_list[[k]] <- merge_sf
}

#3. Create grid blocks for each 100 FOI+covariate dataset (500x500km: 50 blocks each) -------------

for(i in 1:100){
  cell_size_deg <- 500 / 111
  
  # Create a grid of specified block size
  grid_blocks <- st_make_grid(merge_list[[i]], cellsize = c(cell_size_deg, cell_size_deg))
  grid_blocks_sf <- st_as_sf(grid_blocks)
  # Assign a unique block ID to each block
  grid_blocks_sf$blockID <- seq_len(nrow(grid_blocks_sf))
  grid_blocks_list[[i]] <- grid_blocks_sf
  
  # Extract longitude and latitude from geometry
  merge_list[[i]]$long <- st_coordinates(merge_list[[i]])[, 1]
  merge_list[[i]]$lat <- st_coordinates(merge_list[[i]])[, 2]
  
  # Perform a spatial join to identify which blocks contain data points
  blocks_with_data <- st_join(grid_blocks_list[[i]], merge_list[[i]], join = st_intersects)
  blocks_with_data <- blocks_with_data[!is.na(blocks_with_data$foi), ]
  block_data_list[[i]] <- blocks_with_data
}


#4. For each 100 FOIdataset, randomly sample (N-1) train blocks (where N = total N blocks) and leave one block out for test -------- 

for(i in 1:100){
  split <- split(block_data_list[[i]], block_data_list[[i]]$blockID)
  split_list[[i]] <- split
  sampled_indices <- sample(1:length(split_list[[i]]), length(split_list[[i]]) - 1, replace = FALSE)
  non_samp <- setdiff(1:length(split_list[[i]]), sampled_indices)
  train_set <- do.call(rbind, split_list[[i]][sampled_indices])
  train_list[[i]] <- train_set ## seems like almost 99% of data are being used for training
  test_set <- split_list[[i]][non_samp]
  test_list[[i]] <- do.call(rbind, test_set)
}

# ---- NA diagnostics on train_list -------------------------------------------
# train_list[[i]] is an sf object (carries a geometry list-column).
# We must strip geometry before complete.cases() to avoid:
#   "invalid 'type' (list) of argument".
needed_cols <- c("foi", "Tsuit", "PRCP", "GDP", "Albo", "Aegyp", "CHIKRisk")

na_per_set <- sapply(seq_along(train_list), function(i) {
  d <- sf::st_drop_geometry(train_list[[i]])
  sum(!stats::complete.cases(d[, needed_cols, drop = FALSE]))
})
cat(sprintf("Rows with NA in needed cols across %d train sets:\n", length(train_list)))
cat(sprintf("  min=%d  median=%d  max=%d  total=%d\n",
            min(na_per_set), stats::median(na_per_set),
            max(na_per_set), sum(na_per_set)))

# Per-column NA counts in train_list[[1]] (quick column-level diagnosis)
d1 <- sf::st_drop_geometry(train_list[[1]])
cat("Per-column NA counts in train_list[[1]]:\n")
print(sapply(needed_cols, function(c) sum(is.na(d1[[c]]))))

# ---- Train 100 RF models (current covariates only) --------------------------
rf_mod_hyper <- foreach(i = 1:100,
                        .packages = c("randomForest", "data.table", "dplyr", "sf")) %dopar% {
  best_params <- model_results[[i]]$best_params

  d <- sf::st_drop_geometry(train_list[[i]])
  d <- d[stats::complete.cases(d[, needed_cols, drop = FALSE]), , drop = FALSE]

  randomForest(
    foi ~ Tsuit + PRCP + GDP + Albo + Aegyp + CHIKRisk,
    data        = d,
    importance  = TRUE,
    mtry        = best_params$mtry,
    nodesize    = best_params$min.node.size,
    replace     = best_params$replace,
    sampsize    = floor(best_params$sample.fraction * nrow(d)),
    na.action   = na.omit
  )
}

# ---- Stop cluster (also handled by on.exit, but explicit is clearer) --------
parallel::stopCluster(cl)
