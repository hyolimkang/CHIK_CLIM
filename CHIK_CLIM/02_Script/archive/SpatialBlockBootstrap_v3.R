
FOI <- list()
FOI_sf_list <- list()
grid_blocks_list <- list()
block_data_list <- list()
split_list <- list()
train_list <- list()
test_list <- list()
merge_list <- list()

load("MainData/mcmc_tot.RData")
load("MainData/hyperparam_results_raw_0625.RData") # with CHIK risk 
load("MainData/hyperparam_results_log_1010.RData") # with CHIK risk 

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

## plot
p_covs_sf <- st_as_sf(p_covs, coords = c("Longitude", "Latitude"))

ggplot(data = block_data_list[[1]]) +
  geom_sf(aes(fill = foi)) +  # Use meanRMSE for the color fill
  scale_fill_viridis_c() +
  geom_sf(data = p_covs_sf, color = 'red', size = 0.4) +  # Overlay other spatial data if needed
  labs(fill = "FOI") +
  theme_minimal() +
  theme(legend.position = "right")

# hyper-parameter tuning--------------------------------------------------------
# create hyperparameter grid
n_features <- 6

hyper_grid <- expand.grid(
  mtry = floor(n_features * c(.25, .4, 0.5, 0.7, 0.9)),
  min.node.size = c(1, 3, 5, 10), 
  replace = c(TRUE, FALSE),                               
  sample.fraction = seq(0.1, 0.8, by = 0.1),                       
  rmse = NA                                               
)


#5. Run RF model using logit FOI for each 100 FOI dataset-----------------------
### log FOI hyper param tuned
ncpus = 12
cl <- makeCluster(ncpus) 
registerDoParallel(cl) 
start_time <- Sys.time()
rf_mod_hyper <- foreach(i = 1:100, .packages = c("randomForest", "data.table", "dplyr")) %dopar% {
  randomForest(logfoi ~ Tsuit + PRCP + CHIKRisk + Albo + Aegyp + GDP,
                                 data = train_list[[i]],
                                 importance = TRUE,
               mtry = 4,
               nodesize = 1, # nodesize corresponds to min.node.size
               replace = FALSE, # sampling without replacement
               sampsize = floor(0.8 * nrow(train_list[[i]])) 
  )
}

## all model tuning
ncpus = 12
cl <- makeCluster(ncpus-1) 
registerDoParallel(cl) 
start_time <- Sys.time()
rf_mod_hyper <- foreach(i = 1:100, .packages = c("randomForest", "data.table", "dplyr")) %dopar% {
  # Extract the best hyperparameters for the ith model
  best_params <- model_results[[i]]$best_params
  
  # Train the random forest model using the best hyperparameters
  randomForest(
    foi ~ Tsuit + PRCP + GDP + Albo + Aegyp + CHIKRisk,
    data = train_list[[i]],
    importance = TRUE,
    mtry = best_params$mtry,
    nodesize = best_params$min.node.size,  # Ensure this name matches your hyper_grid
    replace = best_params$replace,
    sampsize = floor(best_params$sample.fraction * nrow(train_list[[i]]))  # Convert fraction to actual size
  )
}
stopCluster(cl)
end_time <- Sys.time()
time_taken <- end_time - start_time
print(time_taken) #  4.877229  secs

## raw processing w/o parallel
start_time <- Sys.time()

rf_mod_hyper <- lapply(1:100, function(i) {
  # Extract the best hyperparameters for the ith model
  best_params <- model_results[[i]]$best_params
  
  # Train the random forest model using the best hyperparameters
  randomForest(
    foi ~ Tsuit + PRCP + GDP + Albo + Aegyp + NDVI,
    data = train_list[[i]],
    importance = TRUE,
    mtry = best_params$mtry,
    nodesize = best_params$min.node.size,
    replace = best_params$replace,
    sampsize = floor(best_params$sample.fraction * nrow(train_list[[i]]))
  )
})
end_time <- Sys.time()
time_taken <- end_time - start_time
print(time_taken)

## log hypertune
ncpus = 12
cl <- makeCluster(ncpus) 
registerDoParallel(cl) 
start_time <- Sys.time()
rf_mod_hyper_log <- foreach(i = 1:100, .packages = c("randomForest", "data.table", "dplyr")) %dopar% {
  # Extract the best hyperparameters for the ith model
  best_params <- model_results_log[[i]]$best_params
  
  # Train the random forest model using the best hyperparameters
  randomForest(
    logfoi ~ Tsuit + PRCP + CHIKRisk + Albo + Aegyp + GDP,
    data = train_list[[i]],
    importance = TRUE,
    mtry = best_params$mtry,
    nodesize = best_params$min.node.size,  # Ensure this name matches your hyper_grid
    replace = best_params$replace,
    sampsize = floor(best_params$sample.fraction * nrow(train_list[[i]]))  # Convert fraction to actual size
  )
}
stopCluster(cl)
end_time <- Sys.time()
time_taken <- end_time - start_time
print(time_taken) #  4.877229  secs


### logit FOI 
ncpus = 12
cl <- makeCluster(ncpus) 
registerDoParallel(cl) 
start_time <- Sys.time()
rf_mod_logit <- foreach(i = 1:100, .packages = c("randomForest", "data.table", "dplyr")) %dopar% {
  randomForest(logitfoi ~ Tsuit + PRCP + CHIKRisk + Albo + Aegyp + GDP,
               data = train_list[[i]],
               importance = TRUE) 
}
stopCluster(cl)
end_time <- Sys.time()
time_taken <- end_time - start_time
print(time_taken) #  4.877229  secs

### log FOI (w/o tuning) 
ncpus = 12
cl <- makeCluster(ncpus) 
registerDoParallel(cl) 
start_time <- Sys.time()
rf_mod_log <- foreach(i = 1:100, .packages = c("randomForest", "data.table", "dplyr")) %dopar% {
  randomForest(logfoi ~ Tsuit + PRCP + CHIKRisk + Albo + Aegyp + GDP,
               data = train_list[[i]],
               importance = TRUE) 
}
stopCluster(cl)
end_time <- Sys.time()
time_taken <- end_time - start_time
print(time_taken) #  4.877229  secs


####### raw foi 
ncpus = 15
cl <- makeCluster(ncpus) 
registerDoParallel(cl) 

start_time <- Sys.time()
rf_mod <- foreach(i = 1:100, .packages = c("randomForest", "data.table", "dplyr")) %dopar% {
  randomForest(foi ~ Tsuit + PRCP + CHIKRisk + Albo + Aegyp + GDP,
               data = train_list[[i]],
               importance = TRUE)
}
stopCluster(cl)
end_time <- Sys.time()
time_taken <- end_time - start_time
print(time_taken) #  4.877229  secs


### GAM model 
ncpus = 26
cl <- makeCluster(ncpus) 
registerDoParallel(cl) 

start_time <- Sys.time()
gam_mod <- foreach(i = 1:100, .packages = c("mgcv", "data.table", "dplyr")) %dopar% {
  gam(foi ~ s(Tsuit) + s(PRCP) + s(GDP) + s(Albo) + s(Aegyp) + s(CHIKRisk),
      data = train_list[[i]],
      family = betar(link="logit"))
}
stopCluster(cl)
end_time <- Sys.time()
time_taken <- end_time - start_time
print(time_taken) #  4.877229  secs


# 6. Predict FOI for 1 held-out block in each of 100 FOI dataset ---------------
ncpus = 12
cl <- makeCluster(ncpus) 
registerDoParallel(cl) 

start_time <- Sys.time()
prediction <- foreach(i = 1:100, .packages = c("randomForest", "data.table", "dplyr")) %dopar% {
  predict(rf_mod[[i]], newdata  = test_list[[i]], type = "response")
} # 2.752589 mins
stopCluster(cl)
end_time <- Sys.time()
time_taken <- end_time - start_time


ncpus = 12
#Find out how many cores are available (if you don't already know)
cores<-detectCores()
#Create cluster with desired number of cores, leave one open for the machine         
#core processes
cl <- makeCluster(27)
registerDoParallel(cl) 

start_time <- Sys.time()
prediction_test <- foreach(i = 1:100, .packages = c("randomForest", "data.table", "dplyr")) %dopar% {
  predict(rf_mod_hyper[[i]], newdata  = test_list[[i]], type = "response")
} # 2.752589 mins
stopCluster(cl)
end_time <- Sys.time()
time_taken <- end_time - start_time
print(time_taken) # 2.884095 mins

prediction_test <- lapply(1:100, function(i) {
  predict(rf_mod_hyper[[i]], newdata = test_list[[i]], type = "response")
})

start_time <- Sys.time()
prediction_log_test_hyper <- foreach(i = 1:100, .packages = c("randomForest", "data.table", "dplyr")) %dopar% {
  exp(predict(rf_mod_hyper[[i]], newdata  = test_list[[i]], type = "response"))
} # 2.752589 mins

prediction_log_test_hyper <-lapply(1:100, function(i) {
  exp(predict(rf_mod_hyper[[i]], newdata  = test_list[[i]], type = "response"))
})

stopCluster(cl)
end_time <- Sys.time()
time_taken <- end_time - start_time
print(time_taken) # 2.884095 mins

start_time <- Sys.time()
prediction_log_test <- foreach(i = 1:100, .packages = c("randomForest", "data.table", "dplyr")) %dopar% {
  exp(predict(rf_mod_log[[i]], newdata  = test_list[[i]], type = "response"))
} # 2.752589 mins
stopCluster(cl)
end_time <- Sys.time()
time_taken <- end_time - start_time
print(time_taken) # 2.884095 mins


ncpus = 12
cl <- makeCluster(8) 
registerDoParallel(cl) 

start_time <- Sys.time()
prediction_train <- foreach(i = 1:100, .packages = c("randomForest", "data.table", "dplyr")) %dopar% {
  predict(rf_mod_hyper[[i]], newdata  = train_list[[i]], type = "response")
} # 2.752589 mins
stopCluster(cl)
end_time <- Sys.time()
time_taken <- end_time - start_time
print(time_taken) # 2.884095 mins

prediction_train <- lapply(1:100, function(i) {
  predict(rf_mod_hyper[[i]], newdata = train_list[[i]], type = "response")
})


start_time <- Sys.time()
prediction_log_train_hyper <- foreach(i = 1:100, .packages = c("randomForest", "data.table", "dplyr")) %dopar% {
  exp(predict(rf_mod_hyper[[i]], newdata  = train_list[[i]], type = "response"))
} # 2.752589 mins

prediction_log_train_hyper <- lapply(1:100, function(i) {
  exp(predict(rf_mod_hyper[[i]], newdata  = train_list[[i]], type = "response"))
})
stopCluster(cl)
end_time <- Sys.time()
time_taken <- end_time - start_time
print(time_taken) # 2.884095 mins

start_time <- Sys.time()
prediction_log_train <- foreach(i = 1:100, .packages = c("randomForest", "data.table", "dplyr")) %dopar% {
  exp(predict(rf_mod_log[[i]], newdata  = train_list[[i]], type = "response"))
} # 2.752589 mins
prediction_log_train <- lapply(1:100, function(i) {
  exp(predict(rf_mod_log[[i]], newdata  = train_list[[i]], type = "response"))
})

stopCluster(cl)
end_time <- Sys.time()
time_taken <- end_time - start_time
print(time_taken) # 2.884095 mins

## gam model prediction
ncpus = 18
cl <- makeCluster(ncpus) 
registerDoParallel(cl) 

start_time <- Sys.time()
prediction_gam_train <- foreach(i = 1:100, .packages = c("data.table", "dplyr")) %dopar% {
  predict(gam_mod[[i]], newdata  = train_list[[i]], type = "response")
} # 2.752589 mins
stopCluster(cl)
end_time <- Sys.time()
time_taken <- end_time - start_time
print(time_taken) # 2.884095 mins

ncpus = 18
cl <- makeCluster(ncpus) 
registerDoParallel(cl) 

start_time <- Sys.time()
prediction_gam_test <- foreach(i = 1:100, .packages = c("data.table", "dplyr")) %dopar% {
  predict(gam_mod[[i]], newdata  = test_list[[i]], type = "response")
} # 2.752589 mins
stopCluster(cl)
end_time <- Sys.time()
time_taken <- end_time - start_time
print(time_taken) # 2.884095 mins

#########
ncpus = 12
cl <- makeCluster(ncpus) 
registerDoParallel(cl) 

start_time <- Sys.time()
prediction_log <- foreach(i = 1:100, .packages = c("randomForest", "data.table", "dplyr")) %dopar% {
  exp(predict(rf_mod_log[[i]], newdata  = test_list[[i]], type = "response"))
} # 2.752589 mins
stopCluster(cl)
end_time <- Sys.time()
time_taken <- end_time - start_time
print(time_taken) # 2.884095 mins

# 7. Make predictions using RF models trained for each of 100 FOI dataset-------
ncpus = 12
cl <- makeCluster(ncpus) 
registerDoParallel(cl) 

start_time <- Sys.time()
ensemble_list_hyper_raw <- foreach(i = 1:100, .packages = c("randomForest", "data.table", "dplyr")) %dopar% {
  predict(rf_mod_hyper[[i]], newdata  = pred.data, type = "response")
}
stopCluster(cl)
end_time <- Sys.time()
time_taken <- end_time - start_time
print(time_taken) # 52.13149 mins

# chunk processing
split_data_into_chunks <- function(data, chunk_size) {
  split(data, ceiling(seq_along(1:nrow(data)) / chunk_size))
}
chunk_size <- 10000
chunks <- split_data_into_chunks(pred.data, chunk_size)

ensemble_list_hyper_raw <- list()
start_time <- Sys.time()

for (chunk_index in seq_along(chunks)) {
  chunk <- chunks[[chunk_index]]
  
  chunk_predictions <- lapply(1:100, function(i) {
    predict(rf_mod_hyper[[i]], newdata = chunk, type = "response")
  })
  
  # Combine predictions from all models for the current chunk
  ensemble_list_hyper_raw[[chunk_index]] <- do.call(cbind, chunk_predictions)
}

end_time <- Sys.time()
time_taken <- end_time - start_time
print(time_taken)

## raw processing
start_time <- Sys.time()

ensemble_list_hyper_raw <- lapply(1:length(rf_mod_hyper), function(i) {
  predict(rf_mod_hyper[[i]], newdata = pred.data, type = "response")
})

end_time <- Sys.time()
time_taken <- end_time - start_time
print(time_taken)


######
ncpus = 12  
cl <- makeCluster(ncpus) 
registerDoParallel(cl) 

start_time <- Sys.time()
ensemble_list_logit <- foreach(i = 1:100, .packages = c("randomForest", "data.table", "dplyr")) %dopar% {
  logistic(predict(rf_mod_logit[[i]], newdata  = pred.data, type = "response"))
}
stopCluster(cl)
end_time <- Sys.time()
time_taken <- end_time - start_time
print(time_taken) # 52.13149 mins


# 8. prediction dataframe ------------------------------------------------------

block_list <- list()

# 1. for log transformed
for(i in 1:length(block_data_list)) {
  # Add the Type and Predicted columns to the test and train data
  test_data <- data.frame(
    BlockID = test_list[[i]]$blockID, 
    FOI = exp(test_list[[i]]$logfoi), # convert back the trained logfoi to 0-1 scale 
    logFOI = test_list[[i]]$logfoi, # originally trained with logfoi
    #Predicted = prediction_test[[i]], # exponentiate back to 0-1 scale
    Predicted_log = prediction_log_test[[i]], # predict with log scale
    Type = "Test",
    x  = test_list[[i]]$x
  )
  
  train_data <- data.frame(
    BlockID = train_list[[i]]$blockID, 
    FOI = exp(train_list[[i]]$logfoi),
    logFOI = train_list[[i]]$logfoi,
    #Predicted = prediction_train[[i]],
    Predicted_log = prediction_log_train[[i]],
    Type = "Train",
    x = train_list[[i]]$x
  )
  
  # Combine test and train data
  combined_data <- rbind(test_data, train_data)
  combined_data <- combined_data %>%
    mutate(RMSE = sqrt((FOI - Predicted_log)^2))
  combined_data <- st_as_sf(combined_data)
  # Store the combined data in block_list
  block_list[[i]] <- combined_data
  
}

# 2. for linear scale
for(i in 1:length(block_data_list)) {
  # Add the Type and Predicted columns to the test and train data
  test_data <- data.frame(
    BlockID = test_list[[i]]$blockID, 
    FOI = test_list[[i]]$foi,
    #logFOI = test_list[[i]]$logfoi, # originally trained with logfoi
    Predicted = prediction_test[[i]], # exponentiate back to 0-1 scale
    #Predicted_log = prediction_log_test[[i]], # predict with log scale
    Type = "Test",
    x  = test_list[[i]]$x
  )

  train_data <- data.frame(
    BlockID = train_list[[i]]$blockID, 
    FOI = train_list[[i]]$foi,
    #logFOI = train_list[[i]]$logfoi,
    Predicted = prediction_train[[i]],
    #Predicted_log = prediction_log_train[[i]],
    Type = "Train",
    x = train_list[[i]]$x
  )
  
  # Combine test and train data
  combined_data <- rbind(test_data, train_data)
  combined_data <- combined_data %>% group_by(BlockID) %>%
    mutate(RMSE = sqrt(mean((FOI - Predicted)^2))) %>%  # RMSE computed per BlockID
    ungroup()
  combined_data <- st_as_sf(combined_data)
  # Store the combined data in block_list
  block_list[[i]] <- combined_data
  
}

# 3. log scale hyperparam
for(i in 1:length(block_data_list)) {
  # Add the Type and Predicted columns to the test and train data
  test_data <- data.frame(
    BlockID = test_list[[i]]$blockID, 
    FOI = exp(test_list[[i]]$logfoi), # convert back the trained logfoi to 0-1 scale 
    logFOI = test_list[[i]]$logfoi, # originally trained with logfoi
    #Predicted = prediction_test[[i]], # exponentiate back to 0-1 scale
    Predicted_log = prediction_log_test_hyper[[i]], # predict with log scale
    Type = "Test",
    x  = test_list[[i]]$x
  )
  
  train_data <- data.frame(
    BlockID = train_list[[i]]$blockID, 
    FOI = exp(train_list[[i]]$logfoi),
    logFOI = train_list[[i]]$logfoi,
    #Predicted = prediction_train[[i]],
    Predicted_log = prediction_log_train_hyper[[i]],
    Type = "Train",
    x = train_list[[i]]$x
  )
  
  # Combine test and train data
  combined_data <- rbind(test_data, train_data)
  combined_data <- combined_data %>%
    mutate(RMSE = sqrt((FOI - Predicted_log)^2))
  combined_data <- st_as_sf(combined_data)
  # Store the combined data in block_list
  block_list[[i]] <- combined_data
  
}


### gam model block list
block_list_gam <- list()

for(i in 1:length(block_data_list)) {
  # Add the Type and Predicted columns to the test and train data
  test_data <- data.frame(
    BlockID = test_list[[i]]$blockID, 
    FOI = test_list[[i]]$foi,
    Predicted = prediction_gam_test[[i]],
    Type = "Test",
    x  = test_list[[i]]$x
  )
  
  train_data <- data.frame(
    BlockID = train_list[[i]]$blockID, 
    FOI = train_list[[i]]$foi,
    Predicted = prediction_gam_train[[i]],
    Type = "Train",
    x = train_list[[i]]$x
  )
  
  # Combine test and train data
  combined_data <- rbind(test_data, train_data)
  combined_data <- combined_data %>%
    mutate(RMSE = sqrt((FOI - Predicted)^2))
  combined_data <- st_as_sf(combined_data)
  # Store the combined data in block_list
  block_list_gam[[i]] <- combined_data
  
}

# log scale
for(i in 1:100) {
  avg_rmse_list[[i]] <- block_list[[i]] %>%
    group_by(BlockID, Type) %>%
    summarise(
      model_id = i,
      meanRMSE = sqrt(mean((FOI - Predicted_log)^2, na.rm = TRUE)),
      .groups = 'drop'
    )
}

avg_rmse_list <- list()
# for linear scale
for(i in 1:100) {
  avg_rmse_list[[i]] <- block_list[[i]] %>%
    group_by(BlockID, Type) %>%
    summarise(
      model_id = i,
      meanRMSE = sqrt(mean((FOI - Predicted)^2, na.rm = TRUE)),
      .groups = 'drop'
    )
}

avg_rmse_df <- do.call(rbind, avg_rmse_list)
mean_rmse_allblock <- sapply(avg_rmse_list, function(df) mean(df$meanRMSE, na.rm = TRUE))
(q_allblock <- quantile(mean_rmse_allblock, c(0.025, 0.5, 0.975)))
mean_rmse_by_type <- sapply(avg_rmse_list, function(df) {
  tapply(df$meanRMSE, df$Type, mean, na.rm = TRUE)
})

# out sample rmse
(q_oos <- quantile(mean_rmse_by_type[1,], c(0.025, 0.5, 0.975)))
# in sample rmse
(q_is <- quantile(mean_rmse_by_type[2,], c(0.025, 0.5, 0.975)))

med <- quantile(mean_rmse_allblock, 0.5)
lower_bound <- quantile(mean_rmse_allblock, 0.025)
upper_bound <- quantile(mean_rmse_allblock, 0.975)
(ui_95 <- c(med, lower_bound, upper_bound))

# average rmse by block for all 100 models 
avg_by_block <- avg_rmse_df %>% group_by(BlockID) %>% 
  summarise(
  meanRMSE = mean(meanRMSE, na.rm = TRUE),
  .groups = 'drop'
)


## RMSE summary graph
# Combine all the RMSE data from block_list into one data frame
all_rmse_data <- bind_rows(lapply(1:length(block_list), function(i) {
  block_list[[i]] %>% 
    group_by(BlockID) %>% 
    summarise(meanRMSE = mean(RMSE),  # Summarize RMSE (using mean)
              geometry = first(geometry))%>%
    select(BlockID, RMSE) %>%  # Only keep BlockID and RMSE columns
    mutate(Model = i)  # Add a column to identify which model the data came from
}), .id = "model")

# Summarize RMSE for each BlockID across all models
rmse_summary <- all_rmse_data %>%
  group_by(BlockID) %>%
  summarise(meanRMSE = mean(RMSE),
            medianRMSE = median(RMSE),
            lowerCI = quantile(RMSE, probs = 0.025),
            upperCI = quantile(RMSE, probs = 0.975),
            Uncertainty_RMSE = sd(RMSE))

# Merge `rmse_summary` with your spatial data for plotting (if needed)
# Assuming you have a spatial object `block_sf` that contains geometry for each block
block_sf <- st_as_sf(avg_rmse_list[[1]], coords = c("Longitude", "Latitude"))

# Merge the spatial data with the RMSE summary statistics
plot_data <- st_join(avg_rmse_list[[1]], rmse_summary, by = "BlockID")

# Now, create a plot summarizing the RMSE for each block
p <- ggplot(data = plot_data) +
  geom_sf(aes(fill = Uncertainty_RMSE)) +  # Use meanRMSE for the color fill
  scale_fill_viridis_c() +
  geom_sf(data = p_covs_sf, color = 'red', size = 0.4) +  # Overlay other spatial data if needed
  labs(fill = "Standard Deviation of RMSE ") +
  theme_minimal() +
  theme(legend.position = "right")
ggsave(filename = paste0("Results_figs/blockCV_",
                         gsub("-", "_", Sys.Date()),
                         ".jpg"), p, height=6, width=12, dpi=900)

# 9. plots ---------------------------------------------------------------------
p_covs_sf <- st_as_sf(p_covs, coords = c("Longitude", "Latitude"))

ggplot(data = avg_by_block)+
  geom_sf(aes(fill = meanRMSE)) +
  scale_fill_viridis_c()+
  geom_sf(data = p_covs_sf, color = 'red', size = 2)+
  geom_sf_text(data = blocks_with_data, aes(label = blockID), size = 3, check_overlap = TRUE)+
  labs(title = "Spatial Blocks with Mean RMSE; RF", fill = "Mean RMSE") +
  theme_minimal() +
  theme(legend.position = "right")


pdf("combined_rmse_avg.pdf", width = 20, height = 8)

for(j in seq(1, length(rf_mod_hyper), by = 9)) {
  plot_list <- list()
  
  # Loop through each chunk of 9 models and create ggplot objects
  for(i in j:min(j+8, length(rf_mod_hyper))) {
    p <- ggplot(data = avg_rmse_list[[i]])+
      geom_sf(aes(fill = meanRMSE)) +
      scale_fill_viridis_c()+
      geom_sf(data = p_covs_sf, color = 'red', size = 0.3)+
      geom_sf_text(data = block_list[[i]], aes(label = ifelse(Type == "Test", Type, NA)), size = 3, check_overlap = TRUE)+
      labs(title =  paste("Spatial Blocks with RMSE per block: Graph", i), fill = "Mean RMSE") +
      theme_minimal() +
      theme(legend.position = "right")
    
    
    plot_list[[i - j + 1]] <- p
  }
  
  # Arrange the plots in a 3x3 grid
  grid.arrange(grobs = plot_list, ncol = 3, nrow = 3)
}

dev.off()


x_lim <- c(0, 0.15)
y_lim <- c(0, 0.15)

# logitFOI -> original transformed
ggplot(combined_df_rf, aes(x = foi, y = Predicted_FOI, color = factor(test_ID))) +
  geom_point() +  # Add points
  geom_point(alpha = 0.5) +
  labs(x = "Actual", y = "Fitted", title = "FOI, RF") +  # Labels and title
  geom_abline(slope = 1, linetype = "dashed", color = "red") +
  coord_fixed(ratio = 1) +  # Equal aspect ratio
  xlim(y_lim) +
  ylim(y_lim) +
  theme_bw()

plot_list <- list()

for (i in 1:length(combined_results_rf)) {
  # Create each ggplot object and store in the list
  p <- ggplot(combined_results_rf[[i]], aes(x = FOI.x, y = Predicted_FOI)) +
    geom_point() +
    labs(x = "Actual", y = "Fitted", title = paste("FOI, RF", i)) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    coord_fixed(ratio = 1) +
    xlim(y_lim) +
    ylim(y_lim) +
    theme_bw()
  plot_list[[i]] <- p
}

pdf("combined_plots.pdf", width = 8, height = 8)
# Combine and print or save the plots in batches of 9
for (j in seq(1, length(plot_list), by = 9)) {
  # Combine 9 plots or the remainder if less than 9
  combined_plot <- wrap_plots(plot_list[j:min(j+8, length(plot_list))], ncol = 3)
  
  # If you want to display the plot in an R environment
  print(combined_plot)
}

dev.off()


# each of 100 logit FOI model fitted to the entire data (logit scale)
pdf("combined_rf_alldata_raw.pdf", width = 8, height = 8)
par(mar = c(2, 2, 2, 2))

for(j in seq(1, length(rf_mod), by = 10)) {
  
  # Set up a multi-plot layout: 3 rows, 3 column
  par(mfrow = c(3, 3))
  
  # Loop through each chunk of 10 models
  for(i in j:min(j+9, length(rf_mod))) {
    
    # Predict and plot
    
    # Identify training data points (you will need to adjust this part)
    in_training <- merge_list[[i]]$ID %in% train_list[[i]]$ID
    
    p.rf <- predict(rf_mod[[i]], newdata=merge_list[[i]], type = "response")
    plot(merge_list[[i]]$foi[in_training], p.rf[in_training], asp=1, pch=20, xlab="actual", ylab="fitted", 
         main=paste("FOI, Random Forest", i),
         xaxs="i", yaxs="i", xlim=c(0,0.08),  ylim=c(0,0.08))
    points(merge_list[[i]]$foi[!in_training], p.rf[!in_training], col='red', pch=20)
    grid()
    abline(h=0, v=0, col="gray")
    abline(a=0, b=1)
    
    #text(p.rf, train_list[[i]]$FOI.x, labels=p_covs$ID, cex=0.7, pos=4)
  }
}

dev.off()

ggplot(data = block_list[[12]])+
  geom_point(aes(x = FOI, y = Predicted, color = Type, alpha = Type))+
  scale_color_manual(values = c("Test" = "Red", "Train" = "Black"))+
  scale_alpha_manual(values = c("Train" = 0.1, "Test" = 1)) +
  geom_abline(slope = 1)+
  xlim(c(0,0.3))+
  ylim(c(0,0.3))+
  theme_bw()

## in-out sample fitting
plot_list <- list()

for (i in 1:length(block_list)) {
  # Create each ggplot object and store in the list
  p <- ggplot(data = block_list[[i]])+
    geom_point(aes(x = FOI, y = Predicted, color = Type, alpha = Type))+
    scale_color_manual(values = c("Test" = "Red", "Train" = "Black"))+
    scale_alpha_manual(values = c("Train" = 0.2, "Test" = 1)) +
    geom_abline(slope = 1)+
    #xlim(c(0,0.25))+
    #ylim(c(0,0.25))+
    theme_bw()
  
  plot_list[[i]] <- p
}

pdf("Results_figs/combined_fitting_hyper_raw_1010.pdf", width = 12, height = 8)
# Combine and print or save the plots in batches of 9
for (j in seq(1, length(plot_list), by = 9)) {
  # Combine 9 plots or the remainder if less than 9
  combined_plot <- wrap_plots(plot_list[j:min(j+8, length(plot_list))], ncol = 3)
  
  # If you want to display the plot in an R environment
  print(combined_plot)
}

dev.off()

# 10. Ensemble map -------------------------------------------------------------
template = tsuit

mapRasCon <- function(pred_vec) { 
  output <- template
  values(output) = as.numeric(pred_vec)
  return(output)
}

ensemble_raster <- list()
for(i in 1:length(ensemble_list_hyper_raw)){
  mapresult          <- mapRasCon(ensemble_list_hyper_raw[[i]])
  ensemble_raster[[i]] <- mapresult 
}

bootstrap_stack <- stack(ensemble_raster)
avg_boogstrap_map <- calc(bootstrap_stack, fun = median)
lo_bootstrap_map <- calc(bootstrap_stack, fun = function(x) quantile(x, probs = 0.025, na.rm = TRUE))
hi_bootstrap_map <- calc(bootstrap_stack, fun = function(x) quantile(x, probs = 0.975, na.rm = TRUE))

save(avg_boogstrap_map, file = "avg_bootstrap_unmask_raw_0415.RData")
save(lo_boogstrap_map, file = "lo_bootstrap_unmask_raw_0415.RData")
save(hi_boogstrap_map, file = "hi_bootstrap_unmask_raw_0415.RData")

boot_CHIK_rf <- plotRaster(avg_boogstrap_map)+scale_fill_viridis(option = "rocket", direction = -1)
ggsave(filename = paste0("Map_figs/CHIK_boot_rf_raw",
                         gsub("-", "_", Sys.Date()),
                         ".jpg"), boot_CHIK_rf, height=6, width=12, dpi=900)


bootstrapRasterMid <- as.data.frame(rasterToPoints(avg_boogstrap_map, spatial = T))
bootstrapRasterLo <- as.data.frame(rasterToPoints(lo_boogstrap_map, spatial = T))
bootstrapRasterHi <- as.data.frame(rasterToPoints(hi_boogstrap_map, spatial = T))

save(bootstrapRasterRawMid, file = "bootRaster_raw_0415.RData")
save(bootstrapRasterRawLo, file = "bootRaster_raw_0415_lo.RData")
save(bootstrapRasterRawHi, file = "bootRaster_raw_0415_hi.RData")

save(ensemble_list_hyper_raw, file = "ensemble_list_hyper_raw_0415.RData")
save(ensemble_list_hyper_raw, file = "ensemble_list_hyper_raw_0623.RData") # without CHIK Risk
save(ensemble_list_hyper_raw, file = "ensemble_list_hyper_raw_0625.RData") # updated CHIK Risk

## masking 
chik_binary <- raster("CHIK_binmap_2024_04_05.tif")

chik_mask <- chik_binary * avg_boogstrap_map 

boot_CHIK_rf_mask <- plotRaster(chik_mask)+scale_fill_viridis(option = "rocket", direction = -1)
ggsave(filename = paste0("Map_figs/CHIK_boot_rf_mask_raw",
                         gsub("-", "_", Sys.Date()),
                         ".jpg"), boot_CHIK_rf_mask, height=6, width=12, dpi=900)

bootstrapRasterMask <- as.data.frame(rasterToPoints(chik_mask, spatial = T))

save(bootstrapRasterMask, file = "bootRasterMask_raw_0415.RData")
save(chik_mask, file = "chik_mask_raw_0415.RData")


# 11. var Imp plots ------------------------------------------------------------
# Open a PDF device to store the plots
pdf("combined_rf_varImp.pdf", width = 8, height = 8)

# Loop through the models in chunks
for(j in seq(1, length(rf_mod_hyper), by = 9)) {
  
  # Set up a 3x3 multi-plot layout and adjust margins
  par(mfrow = c(3, 3), mar = c(2, 2, 2, 2))
  
  # Inner loop for plotting each model in the current chunk
  for(i in j:min(j+8, length(rf_mod_hyper))) {
    
    # Generate and plot the variable importance for each model
    varImpPlot(rf_mod_hyper[[i]])
    
  }
}

# Close the PDF device after all plots are done
dev.off()


# pdps
# Assuming plotList is populated with PDPs as per your existing code

# Open a PDF device to store the plots
# if error persists in importance: detach("package:randomForest", unload=TRUE) and then reload the randomforest library

imp <- list()
impvar <- list()
for(i in 1:length(rf_mod_hyper)){
  result <- importance(rf_mod_hyper[[i]])
  imp[[i]] <- result
  impvar[[i]] <-rownames(result)[order(result[,1],decreasing=TRUE)]
  
}

merge_dat <- as.data.frame(merge_list[[1]])
#correct_levels <- levels(train_list[[1]]$DHI)

# Align the levels in the new data
#merge_dat$DHI <- factor(merge_dat$DHI, levels = correct_levels)

plotList <- list()
for(i in seq_along(rf_mod)){
  plotList[[i]] <- list()
  for(k in seq_along(impvar[[i]])){ 
    pdp <- partialPlot(rf_mod_hyper[[i]], merge_dat, impvar[[i]][k],
                                 xlab=impvar[[i]][k], 
                                 main=paste("Partial Dependence on", impvar[[i]][k]),
                                 plot = FALSE)
    pdp$y <- exp(pdp$y)
    plotList[[i]][[k]] <- pdp
  }
}

plotList <- list()

for(i in seq_along(rf_mod_hyper)){
  # Copy merge_dat for this iteration to avoid altering original data
  #aligned_merge_dat <- merge_dat
  
  # Align the factor levels for each factor variable used in the model
  # Assuming you know the factor variables and have the correct levels for each
  #factor_variables <- c("DHI")  # Update this list as needed
  #for(fv in factor_variables) {
  #  correct_levels <- levels(train_list[[i]][[fv]])  # Assuming train_list[[i]] aligns with rf_mod[[i]]
  #  aligned_merge_dat[[fv]] <- factor(aligned_merge_dat[[fv]], levels = correct_levels)
  #}
  
  plotList[[i]] <- list()
  for(k in seq_along(impvar[[i]])){ 
    plotList[[i]][[k]] <- partialPlot(rf_mod_hyper[[i]], merge_dat, impvar[[i]][k],
                                      xlab=impvar[[i]][k], 
                                      main=paste("Partial Dependence on", impvar[[i]][k]),
                                      plot = FALSE)
  }
}

pdf("Results_figs/combined_rf_pdp_1010.pdf", width = 8, height = 8)

plotNumber <- 1  # Initialize a counter for plot numbers

# Loop through all models
for (k in seq_along(rf_mod_hyper)) {
  par(mfrow = c(3, 3), mar = c(4, 4, 4, 4))
  
  # Loop through all variables for each model
  for (i in seq_along(plotList[[k]])) {
    # Plot the current PDP
    plot(plotList[[k]][[i]], type = "l",
         xlab = impvar[[k]][i],
         ylab = "Relative contribution to FOI")
    abline(h = 0, lty = 2)
    
    plotNumber <- plotNumber + 1  # Increment the plot counter
  }
}

# Close the PDF device after all plots are done
dev.off()

# varimp dataframe -------------------------------------------------------------
# if error persists in importance: detach("package:randomForest", unload=TRUE)

varimp <- list()
varexp <- list()
for(i in 1:length(rf_mod_hyper)){
  varimp[[i]] <- importance(rf_mod_hyper[[i]])
  
  # Extract % var explained
  varexp[[i]] <- tail(rf_mod_hyper[[i]]$rsq, 1) * 100
}

varimp <- do.call(rbind, varimp)
varexp <- do.call(rbind, varexp)

covariate <- rep(c("Tsuit", "PRCP","GDP", "Albo", "Aegyp", "CHIKRisk"), times = 100)
model <- rep(c(1:100), each = 6)
combined_imp <- cbind(varimp, covariate, model)
colnames(combined_imp) <- c("incMSE", "incnodepurity", "cov", "model")
combined_imp <- as.data.frame(combined_imp)
combined_imp$incMSE <- as.numeric(combined_imp$incMSE)
combined_imp$incnodepurity <- as.numeric(combined_imp$incnodepurity)
combined_imp <- combined_imp %>% group_by(model) %>% 
  mutate(normInc = incMSE/sum(incMSE))
combined_imp <- combined_imp %>% group_by(model) %>% 
  mutate(normpurity = incnodepurity/sum(incnodepurity))
combined_imp <- combined_imp %>% 
  group_by(model) %>% 
  mutate(normInc = (normInc - min(normInc)) / (max(normInc) - min(normInc)))

avg_imp <- combined_imp %>% group_by(cov)%>%
  summarise(mean_normIncMSE = mean(normInc),
            mean_normIncPurity = mean(normpurity))


write.csv(avg_imp, file = "MainData/avg_imp.csv", row.names = FALSE)
avg_imp$cov <- fct_reorder(avg_imp$cov, avg_imp$mean_normIncMSE, .desc = F)
avg_imp$cov <- fct_reorder(avg_imp$cov, avg_imp$mean_normIncPurity, .desc = F)

varimp1 <- ggplot(avg_imp) +
  geom_col(aes(x = cov, y = mean_normIncMSE), fill = "steelblue3", width = 0.6) +
  coord_flip()+
  labs(y = "Mean normalized % Increase in MSE", x = "Covariates")+
  geom_text(aes(x = cov, y = mean_normIncMSE, label = round(mean_normIncMSE, 2)), 
            color = "black", size = 6)+
  theme_bw()

ggsave(filename = paste0("Results_figs/varImp1_1010",
                         gsub("-", "_", Sys.Date()),
                         ".jpg"), varimp1, height=6, width=8, dpi=900)


varimp2 <- ggplot(avg_imp) +
  geom_col(aes(x = cov, y = mean_normIncPurity), fill = "steelblue3", width = 0.6) +
  coord_flip()+
  labs(y = "Mean normalized % Increase in Purity", x = "Covariates")+
  geom_text(aes(x = cov, y = mean_normIncPurity, label = round(mean_normIncPurity, 2)), 
            color = "black", size = 6)+
  theme_bw()

ggsave(filename = paste0("Results_figs/varImp2_1010",
                         gsub("-", "_", Sys.Date()),
                         ".jpg"), varimp2, height=6, width=8, dpi=900)

### 95%CI pdp
covariates <- c("Tsuit", "PRCP", "GDP", "Albo", "Aegyp", "CHIKRisk")

get_pdp <- function(model, data, var) {
  pd <- partial(model, pred.var = var, train = data, grid.resolution = 100, plot = FALSE)
  pd$covariate <- var  # Add covariate name to the result
  return(pd)
}

# Iterate over models and their corresponding datasets
pdp_list <- map2(rf_mod_hyper, train_data, ~ map(covariates, get_pdp, model = .x, data = .y))






