# Storing all the functions necessary for handling the data post-obtaining
# the chemical descriptors (in this case, post-PaDEL). 
# Contains the functions necessary for pre-processing and feature selection. 

# Packages ----------------------------------------------------------------

if(!require("pacman")) { 
  install.packages("pacman")
  library(pacman)
} else {
  library(pacman)
}

p_load(
  caret, 
  data.table, 
  Matrix, 
  randomForest,
  stringr, 
  tidyverse
  )

# Pre-processing ----------------------------------------------------------

# Finding activity outliers
# Because there's no "good way" to do this, and some QSARs are actually 
# fairly good at managing to capture outliers (> 2 SD away from mean),
# the threshold is allowed to be flexible
# data assumes DelG is in the second column
find_outlier <- function(data, thresh) {
  sd <- find_sd(data[ , 2])
  avg <- mean(data[ , 2])
  # There's definitely a better way to do this, but I'll deal
  # with a for loop for now
  outliers <- c()
  for (i in 1:nrow(data)) {
    if ((data[i, 2] > avg + thresh*sd || data[i, 2] < avg - thresh*sd))
      outliers <- c(outliers, i)
  }
  outliers
}

# centers and scales data
# assume the data only contains predictors
# also removes NAs and replaces them with 0 (this is equivalent to 
# replacing the values with mean after centering)
# returns list of centered and scaled data as well as settings
center_scale <- function(predictors_only) {
  pp_settings <- preProcess(
    predictors_only, 
    na.remove = T, 
    method = c("center", "scale"), 
    verbose = F
  )
  predictors_pp <- predict(pp_settings, predictors_only)
  predictors_pp[is.na(predictors_pp)] <- 0
  list(predictors_pp, pp_settings)
}

# Creating data partitions. Unlike folds, this allows data points
# to appear multiple times, allowing for testing on multiple sets.
# Partitioned on dG
split_train_test <- function(times, data) {
  partitions <- createDataPartition(
    data[, 2], 
    times = times, 
    p = 0.75
  )
  splits <- list()
  for(i in 1:times) {
    splits[[i]] <- list(
      tst = data[-partitions[[i]], ], 
      trn = data[partitions[[i]], ]
    )
  }
  splits
}

# this saves the resulting splits, unlike split_train_test
save_splits <- function(times, data, path) {
  if(!dir.exists(path)) dir.create(path)
  if(!substr(path, nchar(path), nchar(path)) == "/")
    path <- paste0(path, "/")
  partitions <- createDataPartition(
    data[ , 2], 
    times = times,
    p = 0.75
    )
  for(i in 1:times) {
    trn_index <- partitions[[i]]
    tst <- data[-trn_index, ]
    trn <- data[trn_index, ]
    
    tst_path <- paste0(path, i - 1, "-tst", ".RDS")
    trn_path <- paste0(path, i - 1, "-trn", ".RDS")
    saveRDS(tst, tst_path)
    saveRDS(trn, trn_path)
  }
}


# Feature selection -------------------------------------------------------

# trntst: list of trn and tst data
use_rfe <- function(trntst) {
  pred <- dplyr::select(trntst[["trn"]], -guest, -dG)
  obs  <- trntst[["trn"]]$dG
  ctrl <- rfeControl(
    functions = rfFuncs,
    method = "repeatedcv", 
    repeats = 5
    )
  # trying RFE with differing numbers of variables
  subsets <- c(1:5, 10, 15, 20, 25, 50)
  # returns feature selection profile
  rfe(
    x = pred, 
    y = obs,
    sizes = subsets, 
    rfeControl = ctrl
    )
}