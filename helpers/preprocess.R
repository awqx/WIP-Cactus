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

p_load(caret, data.table, Matrix, stringr, tidyverse)

# Pre-processing ----------------------------------------------------------

# Creating data partitions. Unlike folds, this allows data points
# to appear multiple times, allowing for testing on multiple sets.
split_train_test <- function(times, data, path) {
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

# split_train_test modified for dG as col1
tt_split <- function(times, data, path) {
  partitions <- createDataPartition(
    data[ , 1], 
    times = times, 
    p = 0.75
    )
  for(i in 1:times) {
    trn_index <- partitions[[i]]
    tst <- data[-trn_index, ]
    trn <- data[trn_index, ]
    
    tst_path <- paste0(path, "tst", i, ".RDS")
    trn_path <- paste0(path, "trn", i, ".RDS")
    saveRDS(tst, tst_path)
    saveRDS(trn, trn_path)
  }
}

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
# Performing pre-processing on the files created by split_train_test
preprocess_splits <- function(filepath, writepath) {
  files <- list.files(path = filepath) %>% 
    str_extract(., "trn[[:print:]]+") %>% 
    .[!is.na(.)]
  for(n in 1:length(files)) {
    dirpath <- paste0(writepath, n)
    dir.create(dirpath)
    data <- readRDS(paste0(filepath, files[n]))
    info <- data %>% dplyr::select(guest, DelG)
    desc <- data %>% dplyr::select(-guest, -DelG)
    
    # Replacing inconvenient values
    desc <- do.call(
      data.frame, 
      lapply(
        desc, 
        function(x) replace(x, is.infinite(x), NA)
        )
      )
    
    # Initial pre-processing
    pp_settings <- preProcess(
      desc, 
      na.remove = T,
      method = c("knnImpute", "center", "scale"), 
      verbose = F)
    desc_pp <- predict(pp_settings, desc)
    
    # Removing predictors with near-zero variance
    zero_pred <- nearZeroVar(desc_pp)
    desc_pp <- desc_pp[ , -zero_pred]
    
    # Removing highly correlated predictors
    desc_pp <- desc_pp[ , sapply(desc_pp, is_numeric)]
    high_cor <- findCorrelation(
      cor(
        desc_pp, 
        use = "pairwise.complete.obs"
        ),
      0.95 # 0.95 mostly arbitrary
      ) 
    desc_pp <- desc_pp[ , -high_cor]
    
    desc_pp <- cbind(info, desc_pp)
    saveRDS(pp_settings, paste0(dirpath, "/pp-settings.RDS"))
    saveRDS(desc_pp, paste0(dirpath, "/pp.RDS"))
    saveRDS(zero_pred, paste0(dirpath, "/zero-pred.RDS"))
    saveRDS(too_high, paste0(dirpath, "/high-cor.RDS"))
    
    message("Pre-processing of split ", n, " completed.")
  }
}

# preprocess_splits without the need for guests (assuming DelG is col1)
pp_split <- function(filepath, writepath) {
  files <- list.files(path = filepath) %>% 
    str_extract(., "trn[[:print:]]+") %>% 
    .[!is.na(.)]
  for(n in 1:length(files)) {
    dirpath <- paste0(writepath, n)
    dir.create(dirpath)
    data <- readRDS(paste0(filepath, files[n]))
    info <- data[ , 1]
    desc <- data[ , -1]
    
    # Replacing inconvenient values
    desc <- do.call(data.frame, lapply(desc, 
                                       function(x)
                                         replace(x, is.infinite(x), NA)))
    
    # Initial pre-processing
    pp.settings <- preProcess(desc, na.remove = T, 
                              method = c("knnImpute", "center", "scale"), 
                              verbose = F)
    desc.pp <- predict(pp.settings, desc)
    
    # Removing predictors with near-zero variance
    zero.pred <- nearZeroVar(desc.pp)
    desc.pp <- desc.pp[ , -zero.pred]
    
    # Removing highly correlated predictors
    desc.pp <- desc.pp[ , sapply(desc.pp, is.numeric)]
    too.high <- findCorrelation(cor(desc.pp, use = "pairwise.complete.obs"), 0.95) # 0.95 mostly arbitrary
    desc.pp <- desc.pp[ , -too.high]
    
    desc.pp <- cbind(info, desc.pp)
    saveRDS(pp.settings, paste0(dirpath, "/pp.settings.RDS"))
    saveRDS(desc.pp, paste0(dirpath, "/pp.RDS"))
    saveRDS(zero.pred, paste0(dirpath, "/zero.pred.RDS"))
    saveRDS(too.high, paste0(dirpath, "/high.cor.RDS"))
    
    message("Pre-processing of split ", n, " completed.")
  }
}

# So far, only used for y-random
ev.split <- function(data, writepath) {
  ev <- sample_frac(data, size = 0.15)
  md <- data[!row.names(data) %in% row.names(ev), ]
  saveRDS(ev, paste0(writepath, '/extval.RDS'))
  saveRDS(md, paste0(writepath, '/model.data.RDS'))
  return(md)
}

# Feature selection -------------------------------------------------------

# path: location of the pre-processed training dataset
use.rfe <- function(path) {
  trn <- readRDS(path)
  pred <- trn %>% dplyr::select(., -guest, -DelG)
  obs <- trn$DelG
  
  ctrl <- rfeControl(functions = rfFuncs, 
                     method = "repeatedcv", 
                     # verbose = T, # uncomment to monitor
                     repeats = 5)
  # subsets <- c(5, 10, 15, 20, 35, 50, 75, 100)
  subsets <- c(1:5, 10, 15, 20, 25, 50)
  
  rfe.profile <- rfe(x = pred, y = obs, 
                     sizes = subsets, rfeControl = ctrl)
  
  return(rfe.profile)
}

# An RFE for yrand data
use.rfe.yrand <- function(path) {
  trn <- readRDS(path)
  # For some reason, instead of DelG, dG info is listed as 'info'
  pred <- trn %>% dplyr::select(., -info)
  obs <- trn$info
  
  ctrl <- rfeControl(functions = rfFuncs, 
                     method = "repeatedcv", 
                     # verbose = T, # uncomment to monitor
                     repeats = 5)
  # Slightly smaller for time efficency
  subsets <- c(1, 5, 10, 20, 50)
  
  rfe.profile <- rfe(x = pred, y = obs, 
                     sizes = subsets, rfeControl = ctrl)
  
  return(rfe.profile)
}

# pp.path: folder where pre-processed files reside
#     Here, it should be "./pre-process
# write.path: target folder, no backslash at the end
#     "./feature.selection
# cd.type: string ("alpha" or "beta")
# num: integer referring to the trn-tst split
save.rfe <- function(pp.dir, write.dir, cd.type, num) {
  
  info.path <- paste0(pp.dir, "/", cd.type, "/", num, "/pp.RDS")
  rfe.profile <- use.rfe(info.path)
  rfe.name <- paste0("rfe", num, ".RDS")
  target.path <- paste0(write.dir, "/", cd.type, "/", rfe.name)
  saveRDS(rfe.profile, target.path)
  message("RFE of test-train split ", num, " completed.")
  
}

# var.name: strings for variable names
# file.name: strings for filepaths
read.rfe <- function(var.names, file.names) {
  for(n in 1:length(var.names)) {
    temp <- readRDS(file.names[n])
    assign(var.names[n], temp, envir = .GlobalEnv)
  }
}

# A second version of save.rfe for yrand analysis
# Only saves a copy of the important variables
# host.dir: the directory that contains all the files for the cyclodextrin host
    # should end in backslash
# num.trial: the number of repetitions of yrand
# num.split: the number of test splits
# write.dir
rfe.yrand <- function(host.dir, num.trial, num.split) {
  # Create the RFE profiles
  # This has to be looped through the number of yrand trials
  for (i in 1:num.trial) {
    message('Starting RFE of yrand trial ', i)
    # Creates the combinations of filepaths where pre-processed data can be found
    pp.paths <- paste0(host.dir, i, '/pp/', 1:num.split, '/pp.RDS')
    # List of the RFE objects
    rfe.list <- lapply(pp.paths, use.rfe.yrand)
    # Finding the most used variables (same as in 06.feature.selection)
    pred <- unlist(lapply(FUN = predictors, rfe.list))
    pred.uq <- unique(pred)
    pred.pattern <- paste0("^", pred.uq) 
    pred.pattern <- paste0(pred.pattern, "$")
    count <- lapply(FUN = str_count, X = pred.pattern, string = pred) %>%
      lapply(FUN = sum, X = .) %>% unlist()
    varimp <- data.frame(pred.uq, count) %>%
      rename(predictor = pred.uq, frequency = count) %>%
      mutate(predictor = as.character(predictor)) %>%
      .[order(.$frequency, decreasing = T), ]
    # print(head(varimp))
    vars <- varimp %>% filter(frequency == max(varimp$frequency)) %>%
      .$predictor
    # print(vars)
    if(length(vars) < 2) 
      vars <- varimp %>% filter(frequency >= (max(varimp$frequency) - 1)) %>%
      .$predictor
    saveRDS(vars, paste0(host.dir, i, '/vars.RDS'))
    message('RFE of yrand trial ', i, ' completed.')
  }
}



