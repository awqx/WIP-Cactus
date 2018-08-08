# Storing all the functions necessary for handling the data post-obtaining
# the chemical descriptors (in this case, post-PaDEL). 
# Contains the functions necessary for pre-processing and feature selection. 


# Packages ----------------------------------------------------------------

if(!require("pacman")) { 
  install.packages("pacman")
  library(pacman)
} else
  library(pacman)

p_load(caret, data.table, Matrix, stringr, tidyverse)

# Pre-processing ----------------------------------------------------------

split.train.test <- function(times, data, info, path) {
  for(n in 1:times) {
    # Initialize a sample for reference for maxDissim
    data.init <- sample_n(data, 20)
    
    # Generate representative training set
    trn.ind <- maxDissim(a = data.init, b = data, 
                         n = round(nrow(data)*0.75), 
                         na.rm = T)
    tst.data <- data[-trn.ind, ]
    trn.data <- data[trn.ind, ]
    
    # Adding the guest molecule name and correct DelG
    tst <- cbind(info[-trn.ind, ], tst.data)
    trn <- cbind(info[trn.ind, ], trn.data)
    
    # Create logical file path name
    tst.path <- paste0(path, "tst", n, ".RDS")
    trn.path <- paste0(path, "trn", n, ".RDS")
    saveRDS(tst, tst.path)
    saveRDS(trn, trn.path)
    
    message("Split ", n, " completed")
  }
}

# Splitting data by folds instead of maxDissim

# data: data.frame with first column being guest, second column being DelG or
# the response 
# path: string of the target directory, ending with a backslash
split.train.test <- function(times, data, path) {
  folds <- createFolds(data[ , 2], k = times)
  for(i in 1:times) {
    tst.ind <- folds[[i]]
    tst <- data[tst.ind, ]
    trn <- data[-tst.ind, ]
    
    tst.path <- paste0(path, "tst", i, ".RDS")
    trn.path <- paste0(path, "trn", i, ".RDS")
    saveRDS(tst, tst.path)
    saveRDS(trn, trn.path)
  }
}

# Creating data partitions instead of folds. This allows data points
# to appear multiple times, allowing for testing/training on multiple
# sets
split.train.test <- function(times, data, path) {
  parts <- createDataPartition(data[ , 2], times = times, 
                               p = 0.75)
  for(i in 1:times) {
    trn.ind <- parts[[i]]
    tst <- data[-trn.ind, ]
    trn <- data[trn.ind, ]
    
    tst.path <- paste0(path, "tst", i, ".RDS")
    trn.path <- paste0(path, "trn", i, ".RDS")
    saveRDS(tst, tst.path)
    saveRDS(trn, trn.path)
  }
}

# Finding activity outliers
# Because there's no "good way" to do this, and some QSARs are actually 
# fairly good at managing to capture outliers (> 2 SD away from mean)
# I'm using a fairly generous 2.5 standard deviations away from the mean
# (or whatever you set the threshold to be)
# data assumes DelG is in the second column
find.activity.outlier <- function(data, thresh) {
  sd <- find.sd.desc(data[ , 2])
  avg <- mean(data[ , 2])
  # There's definitely a better way to do this, but I'll deal
  # with a for loop for now
  outliers <- c()
  for (i in 1:nrow(data)) {
    if ((data[i, 2] > avg + thresh*sd || data[i, 2] < avg - thresh*sd))
      outliers <- c(outliers, i)
  }
  return(outliers)
}
# Performing pre-processing on the files created by split.train.test
preprocess.splits <- function(filepath, writepath) {
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

