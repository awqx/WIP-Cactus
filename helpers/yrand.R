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