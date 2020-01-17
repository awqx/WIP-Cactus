# Libaries and packages ---------------------------------------------------

if(!require("pacman")) { 
  install.packages("pacman")
  library(pacman)
} else {
  library(pacman)
}

p_load(
  caret, 
  e1071, kernlab, # svm
  randomForest,
  stringr, 
  tidyverse
)

# Data handling -----------------------------------------------------------

# trn: df of train data. contains guest, dG, and descriptors (in order)
# nfolds: number of folds
# seed: random seed
# returns a list of lists. each contains "trn" and "tst"
split_train <- function(trn, nfolds, seed) {
  set.seed(seed)
  fold_list <- createFolds(trn$dG, k = nfolds)
  lapply(
    fold_list, 
    function(x, dat) list("trn" = dat[-x, ], "tst" = dat[x, ]), 
    dat = trn
  )
}

# SVM ---------------------------------------------------------------------

# folds: a list of lists. each list consists of training and testing data,
# labeled as "trn" and "tst", respectively
# param: a list of the parameters passed to the function. requires that
# the list names are the arguments for the model
# fold_avg: return the average rmse and r2 for the set of folds 
  # T for tuning, F for training on the splits
train_svm <- function(folds, kernel, param, fold_avg = T) {
  svm_param <- append(list("kernel" = kernel), param)
  if (fold_avg) {
    fold_results <- do.call(
      rbind, 
      mapply(
        train_svm_fold, 
        fold = folds, 
        name = names(folds),
        MoreArgs = list(parameters = svm_param),
        SIMPLIFY = F
      )
    )
    data.frame(
      svm_param, 
      "rmse" = mean(fold_results$rmse),
      "r2" = mean(fold_results$r2)
    )
  } else {
    do.call(
      rbind, 
      mapply(
        train_svm_fold, 
        fold = folds, 
        name = names(folds),
        MoreArgs = list(parameters = svm_param),
        SIMPLIFY = F
      )
    )
  }
} 

# a helper for train_svm. used in a lapply given a list of folds
# takes a fold, extracts x and y; builds svm model
# name is the fold number
train_svm_fold <- function(fold, name, parameters) {
  svm_param_fold <- append(
    list(
      "x" = fold$trn[, -1:-2],
      "y" = fold$trn$dG     
    ),
    parameters
  )
  tst_x <- fold$tst[, -1:-2]
  tst_y <- fold$tst$dG
  svm_model <- do.call(svm, svm_param_fold)
  svm_df <- data.frame(
    "pred" = predict(svm_model, tst_x), 
    "obs" = tst_y
  )
  data.frame(
    parameters, 
    "fold" = name,
    "rmse" = defaultSummary(svm_df)["RMSE"], 
    "r2" = defaultSummary(svm_df)["Rsquared"]
  )
}

# Random forest -----------------------------------------------------------

# important parameters are ntree, nodesize, and mtry
train_rf <- function(folds, param, fold_avg = T) {
  if (fold_avg) {
    fold_results <- do.call(
      rbind, 
      mapply(
        train_rf_fold, 
        fold = folds, 
        name = names(folds),
        MoreArgs = list(parameters = param),
        SIMPLIFY = F
      )
    )
    data.frame(
      param, 
      "rmse" = mean(fold_results$rmse),
      "r2" = mean(fold_results$r2)
    )
  } else {
    do.call(
      rbind, 
      mapply(
        train_rf_fold, 
        fold = folds, 
        name = names(folds),
        MoreArgs = list(parameters = param),
        SIMPLIFY = F
      )
    )
  }
}

train_rf_fold <- function(fold, name, parameters) {
  rf_param <- append(
    list(
      "x" = fold$trn[, -1:-2],
      "y" = fold$trn$dG     
    ),
    parameters
  )
  tst_x <- fold$tst[, -1:-2]
  tst_y <- fold$tst$dG
  rf_model <- do.call(randomForest, rf_param)
  rf_df <- data.frame(
    "pred" = predict(rf_model, tst_x), 
    "obs" = tst_y
  )
  data.frame(
    parameters, 
    "fold" = name,
    "rmse" = defaultSummary(rf_df)["RMSE"], 
    "r2" = defaultSummary(rf_df)["Rsquared"]
  )
}
