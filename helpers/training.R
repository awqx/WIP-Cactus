# Libaries and packages ---------------------------------------------------

if(!require("pacman")) { 
  install.packages("pacman")
  library(pacman)
} else {
  library(pacman)
}

p_load(
  caret, 
  Cubist,
  e1071, kernlab, # svm
  glmnet,
  pls,
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

# Cubist ------------------------------------------------------------------

# important parameters are committees, extrapolation, rules, sample
# also include random seed
train_cubist <- function(folds, param, seed, fold_avg = T) {
  if (fold_avg) {
    fold_results <- do.call(
      rbind, 
      mapply(
        train_cubist_fold, 
        fold = folds, 
        name = names(folds),
        MoreArgs = list(parameters = param, seed = seed),
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
        train_cubist_fold, 
        fold = folds, 
        name = names(folds),
        MoreArgs = list(parameters = param, seed = seed),
        SIMPLIFY = F
      )
    )
  }
}

train_cubist_fold <- function(fold, name, parameters, seed) {
  ctrl <- cubistControl(
    seed = seed, 
    extrapolation = parameters$extrapolation,
    sample = parameters$sample
  )
  cubist_param <- list(
    "x" = fold$trn[,-1:-2],
    "y" = fold$trn$dG,
    "committees" = parameters$committees,
    "control" = ctrl
  )
  tst_x <- fold$tst[, -1:-2]
  tst_y <- fold$tst$dG
  
  cubist_model <- do.call(cubist, cubist_param)
  cubist_df <- data.frame(
    "pred" = predict(cubist_model, tst_x), 
    "obs" = tst_y
  )
  data.frame(
    parameters, 
    "fold" = name,
    "rmse" = defaultSummary(cubist_df)["RMSE"], 
    "r2" = defaultSummary(cubist_df)["Rsquared"]
  )
}

# PLS ---------------------------------------------------------------------

# important parameters are `ncomp` and `method`
# plsr requires a formula instead of x and y

train_pls <- function(folds, param, fold_avg = T) {
  if (fold_avg) {
    fold_results <- do.call(
      rbind, 
      mapply(
        train_pls_fold, 
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
        train_pls_fold, 
        fold = folds, 
        name = names(folds),
        MoreArgs = list(parameters = param),
        SIMPLIFY = F
      )
    )
  }
}

train_pls_fold <- function(fold, name, parameters) {
  pls_param_fold <- append(
    list(
      "formula" = dG~.,
      "data" = fold$trn[, -1]     
    ),
    parameters
  )
  tst_x <- fold$tst[, -1:-2]
  tst_y <- fold$tst$dG
  pls_model <- do.call(plsr, pls_param_fold)
  # unlike other models, w/ pls you have to 
  # assign the df names after predictions
  pls_df <- data.frame(
    predict(pls_model, tst_x), 
    tst_y
  )
  colnames(pls_df) <- c("pred", "obs")
  data.frame(
    parameters, 
    "fold" = name,
    "rmse" = defaultSummary(pls_df)["RMSE"], 
    "r2" = defaultSummary(pls_df)["Rsquared"]
  )
}

# GLMNet ------------------------------------------------------------------


