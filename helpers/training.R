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

tune_svm <- function(x, y, folds = 0, kernel, combo) {
  # create a list of the parameters
  svm_param <- append(
    list("x" = x, "y" = y, "kernel" = kernel),
    combo
  )

  # this `if`` statement allows for tuning the model when holding the
  # parameters constant and changing the folds of data (pass the folds
  # as a list in a `mapply` with a `do.call` and `rbind`)
  # this is implemented in the "Model" section of each QSAR script
  if (class(folds) == "numeric") {
    svm_model <- do.call(svm, svm_param)
    svm_df <- data.frame(
      "pred" = predict(svm_model, tst_x), 
      "obs" = tst_y
    )
    data.frame(
      svm_param[-1:-2], 
      "rmse" = mean(fold_results$rmse),
      "r2" = mean(fold_results$r2)
    )
  } else {
    fold_results <- do.call(
      rbind, 
      lapply(
        folds, 
        tune_svm_fold, 
        param = svm_param
      )
    )
    data.frame(
      svm_param[-1:-2],
      "rmse" = mean(fold_results$rmse),
      "r2" = mean(fold_results$r2)
    )    
  }
}

# take the values of a fold as test
# helper to tune_svm
tune_svm_fold <- function(fold, param) {
  # splitting into test and train sets based on the fold
  tst_x <- param[["x"]][fold, ]
  tst_y <- param[["y"]][fold]
  param[["x"]] <- param[["x"]][-fold, ]
  param[["y"]] <- param[["y"]][-fold]
  svm_model <- do.call(svm, param)
  svm_df <- data.frame(
    "pred" = predict(svm_model, tst_x), 
    "obs" = tst_y
  )
  data.frame(
    "rmse" = defaultSummary(svm_df)[1], 
    "r2" = defaultSummary(svm_df)[2]
  )
}

