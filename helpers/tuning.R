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


# SVM ---------------------------------------------------------------------

tune_svm <- function(x, y, folds, kernel, combo) {
  # create a list of the parameters
  svm_param <- append(
    list("x" = x, "y" = y, "kernel" = kernel),
    combo
  )
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

