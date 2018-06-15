# Libraries and Packages --------------------------------------------------

library(caret)
library(e1071)
library(kernlab)
library(stringr)
library(tidyverse)

# Functions ---------------------------------------------------------------

# The cost parameter is related to the complexity of the SVM
# Large cost = flexible model responsive to outliers
# Small cost = conservative model less likely to overfit
# 
# kerneltype: linear, polynomial, radial basis, and sigmoid
# nfolds: # of CV folds, 10 is recommended
# seed: random seed
# data: dataframe with "guest" removed (only numerics)

tune.svm.cost <- function(data, nfolds, cost, seed, kerneltype) {
  set.seed(seed)
  fold.list <- createFolds(y = data[ , 1], k = nfolds)
  r2.results <- c(rep(0.0, nfolds))
  rmse.results <- c(rep(0.0, nfolds))
  
  # Cross-validation using k-fold cross-validation
  for (i in 1:nfolds)
  {
    fold <- fold.list[[i]]
    train <- data[-fold, ]
    tst <- data[fold, ]
    x <- train[ , -1]
    y <- train[ , 1]
    svm.cv <- svm(x = x, 
                  y = y,
                  cost = cost,
                  kernel = kerneltype)
    svm.df <- predict(svm.cv, tst[ , -1]) %>% 
      cbind(tst[ , 1]) %>%
      data.frame()
    # Necessary to rename columns so that caret:defaultSummary works
    colnames(svm.df)[1] <- "pred"
    colnames(svm.df)[2] <- "obs"
    rmse.results[i] <- defaultSummary(svm.df)[1] 
    r2.results[i] <- defaultSummary(svm.df)[2] 
  }
  return(data.frame( # Useful for records
    seed = seed, 
    nfolds = nfolds,
    kernel = kerneltype,
    cost = cost,
    rsquared = sum(r2.results) / nfolds, 
    rmse = sum(rmse.results)/nfolds))
}

# Gamma determines the influence of any single point
# Does not work for linear kernels

tune.svm.gamma <- function(data, nfolds, g, kerneltype, seed) {
  set.seed(seed)
  fold.list <- createFolds(y = data[ , 1], k = nfolds)
  r2.results <- c(rep(0.0, nfolds))
  rmse.results <- c(rep(0.0, nfolds))
  for (i in 1:nfolds)
  {
    fold <- fold.list[[i]]
    train <- data[-fold, ]
    tst <- data[fold, ]
    x <- train[ , -1]
    y <- train[ , 1]
    svm.cv <- svm(x = x, 
                  y = y,
                  gamma = g,
                  kernel = kerneltype)
    svm.df <- predict(svm.cv, tst[ , -1]) %>% 
      cbind(tst[ , 1]) %>%
      data.frame()
    colnames(svm.df)[1] <- "pred"
    colnames(svm.df)[2] <- "obs"
    rmse.results[i] <- defaultSummary(svm.df)[1] 
    r2.results[i] <- defaultSummary(svm.df)[2] 
  }
  return(data.frame( # Useful for records
    seed = seed, 
    nfolds = nfolds,
    kernel = kerneltype,
    gamma = g,
    rsquared = sum(r2.results) / nfolds, 
    rmse = sum(rmse.results)/nfolds))
}

tune.svm.epsilon <- function(data, nfolds, e, kerneltype, seed) {
  set.seed(seed)
  fold.list <- createFolds(y = data[ , 1], k = nfolds)
  r2.results <- c(rep(0.0, nfolds))
  rmse.results <- c(rep(0.0, nfolds))
  for (i in 1:nfolds)
  {
    fold <- fold.list[[i]]
    train <- data[-fold, ]
    tst <- data[fold, ]
    x <- train[ , -1]
    y <- train[ , 1]
    svm.cv <- svm(x = x, 
                  y = y,
                  epsilon = e,
                  kernel = kerneltype)
    svm.df <- predict(svm.cv, tst[ , -1]) %>% 
      cbind(tst[ , 1]) %>%
      data.frame()
    colnames(svm.df)[1] <- "pred"
    colnames(svm.df)[2] <- "obs"
    R2 <- defaultSummary(svm.df)[2] 
    rmse.results[i] <- defaultSummary(svm.df)[1] 
    r2.results[i] <- defaultSummary(svm.df)[2] 
  }
  return(data.frame( # Useful for records
    seed = seed, 
    nfolds = nfolds,
    kernel = kerneltype,
    epsilon = e,
    rsquared = sum(r2.results) / nfolds, 
    rmse = sum(rmse.results)/nfolds))
}

# Constant coefficient - only for polynomial and sigmoid
tune.svm.coef <- function(data, nfolds, coef, kerneltype, seed) {
  set.seed(seed)
  fold.list <- createFolds(y = data[ , 1], k = nfolds)
  r2.results <- c(rep(0.0, nfolds))
  rmse.results <- c(rep(0.0, nfolds))
  for (i in 1:nfolds)
  {
    fold <- fold.list[[i]]
    train <- data[-fold, ]
    tst <- data[fold, ]
    x <- train[ , -1]
    y <- train[ , 1]
    svm.cv <- svm(x = x, 
                  y = y,
                  coef0 = coef,
                  kernel = kerneltype)
    svm.df <- predict(svm.cv, tst[ , -1]) %>% 
      cbind(tst[ , 1]) %>%
      data.frame()
    colnames(svm.df)[1] <- "pred"
    colnames(svm.df)[2] <- "obs"
    rmse.results[i] <- defaultSummary(svm.df)[1] 
    r2.results[i] <- defaultSummary(svm.df)[2] 
  }
  return(data.frame( # Useful for records
    seed = seed, 
    nfolds = nfolds,
    kernel = kerneltype,
    coef = coef,
    rsquared = sum(r2.results) / nfolds, 
    rmse = sum(rmse.results)/nfolds))
}

# degree: degree of polynomials for polynomial kernel

tune.svm.degree <- function(data, nfolds, deg, seed) {
  set.seed(seed)
  fold.list <- createFolds(y = data[ , 1], k = nfolds)
  r2.results <- c(rep(0.0, nfolds))
  rmse.results <- c(rep(0.0, nfolds))
  
  # Cross-validation using k-fold cross-validation
  for (i in 1:nfolds)
  {
    fold <- fold.list[[i]]
    train <- data[-fold, ]
    tst <- data[fold, ]
    x <- train[ , -1]
    y <- train[ , 1]
    svm.cv <- svm(x = x, 
                  y = y,
                  degree = deg,
                  kernel = "polynomial")
    svm.df <- predict(svm.cv, tst[ , -1]) %>% 
      cbind(tst[ , 1]) %>%
      data.frame()
    # Necessary to rename columns so that caret:defaultSummary works
    colnames(svm.df)[1] <- "pred"
    colnames(svm.df)[2] <- "obs"
    rmse.results[i] <- defaultSummary(svm.df)[1] 
    r2.results[i] <- defaultSummary(svm.df)[2] 
  }
  return(data.frame( # Useful for records
    seed = seed, 
    nfolds = nfolds,
    kernel = "polynomial",
    degree = deg,
    rsquared = sum(r2.results) / nfolds, 
    rmse = sum(rmse.results)/nfolds))
}


# Combination Tuning ------------------------------------------------------

# These are specific to kernels

# These functions, like other tuning functions, should be in a 
# do.call > rbind > lapply framework, though this one uses
# mapply instead (with SIMPLIFY = F). 

# tune.svm.[kernel] does not call set.seed() to save time. Remember to set.seed()
# outside of the function for reproducibility.

tune.svm.poly <- function(data, nfolds, deg, cost, e, g, coef) { 
  fold.list <- createFolds(y = data[ , 1], k = nfolds)
  r2.results <- c(rep(0.0, nfolds))
  rmse.results <- c(rep(0.0, nfolds))
  
  for (i in 1:nfolds)
  {
    fold <- fold.list[[i]]
    train <- data[-fold, ]
    tst <- data[fold, ]
    x <- train[ , -1]
    y <- train[ , 1]
    svm.cv <- svm(x = x, 
                  y = y,
                  cost = cost, 
                  degree = deg,
                  epsilon = e, 
                  gamma = g, 
                  coef0 = coef,
                  kernel = "polynomial")
    svm.df <- predict(svm.cv, tst[ , -1]) %>% 
      cbind(tst[ , 1]) %>%
      data.frame()
    colnames(svm.df)[1] <- "pred"
    colnames(svm.df)[2] <- "obs"
    rmse.results[i] <- defaultSummary(svm.df)[1] 
    r2.results[i] <- defaultSummary(svm.df)[2] 
  }
  
  # message("Cost: ", cost, " || Epsilon: ", e,
  #         " || Gamma: ", g, " || Coef0: ", coef)
  
  return(data.frame(
    nfolds = nfolds, kernel = "polynomial", 
    degree = deg, cost = cost, epsilon = e, 
    gamma = g, coef0 = coef, 
    rsquared = sum(r2.results)/nfolds,
    rmse = sum(rmse.results)/nfolds)
  )
}

# Radial basis function
tune.svm.rbf <- function(data, nfolds, cost, e, g) { 
  fold.list <- createFolds(y = data[ , 1], k = nfolds)
  r2.results <- c(rep(0.0, nfolds))
  rmse.results <- c(rep(0.0, nfolds))
  
  for (i in 1:nfolds)
  {
    fold <- fold.list[[i]]
    train <- data[-fold, ]
    tst <- data[fold, ]
    x <- train[ , -1]
    y <- train[ , 1]
    svm.cv <- svm(x = x, 
                  y = y,
                  cost = cost, 
                  epsilon = e, 
                  gamma = g, 
                  kernel = "radial")
    svm.df <- predict(svm.cv, tst[ , -1]) %>% 
      cbind(tst[ , 1]) %>%
      data.frame()
    colnames(svm.df)[1] <- "pred"
    colnames(svm.df)[2] <- "obs"
    rmse.results[i] <- defaultSummary(svm.df)[1] 
    r2.results[i] <- defaultSummary(svm.df)[2] 
  }
  
  # message("Cost: ", cost, " || Epsilon: ", e,
  #         " || Gamma: ", g, " || Coef0: ", coef)
  
  return(data.frame(
    nfolds = nfolds, kernel = "rbf", 
    cost = cost, epsilon = e, gamma = g, 
    rsquared = sum(r2.results)/nfolds,
    rmse = sum(rmse.results)/nfolds)
  )
}

# Sigmoid
tune.svm.sig <- function(data, nfolds, cost, e, g, coef) { 
  fold.list <- createFolds(y = data[ , 1], k = nfolds)
  r2.results <- c(rep(0.0, nfolds))
  rmse.results <- c(rep(0.0, nfolds))
  
  for (i in 1:nfolds)
  {
    fold <- fold.list[[i]]
    train <- data[-fold, ]
    tst <- data[fold, ]
    x <- train[ , -1]
    y <- train[ , 1]
    svm.cv <- svm(x = x, 
                  y = y,
                  cost = cost, 
                  epsilon = e, 
                  gamma = g, 
                  coef0 = coef,
                  kernel = "sigmoid")
    svm.df <- predict(svm.cv, tst[ , -1]) %>% 
      cbind(tst[ , 1]) %>%
      data.frame()
    colnames(svm.df)[1] <- "pred"
    colnames(svm.df)[2] <- "obs"
    rmse.results[i] <- defaultSummary(svm.df)[1] 
    r2.results[i] <- defaultSummary(svm.df)[2] 
  }
  
  # message("Cost: ", cost, " || Epsilon: ", e,
  #         " || Gamma: ", g, " || Coef0: ", coef)
  
  return(data.frame(
    nfolds = nfolds, kernel = "sigmoid", 
    cost = cost, epsilon = e, 
    gamma = g, coef0 = coef, 
    rsquared = sum(r2.results)/nfolds,
    rmse = sum(rmse.results)/nfolds)
  )
}
