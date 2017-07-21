# Libraries ---------------------------------------------------------------

library(caret)
library(Cubist)
library(e1071)
library(glmnet)
library(tidyverse)

# SVM ---------------------------------------------------------------------

# Creating a set of function to tune each parameter individually

# Cost - works for all kernels
tune.svm.cost <- function(data, nfolds, cost, kerneltype, seed) {
  set.seed(seed)
  fold.list <- createFolds(y = data[ , 2], k = nfolds)
  results <- c(rep(0.0, nfolds))
  # Cross-validation using k-fold cross-validation
  for (i in 1:nfolds)
  {
    fold <- fold.list[[i]]
    train <- data[-fold, ]
    tst <- data[fold, ]
    x <- train[ , -1:-2]
    y <- train[ , 2]
    svm.cv <- svm(x = x, 
                  y = y,
                  cost = cost,
                  kernel = kerneltype)
    svm.df <- predict(svm.cv, tst[ , -1:-2]) %>% 
      cbind(tst[ , 2]) %>%
      data.frame()
    # Necessary to rename columns so that caret:defaultSummary works
    colnames(svm.df)[1] <- "pred"
    colnames(svm.df)[2] <- "obs"
    R2 <- defaultSummary(svm.df)[2] 
    results[i] <- R2
  }
  return(data.frame( # Useful for records
    data = deparse(substitute(data)), # Turns the variable name into char
    nfolds = nfolds,
    seed = seed,
    kernel = kerneltype,
    cost = cost,
    rsquared = sum(results) / nfolds))
}
# Gamma - works for all kernels except linear 
tune.svm.gamma <- function(data, nfolds, g, kerneltype, seed) {
  set.seed(seed)
  fold.list <- createFolds(y = data[ , 2], k = nfolds)
  results <- c(rep(0.0, nfolds))
  for (i in 1:nfolds)
  {
    fold <- fold.list[[i]]
    train <- data[-fold, ]
    tst <- data[fold, ]
    x <- train[ , -1:-2]
    y <- train[ , 2]
    svm.cv <- svm(x = x, 
                  y = y,
                  gamma = g,
                  kernel = kerneltype)
    svm.df <- predict(svm.cv, tst[ , -1:-2]) %>% 
      cbind(tst[ , 2]) %>%
      data.frame()
    colnames(svm.df)[1] <- "pred"
    colnames(svm.df)[2] <- "obs"
    R2 <- defaultSummary(svm.df)[2] 
    results[i] <- R2
  }
  return(data.frame(
    data = deparse(substitute(data)),
    nfolds = nfolds,
    seed = seed,
    kernel = kerneltype,
    gamma = g,
    rsquared = sum(results) / nfolds))
}
# Epsilon - works for all kernels
tune.svm.epsilon <- function(data, nfolds, e, kerneltype, seed) {
  set.seed(seed)
  fold.list <- createFolds(y = data[ , 2], k = nfolds)
  results <- c(rep(0.0, nfolds))
  for (i in 1:nfolds)
  {
    fold <- fold.list[[i]]
    train <- data[-fold, ]
    tst <- data[fold, ]
    x <- train[ , -1:-2]
    y <- train[ , 2]
    svm.cv <- svm(x = x, 
                  y = y,
                  epsilon = e,
                  kernel = kerneltype)
    svm.df <- predict(svm.cv, tst[ , -1:-2]) %>% 
      cbind(tst[ , 2]) %>%
      data.frame()
    colnames(svm.df)[1] <- "pred"
    colnames(svm.df)[2] <- "obs"
    R2 <- defaultSummary(svm.df)[2] 
    results[i] <- R2
  }
  return(data.frame(
    data = deparse(substitute(data)),
    nfolds = nfolds,
    seed = seed,
    kernel = kerneltype,
    epsilon = e,
    rsquared = sum(results) / nfolds))
}
# Constant coefficient - only for polynomial and sigmoid
tune.svm.coef <- function(data, nfolds, coef, kerneltype, seed) {
  set.seed(seed)
  fold.list <- createFolds(y = data[ , 2], k = nfolds)
  results <- c(rep(0.0, nfolds))
  for (i in 1:nfolds)
  {
    fold <- fold.list[[i]]
    train <- data[-fold, ]
    tst <- data[fold, ]
    x <- train[ , -1:-2]
    y <- train[ , 2]
    svm.cv <- svm(x = x, 
                  y = y,
                  coef0 = coef,
                  kernel = kerneltype)
    svm.df <- predict(svm.cv, tst[ , -1:-2]) %>% 
      cbind(tst[ , 2]) %>%
      data.frame()
    colnames(svm.df)[1] <- "pred"
    colnames(svm.df)[2] <- "obs"
    R2 <- defaultSummary(svm.df)[2] 
    results[i] <- R2
  }
  return(data.frame(
    data = deparse(substitute(data)),
    nfolds = nfolds,
    seed = seed,
    kernel = kerneltype,
    coef = coef,
    rsquared = sum(results) / nfolds))
}


# GLMnet ------------------------------------------------------------------

tune.glm.alpha <- function(data, nfolds, alpha, seed) {
  set.seed(seed)
  fold.list <- createFolds(y = data[ , 2], k = nfolds)
  results <- c(rep(0.0, nfolds))
  
  for(i in 1:nfolds) {
    fold <- fold.list[[i]]
    
    trn.x <- data[-fold, -1:-2]
    trn.y <- data[-fold, 2]
    tst.x <- data[fold, -1:-2]
    tst.y <- data[fold, 2]
    
    glm.mod <- glmnet(x = trn.x, y = trn.y, 
                      alpha = alpha, 
                      family = "mgaussian")
    glm.df <- predict.glmnet(glm.mod, tst.x,
                             s = tail(glm.mod$lambda, n = 1)) %>%
      cbind(tst.y) %>%
      data.frame() %>%
      rename(pred = X1, obs = tst.y)
    
    results[i] <- defaultSummary(glm.df)[2]
  }
  
  return(data.frame( # Useful for records
    data = deparse(substitute(data)), # Turns the variable name into char
    nfolds = nfolds,
    seed = seed,
    alpha = alpha,
    rsquared = sum(results) / nfolds))
}

tune.glm.dfmax <- function(data, nfolds, max, seed) {
  set.seed(seed)
  fold.list <- createFolds(y = data[ , 2], k = nfolds)
  results <- c(rep(0.0, nfolds))
  
  for(i in 1:nfolds) {
    fold <- fold.list[[i]]
    
    trn.x <- data[-fold, -1:-2]
    trn.y <- data[-fold, 2]
    tst.x <- data[fold, -1:-2]
    tst.y <- data[fold, 2]
    
    glm.mod <- glmnet(x = trn.x, y = trn.y, 
                      dfmax = max, 
                      family = "mgaussian")
    glm.df <- predict.glmnet(glm.mod, tst.x,
                             s = tail(glm.mod$lambda, n = 1)) %>%
      cbind(tst.y) %>%
      data.frame() %>%
      rename(pred = X1, obs = tst.y)
    
    results[i] <- defaultSummary(glm.df)[2]
  }
  
  return(data.frame( # Useful for records
    data = deparse(substitute(data)), # Turns the variable name into char
    nfolds = nfolds,
    seed = seed,
    dfmax = max,
    rsquared = sum(results) / nfolds))
}


# Cubist ------------------------------------------------------------------

# Cubist should be fed with matrix or data.frame
tune.cubist.cmte <- function(data, nfolds, cmte, seed) {
  set.seed(seed)
  fold.list <- createFolds(y = data[ , 1], k = nfolds)
  results <- c(rep(0.0, nfolds))
  
  ctrl <- cubistControl(
    seed = seed
  )
  
  for(i in 1:nfolds) {
    fold <- fold.list[[i]]
    
    trn.x <- data[-fold, -1]
    trn.y <- data[-fold, 1]
    tst.x <- data[fold, -1]
    tst.y <- data[fold, 1]
    
    cube <- cubist(x = trn.x, y = trn.y, control = ctrl,
                   committees = cmte)
    cube.df <- predict(cube, tst.x) %>%
      cbind(tst.y) %>%
      data.frame() 
    
    colnames(cube.df)[1] <- "pred"
    colnames(cube.df)[2] <- "obs"

    results[i] <- defaultSummary(cube.df)[2]
  }
  
  return(data.frame( # Useful for records
    data = deparse(substitute(data)), # Turns the variable name into char
    nfolds = nfolds,
    seed = seed,
    committees = cmte,
    rsquared = sum(results) / nfolds))
}


