# Functions for evaluating QSARS according to Tropsha and Golbraikh 2007
# Contains timing for models as well as 10-fold CV

# setwd("~/SREP LAB/qsar)

# Libraries and Packages --------------------------------------------------

# install.packages("tictoc")
library(e1071)
library(tictoc)
library(tidyverse)
library(caret)
library(stats)

# Functions for Tropsha----------------------------------------------------

# Requirements: results is a datframe w/ calculated values labeled as pred 
# and observed values listed as obs
find.k <- function(results) {
  top <- sum(results$pred * results$obs)
  bottom <- sum(results$obs * results$obs)
  return (top/bottom)
}

# Equation writeen in Xu2015 unclear
find.kprime <- function(results) {
  top <- sum(results$pred * results$obs)
  bottom <- sum(results$pred ^ 2)
  return (top/bottom)
}

find.q2.f1 <- function(results) {
  top <- (results$obs - results$pred)^2 %>% sum()
  bottom <- (results$obs - mean(results$pred))^2 %>% sum()
  return (1 - (top/bottom))
}

# Uses caret::defaultSummary
find.r2 <- function(results) {
  return(defaultSummary(results)[2])
}

find.r20 <- function(results) {
  ybarr0 <- find.k(results) * results$obs
  top <- (results$pred - ybarr0)^2 %>% sum()
  bottom <- (results$pred - mean(results$pred))^2 %>% sum()
  return(1 - (top/bottom))
}

# May be wrong due to reasons specified in find.kprime
find.r20prime <- function(results) {
  yr0 <- find.kprime(results) * results$pred
  top <- (results$obs - yr0)^2 %>% sum()
  bottom <- (results$obs - mean(results$obs))^2 %>% sum()
  return(1 - (top/bottom))
}

# 10-fold Cross Validation ------------------------------------------------

# Using the parameters determined by tuning
# Only input is folds. Use folds = number observations for LOOCV

# Data 
cv.svm <- function(data, nfolds) {
  set.seed(10)
  fold.list <- createFolds(y = data[ , 2], k = nfolds)
  results <- c(rep(0.0, nfolds))
  times <- c(rep(0.0, nfolds))
  # Cross-validation using k-fold cross-validation
  
  tic.clearlog()
  for (i in 1:nfolds)
  {
    tic()
    fold <- fold.list[[i]]
    train <- data[-fold, ]
    tst <- data[fold, ]
    x <- train[ , -1:-2]
    y <- train[ , 2]
    svm.cv <- svm(x = x, 
                  y = y,
                  cost = 1024,
                  kernel = "polynomial", 
                  degree = 2, 
                  gamma = 0.5, 
                  epsilon = 0.1, 
                  coef0 = 2)
    svm.df <- predict(svm.cv, tst[ , -1:-2]) %>% 
      cbind(tst[ , 2]) %>%
      data.frame()
    # Necessary to rename columns so that caret:defaultSummary works
    colnames(svm.df)[1] <- "pred"
    colnames(svm.df)[2] <- "obs"
    R2 <- defaultSummary(svm.df)[2] 
    results[i] <- R2
    toc(log = T, quiet = T)
    times[i] <- toc.outmsg(tic = 1, toc = 1, "done")
  }
  return(data.frame(
    rsquared = results, 
    time = unlist(tic.log(format = T))))
}

cv.cube <- function(data, nfolds) {
  set.seed(10)
  fold.list <- createFolds(y = data[ , 1], k = nfolds)
  results <- c(rep(0.0, nfolds))
  
  ctrl <- cubistControl(
    seed = 10, 
    sample = 75
  )
  
  for(i in 1:nfolds) {
    fold <- fold.list[[i]]
    
    trn.x <- data[-fold, -1]
    trn.y <- data[-fold, 1]
    tst.x <- data[fold, -1]
    tst.y <- data[fold, 1]
    
    cube <- cubist(x = trn.x, 
                   y = trn.y, 
                   control = ctrl, 
                   committees = 90)
    cube.df <- predict(cube, tst.x) %>%
      cbind(tst.y) %>%
      data.frame() 
    
    colnames(cube.df)[1] <- "pred"
    colnames(cube.df)[2] <- "obs"
    
    results[i] <- defaultSummary(cube.df)[2]
    message("Fold ", i, " completed.")
  }
  
  return(sum(results)/nfolds)
}

cv.svm(df, 3)

