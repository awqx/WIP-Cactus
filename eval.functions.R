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
find.k <- function(df) {
  df <- lapply(df, as.character) %>% sapply(., as.numeric)
  return(sum(df$pred * df$obs) / sum(df$pred^2))
}

# Equation writeen in Xu2015 unclear
find.kprime <- function(df) {
  df <- lapply(df, as.character) %>% sapply(as.numeric) %>% 
    data.frame()
  return(sum(df$pred * df$obs) / sum(df$obs^2))
}

# Since obsolete
# find.q2.f1 <- function(df) {
#   top <- (df$obs - df$pred)^2 %>% sum()
#   bottom <- (df$obs - mean(df$pred))^2 %>% sum()
#   return (1 - (top/bottom))
# }

# Uses caret::defaultSummary
find.r2 <- function(df) {
  return(defaultSummary(df)[2])
}

find.r02 <- function(df) {
  # basically the regression line from the predictions
  yir0 <- find.k(df) * df$pred
  # mean of predictions
  ybar <- mean(df$pred)
  top <- sum((df$obs - yir0)^2) 
  bottom <- sum((df$pred - ybar)^2) 
  return(1 - (top/bottom))
}

# May be wrong due to reasons specified in find.kprime
find.rprime02 <- function(df) {
  # Regression line from observations
  yir0 <- find.k(df) * df$obs
  # mean of observations
  ybar <- mean(df$obs)
  top <- sum((df$pred - yir0)^2) 
  bottom <- sum((df$obs - ybar)^2)
  return(1 - (top/bottom))
}

eval.tropsha <- function(df) {
  # R2 > 0.6
  a <- find.r2(df) 
  # (R2-R02)/R2 < 0.1 and 0.85 <= k <= 1.15
  b1 <- abs(a - find.r02(df))/find.r2(df) 
  b2 <- find.k(df)
  # (R2-R0'2)/R2 < 0.1 and 0.85 <= k' <= 1.15
  c1 <- abs(a - find.rprime02(df))/find.r2(df)
  c2 <- find.kprime(df)
  # iv <- abs(find.r20(df)-find.r20prime(df))
  results <- c(a, b1, b2, c1, c2)
  names(results) <- c('> 0.6', '< 0.1', '0.85 < x < 1.15', 
                      '< 0.1', '0.85 < x < 1.15')
  return(results)
}

# 10-fold Cross Validation ------------------------------------------------

# Using the parameters determined by tuning
# Only input is folds. Use folds = number observations for LOOCV

# Data 
# cv.svm <- function(data, nfolds) {
#   set.seed(10)
#   fold.list <- createFolds(y = data[ , 2], k = nfolds)
#   results <- c(rep(0.0, nfolds))
#   times <- c(rep(0.0, nfolds))
#   # Cross-validation using k-fold cross-validation
#   
#   tic.clearlog()
#   for (i in 1:nfolds)
#   {
#     tic()
#     fold <- fold.list[[i]]
#     train <- data[-fold, ]
#     tst <- data[fold, ]
#     x <- train[ , -1:-2]
#     y <- train[ , 2]
#     svm.cv <- svm(x = x, 
#                   y = y,
#                   cost = 1024,
#                   kernel = "polynomial", 
#                   degree = 2, 
#                   gamma = 0.5, 
#                   epsilon = 0.1, 
#                   coef0 = 2)
#     svm.df <- predict(svm.cv, tst[ , -1:-2]) %>% 
#       cbind(tst[ , 2]) %>%
#       data.frame()
#     # Necessary to rename columns so that caret:defaultSummary works
#     colnames(svm.df)[1] <- "pred"
#     colnames(svm.df)[2] <- "obs"
#     R2 <- defaultSummary(svm.df)[2] 
#     results[i] <- R2
#     toc(log = T, quiet = T)
#     times[i] <- toc.outmsg(tic = 1, toc = 1, "done")
#   }
#   return(data.frame(
#     rsquared = results, 
#     time = unlist(tic.log(format = T))))
# }
# 
# cv.cube <- function(data, nfolds) {
#   set.seed(10)
#   fold.list <- createFolds(y = data[ , 1], k = nfolds)
#   results <- c(rep(0.0, nfolds))
#   
#   ctrl <- cubistControl(
#     seed = 10, 
#     sample = 75
#   )
#   
#   for(i in 1:nfolds) {
#     fold <- fold.list[[i]]
#     
#     trn.x <- data[-fold, -1]
#     trn.y <- data[-fold, 1]
#     tst.x <- data[fold, -1]
#     tst.y <- data[fold, 1]
#     
#     cube <- cubist(x = trn.x, 
#                    y = trn.y, 
#                    control = ctrl, 
#                    committees = 90)
#     cube.df <- predict(cube, tst.x) %>%
#       cbind(tst.y) %>%
#       data.frame() 
#     
#     colnames(cube.df)[1] <- "pred"
#     colnames(cube.df)[2] <- "obs"
#     
#     results[i] <- defaultSummary(cube.df)[2]
#     message("Fold ", i, " completed.")
#   }
#   
#   return(sum(results)/nfolds)
# }
# 
# cv.svm(df, 3)
# 
