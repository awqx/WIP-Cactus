# Libraries and Packages --------------------------------------------------

library(caret)
library(data.table)
library(e1071)
library(kernlab)
library(Matrix)
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

# This function, like other tuning functions, should be in a 
# do.call > rbind > lapply framework, though this one uses
# mapply instead (with SIMPLIFY = F). 

# tune.svm.combos does not call set.seed() to save time. Remember to set.seed()
# outside of the function for reproducibility.

# Variables are the same as all other tuning functions For kernels that don't
# use all the variables, giving mapply a dummy value works. I would use ellipses
# (...), but they don't work well with do.call, for some reason

tune.svm.combos <- function(data, nfolds, kerneltype, cost, e, g, coef) { 
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
                  kernel = kerneltype)
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
    nfolds = nfolds, kernel = kerneltype, 
    cost = cost, epsilon = e, 
    gamma = g, coef0 = coef, 
    rsquared = sum(r2.results)/nfolds,
    rmse = sum(rmse.results)/nfolds)
    )
  }

# Alpha-CD ----------------------------------------------------------------

#     Loading Data --------------------------------------------------------

dir.create("./tuning")
dir.create("./tuning/svm")
# Reading data with all descriptors
trn.all <- readRDS("./model.data/alpha/trn1.RDS") 
colnames(trn.all) <- str_replace(colnames(trn.all), "-", ".")
trn.guest <- trn.all$guest
trn <- select(trn.all, -guest)

rfe1 <- readRDS("./feature.selection/alpha/rfe1.RDS")
trn.pred <- c("DelG", predictors(rfe1))

trn <- trn[ , colnames(trn) %in% trn.pred]
sprse <- sparse.model.matrix(~., trn)

#     Estimation ----------------------------------------------------------

# Working out where an acceptable range of values exists
# Trying to keep ranges to 7 values (arbitrary)

#     Cost ---

cost.range <- 2*(1:7)
results1.cost <- do.call(rbind, lapply(cost.range, FUN = tune.svm.cost, 
                         data = trn, kerneltype = "polynomial",
                         nfolds = 10, seed = 101)) 
results2.cost <- do.call(rbind, lapply(cost.range, FUN = tune.svm.cost,
                         data = trn, kerneltype = "polynomial",
                         nfolds = 10, seed = 102)) 
results3.cost <- do.call(rbind, lapply(cost.range, FUN = tune.svm.cost,
                         data = trn, kerneltype = "polynomial", 
                         nfolds = 10, seed = 103)) 
results.cost <- rbind(results1.cost, results2.cost, results3.cost)
ggplot(results.cost, aes(x = cost, color = seed, group = seed)) + 
  geom_line(aes(y = rsquared)) + 
  theme_bw()

#     Gamma ---

gamma.range <- 2^(-8:-2)
results1.gamma <- do.call(rbind, lapply(gamma.range, FUN = tune.svm.gamma,
                          data = trn, kerneltype = "polynomial", 
                          nfolds = 10, seed = 101))
results2.gamma <- do.call(rbind, lapply(gamma.range, FUN = tune.svm.gamma,
                                        data = trn, kerneltype = "polynomial", 
                                        nfolds = 10, seed = 102))
results3.gamma <- do.call(rbind, lapply(gamma.range, FUN = tune.svm.gamma,
                                        data = trn, kerneltype = "polynomial", 
                                        nfolds = 10, seed = 103))
results.gamma <- rbind(results1.gamma, results2.gamma, results3.gamma) %>%
  mutate(seed = as.factor(seed))
ggplot(results.gamma, aes(x = gamma, color = seed, group = seed)) + 
  geom_line(aes(y = rsquared)) + 
  theme_bw() + 
  scale_x_continuous(trans = "log2")

#     Epsilon ---

epsilon.range <- 2^(-6:0)
results1.epsilon <- do.call(rbind, lapply(epsilon.range, FUN = tune.svm.epsilon,
                            data = trn, kerneltype = "polynomial",
                            nfolds = 10, seed = 101))
results2.epsilon <- do.call(rbind, lapply(epsilon.range, FUN = tune.svm.epsilon,
                                          data = trn, kerneltype = "polynomial",
                                          nfolds = 10, seed = 102))
results3.epsilon <- do.call(rbind, lapply(epsilon.range, FUN = tune.svm.epsilon,
                                          data = trn, kerneltype = "polynomial",
                                          nfolds = 10, seed = 103))
results4.epsilon <- do.call(rbind, lapply(epsilon.range, FUN = tune.svm.epsilon,
                                          data = trn, kerneltype = "polynomial",
                                          nfolds = 10, seed = 104))
results.epsilon <- rbind(results1.epsilon, results2.epsilon, 
                         results3.epsilon, results4.epsilon) %>%
  mutate(seed = as.factor(seed))
ggplot(results.epsilon, aes(x = epsilon, color = seed, group = seed)) + 
  geom_line(aes(y = rsquared)) + 
  # scale_x_continuous(trans = "log2") + 
  theme_bw()

#     Coef ---

coef.range <- 2^(0:6)
results1.coef <- do.call(rbind, lapply(coef.range, FUN = tune.svm.coef,
                         data = trn, kerneltype = "polynomial", 
                         nfolds = 10, seed = 101))
results2.coef <- do.call(rbind, lapply(coef.range, FUN = tune.svm.coef,
                                       data = trn, kerneltype = "polynomial", 
                                       nfolds = 10, seed = 102))
results3.coef <- do.call(rbind, lapply(coef.range, FUN = tune.svm.coef,
                                       data = trn, kerneltype = "polynomial", 
                                       nfolds = 10, seed = 103))
results4.coef <- do.call(rbind, lapply(coef.range, FUN = tune.svm.coef,
                                       data = trn, kerneltype = "polynomial", 
                                       nfolds = 10, seed = 104))
results.coef <- rbind(results1.coef, results2.coef, 
                      results3.coef, results4.coef) %>%
  mutate(seed = as.factor(seed))
ggplot(results.coef, aes(x = coef, color = seed, group = seed)) + 
  geom_line(aes(y = rsquared)) + 
  scale_x_continuous(trans = "log2") + 
  theme_bw()

dir.create("./tuning/svm/alpha")
saveRDS(results.cost, "./tuning/svm/alpha/poly.cost.RDS")
saveRDS(results.gamma, "./tuning/svm/alpha/poly.gamma.RDS")
saveRDS(results.epsilon, "./tuning/svm/alpha/poly.epsilon.RDS")
saveRDS(results.coef, "./tuning/svm/alpha/poly.coef.RDS")

#     Tuning --------------------------------------------------------------

# 7^4 = 2401 tuning combinations
svm.combos <- expand.grid(cost.range, gamma.range, 
                          epsilon.range, coef.range)
colnames(svm.combos) <- c("cost", "gamma", "epsilon", "coef")
cost.combos <- svm.combos$cost
gamma.combos <- svm.combos$gamma
eps.combos <- svm.combos$epsilon
coef.combos <- svm.combos$coef


set.seed(1001)
system.time(
  results.combos <- do.call(
    rbind,
      mapply(
      FUN = tune.svm.combos,
      cost = cost.combos, 
      e = eps.combos, 
      g = gamma.combos,
      coef = coef.combos,
      MoreArgs = 
        list(nfolds = 5, data = trn, kerneltype = "polynomial"), 
      SIMPLIFY = F
    )
  )
)

# system.time output
# user      system    elapsed 
# 429.69    0.97      443.67

saveRDS(results.combos, "./tuning/svm/alpha/poly.tuning.RDS")

# results.combos[order(results.combos$rsquared, decreasing = T), ] %>% head()
# results.combos[order(results.combos$rmse), ] %>% head()

# Best rsquared (0.526)
# cost = 12, eps = 0.03125, gamma = 0.25, coef0 = 4 (rmse = 3.84)
# Best rmse (3.71)
# cost = 6, eps = 0.03125, gamma = 0.00390625, coef0 = 1 (r2 = 0.501)