# Libraries and Packages --------------------------------------------------

library(caret)
library(Cubist)
library(directlabels)
library(e1071)
library(ggplot2)
library(glmnet)
library(pls)
library(randomForest)
library(stringr)
library(tictoc)
library(tidyverse)

# Data Organization -------------------------------------------------------

set.seed(101)
# setwd("~/SREP LAB/qsar")
# Complete dataset
data.raw <- readRDS("./padel.pp.new.RDS")
data <- data.raw %>% select(-guest, -host, -data.source)
data.folds <- createFolds(data$DelG, 10)

# Suzuki data only
suz.raw <- readRDS("./suz.pp.RDS")
suz <- suz.raw %>% select(-guest, -host, -data.source)
suz.folds <- createFolds(suz$DelG, 10)

# Functions ---------------------------------------------------------------

# Returns the value from a single fold from polynomial SVM
svm.fold <- function(data, folds, foldnum, suzuki) {
  
  tst <- data[folds[[foldnum]], ]
  trn <- data[-folds[[foldnum]], ]
  
  if (suzuki) {
    
    a.trn <- trn %>% filter(alpha > 0)
    b.trn <- trn %>% filter(beta > 0)
    a.tst <- tst %>% filter(alpha > 0)
    b.tst <- tst %>% filter(beta > 0)
    
    svm.alpha <- svm(x = a.trn[ , -1], 
                     y = a.trn[ , 1], 
                     coef0 = 2, 
                     cost = 2048, 
                     epsilon = 0.1, 
                     kernel = "polynomial", 
                     gamma = 0.03, 
                     degree = 2)
    svm.beta <- svm(x = b.trn[ , -1], 
                    y = b.trn[ , 1], 
                    coef0 = 2, 
                    cost = 2048, 
                    epsilon = 0.1, 
                    kernel = "polynomial", 
                    gamma = 0.03, 
                    degree = 2)
    
    alpha.df <- predict(svm.alpha, a.tst[ , -1]) %>%
      cbind(a.tst[ , 1]) %>% data.frame() %>%
      dplyr::rename(., pred = `.`, obs = V2)
    beta.df <- predict(svm.beta, b.tst[ , -1]) %>%
      cbind(b.tst[ , 1]) %>% data.frame() %>%
      dplyr::rename(., pred = `.`, obs = V2)
    
    comp <- rbind(alpha.df, beta.df)
  } else {
    a.trn <- trn %>% filter(alpha > 0)
    b.trn <- trn %>% filter(beta > 0)
    c.trn <- trn %>% filter(gamma > 0)
    
    a.tst <- tst %>% filter(alpha > 0)
    b.tst <- tst %>% filter(beta > 0)
    c.tst <- tst %>% filter(gamma > 0)
    
    
    if (nrow(c.trn) < 1) {
      
      svm.alpha <- svm(x = a.trn[ , -1], 
                       y = a.trn[ , 1], 
                       coef0 = 1.5, 
                       cost = 4096, 
                       epsilon = 0.1, 
                       kernel = "polynomial", 
                       gamma = 0.03, 
                       degree = 2)
      svm.beta <- svm(x = b.trn[ , -1], 
                      y = b.trn[ , 1], 
                      coef0 = 1.5, 
                      cost = 4096, 
                      epsilon = 0.1, 
                      kernel = "polynomial", 
                      gamma = 0.03, 
                      degree = 2)
      
      alpha.df <- predict(svm.alpha, a.tst[ , -1]) %>%
        cbind(a.tst[ , 1]) %>% data.frame() %>%
        dplyr::rename(., pred = `.`, obs = V2)
      beta.df <- predict(svm.beta, b.tst[ , -1]) %>%
        cbind(b.tst[ , 1]) %>% data.frame() %>%
        dplyr::rename(., pred = `.`, obs = V2)
      
      comp <- rbind(alpha.df, beta.df)
    } else {
      svm.alpha <- svm(x = a.trn[ , -1], 
                       y = a.trn[ , 1], 
                       coef0 = 1.5, 
                       cost = 4096, 
                       epsilon = 0.1, 
                       kernel = "polynomial", 
                       gamma = 0.03, 
                       degree = 2)
      svm.beta <- svm(x = b.trn[ , -1], 
                      y = b.trn[ , 1], 
                      coef0 = 1.5, 
                      cost = 4096, 
                      epsilon = 0.1, 
                      kernel = "polynomial", 
                      gamma = 0.03, 
                      degree = 2)
      svm.gamma <- svm(x = c.trn[ , -1], 
                       y = c.trn[ , 1], 
                       coef0 = 1.5, 
                       cost = 4096, 
                       epsilon = 0.1, 
                       kernel = "polynomial", 
                       gamma = 0.03, 
                       degree = 2)
      
      alpha.df <- predict(svm.alpha, a.tst[ , -1]) %>%
        cbind(a.tst[ , 1]) %>% data.frame() %>%
        dplyr::rename(., pred = `.`, obs = V2)
      beta.df <- predict(svm.beta, b.tst[ , -1]) %>%
        cbind(b.tst[ , 1]) %>% data.frame() %>%
        dplyr::rename(., pred = `.`, obs = V2)
      gamma.df <- predict(svm.gamma, c.tst[ , -1]) %>%
        cbind(c.tst[ , 1]) %>% data.frame() %>%
        dplyr::rename(., pred = `.`, obs = V2)
      
      comp <- rbind(alpha.df, beta.df, gamma.df)
    }
  }
  return(data.frame(
    data = deparse(substitute(data)), 
    fold = foldnum, 
    model = "SVM", 
    rmse = defaultSummary(comp)[1],
    rsquared = defaultSummary(comp)[2]
  ))
}

# GLMnet folds - doesn't work with regular data (probably b/c gamma)
glm.fold <- function(data, folds, foldnum, suzuki) {
  
  tst <- data[folds[[foldnum]], ]
  trn <- data[-folds[[foldnum]], ]
  
  if (suzuki) {
    
    a.trn <- trn %>% filter(alpha > 0)
    b.trn <- trn %>% filter(beta > 0)
    a.tst <- tst %>% filter(alpha > 0)
    b.tst <- tst %>% filter(beta > 0)
    
    a.trn <- as.matrix(a.trn)
    b.trn <- as.matrix(b.trn)
    a.tst <- as.matrix(a.tst)
    b.tst <- as.matrix(b.tst)
    
    glm.alpha <- glmnet(x = a.trn[ , -1], 
                        y = a.trn[ , 1], 
                        dfmax = 32, 
                        alpha = 1, 
                        family = "mgaussian")
    glm.beta <- glmnet(x = b.trn[ , -1], 
                       y = b.trn[ , 1], 
                       dfmax = 32, 
                       alpha = 1, 
                       family = "mgaussian")
    
    alpha.df <- predict.glmnet(glm.alpha, a.tst[,-1],
                               s = tail(glm.alpha$lambda, n = 1)) %>%
      cbind(a.tst[, 1]) %>% data.frame()
    colnames(alpha.df)[1] <- "pred"
    colnames(alpha.df)[2] <- "obs"
    beta.df <- predict.glmnet(glm.beta, b.tst[,-1],
                              s = tail(glm.beta$lambda, n = 1)) %>%
      cbind(b.tst[, 1]) %>% data.frame()
    colnames(beta.df)[1] <- "pred"
    colnames(beta.df)[2] <- "obs"
    
    comp <- rbind(alpha.df, beta.df)
  } else {
    a.trn <- trn %>% filter(alpha > 0)
    b.trn <- trn %>% filter(beta > 0)
    c.trn <- trn %>% filter(gamma > 0)
    
    a.tst <- tst %>% filter(alpha > 0)
    b.tst <- tst %>% filter(beta > 0)
    c.tst <- tst %>% filter(gamma > 0)
    
    a.trn <- as.matrix(a.trn)
    b.trn <- as.matrix(b.trn)
    c.trn <- as.matrix(c.trn)
    a.tst <- as.matrix(a.tst)
    b.tst <- as.matrix(b.tst)
    c.tst <- as.matrix(c.tst)
    
    if (nrow(c.trn) <= 1 || nrow(c.tst) <= 1) {
      
      glm.alpha <- glmnet(x = a.trn[ , -1], 
                          y = a.trn[ , 1], 
                          dfmax = 32, 
                          alpha = 1, 
                          family = "mgaussian")
      glm.beta <- glmnet(x = b.trn[ , -1], 
                         y = b.trn[ , 1], 
                         dfmax = 32, 
                         alpha = 1, 
                         family = "mgaussian")
      
      alpha.df <- predict.glmnet(glm.alpha, a.tst[,-1],
                                 s = tail(glm.alpha$lambda, n = 1)) %>%
        cbind(a.tst[, 1]) %>% data.frame()
      colnames(alpha.df)[1] <- "pred"
      colnames(alpha.df)[2] <- "obs"
      beta.df <- predict.glmnet(glm.beta, b.tst[,-1],
                                s = tail(glm.beta$lambda, n = 1)) %>%
        cbind(b.tst[, 1]) %>% data.frame()
      colnames(beta.df)[1] <- "pred"
      colnames(beta.df)[2] <- "obs"
      
      comp <- rbind(alpha.df, beta.df)
    } else {
      glm.alpha <- glmnet(x = a.trn[ , -1], 
                          y = a.trn[ , 1], 
                          dfmax = 32, 
                          alpha = 1, 
                          family = "mgaussian")
      glm.beta <- glmnet(x = b.trn[ , -1], 
                         y = b.trn[ , 1], 
                         dfmax = 32, 
                         alpha = 1, 
                         family = "mgaussian")
      
      glm.gamma <- glmnet(x = c.trn[ , -1], 
                          y = c.trn[ , 1], 
                          dfmax = 32, 
                          alpha = 1, 
                          family = "mgaussian")
      
      alpha.df <- predict.glmnet(glm.alpha, a.tst[,-1],
                                 s = tail(glm.alpha$lambda, n = 1)) %>%
        cbind(a.tst[, 1]) %>% data.frame()
      colnames(alpha.df)[1] <- "pred"
      colnames(alpha.df)[2] <- "obs"
      
      beta.df <- predict.glmnet(glm.beta, b.tst[,-1],
                                s = tail(glm.beta$lambda, n = 1)) %>%
        cbind(b.tst[, 1]) %>% data.frame()
      colnames(beta.df)[1] <- "pred"
      colnames(beta.df)[2] <- "obs"
      
      gamma.df <- predict.glmnet(glm.gamma, c.tst[ , -1],
                                 s = tail(glm.gamma$lambda, n = 1)) %>%
        cbind(c.tst[, 1]) %>% data.frame()
      colnames(gamma.df)[1] <- "pred"
      colnames(gamma.df)[2] <- "obs"
      
      comp <- rbind(alpha.df, beta.df, gamma.df)
    }
  }
  return(data.frame(
    data = deparse(substitute(data)), 
    fold = foldnum, 
    model = "GLMnet", 
    rmse = defaultSummary(comp)[1],
    rsquared = defaultSummary(comp)[2]
  ))
}

# Cubist folds 
cube.fold <- function(data, folds, foldnum, suzuki) {
  
  ctrl <- cubistControl(
    seed = 12, 
    sample = 75
  )
  tst <- data[folds[[foldnum]], ]
  trn <- data[-folds[[foldnum]], ]
  
  if (suzuki) {
    
    a.trn <- trn %>% filter(alpha > 0)
    b.trn <- trn %>% filter(beta > 0)
    a.tst <- tst %>% filter(alpha > 0)
    b.tst <- tst %>% filter(beta > 0)
    
    cube.alpha <- cubist(a.trn[ , -1], a.trn[ , 1],
                         control = ctrl,
                         committees = 75)
    cube.beta <- cubist(b.trn[ , -1], b.trn[ , 1],
                        control = ctrl,
                        committees = 75)
    
    alpha.df <- predict(cube.alpha, a.tst[ , -1]) %>%
      cbind(a.tst[ , 1]) %>% data.frame() %>%
      dplyr::rename(., pred = `.`, obs = V2)
    beta.df <- predict(cube.beta, b.tst[ , -1]) %>%
      cbind(b.tst[ , 1]) %>% data.frame() %>%
      dplyr::rename(., pred = `.`, obs = V2)
    
    comp <- rbind(alpha.df, beta.df)
  } else {
    a.trn <- trn %>% filter(alpha > 0)
    b.trn <- trn %>% filter(beta > 0)
    c.trn <- trn %>% filter(gamma > 0)
    
    a.tst <- tst %>% filter(alpha > 0)
    b.tst <- tst %>% filter(beta > 0)
    c.tst <- tst %>% filter(gamma > 0)
    
    
    if (nrow(c.trn) < 1 || nrow(c.tst) < 1) {
      
      cube.alpha <- cubist(a.trn[ , -1], a.trn[ , 1],
                           control = ctrl,
                           committees = 75)
      cube.beta <- cubist(b.trn[ , -1], b.trn[ , 1],
                          control = ctrl,
                          committees = 75)
      
      alpha.df <- predict(cube.alpha, a.tst[ , -1]) %>%
        cbind(a.tst[ , 1]) %>% data.frame() %>%
        dplyr::rename(., pred = `.`, obs = V2)
      beta.df <- predict(cube.beta, b.tst[ , -1]) %>%
        cbind(b.tst[ , 1]) %>% data.frame() %>%
        dplyr::rename(., pred = `.`, obs = V2)
      
      comp <- rbind(alpha.df, beta.df)
      
    } else {
      cube.alpha <- cubist(a.trn[ , -1], a.trn[ , 1],
                           control = ctrl,
                           committees = 75)
      cube.beta <- cubist(b.trn[ , -1], b.trn[ , 1],
                          control = ctrl,
                          committees = 75)
      cube.gamma <- cubist(c.trn[ , -1], c.trn[ , 1],
                           control = ctrl,
                           committees = 75)
      
      alpha.df <- predict(cube.alpha, a.tst[ , -1]) %>%
        cbind(a.tst[ , 1]) %>% data.frame() %>%
        dplyr::rename(., pred = `.`, obs = V2)
      beta.df <- predict(cube.beta, b.tst[ , -1]) %>%
        cbind(b.tst[ , 1]) %>% data.frame() %>%
        dplyr::rename(., pred = `.`, obs = V2)
      gamma.df <- predict(cube.gamma, c.tst[ , -1]) %>%
        cbind(c.tst[ , 1]) %>% data.frame() %>%
        dplyr::rename(., pred = `.`, obs = V2)
      
      comp <- rbind(alpha.df, beta.df, gamma.df)
    }
  }
  return(data.frame(
    data = deparse(substitute(data)), 
    fold = foldnum, 
    model = "Cubist", 
    rmse = defaultSummary(comp)[1],
    rsquared = defaultSummary(comp)[2]
  ))
}

# Random forest
rf.fold <- function(data, folds, foldnum, suzuki) {
  
  tst <- data[folds[[foldnum]], ]
  trn <- data[-folds[[foldnum]], ]
  
  if (suzuki) {
    
    a.trn <- trn %>% filter(alpha > 0)
    b.trn <- trn %>% filter(beta > 0)
    a.tst <- tst %>% filter(alpha > 0)
    b.tst <- tst %>% filter(beta > 0)
    
    rf.alpha <- randomForest(
      a.trn[,-1], a.trn[, 1],
      ntree = 500, na.action = na.omit,
      mtry = 700, nodesize = 1
    )
    rf.beta <- randomForest(
      b.trn[,-1], b.trn[, 1],
      ntree = 500, na.action = na.omit,
      mtry = 700, nodesize = 1
    )
    
    alpha.df <- predict(rf.alpha, a.tst[ , -1]) %>%
      cbind(a.tst[ , 1]) %>% data.frame() %>%
      dplyr::rename(., pred = `.`, obs = V2)
    beta.df <- predict(rf.beta, b.tst[ , -1]) %>%
      cbind(b.tst[ , 1]) %>% data.frame() %>%
      dplyr::rename(., pred = `.`, obs = V2)
    
    comp <- rbind(alpha.df, beta.df)
    
  } else {
    a.trn <- trn %>% filter(alpha > 0)
    b.trn <- trn %>% filter(beta > 0)
    c.trn <- trn %>% filter(gamma > 0)
    
    a.tst <- tst %>% filter(alpha > 0)
    b.tst <- tst %>% filter(beta > 0)
    c.tst <- tst %>% filter(gamma > 0)
    
    
    if (nrow(c.trn) < 1 || nrow(c.tst) < 1) {
      
      rf.alpha <- randomForest(
        a.trn[,-1], a.trn[, 1],
        ntree = 500, na.action = na.omit,
        mtry = 700, nodesize = 1
      )
      rf.beta <- randomForest(
        b.trn[,-1], b.trn[, 1],
        ntree = 500, na.action = na.omit,
        mtry = 700, nodesize = 1
      )
      
      alpha.df <- predict(rf.alpha, a.tst[ , -1]) %>%
        cbind(a.tst[ , 1]) %>% data.frame() %>%
        dplyr::rename(., pred = `.`, obs = V2)
      beta.df <- predict(rf.beta, b.tst[ , -1]) %>%
        cbind(b.tst[ , 1]) %>% data.frame() %>%
        dplyr::rename(., pred = `.`, obs = V2)
      
      comp <- rbind(alpha.df, beta.df)
      
    } else {
      rf.alpha <- randomForest(
        a.trn[,-1], a.trn[, 1],
        ntree = 500, na.action = na.omit,
        mtry = 700, nodesize = 1
      )
      rf.beta <- randomForest(
        b.trn[,-1], b.trn[, 1],
        ntree = 500, na.action = na.omit,
        mtry = 700, nodesize = 1
      )
      rf.gamma <- randomForest(
        c.trn[,-1], c.trn[, 1],
        ntree = 500, na.action = na.omit,
        mtry = 700, nodesize = 1
      )
      
      alpha.df <- predict(rf.alpha, a.tst[ , -1]) %>%
        cbind(a.tst[ , 1]) %>% data.frame() %>%
        dplyr::rename(., pred = `.`, obs = V2)
      beta.df <- predict(rf.beta, b.tst[ , -1]) %>%
        cbind(b.tst[ , 1]) %>% data.frame() %>%
        dplyr::rename(., pred = `.`, obs = V2)
      gamma.df <- predict(rf.gamma, c.tst[ , -1]) %>%
        cbind(c.tst[ , 1]) %>% data.frame() %>%
        dplyr::rename(., pred = `.`, obs = V2)
      
      comp <- rbind(alpha.df, beta.df, gamma.df)
    }
  }
  return(data.frame(
    data = deparse(substitute(data)), 
    fold = foldnum, 
    model = "Random Forest", 
    rmse = defaultSummary(comp)[1],
    rsquared = defaultSummary(comp)[2]
  ))
}

pls.fold <- function(data, folds, foldnum, suzuki) {
  
  tst <- data[folds[[foldnum]], ]
  trn <- data[-folds[[foldnum]], ]
  
  if (suzuki) {
    
    a.trn <- trn %>% filter(alpha > 0)
    b.trn <- trn %>% filter(beta > 0)
    a.tst <- tst %>% filter(alpha > 0)
    b.tst <- tst %>% filter(beta > 0)
    
    pls.alpha <- plsr(DelG ~ ., 
                      ncomp = 11, 
                      data = a.trn, 
                      validation = "LOO", 
                      method = "oscorespls", 
                      seed = 10)
    pls.beta <- plsr(DelG ~ ., 
                     ncomp = 11, 
                     data = b.trn, 
                     validation = "LOO", 
                     method = "oscorespls", 
                     seed = 10)
    
    alpha.df <- predict(pls.alpha, ncomp = 11, newdata = a.tst) %>%
      cbind(a.tst[ , 1]) %>%
      data.frame() %>%
      dplyr::rename(., pred = `.`, obs = V2)
    beta.df <- predict(pls.beta, ncomp = 11, newdata = b.tst) %>%
      cbind(b.tst[ , 1]) %>%
      data.frame() %>%
      dplyr::rename(., pred = `.`, obs = V2)
    
    comp <- rbind(alpha.df, beta.df)
    
  } else {
    a.trn <- trn %>% filter(alpha > 0)
    b.trn <- trn %>% filter(beta > 0)
    c.trn <- trn %>% filter(gamma > 0)
    
    a.tst <- tst %>% filter(alpha > 0)
    b.tst <- tst %>% filter(beta > 0)
    c.tst <- tst %>% filter(gamma > 0)
    
    
    if (nrow(c.trn) < 1 || nrow(c.tst) < 1) {
      
      pls.alpha <- plsr(DelG ~ ., 
                        ncomp = 11, 
                        data = a.trn, 
                        validation = "LOO", 
                        method = "oscorespls", 
                        seed = 10)
      pls.beta <- plsr(DelG ~ ., 
                       ncomp = 11, 
                       data = b.trn, 
                       validation = "LOO", 
                       method = "oscorespls", 
                       seed = 10)
      
      alpha.df <- predict(pls.alpha, ncomp = 11, newdata = a.tst) %>%
        cbind(a.tst[ , 1]) %>%
        data.frame() %>%
        dplyr::rename(., pred = `.`, obs = V2)
      beta.df <- predict(pls.beta, ncomp = 11, newdata = b.tst) %>%
        cbind(b.tst[ , 1]) %>%
        data.frame() %>%
        dplyr::rename(., pred = `.`, obs = V2)
      
      comp <- rbind(alpha.df, beta.df)
      
    } else {
      a.trn <- trn %>% filter(alpha > 0)
      b.trn <- trn %>% filter(beta > 0)
      a.tst <- tst %>% filter(alpha > 0)
      b.tst <- tst %>% filter(beta > 0)
      
      pls.alpha <- plsr(DelG ~ ., 
                        ncomp = 11, 
                        data = a.trn, 
                        validation = "LOO", 
                        method = "oscorespls", 
                        seed = 10)
      pls.beta <- plsr(DelG ~ ., 
                       ncomp = 11, 
                       data = b.trn, 
                       validation = "LOO", 
                       method = "oscorespls", 
                       seed = 10)
      pls.gamma <- plsr(DelG ~ ., 
                        ncomp = 5, 
                        data = c.trn, 
                        validation = "LOO", 
                        method = "oscorespls", 
                        seed = 10)
      
      alpha.df <- predict(pls.alpha, ncomp = 11, newdata = a.tst) %>%
        cbind(a.tst[ , 1]) %>%
        data.frame() %>%
        dplyr::rename(., pred = `.`, obs = V2)
      beta.df <- predict(pls.beta, ncomp = 11, newdata = b.tst) %>%
        cbind(b.tst[ , 1]) %>%
        data.frame() %>%
        dplyr::rename(., pred = `.`, obs = V2)
      gamma.df <- predict(pls.gamma, c.tst) %>%
        cbind(c.tst[ , 1]) %>% data.frame() %>%
        dplyr::rename(., pred = `.`, obs = V2)
      
      comp <- rbind(alpha.df, beta.df, gamma.df)
    }
  }
  return(data.frame(
    data = deparse(substitute(data)), 
    fold = foldnum, 
    model = "PLS", 
    rmse = defaultSummary(comp)[1],
    rsquared = defaultSummary(comp)[2]
  ))
}

# Creates a dataframe of results from all folds
return.folds <- function(method, data, folds, ...) {
  df <- do.call(rbind, 
                lapply(1:10, FUN = method, data = data, 
                       folds = folds, ...))
  df$data <- deparse(substitute(data))
  rownames(df) <- NULL
  return(df)
}

# Data --------------------------------------------------------------------

# Regular data results (excludes GLM)
tic()
cube.data.folds <- return.folds(cube.fold, data, data.folds, suzuki = F)
toc() # 409.99 s, 468.06, 459.68 ... avg = 445.91

tic()
glm.data.folds <- return.folds(glm.fold, data, data.folds, suzuki = F)
toc() # 4.24, 4.39 s ... 4.314

tic()
pls.data.folds <- return.folds(pls.fold, data, data.folds, suzuki = F)
toc() # 167.27, 184.75 s ... 176.01

tic()
rf.data.folds <- return.folds(rf.fold, data, data.folds, suzuki = F)
toc() # 626.25, 645.35 ... 635.8

tic()
svm.data.folds <- return.folds(svm.fold, data, data.folds, suzuki = F)
toc() # 8.61 sec elapsed, 10.8 sec, 9.42 s, 9.51 ... avg = 9.59

data.results <- rbind(cube.data.folds, glm.data.folds, pls.data.folds, 
                      rf.data.folds, svm.data.folds)
qsar.type <- c("Cubist", "GLMNet", "PLS", "Random Forest", "SVM")
qsar.time <- c(445.91, 4.314, 176.01, 635.8, 9.59)
times <- data.frame(qsar.type, qsar.time)
saveRDS(data.results, "./compiled folds.RDS")

# Suzuki Results (stable)
tic()
svm.suz.folds <- return.folds(svm.fold, suz, suz.folds, suzuki = T)
toc() # 5.32

tic()
glm.suz.folds <- return.folds(glm.fold, suz, suz.folds, suzuki = T)
toc() # 2.88 s

tic()
cube.suz.folds <- return.folds(cube.fold, suz, suz.folds, suzuki = T)
toc() # 256.34 s

tic()
rf.suz.folds <- return.folds(rf.fold, suz, suz.folds, suzuki = T)
toc()

tic()
pls.suz.folds <- return.folds(pls.fold, suz, suz.folds, suzuki = T)
toc()

suz.results <- rbind(svm.suz.folds, glm.suz.folds, cube.suz.folds, 
                     rf.suz.folds, pls.suz.folds)

# Saving data
dir.create("./compiled")
saveRDS(data.results, "./compiled/ri and suzuki results.RDS")


# Graphs ------------------------------------------------------------------

#     R-squared -----------------------------------------------------------

# source("./model.code/graph.formatting.R")
ggplot(data.results, aes(x = fold, y = rsquared, color = model)) +
  theme.2018 + 
  geom_line(size = 1) +
  scale_x_continuous(breaks = seq(1,10, by = 1)) + 
  labs(title = "R-squared of QSARs over different data folds", 
       x = "Fold", y = "R-squared", 
       color = "Model") 
ggsave("./fold results.png")

ggplot(data.results, aes(x = fold, y = rsquared, color = model)) + 
  geom_line(size = 1) + 
  scale_x_discrete(expand=c(0, 1)) +
  scale_colour_discrete(guide = 'none') +
  geom_dl(aes(label = model), method = list(dl.trans(x = x + 0.2), "last.points", cex = 0.8)) +
  geom_dl(aes(label = model), method = list(dl.trans(x = x - 0.2), "first.points", cex = 0.8)) + 
  theme_bw() + 
  labs(title = "R-squared by Fold", x = "Fold", y = "R-squared")

ggplot(suz.results, aes(x = fold, y = rsquared, color = model)) + 
  geom_line(size = 1) + 
  scale_x_discrete(expand=c(0, 1)) +
  scale_colour_discrete(guide = 'none') +
  geom_dl(aes(label = model), method = list(dl.trans(x = x + 0.2), "last.points", cex = 0.8)) +
  geom_dl(aes(label = model), method = list(dl.trans(x = x - 0.2), "first.points", cex = 0.8)) + 
  theme_bw() + 
  labs(title = "R-squared by Fold - Suzuki", x = "Fold", y = "R-squared")

results.comp <- rbind(data.results %>% mutate(Source = "Rekharsky and Inoue/Suzuki"), 
                      suz.results %>% mutate(Source = "Suzuki"))
results.comp$fold <- as.factor(results.comp$fold)

ggplot(results.comp, aes(x = fold, y = rsquared, color = model, group = model)) + 
  geom_line(size = 1) + 
  # scale_x_discrete(expand=c(0, 1)) +
  # scale_colour_discrete(guide = 'none') +
  # geom_dl(aes(label = model), method = list(dl.trans(x = x + 0.2), "last.points", cex = 0.8)) +
  # geom_dl(aes(label = model), method = list(dl.trans(x = x - 0.2), "first.points", cex = 0.8)) + 
  theme_bw() + 
  facet_wrap(~Source) + 
  labs(title = "R-Squared by Fold", x = "Fold", y = "R-squared", color = "Model")
ggsave("./graphs/compiled.png")

# Messing with Different Orientations
temp <- suz.results
temp$fold <- as.factor(temp$fold)
ggplot(temp, aes(x = model, y = rsquared, color = fold, group = fold)) + 
  geom_line(size = 1) + 
  theme_bw() + 
  labs(title = "R-squared by Model - Suzuki", x = "Fold", y = "R-squared")

temp <- data.results
temp$fold <- as.factor(temp$fold)
ggplot(temp, aes(x = model, y = rsquared, color = fold, group = fold)) + 
  geom_line(size = 1) + 
  theme_bw() + 
  labs(title = "R-squared by Model - Data", x = "Fold", y = "R-squared")

#     RMSE ----------------------------------------------------------------

ggplot(results.comp, aes(x = fold, y = rmse, color = model, group = model)) + 
  geom_line(size = 1) + 
  theme_bw() + 
  facet_wrap(~Source) + 
  labs(title = "RMSE by Fold", x = "Fold", y = "R-squared", color = "Model")
ggsave("./graphs/compiled rmse.png")


#     Time ----------------------------------------------------------------

ggplot(times, aes(x = qsar.type, y = qsar.time, fill = qsar.type)) + 
  geom_bar(stat = "identity") + 
  theme.isef + 
  labs(title = NULL, y = "Time, seconds", 
       x = "QSAR type") +
  guides(fill = F) + 
  coord_fixed(ratio = 0.005)
ggsave("./graphs/isef qsar calc time.png", dpi = 450, scale = 1.5)

# External Validation Compilation -----------------------------------------

View(ev.glm)
ev.comp <- rbind(ev.glm, ev.svm, ev.rf, ev.cube, ev.pls)
ggplot(ev.comp, aes(x = obs, y = pred, color = cd.type)) + 
  geom_point() + 
  theme_bw() + 
  facet_wrap(~Model) + 
  geom_abline(intercept = 0, slope = 1) + 
  coord_fixed(xlim = c(-45, 5), ylim = c(-45, 5)) + 
  labs(title = "External Validation", 
       x = "Observed DelG, kJ/mol", y = "Predicted DelG, kJ/mol", 
       color = "Cyclodextrin")
