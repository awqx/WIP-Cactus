# Libraries ---------------------------------------------------------------

if(!require("pacman")) { 
  install.packages("pacman")
  library(pacman)
} else
  library(pacman)

p_load(caret,
       Cubist,
       data.table,
       e1071,
       earth,
       gbm,
       glmnet,
       Matrix,
       pls, 
       randomForest,
       stringr,
       tidyverse)

# Cubist ------------------------------------------------------------------

# Cubist should be fed with matrix or data.frame
tune.cubist.cmte <- function(data, nfolds, cmte, seed) {
  set.seed(seed)
  fold.list <- createFolds(y = data[ , 1], k = nfolds)
  r2.results <- c(rep(0.0, nfolds))
  rmse.results <- c(rep(0.0, nfolds))
  
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
    
    rmse.results[i] <- defaultSummary(cube.df)[1]
    r2.results[i] <- defaultSummary(cube.df)[2]
    
  }
  message(nfolds, "-fold cross-validation of ", cmte, " committees completed.")
  return(data.frame( 
    nfolds = nfolds,
    seed = seed,
    committees = cmte,
    rsquared = sum(r2.results)/nfolds, 
    rmse = sum(rmse.results)/nfolds))
}

tune.cubist.samp <- function(data, nfolds, samp, seed) {
  set.seed(seed)
  fold.list <- createFolds(y = data[ , 1], k = nfolds)
  r2.results <- c(rep(0.0, nfolds))
  rmse.results <- c(rep(0.0, nfolds))
  
  ctrl <- cubistControl(
    seed = seed, 
    sample = samp
  )
  
  for(i in 1:nfolds) {
    fold <- fold.list[[i]]
    
    trn.x <- data[-fold, -1]
    trn.y <- data[-fold, 1]
    tst.x <- data[fold, -1]
    tst.y <- data[fold, 1]
    
    cube <- cubist(x = trn.x, y = trn.y, control = ctrl)
    cube.df <- predict(cube, tst.x) %>%
      cbind(tst.y) %>%
      data.frame() %>%
      .[complete.cases(.), ]
    # for some reason, tuning samp causes NAs, so those cases 
    # must be removed before caret::defaultSummary
    
    colnames(cube.df)[1] <- "pred"
    colnames(cube.df)[2] <- "obs"
    
    rmse.results[i] <- defaultSummary(cube.df)[1]
    r2.results[i] <- defaultSummary(cube.df)[2]
  }
  
  message(nfolds, "-fold cross-validation of ", samp, "% sample size completed.")
  return(data.frame(
    nfolds = nfolds,
    seed = seed,
    samp = samp,
    rsquared = sum(r2.results)/nfolds, 
    rmse = sum(rmse.results)/nfolds))
}

tune.cubist.extra <- function(data, nfolds, extra, seed) {
  set.seed(seed)
  fold.list <- createFolds(y = data[ , 1], k = nfolds)
  r2.results <- c(rep(0.0, nfolds))
  rmse.results <- c(rep(0.0, nfolds))
  
  ctrl <- cubistControl(
    seed = seed, 
    extrapolation = extra
  )
  
  for(i in 1:nfolds) {
    fold <- fold.list[[i]]
    
    trn.x <- data[-fold, -1]
    trn.y <- data[-fold, 1]
    tst.x <- data[fold, -1]
    tst.y <- data[fold, 1]
    
    cube <- cubist(x = trn.x, y = trn.y, control = ctrl)
    cube.df <- predict(cube, tst.x) %>%
      cbind(tst.y) %>%
      data.frame() 
    
    colnames(cube.df)[1] <- "pred"
    colnames(cube.df)[2] <- "obs"
    
    rmse.results[i] <- defaultSummary(cube.df)[1]
    r2.results[i] <- defaultSummary(cube.df)[2]
  }
  message(nfolds, "-fold cross-validation of ", extra, " extrapolation completed.")
  return(data.frame(
    nfolds = nfolds,
    seed = seed,
    extrapolation = extra,
    rsquared = sum(r2.results)/nfolds, 
    rmse = sum(rmse.results)/nfolds))
}

# cubist is the only tuning function with a set.seed function
# because the entire thing needs a constant reminder to set.seed
tune.cubist <- function(data, nfolds, cmte, samp, extra, seed) {
  set.seed(seed)
  fold.list <- createFolds(y = data[ , 1], k = nfolds)
  r2.results <- c(rep(0.0, nfolds))
  rmse.results <- c(rep(0.0, nfolds))
  
  ctrl <- cubistControl(
    seed = seed, 
    extrapolation = extra,
    sample = samp
  )
  
  for(i in 1:nfolds) {
    fold <- fold.list[[i]]
    
    trn.x <- data[-fold, -1]
    trn.y <- data[-fold, 1]
    tst.x <- data[fold, -1]
    tst.y <- data[fold, 1]
    
    cube <- cubist(x = trn.x, y = trn.y,  
                   committees = cmte, 
                   control = ctrl)
    cube.df <- predict(cube, tst.x) %>%
      cbind(tst.y) %>%
      data.frame() 
    
    colnames(cube.df)[1] <- "pred"
    colnames(cube.df)[2] <- "obs"
    
    rmse.results[i] <- defaultSummary(cube.df)[1]
    r2.results[i] <- defaultSummary(cube.df)[2]
  }
  # message(nfolds, "-fold cross-validation of ", cmte, " committees, ",
  #         samp, "% sample size, ",
  #         extra, " extrapolation completed.")
  return(data.frame(
    nfolds = nfolds,
    seed = seed,
    committees = cmte, 
    sample = samp, 
    extrapolation = extra,
    rsquared = sum(r2.results)/nfolds, 
    rmse = sum(rmse.results)/nfolds))
}

tune.cubist <- function(data, nfolds, cmte, extra, seed) {
  set.seed(seed)
  fold.list <- createFolds(y = data[ , 1], k = nfolds)
  r2.results <- c(rep(0.0, nfolds))
  rmse.results <- c(rep(0.0, nfolds))
  ctrl <- cubistControl(
    seed = seed, 
    extrapolation = extra
  )
  
  for(i in 1:nfolds) {
    fold <- fold.list[[i]]
    
    trn.x <- data[-fold, -1]
    trn.y <- data[-fold, 1]
    tst.x <- data[fold, -1]
    tst.y <- data[fold, 1]
    
    cube <- cubist(x = trn.x, y = trn.y,  
                   committees = cmte, 
                   control = ctrl)
    cube.df <- predict(cube, tst.x) %>%
      cbind(tst.y) %>%
      data.frame() 
    
    colnames(cube.df)[1] <- "pred"
    colnames(cube.df)[2] <- "obs"
    
    rmse.results[i] <- defaultSummary(cube.df)[1]
    r2.results[i] <- defaultSummary(cube.df)[2]
  }
  return(data.frame(
    nfolds = nfolds,
    seed = seed,
    committees = cmte, 
    extrapolation = extra,
    rsquared = sum(r2.results, na.rm = T)/nfolds, 
    rmse = sum(rmse.results, na.rm = T)/nfolds))
}




# GBM ---------------------------------------------------------------------

tune.gbm.ntree <- function(data, nfolds, num, seed) {
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
    gb <- gbm.fit(x = x, 
                  y = y,
                  n.trees = num,
                  distribution = 'gaussian', 
                  verbose = F
    )
    gbm.df <- predict(gb, tst[ , -1], 
                      n.trees = num) %>% 
      cbind(tst[ , 1]) %>%
      data.frame()
    # Necessary to rename columns so that caret:defaultSummary works
    colnames(gbm.df)[1] <- "pred"
    colnames(gbm.df)[2] <- "obs"
    r2.results[i] <- defaultSummary(gbm.df)[2] 
    rmse.results[i] <- defaultSummary(gbm.df)[1]
  }
  message(nfolds, "-fold cross-validation of ntree = ", num, " completed.")
  return(data.frame( # Useful for records
    data = deparse(substitute(data)), # Turns the variable name into char
    nfolds = nfolds,
    seed = seed,
    ntree = num,
    rsquared = sum(r2.results) / nfolds, 
    rmse = sum(rmse.results)/nfolds))
}

# Using 500 trees
tune.gbm.shrink <- function(data, nfolds, d, seed) {
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
    gb <- gbm.fit(x = x, 
                  y = y,
                  interaction.depth = d,
                  n.trees = 500,
                  distribution = 'gaussian', 
                  verbose = F
    )
    gbm.df <- predict(gb, tst[ , -1], n.trees = 500) %>% 
      cbind(tst[ , 1]) %>%
      data.frame()
    # Necessary to rename columns so that caret:defaultSummary works
    colnames(gbm.df)[1] <- "pred"
    colnames(gbm.df)[2] <- "obs"
    r2.results[i] <- defaultSummary(gbm.df)[2] 
    rmse.results[i] <- defaultSummary(gbm.df)[1]
  }
  message(nfolds, "-fold cross-validation of depth = ", d, " completed.")
  return(data.frame( # Useful for records
    data = deparse(substitute(data)), # Turns the variable name into char
    nfolds = nfolds,
    seed = seed,
    depth = d,
    rsquared = sum(r2.results) / nfolds, 
    rmse = sum(rmse.results)/nfolds))
}

tune.gbm.shrink <- function(data, nfolds, s, seed) {
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
    gb <- gbm.fit(x = x, 
                  y = y,
                  shrinkage = s,
                  n.trees = 500,
                  distribution = 'gaussian', 
                  verbose = F
    )
    gbm.df <- predict(gb, tst[ , -1], n.trees = 500) %>% 
      cbind(tst[ , 1]) %>%
      data.frame()
    # Necessary to rename columns so that caret:defaultSummary works
    colnames(gbm.df)[1] <- "pred"
    colnames(gbm.df)[2] <- "obs"
    r2.results[i] <- defaultSummary(gbm.df)[2] 
    rmse.results[i] <- defaultSummary(gbm.df)[1]
  }
  message(nfolds, "-fold cross-validation of shrinkage = ", s, " completed.")
  return(data.frame( # Useful for records
    data = deparse(substitute(data)), # Turns the variable name into char
    nfolds = nfolds,
    seed = seed,
    shrinkage = s,
    rsquared = sum(r2.results) / nfolds, 
    rmse = sum(rmse.results)/nfolds))
}

tune.gbm.node <- function(data, nfolds, n, seed) {
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
    gb <- gbm.fit(x = x, 
                  y = y,
                  n.minobsinnode = n,
                  n.trees = 500,
                  distribution = 'gaussian', 
                  verbose = F
    )
    gbm.df <- predict(gb, tst[ , -1], n.trees = 500) %>% 
      cbind(tst[ , 1]) %>%
      data.frame()
    # Necessary to rename columns so that caret:defaultSummary works
    colnames(gbm.df)[1] <- "pred"
    colnames(gbm.df)[2] <- "obs"
    r2.results[i] <- defaultSummary(gbm.df)[2] 
    rmse.results[i] <- defaultSummary(gbm.df)[1]
  }
  message(nfolds, "-fold cross-validation of node observations = ", n, " completed.")
  return(data.frame( # Useful for records
    data = deparse(substitute(data)), # Turns the variable name into char
    nfolds = nfolds,
    seed = seed,
    node = n, 
    rsquared = sum(r2.results) / nfolds, 
    rmse = sum(rmse.results)/nfolds))
}

# num: n.trees
# d: interactions.depth
# s: shrinkage
# n: n.minobsinnode
tune.gbm <- function(data, nfolds, num, d, s, n) {
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
    gb <- gbm.fit(x = x, 
                  y = y,
                  n.trees = num, 
                  interaction.depth = d, 
                  shrinkage = s,
                  n.minobsinnode = n,
                  distribution = 'gaussian', 
                  verbose = F
    )
    gbm.df <- predict(gb, tst[ , -1],
                      n.trees = num) %>% 
      cbind(tst[ , 1]) %>%
      data.frame()
    # Necessary to rename columns so that caret:defaultSummary works
    colnames(gbm.df)[1] <- "pred"
    colnames(gbm.df)[2] <- "obs"
    r2.results[i] <- defaultSummary(gbm.df)[2] 
    rmse.results[i] <- defaultSummary(gbm.df)[1]
  }
  return(data.frame(
    nfolds = nfolds,
    ntree = num, 
    depth = d, 
    shrinkage = s,
    node = n, 
    rsquared = sum(r2.results) / nfolds, 
    rmse = sum(rmse.results)/nfolds))
}


# GLMnet ------------------------------------------------------------------

# # For a matrix 
tune.glm.alpha <- function(data, nfolds, alpha, seed) {
  set.seed(seed)
  fold.list <- createFolds(y = data[ , 1], k = nfolds)
  r2.results <- c(rep(0.0, nfolds))
  rmse.results <- c(rep(0.0, nfolds))
  
  for(i in 1:nfolds) {
    fold <- fold.list[[i]]
    
    trn.x <- data[-fold, -1]
    trn.y <- data[-fold, 1]
    tst.x <- data[fold, -1]
    tst.y <- data[fold, 1]
    
    glm.mod <- glmnet(x = trn.x, y = trn.y,
                      alpha = alpha,
                      family = "mgaussian")
    glm.df <- predict.glmnet(glm.mod, tst.x,
                             s = tail(glm.mod$lambda, n = 1)) %>%
      cbind(tst.y) %>%
      data.frame() %>%
      rename(pred = X1, obs = tst.y)
    
    rmse.results[i] <- defaultSummary(glm.df)[1]
    r2.results[i] <- defaultSummary(glm.df)[2]
  }
  
  return(data.frame(
    nfolds = nfolds,
    seed = seed,
    alpha = alpha,
    rsquared = sum(r2.results)/nfolds,
    rmse = sum(rmse.results)/nfolds))
}

# tune.glm.alpha <- function(data, nfolds, alpha, seed) {
#   set.seed(seed)
#   fold.list <- createFolds(y = data[ , 2], k = nfolds)
#   r2.results <- c(rep(0.0, nfolds))
#   rmse.results <- c(rep(0.0, nfolds))
#   
#   for(i in 1:nfolds) {
#     fold <- fold.list[[i]]
#     
#     trn.x <- data[-fold, -1:-2]
#     trn.y <- data[-fold, 2]
#     tst.x <- data[fold, -1:-2]
#     tst.y <- data[fold, 2]
#     
#     glm.mod <- glmnet(x = trn.x, y = trn.y, 
#                       alpha = alpha, 
#                       family = "mgaussian")
#     glm.df <- predict.glmnet(glm.mod, tst.x,
#                              s = tail(glm.mod$lambda, n = 1)) %>%
#       cbind(tst.y) %>%
#       data.frame() %>%
#       rename(pred = X1, obs = tst.y)
#     
#     rmse.results[i] <- defaultSummary(glm.df)[1]
#     r2.results[i] <- defaultSummary(glm.df)[2]
#   }
#   
#   return(data.frame( 
#     nfolds = nfolds,
#     seed = seed,
#     alpha = alpha,
#     rsquared = sum(r2.results)/nfolds, 
#     rmse = sum(rmse.results)/nfolds))
# }


tune.glm.dfmax <- function(data, nfolds, max, seed) {
  set.seed(seed)
  fold.list <- createFolds(y = data[ , 1], k = nfolds)
  r2.results <- c(rep(0.0, nfolds))
  rmse.results <- c(rep(0.0, nfolds))
  
  for(i in 1:nfolds) {
    fold <- fold.list[[i]]
    
    trn.x <- data[-fold, -1]
    trn.y <- data[-fold, 1]
    tst.x <- data[fold, -1]
    tst.y <- data[fold, 1]
    
    glm.mod <- glmnet(x = trn.x, y = trn.y, 
                      dfmax = max, 
                      family = "mgaussian")
    glm.df <- predict.glmnet(glm.mod, tst.x,
                             s = tail(glm.mod$lambda, n = 1)) %>%
      cbind(tst.y) %>%
      data.frame() %>%
      rename(pred = X1, obs = tst.y)
    
    rmse.results[i] <- defaultSummary(glm.df)[1]
    r2.results[i] <- defaultSummary(glm.df)[2]
  }
  
  return(data.frame( 
    nfolds = nfolds,
    seed = seed,
    dfmax = max,
    rsquared = sum(r2.results)/nfolds, 
    rmse = sum(rmse.results)/nfolds))
}

tune.glm <- function(data, nfolds, a, max) {
  fold.list <- createFolds(y = data[ , 1], k = nfolds)
  r2.results <- c(rep(0.0, nfolds))
  rmse.results <- c(rep(0.0, nfolds))
  
  for(i in 1:nfolds) {
    fold <- fold.list[[i]]
    
    trn.x <- data[-fold, -1]
    trn.y <- data[-fold, 1]
    tst.x <- data[fold, -1]
    tst.y <- data[fold, 1]
    
    glm.mod <- glmnet(x = trn.x, y = trn.y, 
                      dfmax = max, alpha = a,
                      pmax = ncol(trn.x), 
                      family = "mgaussian")
    glm.df <- predict.glmnet(glm.mod, tst.x,
                             s = tail(glm.mod$lambda, n = 1)) %>%
      cbind(tst.y) %>%
      data.frame() %>%
      rename(pred = X1, obs = tst.y)
    
    rmse.results[i] <- defaultSummary(glm.df)[1]
    r2.results[i] <- defaultSummary(glm.df)[2]
  }
  
  return(data.frame( 
    nfolds = nfolds,
    alpha = a,
    dfmax = max,
    rsquared = sum(r2.results)/nfolds, 
    rmse = sum(rmse.results)/nfolds))
}

# MARS --------------------------------------------------------------------

tune.mars.deg <- function(data, nfolds, deg, seed) {
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
    mars <- earth(x = x, y = y, 
                  degree = deg)
    mars.df <- predict(mars, tst[ , -1]) %>% 
      cbind(tst[ , 1]) %>%
      data.frame()
    # Necessary to rename columns so that caret:defaultSummary works
    colnames(mars.df) <- c('pred', 'obs')
    r2.results[i] <- defaultSummary(mars.df)[2] 
    rmse.results[i] <- defaultSummary(mars.df)[1]
  }
  return(data.frame( 
    nfolds = nfolds,
    seed = seed,
    degree = deg,
    rsquared = sum(r2.results) / nfolds, 
    rmse = sum(rmse.results)/nfolds))
}

tune.mars.pen <- function(data, nfolds, pen, seed) {
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
    mars <- earth(x = x, y = y, 
                  penalty = pen)
    mars.df <- predict(mars, tst[ , -1]) %>% 
      cbind(tst[ , 1]) %>%
      data.frame()
    # Necessary to rename columns so that caret:defaultSummary works
    colnames(mars.df) <- c('pred', 'obs')
    r2.results[i] <- defaultSummary(mars.df)[2] 
    rmse.results[i] <- defaultSummary(mars.df)[1]
  }
  return(data.frame( 
    nfolds = nfolds,
    seed = seed,
    penalty = pen,
    rsquared = sum(r2.results) / nfolds, 
    rmse = sum(rmse.results)/nfolds))
}

tune.mars.nk <- function(data, nfolds, nk, seed) {
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
    mars <- earth(x = x, y = y, 
                  nk = nk)
    mars.df <- predict(mars, tst[ , -1]) %>% 
      cbind(tst[ , 1]) %>%
      data.frame()
    # Necessary to rename columns so that caret:defaultSummary works
    colnames(mars.df) <- c('pred', 'obs')
    r2.results[i] <- defaultSummary(mars.df)[2] 
    rmse.results[i] <- defaultSummary(mars.df)[1]
  }
  return(data.frame( 
    nfolds = nfolds,
    seed = seed,
    nk = nk,
    rsquared = sum(r2.results) / nfolds, 
    rmse = sum(rmse.results)/nfolds))
}

tune.mars.thresh <- function(data, nfolds, thresh, seed) {
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
    mars <- earth(x = x, y = y, 
                  thresh = thresh)
    mars.df <- predict(mars, tst[ , -1]) %>% 
      cbind(tst[ , 1]) %>%
      data.frame()
    colnames(mars.df) <- c('pred', 'obs')
    r2.results[i] <- defaultSummary(mars.df)[2] 
    rmse.results[i] <- defaultSummary(mars.df)[1]
  }
  return(data.frame( 
    nfolds = nfolds,
    seed = seed,
    thresh = thresh,
    rsquared = sum(r2.results) / nfolds, 
    rmse = sum(rmse.results)/nfolds))
}

tune.mars.minspan <- function(data, nfolds, minspan, seed) {
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
    mars <- earth(x = x, y = y, 
                  minspan = minspan)
    mars.df <- predict(mars, tst[ , -1]) %>% 
      cbind(tst[ , 1]) %>%
      data.frame()
    colnames(mars.df) <- c('pred', 'obs')
    r2.results[i] <- defaultSummary(mars.df)[2] 
    rmse.results[i] <- defaultSummary(mars.df)[1]
  }
  return(data.frame( 
    nfolds = nfolds,
    seed = seed,
    minspan = minspan,
    rsquared = sum(r2.results) / nfolds, 
    rmse = sum(rmse.results)/nfolds))
}

tune.mars.fk <- function(data, nfolds, fk, seed) {
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
    mars <- earth(x = x, y = y, 
                  fast.k = fk)
    mars.df <- predict(mars, tst[ , -1]) %>% 
      cbind(tst[ , 1]) %>%
      data.frame()
    colnames(mars.df) <- c('pred', 'obs')
    r2.results[i] <- defaultSummary(mars.df)[2] 
    rmse.results[i] <- defaultSummary(mars.df)[1]
  }
  return(data.frame( 
    nfolds = nfolds,
    seed = seed,
    fast.k = fk,
    rsquared = sum(r2.results) / nfolds, 
    rmse = sum(rmse.results)/nfolds))
}

tune.mars <- function(data, nfolds, d, p, nk, t, m, fk) {
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
    mars <- earth(x = x, y = y,
                  degree = d, 
                  penalty = p, 
                  nk = nk, 
                  thresh = t, 
                  minspan = m,
                  fast.k = fk)
    mars.df <- predict(mars, tst[ , -1]) %>% 
      cbind(tst[ , 1]) %>%
      data.frame()
    colnames(mars.df) <- c('pred', 'obs')
    r2.results[i] <- defaultSummary(mars.df)[2] 
    rmse.results[i] <- defaultSummary(mars.df)[1]
  }
  return(data.frame( 
    nfolds = nfolds,
    degree = d, 
    penalty = p, 
    nk = nk, 
    thresh = t,
    minspan = m,
    fast.k = fk,
    rsquared = sum(r2.results) / nfolds, 
    rmse = sum(rmse.results)/nfolds))
}

# PLS ---------------------------------------------------------------------

tune.pls.ncomp <- function(data, nfolds, comp, seed) {
  set.seed(seed)
  fold.list <- createFolds(y = data[ , 1], k = nfolds)
  r2.results <- c(rep(0.0, nfolds))
  rmse.results <- c(rep(0.0, nfolds))
  
  for(i in 1:nfolds) {
    fold <- fold.list[[i]]
    
    trn <- data[-fold, ]
    tst <- data[fold, ]
    
    pls <- plsr(DelG ~., data = trn, 
                ncomp = comp)
    pls.df <- predict(pls, tst[ , -1]) %>%
      cbind(tst[ , 1]) %>%
      data.frame() 
    
    colnames(pls.df)[1] <- "pred"
    colnames(pls.df)[2] <- "obs"
    
    r2.results[i] <- defaultSummary(pls.df)[2]
    rmse.results[i] <- defaultSummary(pls.df)[1]
  }
  # message(comp)
  
  return(data.frame( # Useful for records
    data = deparse(substitute(data)), # Turns the variable name into char
    nfolds = nfolds,
    ncomp = comp,
    rsquared = sum(r2.results) / nfolds, 
    rmse = sum(rmse.results)/nfolds,
    seed = seed))
}

# Method = "kernelpls", "widekernelpls", "simpls", or "oscorespls"
tune.pls.method <- function(data, nfolds, met, seed) {
  set.seed(seed)
  fold.list <- createFolds(y = data[ , 1], k = nfolds)
  r2.results <- c(rep(0.0, nfolds))
  rmse.results <- c(rep(0.0), nfolds)
  
  for(i in 1:nfolds) {
    fold <- fold.list[[i]]
    
    trn <- data[-fold, ]
    tst <- data[fold, ]
    
    pls <- plsr(DelG ~., data = trn, 
                method = met #, ncomp = 11
    )
    pls.df <- predict(pls, tst[ , -1]) %>%
      cbind(tst[ , 1]) %>%
      data.frame() 
    
    colnames(pls.df)[1] <- "pred"
    colnames(pls.df)[2] <- "obs"
    
    rmse.results[i] <- defaultSummary(pls.df)[1]
    r2.results[i] <- defaultSummary(pls.df)[2]
  }
  # I have a message for this particular function in order to keep better
  # track, and to not stop the function too early
  message(met, " has been completed")
  return(data.frame( # Useful for records
    data = deparse(substitute(data)), # Turns the variable name into char
    nfolds = nfolds,
    method = met,
    rsquared = sum(r2.results)/nfolds, 
    rmse = sum(rmse.results)/nfolds,
    seed = seed)) 
}

tune.pls <- function(data, nfolds, comp, met) {
  fold.list <- createFolds(y = data[ , 1], k = nfolds)
  r2.results <- c(rep(0.0, nfolds))
  rmse.results <- c(rep(0.0, nfolds))
  
  for(i in 1:nfolds) {
    fold <- fold.list[[i]]
    
    trn <- data[-fold, ]
    tst <- data[fold, ]
    
    pls.formula <- paste0(colnames(data)[1], ' ~ .') %>% 
      as.formula()
    pls <- plsr(pls.formula, data = trn, 
                ncomp = comp, method = met)
    pls.df <- predict(pls, tst[ , -1]) %>%
      cbind(tst[ , 1]) %>%
      data.frame() 
    
    colnames(pls.df)[1] <- "pred"
    colnames(pls.df)[2] <- "obs"
    
    r2.results[i] <- defaultSummary(pls.df)[2]
    rmse.results[i] <- defaultSummary(pls.df)[1]
  }
  
  return(data.frame( 
    nfolds = nfolds,
    ncomp = comp,
    method = met,
    rsquared = sum(r2.results) / nfolds, 
    rmse = sum(rmse.results)/nfolds
  ))
}

# Random forest -----------------------------------------------------------

tune.rf.ntree <- function(data, nfolds, num, seed) {
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
    rf <- randomForest(x = x, 
                       y = y,
                       ntree = num, 
                       seed = seed, 
                       importance = T)
    rf.df <- predict(rf, tst[ , -1]) %>% 
      cbind(tst[ , 1]) %>%
      data.frame()
    # Necessary to rename columns so that caret:defaultSummary works
    colnames(rf.df)[1] <- "pred"
    colnames(rf.df)[2] <- "obs"
    r2.results[i] <- defaultSummary(rf.df)[2] 
    rmse.results[i] <- defaultSummary(rf.df)[1]
  }
  message(nfolds, "-fold cross-validation of ntree = ", num, " completed.")
  return(data.frame( # Useful for records
    data = deparse(substitute(data)), # Turns the variable name into char
    nfolds = nfolds,
    seed = seed,
    ntree = num,
    rsquared = sum(r2.results) / nfolds, 
    rmse = sum(rmse.results)/nfolds))
}

tune.rf.node <- function(data, nfolds, node, seed) {
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
    rf <- randomForest(x = x, 
                       y = y,
                       nodesize = node, 
                       seed = seed, 
                       importance = T)
    rf.df <- predict(rf, tst[ , -1]) %>% 
      cbind(tst[ , 1]) %>%
      data.frame()
    # Necessary to rename columns so that caret:defaultSummary works
    colnames(rf.df)[1] <- "pred"
    colnames(rf.df)[2] <- "obs"
    rmse.results[i] <- defaultSummary(rf.df)[1]
    r2.results[i] <- defaultSummary(rf.df)[2]
  }
  message(nfolds, "-fold cross-validation of nodesize = ", node, " completed.")
  return(data.frame(
    nfolds = nfolds, 
    seed = seed,
    nodesize = node,
    rsquared = sum(r2.results)/nfolds, 
    rmse = sum(rmse.results)/nfolds))
}

tune.rf.mtry <- function(data, nfolds, m, seed) {
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
    rf <- randomForest(x = x, 
                       y = y,
                       mtry = m, 
                       seed = seed, 
                       importance = T)
    rf.df <- predict(rf, tst[ , -1]) %>% 
      cbind(tst[ , 1]) %>%
      data.frame()
    # Necessary to rename columns so that caret:defaultSummary works
    colnames(rf.df)[1] <- "pred"
    colnames(rf.df)[2] <- "obs"
    rmse.results[i] <- defaultSummary(rf.df)[1]
    r2.results[i] <- defaultSummary(rf.df)[2]
  }
  message(nfolds, "-fold cross-validation of mtry = ", m, " completed.")
  return(data.frame( 
    nfolds = nfolds, 
    seed = seed,
    mtry = m,
    rsquared = sum(r2.results)/nfolds, 
    rmse = sum(rmse.results)/nfolds))
}

tune.rf <- function(data, nfolds, ntree, node, m) {
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
    rf <- randomForest(x = x, 
                       y = y,
                       ntree = ntree,
                       nodesize = node,
                       mtry = m, 
                       seed = seed, 
                       importance = T)
    rf.df <- predict(rf, tst[ , -1]) %>% 
      cbind(tst[ , 1]) %>%
      data.frame()
    # Necessary to rename columns so that caret:defaultSummary works
    colnames(rf.df)[1] <- "pred"
    colnames(rf.df)[2] <- "obs"
    rmse.results[i] <- defaultSummary(rf.df)[1]
    r2.results[i] <- defaultSummary(rf.df)[2]
  }
  # message(nfolds, "-fold cross-validation of ntree = ", ntree,
  #         ", nodesize = ", node,
  #         ", mtry = ", m, 
  #         " completed.")
  return(data.frame( 
    nfolds = nfolds, 
    ntree = ntree, 
    nodesize = node,
    mtry = m,
    rsquared = sum(r2.results)/nfolds, 
    rmse = sum(rmse.results)/nfolds))
}

# SVM ---------------------------------------------------------------------

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
  
  # message("Cost: ", cost, " || Gamma: ", g,
  #         " || Epsilon: ", e, " || Coef0: ", coef, 
  #         " || Degree: ", deg)
  
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
  #         " || Gamma: ", g)
  
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





