# Libraries and Packages --------------------------------------------------

library(caret)
library(randomForest)
library(tidyverse)

# Data Organization -------------------------------------------------------

rpt <- readRDS("./rpt.RDS")
mat.dg <- rpt %>% dplyr::select(., -X1:-log.K.Uncertainty, -DelG.Uncertainty:-`bind.aff, kcal/mol`)

set.seed(48)
trn.ind <- sample(x = 1:nrow(mat.dg), size = round(0.8 * nrow(mat.dg)))
trn <- mat.dg[trn.ind, ]
tst <- mat.dg[-trn.ind, ]
mat.imp <- rfImpute(x = mat.dg[ , -1], y = mat.dg[ , 1], ntree = 400, na.action = na.omit)

# Tuning ------------------------------------------------------------------

rf.param <-
  tuneRF(
    x = rf.imp[,-1],
    y = rf.imp[, 1],
    ntreeTry = 400,
    improve = 0.01,
    stepFactor = 1.5
  )
# 741, 494, 555, ... very variable
rf.param2 <-
  tune.randomForest(
    x = rf.imp[,-1],
    y = rf.imp[, 1],
    mtry = 100 * 1.5 ^ (-1:3),
    nodesize = 30 * 2 ^ (-1:2)
  )
# nodesize = 15, mtry = 337
# getting good results with mtry = 555
rf.param3 <-
  tune.randomForest(
    x = rf.imp[,-1],
    y = rf.imp[, 1],
    mtry = 200 * 2 ^ (0:3),
    nodesize = 10 * 2 ^ (-1:2),
    performances = T
  )
# nodesize = 10, mtry = 800 (but 'invalid)
rf.param4 <- tune.randomForest(
  x = mat.imp[ , -1], 
  y = mat.imp[ , 1], 
  mtry = 555, 
  nodesize = 10 ^ (-1:2), 
  ntree = 100 * 2 ^ (0:3)
)

rf.param4 <- tune.randomForest(
  x = mat.imp[ , -1], 
  y = mat.imp[ , 1], 
  mtry = 555, 
  nodesize = 10 ^ (-1:2), 
  ntree = 100 * 2 ^ (0:3)
)
# ntree = 100, nodesize = 0.1

rf.param5 <- tune.randomForest(
  x = mat.imp[ , -1], 
  y = mat.imp[ , 1], 
  mtry = 555, 
  nodesize = 0.1, 
  ntree = 100 * (1:6)
)
# Random Forest Model -----------------------------------------------------

rf.x <- trn[ , -1]
rf.y <- trn[ , 1]
rf.imp <- rfImpute(x = rf.x, y = rf.y, ntree = 400, na.action = na.omit)
rf <- randomForest(DelG ~., data = trn, ntree = 500, na.action = na.omit, mtry = 555, nodesize = 15)
rf.tst.df <- predict(rf, tst[ , -1]) %>%
  cbind(tst[ , 1]) %>% 
  data.frame()
colnames(rf.tst.df)[1] <- "pred"
colnames(rf.tst.df)[2] <- "obs"
rf.tst.df %>% ggplot(., aes(x = obs, y = pred)) +
  geom_point() + 
  geom_abline(slope = 1, intercept = 0)
defaultSummary(rf.tst.df)

rf.trn.df <- predict(rf, trn[ , -1]) %>%
  cbind(trn[ , 1]) %>% 
  data.frame()
colnames(rf.trn.df)[1] <- "pred"
colnames(rf.trn.df)[2] <- "obs"
rf.trn.df %>% ggplot(., aes(x = obs, y = pred)) +
  geom_point() + 
  geom_abline(slope = 1, intercept = 0)
defaultSummary(rf.trn.df)

rf.trn.df <- rf.trn.df %>%
  mutate(trn.resid = obs - pred)

rf.trn.df %>% ggplot(., aes(x = obs, y = trn.resid)) + 
  geom_point() + 
  geom_hline(yintercept = 0)

rf.tst.df <- rf.tst.df %>%
  mutate(tst.resid = obs - pred)

ggplot(rf.tst.df, aes(x = obs, y = tst.resid)) + 
  geom_point() + 
  geom_hline(yintercept = 0)

# Cross-validation --------------------------------------------------------
set.seed(256)
rf.folds <- createFolds(y = mat.imp[ , 1], k = 10, list = T)
find.rf.cv(mat.imp, rf.folds, 1, 555, 400, 0.1) # 0.63
find.rf.cv(mat.imp, rf.folds, 2, 555, 400, 0.1) # 
find.rf.cv(mat.imp, rf.folds, 3, 555, 400, 0.1)
find.rf.cv(mat.imp, rf.folds, 4, 555, 400, 0.1)
find.rf.cv(mat.imp, rf.folds, 5, 555, 400, 0.1)
find.rf.cv(mat.imp, rf.folds, 6, 555, 400, 0.1)
find.rf.cv(mat.imp, rf.folds, 7, 555, 400, 0.1)
find.rf.cv(mat.imp, rf.folds, 8, 555, 400, 0.1)
find.rf.cv(mat.imp, rf.folds, 9, 555, 400, 0.1)
find.rf.cv(mat.imp, rf.folds, 10, 555, 400, 0.1)

find.rf.cv2(mat.imp, rf.folds, 555, 400, 0.1) # 0.483
find.rf.cv2(mat.imp, rf.folds, 555, 400, 10) # 0.467
find.rf.cv2(mat.imp, rf.folds, 350, 500, 0.1) # 0.476
find.rf.cv2(mat.imp, rf.folds, 350, 400, 0.1) # 0.476
find.rf.cv2(mat.imp, rf.folds, 350, 300, 0.1) # 0.482
find.rf.cv2(mat.imp, rf.folds, 350, 200, 0.1) # 0.470
find.rf.cv2(mat.imp, rf.folds, 350, 100, 0.1) # 0.459
find.rf.cv2(mat.imp, rf.folds, 350, 400, 0.1) # 0.468
find.rf.cv2(mat.imp, rf.folds, 350, 400, 0.25) # 0.468
find.rf.cv2(mat.imp, rf.folds, 350, 200, 0.5) # 0.467
find.rf.cv2(mat.imp, rf.folds, 350, 400, 5) # 0.470

find.rf.cv <- function(data, fold.list, fold.num, m, n, node) {
  fold <- fold.list[[fold.num]]
  train = data[-fold, ]
  tst = data[fold, ]
  x = train[ , -1]
  y = train[ , 1]
  rf.cv <- randomForest(x = x, y = y, mtry = m,
                        ntree = n, nodesize = node)
  rf.df <- predict(rf.cv, tst[ , -1]) %>% 
    cbind(tst[ , 1]) %>%
    data.frame()
  colnames(rf.df)[1] <- "pred"
  colnames(rf.df)[2] <- "obs"
  R2 <- defaultSummary(rf.df)[2]
  return(R2)
}
find.rf.cv2 <- function(data, fold.list, m, n, node) {
  results <- c(rep(0.0, length(fold.list)))
  for (i in 1:length(fold.list))
  {
    fold <- fold.list[[i]]
    train = data[-fold, ]
    tst = data[fold, ]
    x = train[ , -1]
    y = train[ , 1]
    rf.cv <- randomForest(x = x, y = y, mtry = m,
                          ntree = n, nodesize = node)
    rf.df <- predict(rf.cv, tst[ , -1]) %>% 
      cbind(tst[ , 1]) %>%
      data.frame()
    colnames(rf.df)[1] <- "pred"
    colnames(rf.df)[2] <- "obs"
    R2 <- defaultSummary(rf.df)[2]
    results[i] <- R2
  }
  return(sum(results) / length(fold.list))
}





matrix.ga <- readRDS("./GAFS.RDS") %>% data.matrix()
m <- rfImpute(x = matrix.ga[ , -1], y = matrix.ga[ , 1], ntree = 400, na.action = na.omit)
saveRDS(m, "./GAFS.impute.RDS")
set.seed(48)
trn.ind <- sample(x = 1:nrow(m), size = round(0.8 * nrow(m)))
trn <- m[trn.ind, ]
tst <- m[-trn.ind, ]

rf.x <- trn[ , -1]
rf.y <- trn[ , 1]
rf <- randomForest(y ~., data = trn, ntree = 1000, na.action = na.omit, mtry = 100, nodesize = .1)
rf.tst.df <- predict(rf, tst[ , -1]) %>%
  cbind(tst[ , 1]) %>% 
  data.frame()
colnames(rf.tst.df)[1] <- "pred"
colnames(rf.tst.df)[2] <- "obs"
rf.tst.df %>% ggplot(., aes(x = obs, y = pred)) +
  geom_point() + 
  geom_abline(slope = 1, intercept = 0)
defaultSummary(rf.tst.df) # 0.432

rf.trn.df <- predict(rf, trn[ , -1]) %>%
  cbind(trn[ , 1]) %>% 
  data.frame()
colnames(rf.trn.df)[1] <- "pred"
colnames(rf.trn.df)[2] <- "obs"
rf.trn.df %>% ggplot(., aes(x = obs, y = pred)) +
  geom_point() + 
  geom_abline(slope = 1, intercept = 0)
defaultSummary(rf.trn.df) # 0.939

find.rf.cv2(m, rf.folds, 200, 1000, 0.1)
