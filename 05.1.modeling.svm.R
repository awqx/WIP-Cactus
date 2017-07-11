# Libraries and Packages --------------------------------------------------

library(caret)
library(e1071)
library(kernlab)
library(Matrix)
library(stringr)
library(tidyverse)

# Loading Data ------------------------------------------------------------

df <- readRDS("./DelG.df.RDS")
set.seed(25)
trn.ind <- sample(x = 1:nrow(df), size = round(0.7 * nrow(df)))
df.trn <- df[trn.ind, ]
df.trn.x <- df.trn[ , -1]
df.trn.y <- df.trn[ , 1]
df.tst <- df[-trn.ind, ]
df.tst.x <- df.tst[ , -1]
df.tst.y <- df.tst[ , 1]
df.ga <- cbind(df[ , 1], df[ , colnames(df) %in% ga.final])
colnames(df.ga)[1] <- "DelG"
df.ga.trn <- df[trn.ind, ]
df.ga.trn.x <- df.ga.trn[ , -1]
df.ga.trn.y <- df.ga.trn[ , 1]
df.ga.tst <- df[-trn.ind, ]
df.ga.tst.x <- df.ga.tst[ , -1]
df.ga.tst.y <- df.tst[ , 1]

sprse <- readRDS("./DelG.sparse.RDS")
trn.ind <- sample(x = 1:nrow(sprse), size = round(0.7 * nrow(sprse)))
sprse.trn <- sprse[trn.ind, ]
sprse.trn.x <- sprse.trn[ , -1:-2]
sprse.trn.y <- sprse.trn[ , 2]
sprse.tst <- sprse[-trn.ind, ]
sprse.tst.x <- sprse.tst[ , -1:-2]
sprse.tst.y <- sprse.tst[ , 2]

sprse.ga <- sparse.model.matrix(~., df.ga)
sprse.ga.trn <- sprse.ga[trn.ind, ]
sprse.ga.trn.x <- sprse.ga.trn[ , -1:-2]
sprse.ga.trn.y <- sprse.ga.trn[ , 2]
sprse.ga.tst <- sprse.ga[-trn.ind, ]
sprse.ga.tst.x <- sprse.ga.tst[ , -1:-2]
sprse.ga.tst.y <- sprse.ga.tst[ , 2]



# SVM Tuning --------------------------------------------------------------

# This tunes the cost value
# If an error is thrown, just re-run this part
svm.tune.rad <- train(x = df.trn.x, y = df.trn.y, 
                      method = "svmRadial", 
                      tuneLength = 16, 
                      trControl = trainControl(method = "cv"))
# Sigma = 0.000867 Cost = 2
svm.tune.poly <- train(x = df.trn.x, y = df.trn.y, 
                       method = "svmPoly", 
                       tuneLength = 4, 
                       trControl = trainControl(method = "cv"))
# Degree = 3, scale = 0.01, C = 2

# Another method of tuning (using tune.svm)
#   Polynomial ----
#     For Gamma ... gamma = 512 ----
tune.poly.gamma <- tune.svm(
  DelG ~ .,
  data = df.trn,
  gamma = 2 ^ (-4:4),
  degree = 2, 
  kernel = "polynomial",
  performances = T
) # All gammas exactly the same
tune.poly.gamma.df <- tune.poly.gamma$performances
tune.poly.gamma2 <- tune.svm(
  DelG ~ .,
  data = df.trn,
  gamma = 2 ^ (4:12),
  cost = 2,
  degree = 2, 
  kernel = "polynomial",
  performances = T
) # gamma = 256
tune.poly.gamma.df <- rbind(tune.poly.gamma.df, 
                            tune.poly.gamma2$performances)

#     For Cost ... cost = 8 ----
tune.poly.cost1 <- tune.svm(
  DelG ~ .,
  data = df.trn,
  gamma = 512,
  cost = 2 ^ (-4:4),
  degree = 2,
  kernel = "polynomial",
  performances = T
)
tune.poly.cost.df <- tune.poly.cost1$performances
#     For Nu (alt to cost) ... nu = tune.poly.ga -----
tune.poly.nu <- tune.svm(
  DelG ~., 
  data = df.trn, 
  gamma = 4, 
  nu = 0.1 * 2 ^ (-3:3), 
  degree = 2, 
  kernel = "polynomial", 
  performances = T
) # no variation
tune.poly.nu.df <- tune.poly.nu$performances
#     For Degree ... degree = 2 ----
tune.poly.deg <- tune.svm(
  DelG ~ .,
  data = df.trn,
  gamma = 512,
  cost = 8,
  degree = 1:4,
  kernel = "polynomial",
  performances = T
) 
tune.poly.deg.df <- tune.poly.deg$performances
# degree tends to vary between 2 and 3

#   Radial ----
#     For Cost ... cost = 128----
tune.rad.cost <- tune.svm(
  DelG ~ .,
  data = df.trn,
  gamma = 2,
  cost = 2 ^ (-4:4),
  kernel = "radial",
  performances = T
) # cost = 16
tune.rad.cost.df <- tune.rad.cost$performances
tune.rad.cost2 <- tune.svm(
  DelG ~ .,
  data = df.trn,
  gamma = 2,
  cost = 2 ^ (4:8),
  kernel = "radial",
  performances = T
) # cost = 128
tune.rad.cost.df <- rbind(tune.rad.cost.df, tune.rad.cost2$performances)

# Manual Tuning -----------------------------------------------------------

validate.svm.cost <- function(data, nfolds, cost, kerneltype, seed) {
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
                  cost = cost,
                  kernel = kerneltype)
    svm.df <- predict(svm.cv, tst[ , -1:-2]) %>% 
      cbind(tst[ , 2]) %>%
      data.frame()
    colnames(svm.df)[1] <- "pred"
    colnames(svm.df)[2] <- "obs"
    R2 <- defaultSummary(svm.df)[2] 
    results[i] <- R2
  }
  return(sum(results) / nfolds)
}
validate.svm.gamma <- function(data, nfolds, g, kerneltype, seed) {
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
  return(sum(results) / length(fold.list))
}
validate.svm.epsilon <- function(data, nfolds, e, kerneltype, seed) {
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
  return(sum(results) / length(fold.list))
}

validate.svm.cost(sprse.trn, 10, 0.125, "polynomial", 5)
validate.svm.cost(sprse.trn, 10, 2, "polynomial", 5)
validate.svm.cost(sprse.trn, 10, 8, "polynomial", 5) # peaks around here
validate.svm.cost(sprse.trn, 10, 64, "polynomial", 5)
validate.svm.cost(sprse.trn, 10, 256, "polynomial", 5)
validate.svm.cost(sprse.trn, 10, 512, "polynomial", 5)

# Barely any change across a range of gamma values
validate.svm.gamma(sprse.trn, 10, 0.125, "polynomial", 5)
validate.svm.gamma(sprse.trn, 10, 2, "polynomial", 5)
validate.svm.gamma(sprse.trn, 10, 8, "polynomial", 5)
validate.svm.gamma(sprse.trn, 10, 64, "polynomial", 5) 
validate.svm.gamma(sprse.trn, 10, 256, "polynomial", 5)

# epsilon tends toward smaller numbers, 0.0625 seems alright
validate.svm.epsilon(sprse.trn, 10, 0.000488, "polynomial", 144) # peaks around here
validate.svm.epsilon(sprse.trn, 10, 0.0625, "polynomial", 144) # peaks around here
validate.svm.epsilon(sprse.trn, 10, 0.125, "polynomial", 144)
validate.svm.epsilon(sprse.trn, 10, 0.25, "polynomial", 144)
validate.svm.epsilon(sprse.trn, 10, 0.5, "polynomial", 144)
validate.svm.epsilon(sprse.trn, 10, 1, "polynomial", 144)
validate.svm.epsilon(sprse.trn, 10, 2, "polynomial", 144) # peaks around here too
validate.svm.epsilon(sprse.trn, 10, 4, "polynomial", 144)

validate.svm.cost(sprse.ga.trn, 10, 0.125, "polynomial", 5)
validate.svm.cost(sprse.ga.trn, 10, 2, "polynomial", 5)
validate.svm.cost(sprse.ga.trn, 10, 8, "polynomial", 5) # peaks around here
validate.svm.cost(sprse.ga.trn, 10, 64, "polynomial", 5)
validate.svm.cost(sprse.ga.trn, 10, 256, "polynomial", 5)
validate.svm.cost(sprse.ga.trn, 10, 512, "polynomial", 5)




# SVM GAFS Model ----------------------------------------------------------
#   Polynomial ----
svm.poly <- ksvm(x = data.matrix(df.trn.x), y = df.trn.y, 
                 kernel = "polydot",
                 C = 0.0001, 
                 degree = 3, 
                 scale = 0.1)
svm.poly.tst <- predict(svm.poly, df.tst.x) %>% 
  cbind(df.tst.y) %>%
  data.frame()
colnames(svm.poly.tst)[1] <- "pred"
colnames(svm.poly.tst)[2] <- "obs"
ggplot(svm.poly.tst, aes(x = obs, y = pred)) + 
  geom_point() + 
  coord_fixed() + 
  theme_bw() + 
  geom_abline(slope = 1, intercept = 0)
defaultSummary(svm.poly.tst) # 0.515

svm.poly.trn <- predict(svm.poly, df.trn.x) %>% 
  cbind(df.trn.y) %>%
  data.frame()
colnames(svm.poly.trn)[1] <- "pred"
colnames(svm.poly.trn)[2] <- "obs"
ggplot(svm.poly.trn, aes(x = obs, y = pred)) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0)
# Results
defaultSummary(svm.poly.tst)[2] # 0.515
defaultSummary(svm.poly.trn)[2] # 0.902

svm.poly2 <- svm(x = sprse.trn.x, y = sprse.trn.y, 
                 nu = 0.25, gamma = 8, kernel = "polynomial", 
                 degree = 3)
svm.poly2.trn <- predict(svm.poly2, sprse.trn.x) %>% 
  cbind(sprse.trn.y) %>%
  data.frame()
colnames(svm.poly2.trn)[1] <- "pred"
colnames(svm.poly2.trn)[2] <- "obs"
ggplot(svm.poly2.trn, aes(x = obs, y = pred)) + 
  geom_point() + 
  coord_fixed() + 
  theme_bw() + 
  geom_abline(slope = 1, intercept = 0)

svm.poly2.trn <- predict(svm.poly2, sprse.trn.x) %>% 
  cbind(sprse.trn.y) %>%
  data.frame()
colnames(svm.poly2.trn)[1] <- "pred"
colnames(svm.poly2.trn)[2] <- "obs"
ggplot(svm.poly2.trn, aes(x = obs, y = pred)) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0)

defaultSummary(svm.poly2.tst)[2] # 0.357
defaultSummary(svm.poly2.trn)[2] # 0.998

#   RBF ----
svm.rad <- ksvm(x = data.matrix(df.trn.x), y = df.trn.y, 
                kernel = "rbfdot",
                cost = 2048)
svm.rad.tst <- predict(svm.rad, df.tst.x) %>% 
  cbind(df.tst.y) %>%
  data.frame()
colnames(svm.rad.tst)[1] <- "pred"
colnames(svm.rad.tst)[2] <- "obs"
ggplot(svm.rad.tst, aes(x = obs, y = pred)) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0)

svm.rad.trn <- predict(svm.rad, df.trn.x) %>% 
  cbind(df.trn.y) %>%
  data.frame()
colnames(svm.rad.trn)[1] <- "pred"
colnames(svm.rad.trn)[2] <- "obs"
ggplot(svm.rad.trn, aes(x = obs, y = pred)) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0)
# OVERFITTED
defaultSummary(svm.rad.tst)[2] # 0.534
defaultSummary(svm.rad.trn)[2] # 0.610

svm.rad2 <- svm(x = sprse.trn.x, y = sprse.trn.y, 
                kernel = "radial",
                cost = 256)
svm.rad2.tst <- predict(svm.rad2, sprse.tst.x) %>% 
  cbind(sprse.tst.y) %>%
  data.frame()
colnames(svm.rad2.tst)[1] <- "pred"
colnames(svm.rad2.tst)[2] <- "obs"
ggplot(svm.rad2.tst, aes(x = obs, y = pred)) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0)

svm.rad2.trn <- predict(svm.rad2, sprse.trn.x) %>% 
  cbind(sprse.trn.y) %>%
  data.frame()
colnames(svm.rad2.trn)[1] <- "pred"
colnames(svm.rad2.trn)[2] <- "obs"
ggplot(svm.rad2.trn, aes(x = obs, y = pred)) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0)
# OVERFITTED
defaultSummary(svm.rad2.tst)[2] # 0.658 # SHRUG
defaultSummary(svm.rad2.trn)[2] # 0.984

# Cross-Validation --------------------------------------------------------

set.seed(256)
folds <- createFolds(y = sparse.data[ , 2], k = 8, list = T)

validate.svm.r2 <- function(data, fold.list, c, g) {
  results <- c(rep(0.0, length(fold.list)))
  for (i in 1:length(fold.list))
  {
    fold <- fold.list[[i]]
    train <- data[-fold, ]
    tst <- data[fold, ]
    x <- train[ , -1:-2]
    y <- train[ , 2]
    svm.cv <- svm(x = x, y = y, 
                  cost = c, 
                  gamma = g, 
                  kernel = "polynomial", 
                  degree = 3)
    svm.df <- predict(svm.cv, tst[ , -1:-2]) %>% 
      cbind(tst[ , 2]) %>%
      data.frame()
    colnames(svm.df)[1] <- "pred"
    colnames(svm.df)[2] <- "obs"
    R2 <- defaultSummary(svm.df)[2] 
    results[i] <- R2
  }
  return(sum(results) / length(fold.list))
}

validate.svm.r2(sparse.data, folds, 0.25, 0.0625)

# External Validation -----------------------------------------------------

suzuki.beta <- read_csv("C:/Users/Wei Xin/Desktop/2017 06 16 Suzuki beta.csv") %>% 
  mutate(alpha = c(rep(0, 185))) %>%
  mutate(beta = c(rep(1, 185))) %>%
  mutate(gamma = c(rep(0, 185)))
suzuki.alpha <- read_csv("C:/Users/Wei Xin/Desktop/2017 06 16 Suzuki alpha comp.csv") %>%
  mutate(alpha = c(rep(1, 55))) %>%
  mutate(beta = c(rep(0, 55))) %>%
  mutate(gamma = c(rep(0, 55)))
suzuki <- rbind(suzuki.alpha, suzuki.beta)
saveRDS(suzuki, "~/suzuki.raw.RDS")
zero.pred.s <- nearZeroVar(suzuki)
sp.no.zero <- suzuki[ , -zero.pred.s]

st <- preProcess(sp.no.zero[ , -1], na.remove = T, 
                 method = c("center", "scale")) %>%
  predict(., sp.no.zero[ , -1])
zero.pred2.s <- nearZeroVar(st)
st <- st[ , -zero.pred2.s]
suzuki.transformed <- cbind(sp.no.zero[ , 1], st)
saveRDS(suzuki.transformed, "~/suzuki.transformed.RDS")

st.pred <- st[ , -1:-17]
too.high <- findCorrelation(cor(st.pred, use = "pairwise.complete.obs"), 0.95)
corr.pred <- names(st.pred)[too.high]
st.pred <- st.pred[ , -too.high]
st <- cbind(st[ , 1:17], st.pred)


# SVM VIP -----------------------------------------------------------------

df.vip <- cbind(df.dg[ , 1], df.dg[ , colnames(df.dg) %in% highest.vip.desc])
colnames(df.vip)[1] <- "DelG"
sparse.vip <- sparse.model.matrix(~., df.vip)
set.seed(222)
trn.ind <- sample(x = 1:nrow(sparse.vip), size = round(0.75 * nrow(sparse.vip)))
sparse.vip.trn <- sparse.vip[trn.ind, ]
sparse.vip.tst <- sparse.vip[-trn.ind, ]

svm.vip <- svm(x = sparse.vip.trn[ , -1:-2], 
               y = sparse.vip.trn[ ,2], 
               cost = 0.625,
               gamma = 50, 
               kernel = "linear")
svm.vip.tst <- predict(svm.vip, sparse.vip.tst[ , -1:-2]) %>% 
  cbind(sparse.vip.tst[ , 2]) %>% data.frame()
colnames(svm.vip.tst)[1] <- "pred"
colnames(svm.vip.tst)[2] <- "obs"
ggplot(svm.vip.tst, aes(x = obs, y = pred)) + 
  geom_point() + 
  coord_fixed() + 
  theme_bw() + 
  geom_abline(slope = 1, intercept = 0)
defaultSummary(svm.vip.tst)

svm.vip.trn <- predict(svm.vip, sparse.vip.trn[ , -1:-2]) %>% 
  cbind(sparse.vip.trn[ , 2]) %>% data.frame()
colnames(svm.vip.trn)[1] <- "pred"
colnames(svm.vip.trn)[2] <- "obs"
ggplot(svm.vip.trn, aes(x = obs, y = pred)) + 
  geom_point() + 
  coord_fixed() + 
  theme_bw() + 
  geom_abline(slope = 1, intercept = 0)
defaultSummary(svm.vip.trn)
# ----
# ----
# KernLab -----------------------------------------------------------------
# KernLab contains more possible kernels
ksvm.tune <- train(x = svm.x, y = svm.y, method = "svmPoly", 
                   tuneLength = 12, trControl = "cv")

# SVM Tuning --------------------------------------------------------------

svm.x <- sparse.dg.trn[ , -1:-2]
svm.y <- sparse.dg.trn[ , 2]
svm.tune.1 <- tune.svm(DelG ~., data = df.dg.trn, gamma = 2^(-2:2), cost = 2^(-2:2), 
                       kernel = "polynomial", performances = T) 
svm.tune <- svm.tune.1$performances
# svm.tune.1$best.parameters: gamma = 4 ... cost = 1
svm.tune.2 <- tune.svm(DelG ~., data = df.dg.trn, gamma = 2^(-4:-2), cost = 2^(-4:-2), 
                       kernel = "polynomial", performances = T) 
svm.tune <- rbind(svm.tune, svm.tune.2$performances)
# svm.tune.2$best.parameters: gamma = 0.25 ... cost = 0.0625
svm.tune.3 <- tune.svm(DelG ~., data = df.dg.trn, gamma = 2^(-6:-4), cost = 2^(-6:-4), 
                       kernel = "polynomial", performances = T) 
svm.tune <- rbind(svm.tune, svm.tune.3$performances)
# svm.tune.3$best.parameters: gamma = 0.015625 ... cost = 0.03125
svm.tune.4 <- tune.svm(DelG ~., data = df.dg.trn, gamma = 2^(-3:-2), cost = 2^(-4:-0), 
                       kernel = "polynomial", performances = T)
svm.tune <- rbind(svm.tune, svm.tune.4$performances)
# svm.tune.4$best.parameters: gamma: 0.25 ... cost: 0.25
svm.tune.5 <- tune.svm(DelG ~., data = df.dg.trn, gamma = 0.25, cost = 0.25, 
                       kernel = "polynomial", degree = 0:4, performances = T)
# degree = 3
svm.tune.6 <- tune.svm(DelG ~., data = df.dg.trn, gamma = 2 ^ (-4:-2), cost = 2 ^ (-4:-2), 
                       kernel = "polynomial", degree = 3, 
                       epsilon = 2^(-4:0), performances = T)

# Let's try to validate
set.seed(256)
svm.folds <- createFolds(y = sparse.dg[ , 2], k = 10, list = T)
svm.folds2 <- createFolds(y = sparse.dg[ , 2], k = 10, list = T)

find.svm.cv.R2 <- function(data, fold.list, fold.num) {
  fold <- fold.list[[fold.num]]
  train = data[-fold, ]
  tst = data[fold, ]
  x = train[ , -1:-2]
  y = train[ , 2]
  svm.cv <- svm(x = x, y = y, cost = 0.25, 
                gamma = 0.0625, epsilon = 0.0625,
                kernel = "polynomial", degree = 3)
  svm.df <- predict(svm.cv, tst[ , -1:-2]) %>% 
    cbind(tst[ , 2]) %>%
    data.frame()
  colnames(svm.df)[1] <- "pred"
  colnames(svm.df)[2] <- "obs"
  R2 <- defaultSummary(svm.df)[2]
  return(R2)
}
find.svm.cv2 <- function(data, fold.list) {
  results <- c(rep(0.0, length(fold.list)))
  for (i in 1:length(fold.list))
  {
    fold <- fold.list[[i]]
    train = data[-fold, ]
    tst = data[fold, ]
    x = train[ , -1:-2]
    y = train[ , 2]
    svm.cv <- svm(x = x, y = y, cost = 0.25, 
                  gamma = 0.0625, epsilon = 0.0625,
                  kernel = "polynomial", degree = 3)
    svm.df <- predict(svm.cv, tst[ , -1:-2]) %>% 
      cbind(tst[ , 2]) %>%
      data.frame()
    colnames(svm.df)[1] <- "pred"
    colnames(svm.df)[2] <- "obs"
    R2 <- defaultSummary(svm.df)[2] 
    results[i] <- R2
  }
  return(sum(results) / length(fold.list))
}
find.svm.cv3 <- function(data, fold.list, c, g, e, d) {
  results <- c(rep(0.0, length(fold.list)))
  for (i in 1:length(fold.list))
  {
    fold <- fold.list[[i]]
    train = data[-fold, ]
    tst = data[fold, ]
    x = train[ , -1:-2]
    y = train[ , 2]
    svm.cv <- svm(x = x, y = y, cost = c, 
                  gamma = g, epsilon = e,
                  kernel = "polynomial", degree = d)
    svm.df <- predict(svm.cv, tst[ , -1:-2]) %>% 
      cbind(tst[ , 2]) %>%
      data.frame()
    colnames(svm.df)[1] <- "pred"
    colnames(svm.df)[2] <- "obs"
    R2 <- defaultSummary(svm.df)[2] 
    results[i] <- R2
  }
  return(sum(results) / length(fold.list))
}

find.svm.cv.R2(sparse.dg, svm.folds, 1) # 0.438
find.svm.cv.R2(sparse.dg, svm.folds, 2) # 0.539
find.svm.cv.R2(sparse.dg, svm.folds, 3) # 0.306
find.svm.cv.R2(sparse.dg, svm.folds, 4) # 0.769
find.svm.cv.R2(sparse.dg, svm.folds, 5) # 0.736
find.svm.cv.R2(sparse.dg, svm.folds, 6) # 0.562
find.svm.cv.R2(sparse.dg, svm.folds, 7) # 0.005
find.svm.cv.R2(sparse.dg, svm.folds, 8) # 0.665
find.svm.cv.R2(sparse.dg, svm.folds, 9) # 0.860
find.svm.cv.R2(sparse.dg, svm.folds, 10) # 0.466
# Average R2 = 0.534
find.svm.cv2(sparse.dg, svm.folds) # 0.534
find.svm.cv3(sparse.dg, svm.folds, 0.25, 0.0625, 0.0625, 3) # 0.535
find.svm.cv3(sparse.dg, svm.folds, 0.0625, 0.0625, 0.0625 , 3) # 0.535
find.svm.cv3(sparse.dg, svm.folds, 0.5, 0.0625, 0.0625 , 3) # 0.535
find.svm.cv3(sparse.dg, svm.folds, 0.25, 0.25, 0.0625 , 3) # 0.535
find.svm.cv3(sparse.dg, svm.folds, 0.25, 0.0625, 0.25 , 3) # 0.535
find.svm.cv3(sparse.dg, svm.folds, 0.25, 0.0625, 1 , 3) # 0.523
find.svm.cv3(sparse.dg, svm.folds, 0.25, 0.0625, 0.5 , 3) # 0.533
find.svm.cv3(sparse.dg, svm.folds, 0.25, 0.0625, 4 , 3) # 0.458
find.svm.cv3(sparse.dg, svm.folds, 0.25, 0.0625, 0.001 , 3) # 0.534
find.svm.cv3(sparse.dg, svm.folds, 0.25, 0.0625, 0.0625 , 2) # 0.570 
find.svm.cv3(sparse.dg, svm.folds, 0.25, 0.0625, 0.0625, 4) # 0.513
find.svm.cv3(sparse.dg, svm.folds, 0.25, 0.0625, 0.0625, 5) # 0.488
find.svm.cv3(sparse.dg, svm.folds, 0.25, 0.0625, 0.0625, 1) # 0.517
find.svm.cv3(sparse.dg, svm.folds, 0.0625, 0.0625, 0.0625 , 2) # 0.522
find.svm.cv3(sparse.dg, svm.folds, 1.25, 0.0625, 0.0625 , 2) # 0.565
find.svm.cv3(sparse.dg, svm.folds, 5, 0.0625, 0.0625 , 2) # 0.565
find.svm.cv3(sparse.dg, svm.folds, 0.25, 1, 0.0625 , 2) # 0.565
find.svm.cv3(sparse.dg, svm.folds, 0.25, 0.0625, 1, 2) # 0.545
find.svm.cv3(sparse.dg, svm.folds, 2, 5, 2 , 2) # 0.481

fold4 <- svm.folds[[1]]
trn.temp <- sparse.dg[-fold4, ]
tst.temp <- sparse.dg[fold4, ]
x.temp = trn.temp[ , -1:-2]
y.temp = trn.temp[ , 2]
svm.temp <- svm(x = x.temp, y = y.temp, cost = 0.25, 
                gamma = 0.25, epsilon = 0.0625, 
                kernel = "polynomial", degree = 3)
svm.df.temp <- predict(svm.temp, tst.temp[ , -1:-2]) %>%
  cbind(tst.temp[ , 2]) %>%
  data.frame()
colnames(svm.df.temp)[1] <- "predicted"
colnames(svm.df.temp)[2] <- "experimental"
ggplot(svm.df.temp, aes(x = experimental, y = predicted)) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0)

# SVM Model ---------------------------------------------------------------

svm.dg <- svm(x = svm.x, y = svm.y, sparse.dg.trn, cost = 0.25, gamma = 0.0625, epsilon = 0.0625, 
              kernel = "polynomial", degree = 3)
svm.dg.tst <- predict(svm.dg, sparse.dg.tst[ , -1:-2]) %>% 
  cbind(sparse.dg.tst[ , 2]) %>%
  data.frame()
colnames(svm.dg.tst)[1] <- "pred"
colnames(svm.dg.tst)[2] <- "obs"
ggplot(svm.dg.tst, aes(x = obs, y = pred)) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0)

svm.dg.trn <- predict(svm.dg, sparse.dg.trn[ , -1:-2]) %>% 
  cbind(sparse.dg.trn[ , 2]) %>%
  data.frame()
colnames(svm.dg.trn)[1] <- "pred"
colnames(svm.dg.trn)[2] <- "obs"
ggplot(svm.dg.trn, aes(x = obs, y = pred)) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0)
# OVERFITTING????
defaultSummary(svm.dg.tst)[2] # 0.637
defaultSummary(svm.dg.trn)[2] # 0.9998
