dir.create("./tuning/glmnet")

# Libraries ---------------------------------------------------------------

library(caret)
library(ggplot2)
library(glmnet)
# library(Matrix)
library(tidyverse)

# Functions ---------------------------------------------------------------

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
                      pmax = 50, 
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

# Alpha -------------------------------------------------------------------
#     Loading Data --------------------------------------------------------

# Reading data with all descriptors
trn.all <- readRDS("./pre-process/alpha/1/pp.RDS") 
colnames(trn.all) <- str_replace(colnames(trn.all), "-", ".")
trn.guest <- trn.all$guest
trn.df <- select(trn.all, -guest)

rfe1 <- readRDS("./feature.selection/alpha/rfe1.RDS")
trn.pred <- c("DelG", predictors(rfe1))

trn.df <- trn.df[ , colnames(trn.df) %in% trn.pred] 
trn <- data.matrix(trn.df)

#     Estimation ----------------------------------------------------------

#     Alpha ---
alpha.range <- seq(0, 1, 0.1)
results1.alpha <- do.call(rbind, lapply(alpha.range, FUN = tune.glm.alpha, 
                                        data = trn, nfolds = 10, seed = 101))
results2.alpha <- do.call(rbind, lapply(alpha.range, FUN = tune.glm.alpha, 
                                        data = trn, nfolds = 10, seed = 102)) 
results3.alpha <- do.call(rbind, lapply(alpha.range, FUN = tune.glm.alpha, 
                                        data = trn, nfolds = 10, seed = 103)) 
results.alpha <- rbind(results1.alpha, results2.alpha, results3.alpha) %>%
  mutate(seed = as.factor(seed))
ggplot(results.alpha, aes(x = alpha, y = rsquared, 
                          group = seed, color = seed)) + 
  geom_line() + 
  theme_bw()

# Maximum Degrees of Freedom ---

df.range <- 2 ^(1:7)
results1.df <- do.call(rbind, lapply(df.range, FUN = tune.glm.dfmax, 
                                        data = trn, nfolds = 10, seed = 101)) 
results2.df <- do.call(rbind, lapply(df.range, FUN = tune.glm.dfmax, 
                                        data = trn, nfolds = 10, seed = 102)) 
results3.df <- do.call(rbind, lapply(df.range, FUN = tune.glm.dfmax, 
                                        data = trn, nfolds = 10, seed = 103)) 
results.df <- rbind(results1.df, results2.df, results3.df) %>%
  mutate(seed = as.factor(seed))
ggplot(results.df, aes(x = dfmax, y = rsquared, 
                          group = seed, color = seed)) + 
  geom_line() + 
  theme_bw()

dir.create("./tuning/glmnet/alpha")
saveRDS(results.alpha, "./tuning/glmnet/alpha/alpha.RDS")
saveRDS(results.df, "./tuning/glmnet/alpha/dfmax.RDS")

#     Tuning --------------------------------------------------------------

# 11 * 7 = 77 combinations

glm.combos <- expand.grid(alpha.range, df.range)
colnames(glm.combos) <- c("alpha", "dfmax")
alpha.combos <- glm.combos$alpha
df.combos <- glm.combos$dfmax

set.seed(1001)
system.time(
  results.combos <- do.call(
    rbind,
    mapply(
      FUN = tune.glm,
      a = alpha.combos,
      max = df.combos,
      MoreArgs = 
        list(nfolds = 10, data = trn), 
      SIMPLIFY = F
    )
  )
)

# system.time output
# user  system elapsed 
# 14.81    0.00   14.93 

# results.combos[order(results.combos$rsquared, decreasing = T), ] %>% head()
# results.combos[order(results.combos$rmse), ] %>% head()

# rsquared of 0.510, rmse = 3.65
# alpha = 0, dfmax = 128

# best rmse: see above

# rsquared of 0.510, rmse = 3.69
# alpha = 1.0, dfmax = 16

saveRDS(results.combos, "./tuning/glmnet/alpha/tune.RDS")
results.combos$alpha <- as.factor(results.combos$alpha)
results.combos$dfmax <- as.factor(results.combos$dfmax)
ggplot(results.combos, aes(x = alpha, y = dfmax, fill = rsquared)) + 
  geom_raster() + 
  scale_fill_gradientn(colours = terrain.colors(20)) + 
  theme_bw() + 
  labs(title = "GLMNet tuning for alpha-CD", x = "Alpha", y = "Maximum degrees of freedom", 
       fill = "R2")
ggsave("./tuning/glmnet/alpha/tune.png", dpi = 450)

# Beta -------------------------------------------------------------------
#     Loading Data --------------------------------------------------------

# Reading data with all descriptors
trn.all <- readRDS("./pre-process/beta/1/pp.RDS") 
colnames(trn.all) <- str_replace(colnames(trn.all), "-", ".")
trn.guest <- trn.all$guest
trn.df <- select(trn.all, -guest)

rfe1 <- readRDS("./feature.selection/beta/rfe1.RDS")
trn.pred <- c("DelG", predictors(rfe1))

trn.df <- trn.df[ , colnames(trn.df) %in% trn.pred] 
trn <- data.matrix(trn.df)

#     Estimation ----------------------------------------------------------

#     Alpha ---
alpha.range <- seq(0, 1, 0.1)
results1.alpha <- do.call(rbind, lapply(alpha.range, FUN = tune.glm.alpha, 
                                        data = trn, nfolds = 10, seed = 101))
results2.alpha <- do.call(rbind, lapply(alpha.range, FUN = tune.glm.alpha, 
                                        data = trn, nfolds = 10, seed = 102)) 
results3.alpha <- do.call(rbind, lapply(alpha.range, FUN = tune.glm.alpha, 
                                        data = trn, nfolds = 10, seed = 103)) 
results.alpha <- rbind(results1.alpha, results2.alpha, results3.alpha) %>%
  mutate(seed = as.factor(seed))
ggplot(results.alpha, aes(x = alpha, y = rsquared, 
                          group = seed, color = seed)) + 
  geom_line() + 
  theme_bw()

# Maximum Degrees of Freedom ---

df.range <- 2 ^(2:8)
results1.df <- do.call(rbind, lapply(df.range, FUN = tune.glm.dfmax, 
                                     data = trn, nfolds = 10, seed = 101)) 
results2.df <- do.call(rbind, lapply(df.range, FUN = tune.glm.dfmax, 
                                     data = trn, nfolds = 10, seed = 102)) 
results3.df <- do.call(rbind, lapply(df.range, FUN = tune.glm.dfmax, 
                                     data = trn, nfolds = 10, seed = 103)) 
results.df <- rbind(results1.df, results2.df, results3.df) %>%
  mutate(seed = as.factor(seed))
ggplot(results.df, aes(x = dfmax, y = rsquared, 
                       group = seed, color = seed)) + 
  geom_line() + 
  theme_bw()

dir.create("./tuning/glmnet/beta")
saveRDS(results.alpha, "./tuning/glmnet/beta/alpha.RDS")
saveRDS(results.df, "./tuning/glmnet/beta/dfmax.RDS")

#     Tuning --------------------------------------------------------------

# 11 * 7 = 77 combinations

glm.combos <- expand.grid(alpha.range, df.range)
colnames(glm.combos) <- c("alpha", "dfmax")
alpha.combos <- glm.combos$alpha
df.combos <- glm.combos$dfmax

set.seed(1001)
system.time(
  results.combos <- do.call(
    rbind,
    mapply(
      FUN = tune.glm,
      a = alpha.combos,
      max = df.combos,
      MoreArgs = 
        list(nfolds = 10, data = trn), 
      SIMPLIFY = F
    )
  )
)

# system.time output
# user  system elapsed 
# 27.02    0.05   27.08 

# results.combos[order(results.combos$rsquared, decreasing = T), ] %>% head()
# results.combos[order(results.combos$rmse), ] %>% head()

# rsquared of 0.628, rmse = 3.21
# alpha = 0.8, dfmax = 32

# best rmse: see above

# rsquared of 0.607, rmse = 3.19
# alpha = 0.5, dfmax = 32

saveRDS(results.combos, "./tuning/glmnet/alpha/tune.RDS")
results.combos$alpha <- as.factor(results.combos$alpha)
results.combos$dfmax <- as.factor(results.combos$dfmax)
ggplot(results.combos, aes(x = alpha, y = dfmax, fill = rsquared)) + 
  geom_raster() + 
  scale_fill_gradientn(colours = terrain.colors(20)) + 
  theme_bw() + 
  labs(title = "GLMNet tuning for beta-CD", x = "Alpha", y = "Maximum degrees of freedom", 
       fill = "R2")
ggsave("./tuning/glmnet/beta/tune.png", dpi = 450)
