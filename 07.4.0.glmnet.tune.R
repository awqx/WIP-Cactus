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

# Alpha -------------------------------------------------------------------
#     Loading Data --------------------------------------------------------

# Reading data with all descriptors
trn.all <- readRDS("./pre-process/alpha/1/pp.RDS") 
colnames(trn.all) <- str_replace(colnames(trn.all), "-", ".")
trn.guest <- trn.all$guest
trn.df <- select(trn.all, -guest)

features <- readRDS("./feature.selection/alpha.vars.RDS")
trn.df <- trn.df[ , colnames(trn.df) %in% c("DelG", features)] 
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

df.range <- c(0, 1, 5, 10, 25, 50, 75, 150)
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

# 11 * 8 = 88 combinations
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
# 10.39    0.00   10.44 
# 13.34    0.00   13.61 

results.combos[order(results.combos$rsquared, decreasing = T), ] %>% head()
# nfolds alpha dfmax  rsquared     rmse
#    10   0.0   150 0.6364452 2.972314
#    10   0.0    50 0.6277415 2.875179
results.combos[order(results.combos$rmse), ] %>% head()
# nfolds alpha dfmax  rsquared     rmse
#    10   0.0    50 0.6277415 2.875179
#    10   0.0   150 0.6364452 2.97231

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

features <- readRDS("./feature.selection/beta.vars.RDS")
trn.df <- trn.df[ , colnames(trn.df) %in% c("DelG", features)] 
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

df.range <- c(0, 1, 2, 5, 10, 20, 50, 100)
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
#  user  system elapsed 
# 10.72    0.00   10.87
# 13.42    0.04   13.58 

results.combos[order(results.combos$rsquared, decreasing = T), ] %>% head()
# nfolds alpha dfmax  rsquared     rmse
#    10   0.7    20 0.6628094 3.367248
#    10   0.5    20 0.6485293 3.316619
results.combos[order(results.combos$rmse), ] %>% head()
# nfolds alpha dfmax  rsquared     rmse
#    10   1.0    20 0.6207187 3.315109
#    10   0.5    20 0.6485293 3.316619

saveRDS(results.combos, "./tuning/glmnet/beta/tune.RDS")
results.combos$alpha <- as.factor(results.combos$alpha)
results.combos$dfmax <- as.factor(results.combos$dfmax)
ggplot(results.combos, aes(x = alpha, y = dfmax, fill = rsquared)) + 
  geom_raster() + 
  scale_fill_gradientn(colours = terrain.colors(20)) + 
  theme_bw() + 
  labs(title = "GLMNet tuning for beta-CD", x = "Alpha", y = "Maximum degrees of freedom", 
       fill = "R2")
ggsave("./tuning/glmnet/beta/tune.png", dpi = 450)

# Gamma -------------------------------------------------------------------
#     Loading Data --------------------------------------------------------

convert.delg.ka <- function(delg) {
  # exp(n) = e^n
  joules <- delg*1000
  return (exp(joules / (-8.314 * 298)))
}
# Reading data with all descriptors
trn.all <- readRDS("./pre-process/gamma/6/pp.RDS") 
colnames(trn.all) <- str_replace(colnames(trn.all), "-", ".")
trn.guest <- trn.all$guest
trn.df <- select(trn.all, -guest)

features <- readRDS("./feature.selection/gamma.vars.RDS")
trn.df <- trn.df[ , colnames(trn.df) %in% c("DelG", features)] 
trn.logk <- trn.df %>%
  mutate(DelG = convert.delg.ka(DelG)) %>%
  mutate(DelG = log10(DelG))
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

df.range <- c(0, 1, 2, 5, 10, 20, 50, 100)
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

dir.create("./tuning/glmnet/gamma")
saveRDS(results.alpha, "./tuning/glmnet/gamma/alpha.RDS")
saveRDS(results.df, "./tuning/glmnet/gamma/dfmax.RDS")

#     Tuning --------------------------------------------------------------

# 11 * 8 = 88 combinations
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
#  user  system elapsed 
# 9.97    0.00   10.16

results.combos[order(results.combos$rsquared, decreasing = T), ] %>% head()
results.combos[order(results.combos$rmse), ] %>% head()

# nfolds alpha dfmax  rsquared     rmse
#    10   0.3    50 0.3675429 1.672413
#    10   0.6    10 0.3010348 1.607742
#    10   0.7    20 0.3522579 1.647089

saveRDS(results.combos, "./tuning/glmnet/gamma/tune.RDS")
results.combos$alpha <- as.factor(results.combos$alpha)
results.combos$dfmax <- as.factor(results.combos$dfmax)
ggplot(results.combos, aes(x = alpha, y = dfmax, fill = rsquared)) + 
  geom_raster() + 
  scale_fill_gradientn(colours = terrain.colors(20)) + 
  theme_bw() + 
  labs(title = "GLMNet tuning for gamma-CD", x = "Alpha", y = "Maximum degrees of freedom", 
       fill = "R2")
ggsave("./tuning/glmnet/gamma/tune.png", dpi = 450)

results.dg <- results.combos

trn <- data.matrix(trn.logk)
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

results.combos.dg <- results.combos
