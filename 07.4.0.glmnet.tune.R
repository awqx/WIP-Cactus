# Libraries ---------------------------------------------------------------

library(caret)
library(ggplot2)
library(glmnet)
library(Matrix)
library(tidyverse)


# Loading Data ------------------------------------------------------------

dir.create("./tuning/glmnet")

df.raw <- readRDS("./data/padel.pp.RDS")
df <- df.raw %>% select(-guest:-host) %>%
  select(-data.source)
mat <- sparse.model.matrix(~., df)

# Functions ---------------------------------------------------------------

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

# Tuning ------------------------------------------------------------------

#     Seed 1 --------------------------------------------------------------

glm.alpha <- do.call(rbind, lapply(seq(0, 1, 0.1),
                                   FUN = tune.glm.alpha,
                                   data = mat, nfolds = 10, seed = 1))

glm.dfmax <- do.call(rbind, lapply(2 ^ (2:9),
                                   FUN = tune.glm.dfmax,
                                   data = mat, nfolds = 10, seed = 1))

#     Seed 2 --------------------------------------------------------------

glm2.alpha <- do.call(rbind, lapply(seq(0, 1, 0.1),
                                   FUN = tune.glm.alpha,
                                   data = mat, nfolds = 10, seed = 2))

glm2.dfmax <- do.call(rbind, lapply(2 ^ (2:9),
                                   FUN = tune.glm.dfmax,
                                   data = mat, nfolds = 10, seed = 2))

#     Seed 3 --------------------------------------------------------------

glm3.alpha <- do.call(rbind, lapply(seq(0, 1, 0.1),
                                    FUN = tune.glm.alpha,
                                    data = mat, nfolds = 10, seed = 3))

glm3.dfmax <- do.call(rbind, lapply(2 ^ (2:9),
                                    FUN = tune.glm.dfmax,
                                    data = mat, nfolds = 10, seed = 3))

#     Seed 4 --------------------------------------------------------------

glm4.alpha <- do.call(rbind, lapply(seq(0, 1, 0.1),
                                    FUN = tune.glm.alpha,
                                    data = mat, nfolds = 10, seed = 4))

glm4.dfmax <- do.call(rbind, lapply(2 ^ (2:9),
                                    FUN = tune.glm.dfmax,
                                    data = mat, nfolds = 10, seed = 4))

#     Seed 5 --------------------------------------------------------------

glm5.alpha <- do.call(rbind, lapply(seq(0, 1, 0.1),
                                    FUN = tune.glm.alpha,
                                    data = mat, nfolds = 10, seed = 5))

glm5.dfmax <- do.call(rbind, lapply(2 ^ (2:9),
                                    FUN = tune.glm.dfmax,
                                    data = mat, nfolds = 10, seed = 5))

#     Seed 6 --------------------------------------------------------------

glm6.alpha <- do.call(rbind, lapply(seq(0, 1, 0.1),
                                    FUN = tune.glm.alpha,
                                    data = mat, nfolds = 10, seed = 6))

glm6.dfmax <- do.call(rbind, lapply(2 ^ (2:9),
                                    FUN = tune.glm.dfmax,
                                    data = mat, nfolds = 10, seed = 6))

#     Compilation ---------------------------------------------------------

alpha.comp <- rbind(glm.alpha, glm2.alpha, glm3.alpha, 
                    glm4.alpha, glm5.alpha, glm6.alpha)

dfmax.comp <- rbind(glm.dfmax, glm2.dfmax, glm3.dfmax, 
                    glm4.dfmax, glm5.dfmax, glm6.dfmax)

saveRDS(alpha.comp, "./tuning/glmnet/alpha.results.RDS")
saveRDS(dfmax.comp, "./tuning/glmnet/dfmax.results.RDS")

# Graphs ------------------------------------------------------------------

# Alpha
temp <- alpha.comp
temp$seed <- as.factor(temp$seed)
ggplot(temp, aes(x = alpha, y = rsquared, group = seed, color = seed)) + 
  geom_line() + 
  theme_bw() + 
  labs(x = "Alpha", y = "R-squared", 
       title = "GLMnet - Tuning Alpha", 
       color = "Random Seed")
ggsave("./tuning/glmnet/2017-07-20 glmnet alpha.png")

# DFMax
temp <- dfmax.comp
temp$seed <- as.factor(temp$seed)
ggplot(temp, aes(x = dfmax, y = rsquared, group = seed, color = seed)) + 
  geom_line() + 
  theme_bw() + 
  scale_x_continuous(trans = "log2") + 
  labs(x = "Max Number of Variables", y = "R-squared", 
       title = "GLMnet - Tuning Number of Variables", 
       color = "Random Seed")
ggsave("./tuning/glmnet/2017-07-20 glmnet dfmax.png")
