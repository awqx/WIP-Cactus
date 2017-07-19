# Libraries and Packages --------------------------------------------------

library(caret)
library(data.table)
library(e1071)
library(kernlab)
library(Matrix)
library(stringr)
library(tidyverse)

# Loading Data ------------------------------------------------------------

setwd("~/SREP LAB/qsar")
dir.create("./tuning")
dir.create("./tuning/svm")
# Reading data with all descriptors
df.raw <- readRDS("./padel.pp.new.RDS") 
df <- df.raw %>%
  select(., -guest:-host) %>%
  select(., -data.source)

# df.ri <- df.raw %>% 
#   select(., -guest:-host) %>%
#   filter(data.source == "rekharsky.inoue") %>%
#   select(., -data.source)
# df.suz <- df.raw %>% 
#   select(., -guest:-host) %>%
#   filter(data.source == "suzuki") %>%
#   select(., -data.source)

# Reading data selected via genetic alg
ga.pred <- readRDS("./feature.select/GAFS names.RDS")
ga <- df.raw[ , colnames(df.raw) %in% ga.pred]
ga <- cbind(df.raw$DelG, ga) %>%
  rename(DelG = `df.raw$DelG`)

# Creating train-test splits on the dataset
set.seed(25)
trn.ind <- sample(x = 1:nrow(df), size = round(0.7 * nrow(df)))
trn <- df[trn.ind, ]
tst <- df[-trn.ind, ]

ga.trn <- ga[trn.ind, ]
ga.tst <- ga[-trn.ind, ]

sprse <- sparse.model.matrix(~., df)
trn.ind <- sample(x = 1:nrow(sprse), size = round(0.7 * nrow(sprse)))
sprse.trn <- sprse[trn.ind, ]
sprse.tst <- sprse[-trn.ind, ]

sprse.ga <- sparse.model.matrix(~., ga)
sprse.ga.trn <- sprse.ga[trn.ind, ]
sprse.ga.tst <- sprse.ga[-trn.ind, ]

# Tuning Functions --------------------------------------------------------

# Creating a set of function to tune each parameter individually

# Cost - works for all kernels
tune.svm.cost <- function(data, nfolds, cost, kerneltype, seed) {
  set.seed(seed)
  fold.list <- createFolds(y = data[ , 2], k = nfolds)
  results <- c(rep(0.0, nfolds))
  # Cross-validation using k-fold cross-validation
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
    # Necessary to rename columns so that caret:defaultSummary works
    colnames(svm.df)[1] <- "pred"
    colnames(svm.df)[2] <- "obs"
    R2 <- defaultSummary(svm.df)[2] 
    results[i] <- R2
  }
  return(data.frame( # Useful for records
    data = deparse(substitute(data)), # Turns the variable name into char
    nfolds = nfolds,
    seed = seed,
    kernel = kerneltype,
    cost = cost,
    rsquared = sum(results) / nfolds))
}
# Gamma - works for all kernels except linear 
tune.svm.gamma <- function(data, nfolds, g, kerneltype, seed) {
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
  return(data.frame(
    data = deparse(substitute(data)),
    nfolds = nfolds,
    seed = seed,
    kernel = kerneltype,
    gamma = g,
    rsquared = sum(results) / nfolds))
}
# Epsilon - works for all kernels
tune.svm.epsilon <- function(data, nfolds, e, kerneltype, seed) {
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
  return(data.frame(
    data = deparse(substitute(data)),
    nfolds = nfolds,
    seed = seed,
    kernel = kerneltype,
    epsilon = e,
    rsquared = sum(results) / nfolds))
}
# Constant coefficient - only for polynomial and sigmoid
tune.svm.coef <- function(data, nfolds, coef, kerneltype, seed) {
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
                  coef0 = coef,
                  kernel = kerneltype)
    svm.df <- predict(svm.cv, tst[ , -1:-2]) %>% 
      cbind(tst[ , 2]) %>%
      data.frame()
    colnames(svm.df)[1] <- "pred"
    colnames(svm.df)[2] <- "obs"
    R2 <- defaultSummary(svm.df)[2] 
    results[i] <- R2
  }
  return(data.frame(
    data = deparse(substitute(data)),
    nfolds = nfolds,
    seed = seed,
    kernel = kerneltype,
    coef = coef,
    rsquared = sum(results) / nfolds))
}

# All Predictors ----------------------------------------------------------

#     Seed 1 ----
results.cost <-
  do.call(
    rbind,
    lapply(
      2 ^ (2:12),
      FUN = tune.svm.cost,
      data = sprse.trn,
      kerneltype = "polynomial",
      nfolds = 10,
      seed = 1
    )
  ) 

results.gamma <-
  do.call(
    rbind,
    lapply(
      2 ^ (-12:1),
      FUN = tune.svm.gamma,
      data = sprse.trn,
      kerneltype = "polynomial",
      nfolds = 10,
      seed = 1
    )
  )

results.epsilon <-
  do.call(
    rbind,
    lapply(
      2 ^ (-5:2),
      FUN = tune.svm.epsilon,
      data = sprse.trn,
      kerneltype = "polynomial",
      nfolds = 10,
      seed = 1
    )
  )

results.coef <-
  do.call(
    rbind,
    lapply(
      2 ^ (-4:4),
      FUN = tune.svm.coef,
      data = sprse.trn,
      kerneltype = "polynomial",
      nfolds = 10,
      seed = 1
    )
  )

#     Seed 2 ----
# Peaks at 8
results2.cost <-
  do.call(
    rbind,
    lapply(
      2 ^ (2:12),
      FUN = tune.svm.cost,
      data = sprse.trn,
      kerneltype = "polynomial",
      nfolds = 10,
      seed = 2
    )
  ) 
results.cost.comp <- rbind(results.cost, results2.cost)

results2.gamma <-
  do.call(
    rbind,
    lapply(
      2 ^ (-12:1),
      FUN = tune.svm.gamma,
      data = sprse.trn,
      kerneltype = "polynomial",
      nfolds = 10,
      seed = 2
    )
  )
results.gamma.comp <- rbind(results.gamma, results2.gamma)

results2.epsilon <-
  do.call(
    rbind,
    lapply(
      2 ^ (-5:2),
      FUN = tune.svm.epsilon,
      data = sprse.trn,
      kerneltype = "polynomial",
      nfolds = 10,
      seed = 2
    )
  )
results.epsilon.comp <- rbind(results.epsilon, results2.epsilon)

results2.coef <-
  do.call(
    rbind,
    lapply(
      2 ^ (-4:4),
      FUN = tune.svm.coef,
      data = sprse.trn,
      kerneltype = "polynomial",
      nfolds = 10,
      seed = 2
    )
  )
results.coef.comp <- rbind(results.coef, results2.coef)

#     Seed 3 ----
results3.cost <-
  do.call(
    rbind,
    lapply(
      2 ^ (2:12),
      FUN = tune.svm.cost,
      data = sprse.trn,
      kerneltype = "polynomial",
      nfolds = 10,
      seed = 3
    )
  ) # 
results.cost.comp <- rbind(results.cost.comp, results3.cost)

results3.gamma <-
  do.call(
    rbind,
    lapply(
      2 ^ (-12:1),
      FUN = tune.svm.gamma,
      data = sprse.trn,
      kerneltype = "polynomial",
      nfolds = 10,
      seed = 3
    )
  )
results.gamma.comp <- rbind(results.gamma.comp, results3.gamma)

results3.epsilon <-
  do.call(
    rbind,
    lapply(
      2 ^ (-5:2),
      FUN = tune.svm.epsilon,
      data = sprse.trn,
      kerneltype = "polynomial",
      nfolds = 10,
      seed = 3
    )
  )
results.epsilon.comp <- rbind(results.epsilon.comp, results3.epsilon)

results3.coef <-
  do.call(
    rbind,
    lapply(
      2 ^ (-4:4),
      FUN = tune.svm.coef,
      data = sprse.trn,
      kerneltype = "polynomial",
      nfolds = 10,
      seed = 3
    )
  )
results.coef.comp <- rbind(results.coef.comp, results3.coef)

#     Seed 4 ----
results4.cost <-
  do.call(
    rbind,
    lapply(
      2 ^ (2:12),
      FUN = tune.svm.cost,
      data = sprse.trn,
      kerneltype = "polynomial",
      nfolds = 10,
      seed = 4
    )
  ) # 
results.cost.comp <- rbind(results.cost.comp, results4.cost)

results4.epsilon <-
  do.call(
    rbind,
    lapply(
      2 ^ (-5:2),
      FUN = tune.svm.epsilon,
      data = sprse.trn,
      kerneltype = "polynomial",
      nfolds = 10,
      seed = 4
    )
  )
results.epsilon.comp <- rbind(results.epsilon.comp, results4.epsilon)

results4.coef <-
  do.call(
    rbind,
    lapply(
      2 ^ (-4:4),
      FUN = tune.svm.coef,
      data = sprse.trn,
      kerneltype = "polynomial",
      nfolds = 10,
      seed = 4
    )
  )
results.coef.comp <- rbind(results.coef.comp, results4.coef)
#     Seed 5 ----
results5.cost <-
  do.call(
    rbind,
    lapply(
      2 ^ (2:12),
      FUN = tune.svm.cost,
      data = sprse.trn,
      kerneltype = "polynomial",
      nfolds = 10,
      seed = 5
    )
  ) # 
results.cost.comp <- rbind(results.cost.comp, results5.cost)

results5.epsilon <-
  do.call(
    rbind,
    lapply(
      2 ^ (-5:2),
      FUN = tune.svm.epsilon,
      data = sprse.trn,
      kerneltype = "polynomial",
      nfolds = 10,
      seed = 5
    )
  )
results.epsilon.comp <- rbind(results.epsilon.comp, results5.epsilon)

results5.coef <-
  do.call(
    rbind,
    lapply(
      2 ^ (-4:4),
      FUN = tune.svm.coef,
      data = sprse.trn,
      kerneltype = "polynomial",
      nfolds = 10,
      seed = 5
    )
  )
results.coef.comp <- rbind(results.coef.comp, results5.coef)
#     Seed 6 ----
results6.cost <-
  do.call(
    rbind,
    lapply(
      2 ^ (2:12),
      FUN = tune.svm.cost,
      data = sprse.trn,
      kerneltype = "polynomial",
      nfolds = 10,
      seed = 7
    )
  ) # 
results.cost.comp <- rbind(results.cost.comp, results6.cost)

results6.epsilon <-
  do.call(
    rbind,
    lapply(
      2 ^ (-5:2),
      FUN = tune.svm.epsilon,
      data = sprse.trn,
      kerneltype = "polynomial",
      nfolds = 10,
      seed = 7
    )
  )
results.epsilon.comp <- rbind(results.epsilon.comp, results6.epsilon)

results6.coef <-
  do.call(
    rbind,
    lapply(
      2 ^ (-4:4),
      FUN = tune.svm.coef,
      data = sprse.trn,
      kerneltype = "polynomial",
      nfolds = 10,
      seed = 7
    )
  )
results.coef.comp <- rbind(results.coef.comp, results6.coef)

#     Compilation -----------------------------------------------------

saveRDS(results.cost.comp, "./tuning/svm/poly.tune.cost.RDS")
saveRDS(results.gamma.comp, "./tuning/svm/poly.tune.gamma.RDS")
saveRDS(results.epsilon.comp, "./tuning/svm/poly.tune.epsilon.RDS")
saveRDS(results.coef.comp, "./tuning/svm/poly.tune.coef.RDS")

#####
# GAFS Predictors ---------------------------------------------------------

#     Seed 1 ----
results.ga.cost <-
  do.call(
    rbind,
    lapply(
      2 ^ (2:12),
      FUN = tune.svm.cost,
      data = sprse.ga.trn,
      kerneltype = "polynomial",
      nfolds = 10,
      seed = 1
    )
  ) 

results.ga.gamma <-
  do.call(
    rbind,
    lapply(
      2 ^ (-12:1),
      FUN = tune.svm.gamma,
      data = sprse.ga.trn,
      kerneltype = "polynomial",
      nfolds = 10,
      seed = 1
    )
  )

results.ga.epsilon <-
  do.call(
    rbind,
    lapply(
      2 ^ (-5:2),
      FUN = tune.svm.epsilon,
      data = sprse.ga.trn,
      kerneltype = "polynomial",
      nfolds = 10,
      seed = 1
    )
  )

results.ga.cost <-
  do.call(
    rbind,
    lapply(
      2 ^ (-4:4),
      FUN = tune.svm.coef,
      data = sprse.ga.trn,
      kerneltype = "polynomial",
      nfolds = 10,
      seed = 1
    )
  )

#     Seed 2 ----
# Peaks at 8
results.ga2.cost <-
  do.call(
    rbind,
    lapply(
      2 ^ (2:12),
      FUN = tune.svm.cost,
      data = sprse.ga.trn,
      kerneltype = "polynomial",
      nfolds = 10,
      seed = 2
    )
  ) 

results.ga2.gamma <-
  do.call(
    rbind,
    lapply(
      2 ^ (-12:1),
      FUN = tune.svm.gamma,
      data = sprse.ga.trn,
      kerneltype = "polynomial",
      nfolds = 10,
      seed = 2
    )
  )

results.ga2.epsilon <-
  do.call(
    rbind,
    lapply(
      2 ^ (-5:2),
      FUN = tune.svm.epsilon,
      data = sprse.ga.trn,
      kerneltype = "polynomial",
      nfolds = 10,
      seed = 2
    )
  )

results.ga2.coef <-
  do.call(
    rbind,
    lapply(
      2 ^ (-4:4),
      FUN = tune.svm.coef,
      data = sprse.ga.trn,
      kerneltype = "polynomial",
      nfolds = 10,
      seed = 2
    )
  )

#     Seed 3 ----
results.ga3.cost <-
  do.call(
    rbind,
    lapply(
      2 ^ (2:12),
      FUN = tune.svm.cost,
      data = sprse.ga.trn,
      kerneltype = "polynomial",
      nfolds = 10,
      seed = 3
    )
  ) 

results.ga3.gamma <-
  do.call(
    rbind,
    lapply(
      2 ^ (-12:1),
      FUN = tune.svm.gamma,
      data = sprse.ga.trn,
      kerneltype = "polynomial",
      nfolds = 10,
      seed = 3
    )
  )

results.ga3.epsilon <-
  do.call(
    rbind,
    lapply(
      2 ^ (-5:2),
      FUN = tune.svm.epsilon,
      data = sprse.ga.trn,
      kerneltype = "polynomial",
      nfolds = 10,
      seed = 3
    )
  )

results.ga3.coef <-
  do.call(
    rbind,
    lapply(
      2 ^ (-4:4),
      FUN = tune.svm.coef,
      data = sprse.ga.trn,
      kerneltype = "polynomial",
      nfolds = 10,
      seed = 3
    )
  )

#     Seed 4 ----
results.ga4.cost <-
  do.call(
    rbind,
    lapply(
      2 ^ (2:12),
      FUN = tune.svm.cost,
      data = sprse.ga.trn,
      kerneltype = "polynomial",
      nfolds = 10,
      seed = 4
    )
  )

results.ga4.epsilon <-
  do.call(
    rbind,
    lapply(
      2 ^ (-5:2),
      FUN = tune.svm.epsilon,
      data = sprse.ga.trn,
      kerneltype = "polynomial",
      nfolds = 10,
      seed = 4
    )
  )

results.ga4.coef <-
  do.call(
    rbind,
    lapply(
      2 ^ (-4:4),
      FUN = tune.svm.coef,
      data = sprse.ga.trn,
      kerneltype = "polynomial",
      nfolds = 10,
      seed = 4
    )
  )

#     Seed 5 ----
results.ga5.cost <-
  do.call(
    rbind,
    lapply(
      2 ^ (2:12),
      FUN = tune.svm.cost,
      data = sprse.ga.trn,
      kerneltype = "polynomial",
      nfolds = 10,
      seed = 5
    )
  ) 

results.ga5.epsilon <-
  do.call(
    rbind,
    lapply(
      2 ^ (-5:2),
      FUN = tune.svm.epsilon,
      data = sprse.ga.trn,
      kerneltype = "polynomial",
      nfolds = 10,
      seed = 5
    )
  )

results.ga5.coef <-
  do.call(
    rbind,
    lapply(
      2 ^ (-4:4),
      FUN = tune.svm.coef,
      data = sprse.ga.trn,
      kerneltype = "polynomial",
      nfolds = 10,
      seed = 5
    )
  )
#     Seed 6 ----
results.ga6.cost <-
  do.call(
    rbind,
    lapply(
      2 ^ (2:12),
      FUN = tune.svm.cost,
      data = sprse.ga.trn,
      kerneltype = "polynomial",
      nfolds = 10,
      seed = 7
    )
  ) # 

results.ga6.epsilon <-
  do.call(
    rbind,
    lapply(
      2 ^ (-5:2),
      FUN = tune.svm.epsilon,
      data = sprse.ga.trn,
      kerneltype = "polynomial",
      nfolds = 10,
      seed = 7
    )
  )

results.ga6.coef <-
  do.call(
    rbind,
    lapply(
      2 ^ (-4:4),
      FUN = tune.svm.coef,
      data = sprse.ga.trn,
      kerneltype = "polynomial",
      nfolds = 10,
      seed = 7
    )
  )

#     Compilation -----------------------------------------------------
ga.coef.comp <- rbind(results.ga.coef, results.ga2.coef, 
                      results.ga3.coef, results.ga4.coef, 
                      results.ga5.coef, results.ga6.coef)
ga.cost.comp <- rbind(results.ga.cost, results.ga2.cost, 
                      results.ga3.cost, results.ga4.cost, 
                      results.ga5.cost, results.ga6.cost)
ga.gamma.comp <- rbind(results.ga.gamma, results.ga2.gamma, 
                      results.ga3.gamma)
ga.epsilon.comp <- rbind(results.ga.epsilon, results.ga2.epsilon, 
                      results.ga3.epsilon, results.ga4.epsilon, 
                      results.ga5.epsilon, results.ga6.epsilon)


saveRDS(ga.cost.comp, "./tuning/svm/ga.tune.cost.RDS")
saveRDS(ga.gamma.comp, "./tuning/svm/ga.tune.gamma.RDS")
saveRDS(ga.epsilon.comp, "./tuning/svm/ga.tune.epsilon.RDS")
saveRDS(ga.coef.comp, "./tuning/svm/ga.tune.coef.RDS")

# Plots -------------------------------------------------------------------
#     All Predictors ---------------------------------------------------
# Creating temporary dataframes in case more data points are added. Coercing the
# random seeds to factors prevents addition of new data, or at least makes it
# more difficult.
temp.data <- results.cost.comp
temp.data$seed <- as.factor(temp.data$seed)
ggplot(data = temp.data, aes(x = cost, y = rsquared, 
                             group = seed, color = seed)) + 
  geom_line(size = 1) + 
  scale_x_continuous(trans = 'log2') + 
  labs(x = "Cost", y = "R-squared", 
       title = "Polynomial Kernel - Cost", 
       color = "Random Seed") + 
  theme_bw()
ggsave(filename = "./tuning/svm/2017-07-18 poly cost tune.png")

temp.data <- results.gamma.comp
temp.data$seed <- as.factor(temp.data$seed)
ggplot(data = temp.data, aes(x = gamma, y = rsquared,
                             group = seed, color = seed)) + 
  geom_line(size = 1) + 
  scale_x_continuous(trans = 'log2') + 
  labs(x = "Gamma", y = "R-squared", 
       title = "Polynomial Kernel - Gamma", 
       color = "Random Seed") + 
  theme_bw()
ggsave(filename = "./tuning/svm/2017-07-18 poly gamma tune.png")

temp.data <- results.epsilon.comp
temp.data$seed <- as.factor(temp.data$seed)
ggplot(data = temp.data, aes(x = epsilon, y = rsquared, 
                             group = seed, color = seed)) + 
  geom_line(size = 1) + 
  scale_x_continuous(trans = 'log2') + 
  labs(x = "Epsilon", y = "R-squared", 
       title = "Polynomial Kernel - Epsilon", 
       color = "Random Seed") + 
  theme_bw()
ggsave("./tuning/svm/2017-07-18 poly epsilon tune.png")

temp.data <- results.coef.comp
temp.data$seed <- as.factor(temp.data$seed)
ggplot(data = temp.data, aes(x = coef, y = rsquared, 
                             group = seed, color = seed)) + 
  geom_line(size = 1) + 
  labs(x = "Constant Coefficient", y = "R-squared", 
       title = "Polynomial Kernel - Coefficient", 
       color = "Random Seed") + 
  scale_x_continuous(trans = 'log2') + 
  theme_bw()
ggsave("./tuning/svm/2017-07-18 poly coef tune.png")


#     GAFS ----------------------------------------------

temp.data <- ga.cost.comp
temp.data$seed <- as.factor(temp.data$seed)
ggplot(data = temp.data, aes(x = cost, y = rsquared, 
                             group = seed, color = seed)) + 
  geom_line(size = 1) + 
  scale_x_continuous(trans = 'log2') + 
  labs(x = "Cost", y = "R-squared", 
       title = "GAFS Polynomial Kernel - Cost", 
       color = "Random Seed") + 
  theme_bw()
ggsave(filename = "./tuning/svm/2017-07-18 ga cost tune.png")

temp.data <- ga.gamma.comp
temp.data$seed <- as.factor(temp.data$seed)
ggplot(data = temp.data, aes(x = gamma, y = rsquared,
                             group = seed, color = seed)) + 
  geom_line(size = 1) + 
  scale_x_continuous(trans = 'log2') + 
  labs(x = "Gamma", y = "R-squared", 
       title = "GAFS Polynomial Kernel - Gamma", 
       color = "Random Seed") + 
  theme_bw()
ggsave(filename = "./tuning/svm/2017-07-18 ga gamma tune.png")

temp.data <- ga.epsilon.comp
temp.data$seed <- as.factor(temp.data$seed)
ggplot(data = temp.data, aes(x = epsilon, y = rsquared, 
                             group = seed, color = seed)) + 
  geom_line(size = 1) + 
  scale_x_continuous(trans = 'log2') + 
  labs(x = "Epsilon", y = "R-squared", 
       title = "GAFS Polynomial Kernel - Epsilon", 
       color = "Random Seed") + 
  theme_bw()
ggsave("./tuning/svm/2017-07-18 ga epsilon tune.png")

temp.data <- ga.coef.comp
temp.data$seed <- as.factor(temp.data$seed)
ggplot(data = temp.data, aes(x = coef, y = rsquared, 
                             group = seed, color = seed)) + 
  geom_line(size = 1) + 
  labs(x = "Constant Coefficient", y = "R-squared", 
       title = "GAFS Polynomial Kernel - Coefficient", 
       color = "Random Seed") + 
  scale_x_continuous(trans = 'log2') + 
  theme_bw()
ggsave("./tuning/svm/2017-07-18 ga coef tune.png")
