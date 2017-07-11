# Libraries and Packages --------------------------------------------------

library(caret)
library(data.table)
library(e1071)
library(kernlab)
library(Matrix)
library(stringr)
library(tidyverse)

# Loading Data ------------------------------------------------------------

df <- readRDS("./DelG.df.RDS") %>%
  select(., -guest)
rf.ga <- readRDS("./genetic alg.RDS")
ga.final <- rf.ga$ga$final
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

sprse <- sparse.model.matrix(~., df)
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

# Manual Tuning -----------------------------------------------------------
# Creating a set of function to tune each parameter individually

# Cost - works for all kernels
validate.svm.cost <- function(data, nfolds, cost, kerneltype, seed) {
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
  return(data.frame(
    data = deparse(substitute(data)),
    nfolds = nfolds,
    seed = seed,
    kernel = kerneltype,
    gamma = g,
    rsquared = sum(results) / nfolds))
}
# Epsilon - works for all kernels
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
  return(data.frame(
    data = deparse(substitute(data)),
    nfolds = nfolds,
    seed = seed,
    kernel = kerneltype,
    epsilon = e,
    rsquared = sum(results) / nfolds))
}
# Constant coefficient - only for polynomial and sigmoid
validate.svm.coef <- function(data, nfolds, coef, kerneltype, seed) {
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

#     Polynomial Kernel ---------------------------------------------------
#         Seed 5 ----
# Peaks at 8
results.cost <-
  do.call(
    rbind,
    lapply(
      2 ^ (2:12),
      FUN = validate.svm.cost,
      data = sprse.trn,
      kerneltype = "polynomial",
      nfolds = 10,
      seed = 5
    )
  ) # Peak around 8

# Keeps reaching "max number of iterations", which may be a factor ^^^
results.gamma <-
  do.call(
    rbind,
    lapply(
      2 ^ (-5:8),
      FUN = validate.svm.gamma,
      data = sprse.trn,
      kerneltype = "polynomial",
      nfolds = 10,
      seed = 5
    )
  )

# Peaks around 1 or 2
results.epsilon <-
  do.call(
    rbind,
    lapply(
      2 ^ (-5:2),
      FUN = validate.svm.epsilon,
      data = sprse.trn,
      kerneltype = "polynomial",
      nfolds = 10,
      seed = 5
    )
  )
# peaks at 2
# probably has the most impact
results.coef <-
  do.call(
    rbind,
    lapply(
      2 ^ (-4:4),
      FUN = validate.svm.coef,
      data = sprse.trn,
      kerneltype = "polynomial",
      nfolds = 10,
      seed = 5
    )
  )

# Very little variation, possibly due to reaching max iterations
results.ga.cost <-
  do.call(
    rbind,
    lapply(
      2 ^ (-5:8),
      FUN = validate.svm.gamma,
      data = sprse.ga.trn,
      kerneltype = "polynomial",
      nfolds = 10,
      seed = 5
    )
  )
#         Seed 6 ----
# Peaks at 8
results2.cost <-
  do.call(
    rbind,
    lapply(
      2 ^ (2:12),
      FUN = validate.svm.cost,
      data = sprse.trn,
      kerneltype = "polynomial",
      nfolds = 10,
      seed = 6
    )
  ) 
results.cost.comp <- rbind(results.cost, results2.cost)

results2.gamma <-
  do.call(
    rbind,
    lapply(
      2 ^ (-5:8),
      FUN = validate.svm.gamma,
      data = sprse.trn,
      kerneltype = "polynomial",
      nfolds = 10,
      seed = 6
    )
  )
results.gamma.comp <- rbind(results.gamma, results2.gamma)

results2.epsilon <-
  do.call(
    rbind,
    lapply(
      2 ^ (-5:2),
      FUN = validate.svm.epsilon,
      data = sprse.trn,
      kerneltype = "polynomial",
      nfolds = 10,
      seed = 6
    )
  )
results.epsilon.comp <- rbind(results.epsilon, results2.epsilon)

results2.coef <-
  do.call(
    rbind,
    lapply(
      2 ^ (-4:4),
      FUN = validate.svm.coef,
      data = sprse.trn,
      kerneltype = "polynomial",
      nfolds = 10,
      seed = 6
    )
  )
results.coef.comp <- rbind(results.coef, results2.coef)

#         Seed 24 ----
results3.cost <-
  do.call(
    rbind,
    lapply(
      2 ^ (2:12),
      FUN = validate.svm.cost,
      data = sprse.trn,
      kerneltype = "polynomial",
      nfolds = 10,
      seed = 24
    )
  ) # 
results.cost.comp <- rbind(results.cost.comp, results3.cost)

results3.gamma <-
  do.call(
    rbind,
    lapply(
      2 ^ (-5:8),
      FUN = validate.svm.gamma,
      data = sprse.trn,
      kerneltype = "polynomial",
      nfolds = 10,
      seed = 24
    )
  )
results.gamma.comp <- rbind(results.gamma.comp, results3.gamma)

results3.epsilon <-
  do.call(
    rbind,
    lapply(
      2 ^ (-5:2),
      FUN = validate.svm.epsilon,
      data = sprse.trn,
      kerneltype = "polynomial",
      nfolds = 10,
      seed = 24
    )
  )
results.epsilon.comp <- rbind(results.epsilon.comp, results3.epsilon)

results3.coef <-
  do.call(
    rbind,
    lapply(
      2 ^ (-4:4),
      FUN = validate.svm.coef,
      data = sprse.trn,
      kerneltype = "polynomial",
      nfolds = 10,
      seed = 24
    )
  )
results.coef.comp <- rbind(results.coef.comp, results3.coef)

#         Seed 512 ----
results4.cost <-
  do.call(
    rbind,
    lapply(
      2 ^ (2:12),
      FUN = validate.svm.cost,
      data = sprse.trn,
      kerneltype = "polynomial",
      nfolds = 10,
      seed = 512
    )
  ) # 
results.cost.comp <- rbind(results.cost.comp, results4.cost)

results4.epsilon <-
  do.call(
    rbind,
    lapply(
      2 ^ (-5:2),
      FUN = validate.svm.epsilon,
      data = sprse.trn,
      kerneltype = "polynomial",
      nfolds = 10,
      seed = 512
    )
  )
results.epsilon.comp <- rbind(results.epsilon.comp, results4.epsilon)

results4.coef <-
  do.call(
    rbind,
    lapply(
      2 ^ (-4:4),
      FUN = validate.svm.coef,
      data = sprse.trn,
      kerneltype = "polynomial",
      nfolds = 10,
      seed = 512
    )
  )
results.coef.comp <- rbind(results.coef.comp, results4.coef)
#         Seed 333 ----
results5.cost <-
  do.call(
    rbind,
    lapply(
      2 ^ (2:12),
      FUN = validate.svm.cost,
      data = sprse.trn,
      kerneltype = "polynomial",
      nfolds = 10,
      seed = 333
    )
  ) # 
results.cost.comp <- rbind(results.cost.comp, results5.cost)

results5.epsilon <-
  do.call(
    rbind,
    lapply(
      2 ^ (-5:2),
      FUN = validate.svm.epsilon,
      data = sprse.trn,
      kerneltype = "polynomial",
      nfolds = 10,
      seed = 333
    )
  )
results.epsilon.comp <- rbind(results.epsilon.comp, results5.epsilon)

results5.coef <-
  do.call(
    rbind,
    lapply(
      2 ^ (-4:4),
      FUN = validate.svm.coef,
      data = sprse.trn,
      kerneltype = "polynomial",
      nfolds = 10,
      seed = 333
    )
  )
results.coef.comp <- rbind(results.coef.comp, results5.coef)
#         Seed 193 ----
results6.cost <-
  do.call(
    rbind,
    lapply(
      2 ^ (2:12),
      FUN = validate.svm.cost,
      data = sprse.trn,
      kerneltype = "polynomial",
      nfolds = 10,
      seed = 193
    )
  ) # 
results.cost.comp <- rbind(results.cost.comp, results6.cost)

results6.epsilon <-
  do.call(
    rbind,
    lapply(
      2 ^ (-5:2),
      FUN = validate.svm.epsilon,
      data = sprse.trn,
      kerneltype = "polynomial",
      nfolds = 10,
      seed = 193
    )
  )
results.epsilon.comp <- rbind(results.epsilon.comp, results6.epsilon)

results6.coef <-
  do.call(
    rbind,
    lapply(
      2 ^ (-4:4),
      FUN = validate.svm.coef,
      data = sprse.trn,
      kerneltype = "polynomial",
      nfolds = 10,
      seed = 193
    )
  )
results.coef.comp <- rbind(results.coef.comp, results6.coef)

#         Compilation -----------------------------------------------------

saveRDS(results.cost.comp, "./tuning/svm/poly.tune.cost.RDS")
saveRDS(results.gamma.comp, "./tuning/svm/poly.tune.gamma.RDS")
saveRDS(results.epsilon.comp, "./tuning/svm/poly.tune.epsilon.RDS")
saveRDS(results.coef.comp, "./tuning/svm/poly.tune.coef.RDS")

#     Radial Kernel -------------------------------------------------------
#         Seed 5 ----------------------------------------------------------
# Peaks at 8
results.rbf.cost <-
  do.call(
    rbind,
    lapply(
      2 ^ (0:12),
      FUN = validate.svm.cost,
      data = sprse.trn,
      kerneltype = "radial",
      nfolds = 10,
      seed = 5
    )
  ) 

# Peaks around 0.0005
results.rbf.gamma <-
  do.call(
    rbind,
    lapply(
      2 ^ (-14:1),
      FUN = validate.svm.gamma,
      data = sprse.trn,
      kerneltype = "radial",
      nfolds = 10,
      seed = 5
    )
  )

# Peaks around 1
results.rbf.epsilon <-
  do.call(
    rbind,
    lapply(
      2 ^ (-8:3),
      FUN = validate.svm.epsilon,
      data = sprse.trn,
      kerneltype = "radial",
      nfolds = 10,
      seed = 5
    )
  )

# Peaks around 0.0009765625 or 0.0001
results.rbf.ga.cost <-
  do.call(
    rbind,
    lapply(
      2 ^ (-14:-4),
      FUN = validate.svm.gamma,
      data = sprse.ga.trn,
      kerneltype = "radial",
      nfolds = 10,
      seed = 5
    )
  )


#         Seed 6 ----------------------------------------------------------
results2.rbf.cost <-
  do.call(
    rbind,
    lapply(
      2 ^ (0:12),
      FUN = validate.svm.cost,
      data = sprse.trn,
      kerneltype = "radial",
      nfolds = 10,
      seed = 6
    )
  ) 
rbf.cost.comp <- rbind(results.rbf.cost, results2.rbf.cost)

results2.rbf.gamma <-
  do.call(
    rbind,
    lapply(
      2 ^ (-14:1),
      FUN = validate.svm.gamma,
      data = sprse.trn,
      kerneltype = "radial",
      nfolds = 10,
      seed = 6
    )
  )
rbf.gamma.comp <- rbind(results.rbf.gamma, results2.rbf.gamma)

results2.rbf.epsilon <-
  do.call(
    rbind,
    lapply(
      2 ^ (-8:3),
      FUN = validate.svm.epsilon,
      data = sprse.trn,
      kerneltype = "radial",
      nfolds = 10,
      seed = 6
    )
  )
rbf.epsilon.comp <- rbind(results.rbf.epsilon, results2.rbf.epsilon)

#         Seed 24  --------------------------------------------------------

results3.rbf.cost <-
  do.call(
    rbind,
    lapply(
      2 ^ (0:12),
      FUN = validate.svm.cost,
      data = sprse.trn,
      kerneltype = "radial",
      nfolds = 10,
      seed = 24
    )
  ) 
rbf.cost.comp <- rbind(rbf.cost.comp, results3.rbf.cost)

results3.rbf.gamma <-
  do.call(
    rbind,
    lapply(
      2 ^ (-14:1),
      FUN = validate.svm.gamma,
      data = sprse.trn,
      kerneltype = "radial",
      nfolds = 10,
      seed = 24
    )
  )
rbf.gamma.comp <- rbind(rbf.gamma.comp, results3.rbf.gamma)

results3.rbf.epsilon <-
  do.call(
    rbind,
    lapply(
      2 ^ (-8:3),
      FUN = validate.svm.epsilon,
      data = sprse.trn,
      kerneltype = "radial",
      nfolds = 10,
      seed = 24
    )
  )
rbf.epsilon.comp <- rbind(rbf.epsilon.comp, results3.rbf.epsilon)

#         Seed 512  -------------------------------------------------------

results4.rbf.cost <-
  do.call(
    rbind,
    lapply(
      2 ^ (0:12),
      FUN = validate.svm.cost,
      data = sprse.trn,
      kerneltype = "radial",
      nfolds = 10,
      seed = 512
    )
  ) 
rbf.cost.comp <- rbind(rbf.cost.comp, results4.rbf.cost)

results4.rbf.gamma <-
  do.call(
    rbind,
    lapply(
      2 ^ (-14:1),
      FUN = validate.svm.gamma,
      data = sprse.trn,
      kerneltype = "radial",
      nfolds = 10,
      seed = 512
    )
  )
rbf.gamma.comp <- rbind(rbf.gamma.comp, results4.rbf.gamma)

results4.rbf.epsilon <-
  do.call(
    rbind,
    lapply(
      2 ^ (-8:3),
      FUN = validate.svm.epsilon,
      data = sprse.trn,
      kerneltype = "radial",
      nfolds = 10,
      seed = 512
    )
  )
rbf.epsilon.comp <- rbind(rbf.epsilon.comp, results4.rbf.epsilon)
#         Seed 333  -------------------------------------------------------

results5.rbf.cost <-
  do.call(
    rbind,
    lapply(
      2 ^ (0:12),
      FUN = validate.svm.cost,
      data = sprse.trn,
      kerneltype = "radial",
      nfolds = 10,
      seed = 333
    )
  ) 
rbf.cost.comp <- rbind(rbf.cost.comp, results5.rbf.cost)

results5.rbf.gamma <-
  do.call(
    rbind,
    lapply(
      2 ^ (-14:1),
      FUN = validate.svm.gamma,
      data = sprse.trn,
      kerneltype = "radial",
      nfolds = 10,
      seed = 333
    )
  )
rbf.gamma.comp <- rbind(rbf.gamma.comp, results5.rbf.gamma)

results5.rbf.epsilon <-
  do.call(
    rbind,
    lapply(
      2 ^ (-8:3),
      FUN = validate.svm.epsilon,
      data = sprse.trn,
      kerneltype = "radial",
      nfolds = 10,
      seed = 333
    )
  )
rbf.epsilon.comp <- rbind(rbf.epsilon.comp, results5.rbf.epsilon)
#         Seed 193  -------------------------------------------------------

results6.rbf.cost <-
  do.call(
    rbind,
    lapply(
      2 ^ (0:12),
      FUN = validate.svm.cost,
      data = sprse.trn,
      kerneltype = "radial",
      nfolds = 10,
      seed = 193
    )
  ) 
rbf.cost.comp <- rbind(rbf.cost.comp, results6.rbf.cost)

results6.rbf.gamma <-
  do.call(
    rbind,
    lapply(
      2 ^ (-14:1),
      FUN = validate.svm.gamma,
      data = sprse.trn,
      kerneltype = "radial",
      nfolds = 10,
      seed = 193
    )
  )
rbf.gamma.comp <- rbind(rbf.gamma.comp, results6.rbf.gamma)

results6.rbf.epsilon <-
  do.call(
    rbind,
    lapply(
      2 ^ (-8:3),
      FUN = validate.svm.epsilon,
      data = sprse.trn,
      kerneltype = "radial",
      nfolds = 10,
      seed = 193
    )
  )
rbf.epsilon.comp <- rbind(rbf.epsilon.comp, results6.rbf.epsilon)

#         Compilation -----------------------------------------------------
saveRDS(rbf.cost.comp, "./tuning/svm/rbf.tune.cost.RDS")
saveRDS(rbf.gamma.comp, "./tuning/svm/rbf.tune.gamma.RDS")
saveRDS(rbf.epsilon.comp, "./tuning/svm/rbf.tune.epsilon.RDS")

#     Sigmoid kernel ----------------------------------------------------------
# Cost = 2
results.sig.cost <-
  do.call(
    rbind,
    lapply(
      2 ^ (-6:6),
      FUN = validate.svm.cost,
      data = sprse.trn,
      kerneltype = "sigmoid",
      nfolds = 10,
      seed = 5
    )
  ) 

# Peaks around 0.0005
results.sig.gamma <-
  do.call(
    rbind,
    lapply(
      2 ^ (-12:1),
      FUN = validate.svm.gamma,
      data = sprse.trn,
      kerneltype = "sigmoid",
      nfolds = 10,
      seed = 5
    )
  )

# Epsilon = 0.25
results.sig.epsilon <-
  do.call(
    rbind,
    lapply(
      2 ^ (-8:3),
      FUN = validate.svm.epsilon,
      data = sprse.trn,
      kerneltype = "sigmoid",
      nfolds = 10,
      seed = 5
    )
  )

# Peaks around 0.004
results.sig.ga.cost <-
  do.call(
    rbind,
    lapply(
      2 ^ (-14:-4),
      FUN = validate.svm.gamma,
      data = sprse.ga.trn,
      kerneltype = "sigmoid",
      nfolds = 10,
      seed = 5
    )
  )


# Plots -------------------------------------------------------------------
#     Polynomial graphs ---------------------------------------------------
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
ggsave(filename = "./tuning/svm/2017-07-11 poly cost tune.png")

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
ggsave(filename = "./tuning/svm/2017-07-11 poly gamma tune.png")

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
ggsave("./tuning/svm/2017-07-11 poly epsilon tune.png")

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
ggsave("./tuning/svm/2017-07-11 poly coef tune.png")

#     RBF graphs ----------------------------------------------------------

temp.data <- rbf.cost.comp
temp.data$seed <- as.factor(temp.data$seed)
ggplot(data = temp.data, aes(x = cost, y = rsquared, 
                             group = seed, color = seed)) + 
  geom_line(size = 1) + 
  scale_x_continuous(trans = 'log2') + 
  labs(x = "Cost", y = "R-squared", 
       title = "Radial Basis Kernel - Cost", 
       color = "Random Seed") + 
  theme_bw()
ggsave(filename = "./tuning/svm/2017-07-11 rbf cost tune.png")

temp.data <- rbf.gamma.comp
temp.data$seed <- as.factor(temp.data$seed)
ggplot(data = temp.data, aes(x = gamma, y = rsquared,
                             group = seed, color = seed)) + 
  geom_line(size = 1) + 
  scale_x_continuous(trans = 'log2') + 
  labs(x = "Gamma", y = "R-squared", 
       title = "Radial Basis Kernel - Gamma", 
       color = "Random Seed") + 
  theme_bw()
ggsave(filename = "./tuning/svm/2017-07-11 rbf gamma tune.png")

temp.data <- rbf.epsilon.comp
temp.data$seed <- as.factor(temp.data$seed)
ggplot(data = temp.data, aes(x = epsilon, y = rsquared, 
                             group = seed, color = seed)) + 
  geom_line(size = 1) + 
  scale_x_continuous(trans = 'log2') + 
  labs(x = "Epsilon", y = "R-squared", 
       title = "Radial Basis Kernel - Epsilon", 
       color = "Random Seed") + 
  theme_bw()
ggsave("./tuning/svm/2017-07-11 rbf epsilon tune.png")
# Final parameters:
# for polynomial:coef0 = 2, gamma = 0.0625, cost = 8, epsilon = 1
# for radial: cost = 8, gamma = 0.0005, epsilon = 1

#     Sigmoid graphs ----------------------------------------------------------
ggplot(data = results.sig.cost, aes(x = cost, y = rsquared)) + 
  geom_line() + 
  scale_x_continuous(trans = 'log2') + 
  theme_bw()

ggplot(data = results.sig.gamma, aes(x = gamma, y = rsquared)) + 
  geom_line() + 
  scale_x_continuous(trans = 'log2') + 
  theme_bw()

ggplot(data = results.sig.epsilon, aes(x = epsilon, y = rsquared)) + 
  geom_line() + 
  scale_x_continuous(trans = 'log2') + 
  theme_bw()

ggplot(data = results.sig.ga.cost, aes(x = gamma, y = rsquared)) + 
  geom_line() + 
  scale_x_continuous(trans = 'log2') + 
  theme_bw()
