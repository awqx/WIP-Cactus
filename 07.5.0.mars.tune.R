# MARS - Multivariate adaptive regression splines
dir.create("./models/mars")

# Libraries ---------------------------------------------------------------

if(!require("pacman"))
  install.packages("pacman")
library(pacman)
p_load(caret, tidyverse, earth)
source("./07.model.functions.R")

# Functions ---------------------------------------------------------------

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

# Alpha -------------------------------------------------------------------

#     Data organization ---------------------------------------------------

trn.all <- readRDS("./pre-process/alpha/1/pp.RDS") 
colnames(trn.all) <- str_replace(colnames(trn.all), "-", ".")
trn.guest <- trn.all$guest
trn <- select(trn.all, -guest)

features <- readRDS("./feature.selection/alpha.vars.RDS")
trn <- trn[ , colnames(trn) %in% c("DelG", features)]

dir.create("./tuning/mars")

#     Estimation ----------------------------------------------------------

# Degree ---
deg.range <- c(1:5)
results1.deg <- do.call(rbind, lapply(deg.range, FUN = tune.mars.deg, 
                                      data = trn, nfolds = 10, seed = 101))
results2.deg <- do.call(rbind, lapply(deg.range, FUN = tune.mars.deg, 
                                      data = trn, nfolds = 10, seed = 102))
results3.deg <- do.call(rbind, lapply(deg.range, FUN = tune.mars.deg, 
                                      data = trn, nfolds = 10, seed = 103))
results.deg <- rbind(results1.deg, results2.deg, results3.deg) %>%
  mutate(seed = as.factor(seed))
ggplot(results.deg, aes(x = degree, y = rsquared, color = seed)) + 
       geom_line()

# Penalty ---
pen.range <- c(0, 2, 4, 6, 8)
results1.pen <- do.call(rbind, lapply(pen.range, FUN = tune.mars.pen, 
                                      data = trn, nfolds = 10, seed = 101))
results2.pen <- do.call(rbind, lapply(pen.range, FUN = tune.mars.pen, 
                                      data = trn, nfolds = 10, seed = 102))
results3.pen <- do.call(rbind, lapply(pen.range, FUN = tune.mars.pen, 
                                      data = trn, nfolds = 10, seed = 103))
results.pen <- rbind(results1.pen, results2.pen, results3.pen) %>%
  mutate(seed = as.factor(seed))
ggplot(results.pen, aes(x = penalty, y = rsquared, color = seed)) + 
  geom_line() + 
  theme_bw()

# nk: Model terms
nk.range <- c(10, 15, 25, 35, 50, 60)
results1.nk <- do.call(rbind, lapply(nk.range, FUN = tune.mars.nk, 
                                      data = trn, nfolds = 10, seed = 101))
results2.nk <- do.call(rbind, lapply(nk.range, FUN = tune.mars.nk, 
                                      data = trn, nfolds = 10, seed = 102))
results3.nk <- do.call(rbind, lapply(nk.range, FUN = tune.mars.nk, 
                                      data = trn, nfolds = 10, seed = 103))
results.nk <- rbind(results1.nk, results2.nk, results3.nk) %>%
  mutate(seed = as.factor(seed))
ggplot(results.nk, aes(x = nk, y = rsquared, color = seed)) + 
  geom_line() + 
  theme_bw()

# Threshold
thresh.range <- c(0, 0.01, 0.025, 0.05, 0.1, 0.2)
results1.thresh <- do.call(rbind, lapply(thresh.range, 
                                         FUN = tune.mars.thresh, 
                                         data = trn, nfolds = 10, 
                                         seed = 101))
results2.thresh <- do.call(rbind, lapply(thresh.range, 
                                         FUN = tune.mars.thresh, 
                                         data = trn, nfolds = 10, 
                                         seed = 102))
results3.thresh <- do.call(rbind, lapply(thresh.range, 
                                         FUN = tune.mars.thresh,
                                         data = trn, nfolds = 10, 
                                         seed = 103))
results.thresh <- rbind(results1.thresh, 
                        results2.thresh, 
                        results3.thresh) %>%
  mutate(seed = as.factor(seed))
ggplot(results.thresh, aes(x = thresh, y = rsquared, color = seed)) + 
  geom_line() + 
  theme_bw()

# minspan: minimum number of observations between knots
minspan.range <- c(0, 10, 25, 35, 60, 85, 110)
results1.minspan <- do.call(rbind, lapply(minspan.range, 
                                          FUN = tune.mars.minspan, 
                                          data = trn, nfolds = 10, 
                                          seed = 101))
results2.minspan <- do.call(rbind, lapply(minspan.range, 
                                          FUN = tune.mars.minspan, 
                                          data = trn, nfolds = 10, 
                                          seed = 102))
results3.minspan <- do.call(rbind, lapply(minspan.range, 
                                          FUN = tune.mars.minspan, 
                                          data = trn, nfolds = 10, 
                                          seed = 103))
results.minspan <- rbind(results1.minspan, 
                         results2.minspan, 
                         results3.minspan) %>%
  mutate(seed = as.factor(seed))
ggplot(results.minspan, aes(x = minspan, y = rsquared, color = seed)) + 
  geom_line() + 
  theme_bw()

# fast.k: max number of parent terms
fk.range <- c(0, 10, 20, 30, 40)
results1.fk <- do.call(rbind, lapply(fk.range, FUN = tune.mars.fk, 
                                     data = trn, nfolds = 10, seed = 101))
results2.fk <- do.call(rbind, lapply(fk.range, FUN = tune.mars.fk, 
                                     data = trn, nfolds = 10, seed = 102))
results3.fk <- do.call(rbind, lapply(fk.range, FUN = tune.mars.fk, 
                                     data = trn, nfolds = 10, seed = 103))
results.fk <- rbind(results1.fk, results2.fk, results3.fk) %>%
  mutate(seed = as.factor(seed))
ggplot(results.fk, aes(x = fast.k, y = rsquared, color = seed)) + 
  geom_line() + 
  theme_bw()

#     Tuning --------------------------------------------------------------

# This is a massive undertaking due to the number of tuning params
deg.range <- 1:4
pen.range <- c(0, 2, 4, 6)
nk.range <- c(15, 25, 35, 50)
thresh.range <- c(0, 0.015, 0.5, 0.1)
minspan.range <- c(0, 10, 25, 40, 70)
fk.range <- c(0, 20, 30, 40)
mars.combos <- expand.grid(deg.range, pen.range, nk.range, 
                           thresh.range, minspan.range, fk.range)
d.combos <- mars.combos$Var1
p.combos <- mars.combos$Var2
nk.combos <- mars.combos$Var3
t.combos <- mars.combos$Var4
m.combos <- mars.combos$Var5
fk.combos <- mars.combos$Var6

set.seed(1001)
system.time(
  results.combos <- do.call(
    rbind, 
    mapply(
      FUN = tune.mars, 
      d = d.combos, 
      p = p.combos, 
      nk = nk.combos, 
      t = t.combos, 
      m = m.combos, 
      fk = fk.combos, 
      MoreArgs = 
        list(nfolds = 10, data = trn), 
      SIMPLIFY = F
    )
  )
)

# System time
# user  system elapsed 
# 803.97    1.98  888.59 

results.combos[order(results.combos$rsquared, decreasing = T), ] %>%
  head()
# nfolds degree penalty nk thresh minspan fast.k  rsquared     rmse
#    10      1       0 35  0.000       0     30 0.6784813 2.770432
#    10      1       4 35  0.000       0     30 0.6782463 2.716677
#    10      1       4 50  0.000      10     40 0.6702015 2.686886
results.combos[order(results.combos$rmse), ] %>%
  head()
# nfolds degree penalty nk thresh minspan fast.k  rsquared     rmse
#    10      1       4 15  0.000       0     20 0.6651649 2.593427
#    10      1       6 25  0.000       0     20 0.6493405 2.666962
#    10      1       4 50  0.000      10     40 0.6702015 2.686886
dir.create('tuning/mars/alpha')
saveRDS(results.combos, 'tuning/mars/alpha/tune.RDS')
saveRDS(results.deg, 'tuning/mars/alpha/degree.RDS')
saveRDS(results.pen, 'tuning/mars/alpha/penalty.RDS')
saveRDS(results.nk, 'tuning/mars/alpha/nk.RDS')
saveRDS(results.thresh, 'tuning/mars/alpha/thresh.RDS')
saveRDS(results.minspan, 'tuning/mars/alpha/minspan.RDS')
saveRDS(results.fk, 'tuning/mars/alpha/fast.k.RDS')

# Beta -------------------------------------------------------------------

#     Data organization ---------------------------------------------------

trn.all <- readRDS("./pre-process/beta/1/pp.RDS") 
colnames(trn.all) <- str_replace(colnames(trn.all), "-", ".")
trn.guest <- trn.all$guest
trn <- select(trn.all, -guest)

features <- readRDS("./feature.selection/beta.vars.RDS")
trn <- trn[ , colnames(trn) %in% c("DelG", features)]

#     Estimation ----------------------------------------------------------

# Degree ---
deg.range <- c(1:5)
results1.deg <- do.call(rbind, lapply(deg.range, FUN = tune.mars.deg, 
                                      data = trn, nfolds = 10, seed = 101))
results2.deg <- do.call(rbind, lapply(deg.range, FUN = tune.mars.deg, 
                                      data = trn, nfolds = 10, seed = 102))
results3.deg <- do.call(rbind, lapply(deg.range, FUN = tune.mars.deg, 
                                      data = trn, nfolds = 10, seed = 103))
results.deg <- rbind(results1.deg, results2.deg, results3.deg) %>%
  mutate(seed = as.factor(seed))
ggplot(results.deg, aes(x = degree, y = rsquared, color = seed)) + 
  geom_line() + 
  theme_bw()

# Penalty ---
pen.range <- c(0, 2, 4, 6, 8)
results1.pen <- do.call(rbind, lapply(pen.range, FUN = tune.mars.pen, 
                                      data = trn, nfolds = 10, seed = 101))
results2.pen <- do.call(rbind, lapply(pen.range, FUN = tune.mars.pen, 
                                      data = trn, nfolds = 10, seed = 102))
results3.pen <- do.call(rbind, lapply(pen.range, FUN = tune.mars.pen, 
                                      data = trn, nfolds = 10, seed = 103))
results.pen <- rbind(results1.pen, results2.pen, results3.pen) %>%
  mutate(seed = as.factor(seed))
ggplot(results.pen, aes(x = penalty, y = rsquared, color = seed)) + 
  geom_line() + 
  theme_bw()

# nk: Model terms
nk.range <- c(10, 15, 25, 35, 50, 60)
results1.nk <- do.call(rbind, lapply(nk.range, FUN = tune.mars.nk, 
                                     data = trn, nfolds = 10, seed = 101))
results2.nk <- do.call(rbind, lapply(nk.range, FUN = tune.mars.nk, 
                                     data = trn, nfolds = 10, seed = 102))
results3.nk <- do.call(rbind, lapply(nk.range, FUN = tune.mars.nk, 
                                     data = trn, nfolds = 10, seed = 103))
results.nk <- rbind(results1.nk, results2.nk, results3.nk) %>%
  mutate(seed = as.factor(seed))
ggplot(results.nk, aes(x = nk, y = rsquared, color = seed)) + 
  geom_line() + 
  theme_bw()

# Threshold
thresh.range <- c(0, 0.01, 0.025, 0.05, 0.1, 0.2)
results1.thresh <- do.call(rbind, lapply(thresh.range, 
                                         FUN = tune.mars.thresh, 
                                         data = trn, nfolds = 10, 
                                         seed = 101))
results2.thresh <- do.call(rbind, lapply(thresh.range, 
                                         FUN = tune.mars.thresh, 
                                         data = trn, nfolds = 10, 
                                         seed = 102))
results3.thresh <- do.call(rbind, lapply(thresh.range, 
                                         FUN = tune.mars.thresh,
                                         data = trn, nfolds = 10, 
                                         seed = 103))
results.thresh <- rbind(results1.thresh, 
                        results2.thresh, 
                        results3.thresh) %>%
  mutate(seed = as.factor(seed))
ggplot(results.thresh, aes(x = thresh, y = rsquared, color = seed)) + 
  geom_line() + 
  theme_bw()

# minspan: minimum number of observations between knots
minspan.range <- c(0, 10, 25, 35, 60, 85, 110)
results1.minspan <- do.call(rbind, lapply(minspan.range, 
                                          FUN = tune.mars.minspan, 
                                          data = trn, nfolds = 10, 
                                          seed = 101))
results2.minspan <- do.call(rbind, lapply(minspan.range, 
                                          FUN = tune.mars.minspan, 
                                          data = trn, nfolds = 10, 
                                          seed = 102))
results3.minspan <- do.call(rbind, lapply(minspan.range, 
                                          FUN = tune.mars.minspan, 
                                          data = trn, nfolds = 10, 
                                          seed = 103))
results.minspan <- rbind(results1.minspan, 
                         results2.minspan, 
                         results3.minspan) %>%
  mutate(seed = as.factor(seed))
ggplot(results.minspan, aes(x = minspan, y = rsquared, color = seed)) + 
  geom_line() + 
  theme_bw()

# fast.k: max number of parent terms
fk.range <- c(0, 10, 20, 30, 40)
results1.fk <- do.call(rbind, lapply(fk.range, FUN = tune.mars.fk, 
                                     data = trn, nfolds = 10, seed = 101))
results2.fk <- do.call(rbind, lapply(fk.range, FUN = tune.mars.fk, 
                                     data = trn, nfolds = 10, seed = 102))
results3.fk <- do.call(rbind, lapply(fk.range, FUN = tune.mars.fk, 
                                     data = trn, nfolds = 10, seed = 103))
results.fk <- rbind(results1.fk, results2.fk, results3.fk) %>%
  mutate(seed = as.factor(seed))
ggplot(results.fk, aes(x = fast.k, y = rsquared, color = seed)) + 
  geom_line() + 
  theme_bw()

#     Tuning --------------------------------------------------------------

# This is a massive undertaking due to the number of tuning params
deg.range <- 1:4
pen.range <- c(0, 2, 4, 6)
nk.range <- c(15, 25, 35, 50)
thresh.range <- c(0, 0.015, 0.5, 0.1)
minspan.range <- c(0, 10, 25, 40, 70)
fk.range <- c(10, 20, 30, 40)
mars.combos <- expand.grid(deg.range, pen.range, nk.range, 
                           thresh.range, minspan.range, fk.range)
d.combos <- mars.combos$Var1
p.combos <- mars.combos$Var2
nk.combos <- mars.combos$Var3
t.combos <- mars.combos$Var4
m.combos <- mars.combos$Var5
fk.combos <- mars.combos$Var6

set.seed(1001)
system.time(
  results.combos <- do.call(
    rbind, 
    mapply(
      FUN = tune.mars, 
      d = d.combos, 
      p = p.combos, 
      nk = nk.combos, 
      t = t.combos, 
      m = m.combos, 
      fk = fk.combos, 
      MoreArgs = 
        list(nfolds = 10, data = trn), 
      SIMPLIFY = F
    )
  )
)

# System time
# user  system elapsed 
# 803.97    1.98  888.59 

results.combos[order(results.combos$rsquared, decreasing = T), ] %>%
  head()
# nfolds degree penalty nk thresh minspan fast.k  rsquared     rmse
#    10      2       4 25  0.015      70     20 0.7313526 2.470377
#    10      1       6 15  0.000      25     30 0.7238435 2.423611
#    10      2       4 35  0.015      40     20 0.7105351 2.420387
results.combos[order(results.combos$rmse), ] %>%
  head()
# nfolds degree penalty nk thresh minspan fast.k  rsquared     rmse
#    10      2       4 35  0.015      40     20 0.7105351 2.420387
#    10      1       6 15  0.000      25     30 0.7238435 2.423611
#    10      2       4 25  0.015      70     20 0.7313526 2.470377
dir.create('tuning/mars/beta')
saveRDS(results.combos, 'tuning/mars/beta/tune.RDS')
saveRDS(results.deg, 'tuning/mars/beta/degree.RDS')
saveRDS(results.pen, 'tuning/mars/beta/penalty.RDS')
saveRDS(results.nk, 'tuning/mars/beta/nk.RDS')
saveRDS(results.thresh, 'tuning/mars/beta/thresh.RDS')
saveRDS(results.minspan, 'tuning/mars/beta/minspan.RDS')
saveRDS(results.fk, 'tuning/mars/beta/fast.k.RDS')

# Gamma -------------------------------------------------------------------
#     Data organization ---------------------------------------------------

trn.all <- readRDS("./pre-process/gamma/1/pp.RDS") 
colnames(trn.all) <- str_replace(colnames(trn.all), "-", ".")
trn.guest <- trn.all$guest
trn <- select(trn.all, -guest)

features <- readRDS("./feature.selection/gamma.vars.RDS")
trn <- trn[ , colnames(trn) %in% c("DelG", features)]

#     Estimation ----------------------------------------------------------

# Degree ---
deg.range <- c(1:5)
results1.deg <- do.call(rbind, lapply(deg.range, FUN = tune.mars.deg, 
                                      data = trn, nfolds = 10, seed = 101))
results2.deg <- do.call(rbind, lapply(deg.range, FUN = tune.mars.deg, 
                                      data = trn, nfolds = 10, seed = 102))
results3.deg <- do.call(rbind, lapply(deg.range, FUN = tune.mars.deg, 
                                      data = trn, nfolds = 10, seed = 103))
results.deg <- rbind(results1.deg, results2.deg, results3.deg) %>%
  mutate(seed = as.factor(seed))
ggplot(results.deg, aes(x = degree, y = rsquared, color = seed)) + 
  geom_line() + 
  theme_bw()

# Penalty ---
pen.range <- c(0:4)
results1.pen <- do.call(rbind, lapply(pen.range, FUN = tune.mars.pen, 
                                      data = trn, nfolds = 10, seed = 101))
results2.pen <- do.call(rbind, lapply(pen.range, FUN = tune.mars.pen, 
                                      data = trn, nfolds = 10, seed = 102))
results3.pen <- do.call(rbind, lapply(pen.range, FUN = tune.mars.pen, 
                                      data = trn, nfolds = 10, seed = 103))
results.pen <- rbind(results1.pen, results2.pen, results3.pen) %>%
  mutate(seed = as.factor(seed))
ggplot(results.pen, aes(x = penalty, y = rsquared, color = seed)) + 
  geom_line() + 
  theme_bw()

# nk: Model terms
nk.range <- c(10, 15, 25, 35, 50, 60)
results1.nk <- do.call(rbind, lapply(nk.range, FUN = tune.mars.nk, 
                                     data = trn, nfolds = 10, seed = 101))
results2.nk <- do.call(rbind, lapply(nk.range, FUN = tune.mars.nk, 
                                     data = trn, nfolds = 10, seed = 102))
results3.nk <- do.call(rbind, lapply(nk.range, FUN = tune.mars.nk, 
                                     data = trn, nfolds = 10, seed = 103))
results.nk <- rbind(results1.nk, results2.nk, results3.nk) %>%
  mutate(seed = as.factor(seed))
ggplot(results.nk, aes(x = nk, y = rsquared, color = seed)) + 
  geom_line() + 
  theme_bw()

# Threshold
thresh.range <- c(0, 0.01, 0.025, 0.05, 0.1)
results1.thresh <- do.call(rbind, lapply(thresh.range, 
                                         FUN = tune.mars.thresh, 
                                         data = trn, nfolds = 10, 
                                         seed = 101))
results2.thresh <- do.call(rbind, lapply(thresh.range, 
                                         FUN = tune.mars.thresh, 
                                         data = trn, nfolds = 10, 
                                         seed = 102))
results3.thresh <- do.call(rbind, lapply(thresh.range, 
                                         FUN = tune.mars.thresh,
                                         data = trn, nfolds = 10, 
                                         seed = 103))
results.thresh <- rbind(results1.thresh, 
                        results2.thresh, 
                        results3.thresh) %>%
  mutate(seed = as.factor(seed))
ggplot(results.thresh, aes(x = thresh, y = rsquared, color = seed)) + 
  geom_line() + 
  theme_bw()

# minspan: minimum number of observations between knots
minspan.range <- c(0, 10, 25, 35, 60)
results1.minspan <- do.call(rbind, lapply(minspan.range, 
                                          FUN = tune.mars.minspan, 
                                          data = trn, nfolds = 10, 
                                          seed = 101))
results2.minspan <- do.call(rbind, lapply(minspan.range, 
                                          FUN = tune.mars.minspan, 
                                          data = trn, nfolds = 10, 
                                          seed = 102))
results3.minspan <- do.call(rbind, lapply(minspan.range, 
                                          FUN = tune.mars.minspan, 
                                          data = trn, nfolds = 10, 
                                          seed = 103))
results.minspan <- rbind(results1.minspan, 
                         results2.minspan, 
                         results3.minspan) %>%
  mutate(seed = as.factor(seed))
ggplot(results.minspan, aes(x = minspan, y = rsquared, color = seed)) + 
  geom_line() + 
  theme_bw()

# fast.k: max number of parent terms
fk.range <- c(0, 10, 20, 30, 40)
results1.fk <- do.call(rbind, lapply(fk.range, FUN = tune.mars.fk, 
                                     data = trn, nfolds = 10, seed = 101))
results2.fk <- do.call(rbind, lapply(fk.range, FUN = tune.mars.fk, 
                                     data = trn, nfolds = 10, seed = 102))
results3.fk <- do.call(rbind, lapply(fk.range, FUN = tune.mars.fk, 
                                     data = trn, nfolds = 10, seed = 103))
results.fk <- rbind(results1.fk, results2.fk, results3.fk) %>%
  mutate(seed = as.factor(seed))
ggplot(results.fk, aes(x = fast.k, y = rsquared, color = seed)) + 
  geom_line() + 
  theme_bw()

#     Tuning --------------------------------------------------------------

# This is a massive undertaking due to the number of tuning params
deg.range <- 1:4
pen.range <- 0:4
nk.range <- c(1, 15, 25, 40)
thresh.range <- c(0, 0.015, 0.5, 0.1)
minspan.range <- c(0, 10, 25, 45)
fk.range <- c(0, 10, 20, 30)
mars.combos <- expand.grid(deg.range, pen.range, nk.range, 
                           thresh.range, minspan.range, fk.range)
d.combos <- mars.combos$Var1
p.combos <- mars.combos$Var2
nk.combos <- mars.combos$Var3
t.combos <- mars.combos$Var4
m.combos <- mars.combos$Var5
fk.combos <- mars.combos$Var6

set.seed(1001)
system.time(
  results.combos <- do.call(
    rbind, 
    mapply(
      FUN = tune.mars, 
      d = d.combos, 
      p = p.combos, 
      nk = nk.combos, 
      t = t.combos, 
      m = m.combos, 
      fk = fk.combos, 
      MoreArgs = 
        list(nfolds = 10, data = trn), 
      SIMPLIFY = F
    )
  )
)

# System time
# user  system elapsed 
# 803.97    1.98  888.59 

results.combos[order(results.combos$rsquared, decreasing = T), ] %>%
  head()
# nfolds degree penalty nk thresh minspan fast.k  rsquared     rmse
#    10      2       4 25  0.015      70     20 0.7313526 2.470377
#    10      1       6 15  0.000      25     30 0.7238435 2.423611
#    10      2       4 35  0.015      40     20 0.7105351 2.420387
results.combos[order(results.combos$rmse), ] %>%
  head()
# nfolds degree penalty nk thresh minspan fast.k  rsquared     rmse
#    10      2       4 35  0.015      40     20 0.7105351 2.420387
#    10      1       6 15  0.000      25     30 0.7238435 2.423611
#    10      2       4 25  0.015      70     20 0.7313526 2.470377
dir.create('tuning/mars/gamma')
saveRDS(results.combos, 'tuning/mars/gamma/tune.RDS')
saveRDS(results.deg, 'tuning/mars/gamma/degree.RDS')
saveRDS(results.pen, 'tuning/mars/gamma/penalty.RDS')
saveRDS(results.nk, 'tuning/mars/gamma/nk.RDS')
saveRDS(results.thresh, 'tuning/mars/gamma/thresh.RDS')
saveRDS(results.minspan, 'tuning/mars/gamma/minspan.RDS')
saveRDS(results.fk, 'tuning/mars/gamma/fast.k.RDS')
