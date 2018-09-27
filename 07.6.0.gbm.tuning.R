source('./07.model.functions.R')
dir.create('./tuning/gbm')
p_load(gbm)

# Functions ---------------------------------------------------------------

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

# Alpha -------------------------------------------------------------------

#     Data Organization ---------------------------------------------------

trn.all <- readRDS("./pre-process/alpha/1/pp.RDS") 
colnames(trn.all) <- str_replace(colnames(trn.all), "-", ".")
trn.guest <- trn.all$guest
trn <- select(trn.all, -guest)

features <- readRDS("./feature.selection/alpha.vars.RDS")
trn <- trn[ , colnames(trn) %in% c("DelG", features)]

#     Estimation ----------------------------------------------------------

# Number of trees ---
# For some reason, ntree gets up to frankly ridiculous levels
# Probably should cap it at 5000 for sanity's sake, though it can go 
# higher (it increases through 10000 trees)
ntree.range <- c(250, 500, 750, 1250, 2000, 3000, 5000)
results1.ntree <- do.call(rbind, lapply(ntree.range, 
                                        FUN = tune.gbm.ntree, 
                                        data = trn, 
                                        nfolds =10, seed = 101))
results2.ntree <- do.call(rbind, lapply(ntree.range, 
                                        FUN = tune.gbm.ntree, 
                                        data = trn, 
                                        nfolds =10, seed = 102))
results3.ntree <- do.call(rbind, lapply(ntree.range, 
                                        FUN = tune.gbm.ntree, 
                                        data = trn, 
                                        nfolds =10, seed = 103))
results.ntree <- rbind(results1.ntree, results2.ntree, results3.ntree) %>%
  mutate(seed = as.factor(seed))
ggplot(results.ntree, aes(x = ntree, y = rsquared, color = seed)) + 
  geom_line() + 
  theme_bw()

# Interaction depth ---
depth.range <- 1:5
results1.depth <- do.call(rbind, lapply(depth.range, 
                                        FUN = tune.gbm.depth, 
                                        data = trn, 
                                        nfolds =10, seed = 101))
results2.depth <- do.call(rbind, lapply(depth.range, 
                                        FUN = tune.gbm.depth, 
                                        data = trn, 
                                        nfolds =10, seed = 102))
results3.depth <- do.call(rbind, lapply(depth.range, 
                                        FUN = tune.gbm.depth, 
                                        data = trn, 
                                        nfolds =10, seed = 103))
results.depth <- rbind(results1.depth, results2.depth, results3.depth) %>%
  mutate(seed = as.factor(seed))
ggplot(results.depth, aes(x = depth, y = rsquared, color = seed)) + 
  geom_line() + 
  theme_bw()

# Shrinkage (learning rate)---
shrink.range <- c(0.01, 0.05, 0.1, 0.15, 0.25, 0.5)
results1.shrink <- do.call(rbind, lapply(shrink.range, 
                                         FUN = tune.gbm.shrink, 
                                         data = trn, 
                                         nfolds =10, seed = 101))
results2.shrink <- do.call(rbind, lapply(shrink.range, 
                                         FUN = tune.gbm.shrink, 
                                         data = trn, 
                                         nfolds =10, seed = 102))
results3.shrink <- do.call(rbind, lapply(shrink.range, 
                                         FUN = tune.gbm.shrink, 
                                         data = trn, 
                                         nfolds =10, seed = 103))
results.shrink <- rbind(results1.shrink, results2.shrink, results3.shrink) %>%
  mutate(seed = as.factor(seed))
ggplot(results.shrink, aes(x = shrinkage, y = rsquared, color = seed)) + 
  geom_line() + 
  theme_bw()

# Minimum observations per node
node.range <- c(1, 5, 10, 15, 20, 25)
results1.node <- do.call(rbind, lapply(node.range, 
                                         FUN = tune.gbm.node, 
                                         data = trn, 
                                         nfolds =10, seed = 101))
results2.node <- do.call(rbind, lapply(node.range, 
                                         FUN = tune.gbm.node, 
                                         data = trn, 
                                         nfolds =10, seed = 102))
results3.node <- do.call(rbind, lapply(node.range, 
                                         FUN = tune.gbm.node, 
                                         data = trn, 
                                         nfolds =10, seed = 103))
results.node <- rbind(results1.node, results2.node, results3.node) %>%
  mutate(seed = as.factor(seed))
ggplot(results.node, aes(x = node, y = rsquared, color = seed)) + 
  geom_line() + 
  theme_bw()

#     Tuning --------------------------------------------------------------

gbm.combos <- expand.grid(ntree.range, 
                          depth.range, 
                          shrink.range, 
                          node.range)
ntree.combos <- gbm.combos$Var1
depth.combos <- gbm.combos$Var2
shrink.combos <- gbm.combos$Var3
node.combos <- gbm.combos$Var4

set.seed(1001)
# I reduce the number of folds because gbmorest takes a while
system.time(
  alpha.tune <- do.call(
    rbind,
    mapply(
      FUN = tune.gbm,
      num = ntree.combos, 
      d = depth.combos, 
      s = shrink.combos,
      n = node.combos, 
      MoreArgs = 
        list(nfolds = 5, data = trn), 
      SIMPLIFY = F
    )
  )
)
# user  system elapsed 
# 957.60    4.06 1029.32 
alpha.tune$data <- 'trn'
alpha.tune[order(alpha.tune$rsquared, decreasing = T), ] %>% head()
# data nfolds ntree depth shrinkage node  rsquared    rmse
#  trn      5  2000     4      0.01    1 0.7509958 2.23541
#  trn      5  2000     4      0.01    1 0.7509958 2.23541
alpha.tune[order(alpha.tune$rmse), ] %>% head()
# data nfolds ntree depth shrinkage node  rsquared    rmse
# trn      5  2000     4      0.01    1 0.7509958 2.23541
# trn      5  2000     4      0.01    1 0.7509958 2.23541

dir.create('tuning/gbm/alpha')
saveRDS(alpha.tune, 'tuning/gbm/alpha/tune.RDS')
saveRDS(results.ntree, 'tuning/gbm/alpha/ntree.RDS')
saveRDS(results.depth, 'tuning/gbm/alpha/depth.RDS')
saveRDS(results.shrink, 'tuning/gbm/alpha/shrink.RDS')
saveRDS(results.node, 'tuning/gbm/alpha/node.RDS')

# Beta -------------------------------------------------------------------

#     Data Organization ---------------------------------------------------

trn.all <- readRDS("./pre-process/beta/1/pp.RDS") 
colnames(trn.all) <- str_replace(colnames(trn.all), "-", ".")
trn.guest <- trn.all$guest
trn <- select(trn.all, -guest)

features <- readRDS("./feature.selection/beta.vars.RDS")
trn <- trn[ , colnames(trn) %in% c("DelG", features)]

#     Estimation ----------------------------------------------------------

# Number of trees ---
# For some reason, ntree gets up to frankly ridiculous levels
# Probably should cap it at 5000 for sanity's sake, though it can go 
# higher (it increases through 10000 trees)
ntree.range <- c(250, 500, 750, 1250, 2000, 3000, 5000)
results1.ntree <- do.call(rbind, lapply(ntree.range, 
                                        FUN = tune.gbm.ntree, 
                                        data = trn, 
                                        nfolds =10, seed = 101)) %>% print()
results2.ntree <- do.call(rbind, lapply(ntree.range, 
                                        FUN = tune.gbm.ntree, 
                                        data = trn, 
                                        nfolds =10, seed = 102))
results3.ntree <- do.call(rbind, lapply(ntree.range, 
                                        FUN = tune.gbm.ntree, 
                                        data = trn, 
                                        nfolds =10, seed = 103))
results.ntree <- rbind(results1.ntree, results2.ntree, results3.ntree) %>%
  mutate(seed = as.factor(seed))
ggplot(results.ntree, aes(x = ntree, y = rsquared, color = seed)) + 
  geom_line() + 
  theme_bw()

# Interaction depth ---
depth.range <- 2:6
results1.depth <- do.call(rbind, lapply(depth.range, 
                                        FUN = tune.gbm.depth, 
                                        data = trn, 
                                        nfolds =10, seed = 101))
results2.depth <- do.call(rbind, lapply(depth.range, 
                                        FUN = tune.gbm.depth, 
                                        data = trn, 
                                        nfolds =10, seed = 102))
results3.depth <- do.call(rbind, lapply(depth.range, 
                                        FUN = tune.gbm.depth, 
                                        data = trn, 
                                        nfolds =10, seed = 103))
results.depth <- rbind(results1.depth, results2.depth, results3.depth) %>%
  mutate(seed = as.factor(seed))
ggplot(results.depth, aes(x = depth, y = rsquared, color = seed)) + 
  geom_line() + 
  theme_bw()

# Shrinkage (learning rate)---
shrink.range <- c(0.01, 0.05, 0.1, 0.15, 0.25)
results1.shrink <- do.call(rbind, lapply(shrink.range, 
                                         FUN = tune.gbm.shrink, 
                                         data = trn, 
                                         nfolds =10, seed = 101))
results2.shrink <- do.call(rbind, lapply(shrink.range, 
                                         FUN = tune.gbm.shrink, 
                                         data = trn, 
                                         nfolds =10, seed = 102))
results3.shrink <- do.call(rbind, lapply(shrink.range, 
                                         FUN = tune.gbm.shrink, 
                                         data = trn, 
                                         nfolds =10, seed = 103))
results.shrink <- rbind(results1.shrink, results2.shrink, results3.shrink) %>%
  mutate(seed = as.factor(seed))
ggplot(results.shrink, aes(x = shrinkage, y = rsquared, color = seed)) + 
  geom_line() + 
  theme_bw()

# Minimum observations per node
node.range <- c(1, 5, 10, 15, 20, 25)
results1.node <- do.call(rbind, lapply(node.range, 
                                       FUN = tune.gbm.node, 
                                       data = trn, 
                                       nfolds =10, seed = 101))
results2.node <- do.call(rbind, lapply(node.range, 
                                       FUN = tune.gbm.node, 
                                       data = trn, 
                                       nfolds =10, seed = 102))
results3.node <- do.call(rbind, lapply(node.range, 
                                       FUN = tune.gbm.node, 
                                       data = trn, 
                                       nfolds =10, seed = 103))
results.node <- rbind(results1.node, results2.node, results3.node) %>%
  mutate(seed = as.factor(seed))
ggplot(results.node, aes(x = node, y = rsquared, color = seed)) + 
  geom_line() + 
  theme_bw()

#     Tuning --------------------------------------------------------------

gbm.combos <- expand.grid(ntree.range, 
                          depth.range, 
                          shrink.range, 
                          node.range)
ntree.combos <- gbm.combos$Var1
depth.combos <- gbm.combos$Var2
shrink.combos <- gbm.combos$Var3
node.combos <- gbm.combos$Var4

set.seed(1001)
# I reduce the number of folds because gbmorest takes a while
system.time(
  beta.tune <- do.call(
    rbind,
    mapply(
      FUN = tune.gbm,
      num = ntree.combos, 
      d = depth.combos, 
      s = shrink.combos,
      n = node.combos, 
      MoreArgs = 
        list(nfolds = 5, data = trn), 
      SIMPLIFY = F
    )
  )
)
# user  system elapsed 
# 1449.98    1.99 1493.29 
beta.tune[order(beta.tune$rsquared, decreasing = T), ] %>% head()
# nfolds ntree depth shrinkage node  rsquared     rmse
#    5   250     2      0.10    1 0.7569815 2.207881
#    5   500     5      0.01    1 0.7530211 2.262719
beta.tune[order(beta.tune$rmse), ] %>% head()
# nfolds ntree depth shrinkage node  rsquared     rmse
#      5   250     2      0.10    1 0.7569815 2.207881
#      5  2000     5      0.01    1 0.7502068 2.224097

dir.create('tuning/gbm/beta')
saveRDS(beta.tune, 'tuning/gbm/beta/tune.RDS')
saveRDS(results.ntree, 'tuning/gbm/beta/ntree.RDS')
saveRDS(results.depth, 'tuning/gbm/beta/depth.RDS')
saveRDS(results.shrink, 'tuning/gbm/beta/shrink.RDS')
saveRDS(results.node, 'tuning/gbm/beta/node.RDS')