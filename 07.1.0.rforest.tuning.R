# Libraries and Packages --------------------------------------------------

library(caret)
library(randomForest)
library(tidyverse)

# Tuning Functions -----------------------------------------------------

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
  message(nfolds, "-fold cross-validation of ntree = ", ntree,
          ", nodesize = ", node,
          ", mtry = ", m, 
          " completed.")
  return(data.frame( 
    nfolds = nfolds, 
    ntree = ntree, 
    nodesize = node,
    mtry = m,
    rsquared = sum(r2.results)/nfolds, 
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

dir.create("./tuning/rforest")

#     Estimation ----------------------------------------------------------

# Number of trees ---
ntree.range <- c(25, 50, 75, 100, 250, 500)
results1.ntree <- do.call(rbind, lapply(ntree.range, 
                                        FUN = tune.rf.ntree, 
                                        data = trn, 
                                        nfolds =10, seed = 101))
results2.ntree <- do.call(rbind, lapply(ntree.range, 
                                        FUN = tune.rf.ntree, 
                                        data = trn, 
                                        nfolds =10, seed = 102))
results3.ntree <- do.call(rbind, lapply(ntree.range, 
                                        FUN = tune.rf.ntree, 
                                        data = trn, 
                                        nfolds =10, seed = 103))
results.ntree <- rbind(results1.ntree, results2.ntree, results3.ntree) %>%
  mutate(seed = as.factor(seed))
ggplot(results.ntree, aes(x = ntree, y = rsquared, color = seed)) + 
  geom_line() + 
  theme_bw()

#     Node size ---
node.range <- c(1:5, 10, 25)
results1.node <- do.call(rbind, lapply(node.range, 
                                       FUN = tune.rf.node, 
                                       data = trn, 
                                       nfolds =10, seed = 101))
results2.node <- do.call(rbind, lapply(node.range, 
                                        FUN = tune.rf.node, 
                                        data = trn, 
                                        nfolds =10, seed = 102))
results3.node <- do.call(rbind, lapply(node.range, 
                                        FUN = tune.rf.node, 
                                        data = trn, 
                                        nfolds =10, seed = 103))
results.node <- rbind(results1.node, results2.node, results3.node) %>%
  mutate(seed = as.factor(seed))
ggplot(results.node, aes(x = nodesize, y = rsquared, color = seed)) + 
  geom_line() + 
  theme_bw()

#     Number of variables (mtry) ---
mtry.range <- c(1:5, 8, 15)
results1.mtry <- do.call(rbind, lapply(mtry.range, 
                                       FUN = tune.rf.mtry, 
                                       data = trn, 
                                       nfolds =10, seed = 101))
results2.mtry <- do.call(rbind, lapply(mtry.range, 
                                       FUN = tune.rf.mtry, 
                                       data = trn, 
                                       nfolds =10, seed = 102))
results3.mtry <- do.call(rbind, lapply(mtry.range, 
                                       FUN = tune.rf.mtry, 
                                       data = trn, 
                                       nfolds =10, seed = 103))
results.mtry <- rbind(results1.mtry, results2.mtry, results3.mtry) %>%
  mutate(seed = as.factor(seed))
ggplot(results.mtry, aes(x = mtry, y = rsquared, color = seed)) + 
  geom_line() + 
  theme_bw()

#     Tuning --------------------------------------------------------------

# 6*6*7 = 252 combinations
rf.combos <- expand.grid(ntree.range, node.range, mtry.range)
colnames(rf.combos) <- c("ntree", "node", "mtry")
ntree.combos <- rf.combos$ntree
node.combos <- rf.combos$node
mtry.combos <- rf.combos$mtry

set.seed(1001)
# I reduce the number of folds because rforest takes a while
system.time(
  results.combos <- do.call(
    rbind,
    mapply(
      FUN = tune.rf,
      ntree = ntree.combos, 
      node = node.combos, 
      m = mtry.combos, 
      MoreArgs = 
        list(nfolds = 5, data = trn), 
      SIMPLIFY = F
    )
  )
)
#  user  system elapsed 
# 122.26    0.79  129.89 
results.combos[order(results.combos$rsquared, decreasing = T), ] %>% head()
# nfolds ntree nodesize mtry  rsquared     rmse
#      5    75        3   15 0.6963144 2.606530
#      5   250        1    4 0.6946227 2.624331
results.combos[order(results.combos$rmse), ] %>% head()
# nfolds ntree nodesize mtry  rsquared     rmse
#      5    50        3   15 0.6905937 2.534911
#      5   100        1    4 0.6808741 2.537100


dir.create("./tuning/rforest/alpha")
saveRDS(results.combos, "./tuning/rforest/alpha/tune.RDS")
saveRDS(results.ntree, "./tuning/rforest/alpha/ntree.RDS")
saveRDS(results.node, "./tuning/rforest/alpha/nodesize.RDS")
saveRDS(results.mtry, "./tuning/rforest/alpha/mtry.RDS")


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
ntree.range <- c(25, 50, 75, 100, 250, 500)
results1.ntree <- do.call(rbind, lapply(ntree.range, 
                                        FUN = tune.rf.ntree, 
                                        data = trn, 
                                        nfolds =10, seed = 101)) 
results2.ntree <- do.call(rbind, lapply(ntree.range, 
                                        FUN = tune.rf.ntree, 
                                        data = trn, 
                                        nfolds =10, seed = 102))
results3.ntree <- do.call(rbind, lapply(ntree.range, 
                                        FUN = tune.rf.ntree, 
                                        data = trn, 
                                        nfolds =10, seed = 103))
results.ntree <- rbind(results1.ntree, results2.ntree, results3.ntree) %>%
  mutate(seed = as.factor(seed))
ggplot(results.ntree, aes(x = ntree, y = rsquared, color = seed)) + 
  geom_line() + 
  theme_bw()

#     Node size ---
node.range <- c(1, 2, 5, 10, 25, 50)
results1.node <- do.call(rbind, lapply(node.range, 
                                       FUN = tune.rf.node, 
                                       data = trn, 
                                       nfolds =10, seed = 101)) 
results2.node <- do.call(rbind, lapply(node.range, 
                                       FUN = tune.rf.node, 
                                       data = trn, 
                                       nfolds =10, seed = 102))
results3.node <- do.call(rbind, lapply(node.range, 
                                       FUN = tune.rf.node, 
                                       data = trn, 
                                       nfolds =10, seed = 103))
results.node <- rbind(results1.node, results2.node, results3.node) %>%
  mutate(seed = as.factor(seed))
ggplot(results.node, aes(x = nodesize, y = rsquared, color = seed)) + 
  geom_line() + 
  theme_bw() 

#     Number of variables (mtry) ---
mtry.range <- c(1:5, 8, 15)
results1.mtry <- do.call(rbind, lapply(mtry.range, 
                                       FUN = tune.rf.mtry, 
                                       data = trn,  
                                       nfolds =10, seed = 101)) 
results2.mtry <- do.call(rbind, lapply(mtry.range, 
                                       FUN = tune.rf.mtry, 
                                       data = trn, 
                                       nfolds =10, seed = 102))
results3.mtry <- do.call(rbind, lapply(mtry.range, 
                                       FUN = tune.rf.mtry, 
                                       data = trn, 
                                       nfolds =10, seed = 103))
results.mtry <- rbind(results1.mtry, results2.mtry, results3.mtry) %>%
  mutate(seed = as.factor(seed))
ggplot(results.mtry, aes(x = mtry, y = rsquared, color = seed)) + 
  geom_line() + 
  theme_bw()

#     Tuning --------------------------------------------------------------

# 5*6*7 = 210 combinations
rf.combos <- expand.grid(ntree.range, node.range, mtry.range)
colnames(rf.combos) <- c("ntree", "node", "mtry")
ntree.combos <- rf.combos$ntree
node.combos <- rf.combos$node
mtry.combos <- rf.combos$mtry

set.seed(1001)
system.time(
  results.combos <- do.call(
    rbind,
    mapply(
      FUN = tune.rf,
      ntree = ntree.combos, 
      node = node.combos, 
      m = mtry.combos, 
      MoreArgs = 
        list(nfolds = 10, data = trn), 
      SIMPLIFY = F
    )
  )
)
#  user  system elapsed 
# 452.75    1.22  466.45 
results.combos[order(results.combos$rsquared, decreasing = T), ] %>% head()
# nfolds ntree nodesize mtry  rsquared     rmse
#     10    50        5    4 0.7761236 2.225618
#     10    75        1    1 0.7672687 2.280025
results.combos[order(results.combos$rmse), ] %>% head()
# nfolds ntree nodesize mtry  rsquared     rmse
#    10    50        5    4 0.7761236 2.225618
#    10    75        1    5 0.7647129 2.234098

dir.create("./tuning/rforest/beta")
saveRDS(results.combos, "./tuning/rforest/beta/tune.RDS")
saveRDS(results.ntree, "./tuning/rforest/beta/ntree.RDS")
saveRDS(results.node, "./tuning/rforest/beta/nodesize.RDS")
saveRDS(results.mtry, "./tuning/rforest/beta/mtry.RDS")

# Gamma -------------------------------------------------------------------

#     Data Organization ---------------------------------------------------

trn.all <- readRDS("./pre-process/gamma/1/pp.RDS") 
colnames(trn.all) <- str_replace(colnames(trn.all), "-", ".")
trn.guest <- trn.all$guest
trn <- select(trn.all, -guest)

features <- readRDS("./feature.selection/gamma.vars.RDS")
trn <- trn[ , colnames(trn) %in% c("DelG", features)]

#     Estimation ----------------------------------------------------------

# Number of trees ---
ntree.range <- c(50, 75, 100, 250, 400, 700)
results1.ntree <- do.call(rbind, lapply(ntree.range, 
                                        FUN = tune.rf.ntree, 
                                        data = trn, 
                                        nfolds =10, seed = 101)) 
results2.ntree <- do.call(rbind, lapply(ntree.range, 
                                        FUN = tune.rf.ntree, 
                                        data = trn, 
                                        nfolds =10, seed = 102))
results3.ntree <- do.call(rbind, lapply(ntree.range, 
                                        FUN = tune.rf.ntree, 
                                        data = trn, 
                                        nfolds =10, seed = 103))
results.ntree <- rbind(results1.ntree, results2.ntree, results3.ntree) %>%
  mutate(seed = as.factor(seed))
ggplot(results.ntree, aes(x = ntree, y = rsquared, color = seed)) + 
  geom_line() + 
  theme_bw()

#     Node size ---
node.range <- c(1, 2, 5, 10, 25, 50)
results1.node <- do.call(rbind, lapply(node.range, 
                                       FUN = tune.rf.node, 
                                       data = trn, 
                                       nfolds =10, seed = 101)) 
results2.node <- do.call(rbind, lapply(node.range, 
                                       FUN = tune.rf.node, 
                                       data = trn, 
                                       nfolds =10, seed = 102))
results3.node <- do.call(rbind, lapply(node.range, 
                                       FUN = tune.rf.node, 
                                       data = trn, 
                                       nfolds =10, seed = 103))
results.node <- rbind(results1.node, results2.node, results3.node) %>%
  mutate(seed = as.factor(seed))
ggplot(results.node, aes(x = nodesize, y = rsquared, color = seed)) + 
  geom_line() + 
  theme_bw() 

#     Number of variables (mtry) ---
mtry.range <- c(1, 2, 4, 8, 15, 20, 25)
results1.mtry <- do.call(rbind, lapply(mtry.range, 
                                       FUN = tune.rf.mtry, 
                                       data = trn,  
                                       nfolds =10, seed = 101)) 
results2.mtry <- do.call(rbind, lapply(mtry.range, 
                                       FUN = tune.rf.mtry, 
                                       data = trn, 
                                       nfolds =10, seed = 102))
results3.mtry <- do.call(rbind, lapply(mtry.range, 
                                       FUN = tune.rf.mtry, 
                                       data = trn, 
                                       nfolds =10, seed = 103))
results.mtry <- rbind(results1.mtry, results2.mtry, results3.mtry) %>%
  mutate(seed = as.factor(seed))
ggplot(results.mtry, aes(x = mtry, y = rsquared, color = seed)) + 
  geom_line() + 
  theme_bw()

#     Tuning --------------------------------------------------------------

# 5*6*7 = 210 combinations
rf.combos <- expand.grid(ntree.range, node.range, mtry.range)
colnames(rf.combos) <- c("ntree", "node", "mtry")
ntree.combos <- rf.combos$ntree
node.combos <- rf.combos$node
mtry.combos <- rf.combos$mtry

set.seed(1001)
system.time(
  results.combos <- do.call(
    rbind,
    mapply(
      FUN = tune.rf,
      ntree = ntree.combos, 
      node = node.combos, 
      m = mtry.combos, 
      MoreArgs = 
        list(nfolds = 10, data = trn), 
      SIMPLIFY = F
    )
  )
)
# user  system elapsed 
# 200.83    0.30  201.37 
results.combos[order(results.combos$rsquared, decreasing = T), ] %>% head()
# nfolds ntree nodesize mtry  rsquared     rmse
# 10   700        1   15 0.3902417 1.656713
# 10   250        1    1 0.3743013 1.582166
results.combos[order(results.combos$rmse), ] %>% head()
# nfolds ntree nodesize mtry  rsquared     rmse
# 10   400        5   15 0.3238978 1.536768
# 10    75       10    4 0.3225326 1.540947

dir.create("./tuning/rforest/gamma")
saveRDS(results.combos, "./tuning/rforest/gamma/tune.RDS")
saveRDS(results.ntree, "./tuning/rforest/gamma/ntree.RDS")
saveRDS(results.node, "./tuning/rforest/gamma/nodesize.RDS")
saveRDS(results.mtry, "./tuning/rforest/gamma/mtry.RDS")
