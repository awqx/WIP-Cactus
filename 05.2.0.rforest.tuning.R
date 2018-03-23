# Libraries and Packages --------------------------------------------------

library(caret)
library(randomForest)
library(tidyverse)

# Data Organization -------------------------------------------------------

rpt <- readRDS("./data/padel.pp.RDS")
mat.dg <- rpt %>% dplyr::select(., -guest:-host, -data.source)

set.seed(1)
trn.ind <- sample(x = 1:nrow(mat.dg), size = round(0.8 * nrow(mat.dg)))
trn <- mat.dg[trn.ind, ]
tst <- mat.dg[-trn.ind, ]

# dir.create("./tuning")
dir.create("./tuning/rforest")

# Tuning Functions ---------------------------------------------------------

tune.rf.ntree <- function(data, nfolds, num, seed) {
  set.seed(seed)
  fold.list <- createFolds(y = data[ , 1], k = nfolds)
  results <- c(rep(0.0, nfolds))
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
    R2 <- defaultSummary(rf.df)[2] 
    results[i] <- R2
    message(paste0("Fold ", i, ", ", num, " trees, has been completed."))
  }
  return(data.frame( # Useful for records
    data = deparse(substitute(data)), # Turns the variable name into char
    nfolds = nfolds,
    seed = seed,
    ntree = num,
    rsquared = sum(results) / nfolds))
}

tune.rf.node <- function(data, nfolds, node, seed) {
  set.seed(seed)
  fold.list <- createFolds(y = data[ , 1], k = nfolds)
  results <- c(rep(0.0, nfolds))
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
    R2 <- defaultSummary(rf.df)[2] 
    results[i] <- R2
    message(paste0("Fold ", i, ", nodesize of ", node, ", has been completed."))
  }
  return(data.frame( # Useful for records
    data = deparse(substitute(data)), # Turns the variable name into char
    nfolds = nfolds, 
    seed = seed,
    nodesize = node,
    rsquared = sum(results) / nfolds))
}

tune.rf.mtry <- function(data, nfolds, m, seed) {
  set.seed(seed)
  fold.list <- createFolds(y = data[ , 1], k = nfolds)
  results <- c(rep(0.0, nfolds))
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
    R2 <- defaultSummary(rf.df)[2] 
    results[i] <- R2
    message(paste0("Fold ", i, ", mtry of ", m, ", has been completed."))
  }
  return(data.frame( # Useful for records
    data = deparse(substitute(data)), # Turns the variable name into char
    nfolds = nfolds, 
    seed = seed,
    mtry = m,
    rsquared = sum(results) / nfolds))
}

# Tuning ------------------------------------------------------------------

#     Seed 1 --------------------------------------------------------------

ntree1 <-
  do.call(
    rbind,
    lapply(
      100 * seq(3, 19, 2),
      FUN = tune.rf.ntree,
      data = mat.dg,
      nfolds = 10,
      seed = 1
    )
  )

node1 <-
  do.call(
    rbind,
    lapply(
      2 ^ (-1:6),
      FUN = tune.rf.node,
      data = mat.dg,
      nfolds = 10,
      seed = 1
    )
  )

mtry1 <-
  do.call(
    rbind,
    lapply(
      100 * seq(1, 7, 2),
      FUN = tune.rf.mtry,
      data = mat.dg,
      nfolds = 10,
      seed = 1
    )
  )

#     Seed 2 --------------------------------------------------------------

ntree2 <-
  do.call(
    rbind,
    lapply(
      100 * seq(3, 19, 2),
      FUN = tune.rf.ntree,
      data = mat.dg,
      nfolds = 10,
      seed = 2
    )
  )

node2 <-
  do.call(
    rbind,
    lapply(
      2 ^ (-1:6),
      FUN = tune.rf.node,
      data = mat.dg,
      nfolds = 10,
      seed = 2
    )
  )

mtry2 <-
  do.call(
    rbind,
    lapply(
      100 * seq(1, 7, 2),
      FUN = tune.rf.mtry,
      data = mat.dg,
      nfolds = 10,
      seed = 2
    )
  )

# Graphs ------------------------------------------------------------------

dir.create("./tuning/rforest/graphs")
# Number of Trees
temp.data <- ntree.comp
temp.data$seed <- as.factor(temp.data$seed)
ggplot(temp.data, aes(x = ntree, y = rsquared,
                       group = seed, color = seed)) + 
  geom_line(size = 1) + 
  theme_bw() + 
  labs(x = "Number of Trees", y = "R-squared", 
       title = "Random Forest Tuning - Number of Trees")
ggsave("./tuning/rforest/graphs/2017-07-13 ntree.png")

# Size of the nodes
temp.data <- node.comp
temp.data$seed <- as.factor(temp.data$seed)
ggplot(temp.data, aes(x = nodesize, y = rsquared,
                      group = seed, color = seed)) + 
  geom_line(size = 1) + 
  theme_bw() + 
  labs(x = "Size of Node", y = "R-squared", 
       title = "Random Forest Tuning - Node Size")
ggsave("./tuning/rforest/graphs/2017-07-13 nodesize.png")

# Mtry
temp.data <- mtry.comp
temp.data$seed <- as.factor(temp.data$seed)
ggplot(temp.data, aes(x = mtry, y = rsquared,
                      group = seed, color = seed)) + 
  geom_line(size = 1) + 
  theme_bw() + 
  labs(x = "Samples per Tree", y = "R-squared", 
       title = "Random Forest Tuning - Samples per Tree")
ggsave("./tuning/rforest/graphs/2017-07-13 mtry.png")
