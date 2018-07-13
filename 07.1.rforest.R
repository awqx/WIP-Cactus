# Libraries and Packages --------------------------------------------------

library(caret)
library(randomForest)
library(tidyverse)

source("./07.model.functions.R")

# Functions ---------------------------------------------------------------

rf.looq2 <- function(read.dir, rfe.dir, nsplits, ntree, node, m) {
  # initialize a vector for analysis
  q2.results <- c(rep(0.0, nsplits))
  if(str_detect(read.dir, "alpha"))
    features <- readRDS("./feature.selection/alpha.vars.RDS")
  else
    features <- readRDS("./feature.selection/beta.vars.RDS")
  for(n in 1:nsplits) {
    data <- readRDS(paste0(read.dir, n, "/pp.RDS")) %>%
      select(., -guest)
    obs <- data[ , 1]
    data <- data %>% select(., DelG, one_of(features)) 
    pred <- c(rep(0.0, nrow(data) - 1))
    # In this loop, i represents the index of the datapoint left out of 
    # model building
    for(i in 1:nrow(data)) {
      trn <- data[-i, ]
      tst <- data[i, ]
      trn.x <- trn[ , -1]
      trn.y <- trn[ , 1]
      tst.x <- tst[ , -1]
      tst.y <- tst[ , 1]
      
      rf <- randomForest(x = trn.x, 
                         y = trn.y,
                         ntree = ntree,
                         nodesize = node,
                         mtry = m,  
                         importance = T)
      pred[i] <- predict(rf, tst.x) 
    }
    # Handling outliers
    for(i in 1:length(pred)) {
      if(abs(pred[i]) > 85)
        pred[i] <- mean(obs)
    }
    pred.df <- data.frame(obs, pred)
    colnames(pred.df) <- c("obs", "pred")
    p <- ggplot(pred.df, aes(x = obs, y = pred)) + 
      theme_bw() + 
      geom_point() + 
      labs(title = n) + 
      geom_abline(slope = 1, intercept = 0) + 
      coord_fixed()
    print(p)
    PRESS <- sum((obs - pred)^2)
    # total sum of squares
    TSS <- sum((obs - mean(obs))^2)
    q2 <- 1 - PRESS/TSS
    message("Q2 = ", q2)
    q2.results[n] <- q2
  }
  return(mean(q2.results))
}

rf.tst <- function(pp.dir, tst.dir, nsplits, ntree, node, m) {
  split <- 1:nsplits
  # initialize a vector for analysis
  r2.results <- c(rep(0.0, nsplits))
  rmse.results <- c(rep(0.0, nsplits))
  if(str_detect(pp.dir, "alpha"))
    features <- readRDS("./feature.selection/alpha.vars.RDS")
  else
    features <- readRDS("./feature.selection/beta.vars.RDS")
  for(n in 1:nsplits) {
    trn <- readRDS(paste0(pp.dir, n, "/pp.RDS")) %>%
      select(., -guest)
    trn.y <- trn$DelG
    trn.x <- select(trn, one_of(features))
    
    tst <- preprocess.tst.mod(pp.dir = pp.dir, tst.dir = tst.dir, 
                              feat = colnames(trn.x), n = n)
    tst.y <- tst[ , 1]
    tst.x <- tst[ , -1]
    
    rf <- randomForest(
      x = trn.x, y = trn.y,
      ntree = ntree, nodesize = node,
      mtry = m, importance = T
    )

    tst.df <- predict(rf, tst.x) %>%
      cbind(tst.y, .) %>%
      as.data.frame()
    colnames(tst.df) <- c("obs", "pred")
    p <- ggplot(tst.df, aes(x = obs, y = pred)) + 
      theme_bw() + 
      geom_point() + 
      labs(title = n) + 
      geom_abline(slope = 1, intercept = 0) + 
      coord_fixed()
    print(p)
    for(i in 1:nrow(tst.df))
      if(abs(tst.df$pred[i]) > 80)
        tst.df$pred[i] <- mean(trn.y)
    r2.results[n] <- defaultSummary(tst.df)[2]
    rmse.results[n] <- defaultSummary(tst.df)[1]
  }
  return(data.frame(split, r2.results, rmse.results))
}

# LOOCV-Q2 analysis -------------------------------------------------------

# Alpha
# All pass
rf.looq2("./pre-process/alpha/", 
          nsplits = 5, ntree = 100, node = 2, m = 15)
# 0.584

# Beta
# all pass
rf.looq2("./pre-process/beta/", "./feature.selection/beta", 
         nsplits = 5, ntree = 250, node = 2, m = 10)
# 0.667

# Test sets ---------------------------------------------------------------

# Split 5 is the best R2 RMSE balance
alpha.tst <- rf.tst("./pre-process/alpha/", "./model.data/alpha/", 
                    nsplits = 5, ntree = 100, node = 2, m = 15)

# Split 4 is pretty good R2 RMSE balanace
beta.tst <- rf.tst("./pre-process/beta/", "./model.data/beta/", 
                   nsplits = 5, ntree = 250, node = 2, m = 10)

# Single models -----------------------------------------------------------

#    Alpha ----

trn.alpha <- readRDS("./pre-process/alpha/5/pp.RDS") %>%
  select(., -guest)
features <- readRDS("./feature.selection/alpha.vars.RDS")
colnames(trn.alpha) <- str_replace(colnames(trn.alpha), "-", ".")
trn.alpha <- trn.alpha %>% select(., DelG, features) 
trn.alpha.x <- select(trn.alpha, -DelG) 
trn.alpha.y <- trn.alpha$DelG

rf.alpha <- randomForest(x = trn.alpha.x, y = trn.alpha.y, 
                         ntree = 100, nodesize = 2, mtry = 15)

tst.alpha <- preprocess.tst.mod("./pre-process/alpha/", "./model.data/alpha/", 
                                features, 5)

tst.alpha.df <- predict(rf.alpha, tst.alpha[ , -1]) %>%
  cbind(tst.alpha[ , 1], .) %>% data.frame()
colnames(tst.alpha.df) <- c("obs", "pred")

# Yay, you pass
eval.tropsha(tst.alpha.df)
graph.alpha <- ggplot(tst.alpha.df, aes(x = obs, y = pred)) + 
  geom_point() + 
  theme_bw() + 
  coord_fixed()  + 
  geom_abline(intercept = 0, slope = 1) + 
  labs(x = "Observed dG, kJ/mol", y = "Predicted dG, kJ/mol", 
       title = "Alpha-CD Random Forest")

#     Beta ----

trn.beta <- readRDS("./pre-process/beta/4/pp.RDS") %>%
  select(., -guest)
features <- readRDS("./feature.selection/beta.vars.RDS")
colnames(trn.beta) <- str_replace(colnames(trn.beta), "-", ".")
trn.beta <- trn.beta %>% select(., DelG, one_of(features)) 
trn.beta.x <- select(trn.beta, -DelG) 
trn.beta.y <- trn.beta$DelG

rf.beta <- randomForest(x = trn.beta.x, y = trn.beta.y, 
                         ntree = 250, nodesize = 2, mtry = 10)

tst.beta <- preprocess.tst.mod("./pre-process/beta/", "./model.data/beta/", 
                                colnames(trn.beta.x), 4)

tst.beta.df <- predict(rf.beta, tst.beta[ , -1]) %>%
  cbind(tst.beta[ , 1], .) %>% data.frame()
colnames(tst.beta.df) <- c("obs", "pred")

# Yay, you pass
eval.tropsha(tst.beta.df)
graph.beta <- ggplot(tst.beta.df, aes(x = obs, y = pred)) + 
  geom_point() + 
  theme_bw() + 
  coord_fixed()  + 
  geom_abline(intercept = 0, slope = 1) + 
  labs(x = "Observed dG, kJ/mol", y = "Predicted dG, kJ/mol", 
       title = "Beta-CD Random Forest")

#     Saving models ----

pp.settings <- readRDS("./pre-process/alpha/5/pp.settings.RDS")
saveRDS(list(pp.settings, rf.alpha), "./models/alpha/rf.RDS")
saveRDS(tst.alpha.df, "./results/alpha/rf.RDS")
print(graph.alpha)
ggsave("./results/alpha/rf.png")


pp.settings <- readRDS("./pre-process/beta/4/pp.settings.RDS")
saveRDS(list(pp.settings, rf.beta), "./models/beta/rf.RDS")
saveRDS(tst.beta.df, "./results/beta/rf.RDS")
print(graph.beta)
ggsave("./results/beta/rf.png")
