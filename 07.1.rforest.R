# Libraries and Packages --------------------------------------------------

library(caret)
library(randomForest)
library(tidyverse)

source("./07.model.functions.R")

# Functions ---------------------------------------------------------------

rf.looq2 <- function(read.dir, nsplits, ntree, node, m) {
  trn.split <- 1:nsplits
  q2.results <- c(rep(0.0, nsplits))
  
  features <- find.features(read.dir)
  
  for(n in 1:nsplits) {
    data <- readRDS(paste0(read.dir, n, "/pp.RDS")) %>%
      select(., -guest)
    obs <- data[ , 1]
    data <- data %>% select(., DelG, one_of(features)) 
    pred <- c(rep(0.0, nrow(data) - 1))
    
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
    # # Handling outliers
    # for(i in 1:length(pred)) {
    #   if(abs(pred[i]) > 85)
    #     pred[i] <- mean(obs)
    # }
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
    TSS <- sum((obs - mean(obs))^2)
    q2 <- 1 - PRESS/TSS
    message("Q2 = ", q2)
    q2.results[n] <- q2
  }
  return(data.frame(trn.split, q2.results))
}

rf.tst <- function(pp.dir, tst.dir, nsplits, ntree, node, m) {
  split <- 1:nsplits
  # r2.results <- c(rep(0.0, nsplits))
  # rmse.results <- c(rep(0.0, nsplits))
  
  if(str_detect(pp.dir, "alpha"))
    features <- readRDS("./feature.selection/alpha.vars.RDS")
  else if(str_detect(pp.dir, 'beta'))
    features <- readRDS("./feature.selection/beta.vars.RDS")
  else
    features <- readRDS('feature.selection/gamma.vars.RDS')
  
  results.all <- data.frame() 
  
  for(i in 1:nsplits) {
    trn <- readRDS(paste0(pp.dir, i, "/pp.RDS")) %>%
      select(., -guest)
    trn.y <- trn[ , 1]
    trn.x <- trn %>% select(., one_of(features))
    
    # Refer to 07.0.1.svm.tune.poly.R for values
    rf <- randomForest(
      x = trn.x, y = trn.y,
      ntree = ntree, nodesize = node,
      mtry = m, importance = T
    )
    
    tst.all <- data.frame()
    for(j in 1:nsplits) {
      tst <- preprocess.tst.mod(pp.dir = pp.dir, tst.dir = tst.dir, 
                                feat = colnames(trn.x), n = j)
      tst.y <- tst[ , 1]
      tst.x <- tst[ , -1]
      
      tst.df <- predict(rf, tst.x) %>%
        cbind(tst.y, .) %>% as.data.frame() 
      colnames(tst.df) <- c("obs", "pred")
      tst.df <- tst.df %>% mutate(split = j)
      # Controlling for extreme outliers
      for(k in 1:nrow(tst.df))
        if(abs(tst.df$pred[k]) > 100)
          tst.df$pred[k] <- mean(trn.y)
      
      results.tst <- data.frame(trn.split = i, tst.split = j, 
                                r2 = defaultSummary(tst.df)[2], 
                                rmse = defaultSummary(tst.df)[1])
      results.all <- rbind(results.all, results.tst)
      tst.all <- rbind(tst.all, tst.df)
    }
    
    tst.all$split <- as.factor(tst.all$split)
    p <- ggplot(tst.all, aes(x = obs, y = pred, color = split)) + 
      theme_bw() + 
      geom_point() + 
      labs(title = i) + 
      geom_abline(slope = 1, intercept = 0) + 
      coord_fixed()
    print(p)
  }
  row.names(results.all) <- NULL
  return(results.all)
}

rf.tst.splits <- function(pp.dir, tst.dir, feat, nsplits, model) {
  tst.all <- data.frame()
  for (j in 1:nsplits) {
    tst <- preprocess.tst.mod(pp.dir = pp.dir, tst.dir = tst.dir, 
                              feat = feat, n = j)
    tst.y <- tst[ , 1]
    tst.x <- tst[ , -1]
    
    tst.df <- predict(model, tst.x) %>%
      cbind(tst.y, .) %>% as.data.frame() 
    colnames(tst.df) <- c("obs", "pred")
    tst.df <- tst.df %>% mutate(split = j)
    
    for(k in 1:nrow(tst.df))
      if(abs(tst.df$pred[k]) > 100)
        tst.df$pred[k] <- mean(trn.y)
    
    tst.all <- rbind(tst.all, tst.df)
  }
  tst.all$split <- as.factor(tst.all$split)
  return(tst.all)
}
# LOOCV-Q2 analysis -------------------------------------------------------

# Alpha
# Everything except 2
rf.alpha.q2 <- rf.looq2("./pre-process/alpha/", nsplits = 10,
                        ntree = 75, node = 3, m = 13)


# Beta
# All pass
rf.beta.q2 <-  rf.looq2("./pre-process/beta/", nsplits = 10,
                         ntree = 50, node = 5, m = 4)

# Gamma - none pass
rf.gamma.q2 <-  rf.looq2("./pre-process/gamma/", nsplits = 10,
                        ntree = 250, node = 5, m = 10)

# Test sets ---------------------------------------------------------------

# ALPHA - 6
# Q2 = 0.6321
alpha.tst <- rf.tst("./pre-process/alpha/", "./model.data/alpha/", 
                    nsplits = 10, ntree = 75, node = 3, m = 13)
alpha.1to1 <- alpha.tst %>% filter(trn.split == tst.split)
alpha.avg <- data.table(alpha.tst, key = 'trn.split')
alpha.avg <- alpha.avg[ , list(r2 = mean(r2), 
                               rmse = mean(rmse)), 
                        by = trn.split] %>% print()
saveRDS(alpha.avg, 'results/alpha/rf.avg.RDS')
# BETA - 3
# Q2 for this is 0.6733
# All of them are pretty good
beta.tst <- rf.tst("./pre-process/beta/", "./model.data/beta/", 
                   nsplits = 10, ntree = 50, node = 5, m = 4)
beta.1to1 <- beta.tst %>% filter(trn.split == tst.split)
beta.avg <- data.table(beta.tst, key = 'trn.split')
beta.avg <- beta.avg[ , list(r2 = mean(r2), 
                               rmse = mean(rmse)), 
                        by = trn.split] %>% print()
saveRDS(beta.avg, 'results/beta/rf.avg.RDS')

# GAMMA - 4
# Q2 is 0.2849
gamma.tst <- rf.tst("./pre-process/gamma/", "./model.data/gamma/", 
                    nsplits = 10, ntree = 250, node = 5, m = 10)
gamma.1to1 <- gamma.tst %>% filter(trn.split == tst.split)
gamma.avg <- data.table(gamma.tst, key = 'trn.split')
gamma.avg <- gamma.avg[ , list(r2 = mean(r2), 
                             rmse = mean(rmse)), 
                      by = trn.split] %>% print()
saveRDS(gamma.avg, 'results/gamma/rf.avg.RDS')


# Single models -----------------------------------------------------------

#    Alpha ----
trn.alpha <- readRDS("./pre-process/alpha/6/pp.RDS") %>%
  select(., -guest)
features <- readRDS("./feature.selection/alpha.vars.RDS")
colnames(trn.alpha) <- str_replace(colnames(trn.alpha), "-", ".")
trn.alpha <- trn.alpha %>% select(., DelG, features) 
trn.alpha.x <- select(trn.alpha, -DelG) 
trn.alpha.y <- trn.alpha$DelG

rf.alpha <- randomForest(x = trn.alpha.x, y = trn.alpha.y, 
                         ntree = 100, nodesize = 3, mtry = 13)
tst.alpha <- rf.tst.splits("pre-process/alpha/", "model.data/alpha/", 
                            features, 10, rf.alpha) %>% print()

# Yay, you pass
eval.tropsha(tst.alpha)
graph.alpha <- ggplot(tst.alpha, aes(x = obs, y = pred, color = split)) + 
  geom_point() + 
  theme_bw() + 
  coord_fixed()  + 
  geom_abline(intercept = 0, slope = 1) + 
  labs(x = "Observed dG, kJ/mol", y = "Predicted dG, kJ/mol", 
       title = "Alpha-CD Random Forest", color = "Test split")
print(graph.alpha)

#     Beta ----

trn.beta <- readRDS("./pre-process/beta/3/pp.RDS") %>%
  select(., -guest)
features <- readRDS("./feature.selection/beta.vars.RDS")
colnames(trn.beta) <- str_replace(colnames(trn.beta), "-", ".")
trn.beta <- trn.beta %>% select(., DelG, features) 
trn.beta.x <- select(trn.beta, -DelG) 
trn.beta.y <- trn.beta$DelG

rf.beta <- randomForest(x = trn.beta.x, y = trn.beta.y, 
                         ntree = 50, nodesize = 5, mtry = 4)
tst.beta <- rf.tst.splits("pre-process/beta/", "model.data/beta/", 
                           features, 10, rf.beta) 

# Yay, you pass
eval.tropsha(tst.beta)
graph.beta <- ggplot(tst.beta, aes(x = obs, y = pred, color = split)) + 
  geom_point() + 
  theme_bw() + 
  coord_fixed()  + 
  geom_abline(intercept = 0, slope = 1) + 
  labs(x = "Observed dG, kJ/mol", y = "Predicted dG, kJ/mol", 
       title = "Beta-CD Random Forest", color = "Test split")
print(graph.beta)

#    Gamma ----

trn.gamma <- readRDS("./pre-process/gamma/4/pp.RDS") %>%
  select(., -guest)
features <- readRDS("./feature.selection/gamma.vars.RDS")
colnames(trn.gamma) <- str_replace(colnames(trn.gamma), "-", ".")
trn.gamma <- trn.gamma %>% select(., DelG, features) 
trn.gamma.x <- select(trn.gamma, -DelG) 
trn.gamma.y <- trn.gamma$DelG

rf.gamma <- randomForest(x = trn.gamma.x, y = trn.gamma.y, 
                         ntree = 250, nodesize = 5, mtry = 10)
tst.gamma <- rf.tst.splits("pre-process/gamma/", "model.data/gamma/", 
                           features, 10, rf.gamma)
# Oh no, you fail
eval.tropsha(tst.gamma)
graph.gamma <- ggplot(tst.gamma, aes(x = obs, y = pred, color = split)) + 
  geom_point() + 
  theme_bw() + 
  coord_fixed()  + 
  geom_abline(intercept = 0, slope = 1) + 
  labs(x = "Observed dG, kJ/mol", y = "Predicted dG, kJ/mol", 
       title = "Gamma-CD Random Forest", color = "Test split")
print(graph.gamma)

#     Saving models ----

pp.settings <- readRDS("./pre-process/alpha/6/pp.settings.RDS")
saveRDS(list(pp.settings, rf.alpha), "./models/alpha/rf.RDS")
saveRDS(tst.alpha.df, "./results/alpha/rf.RDS")
print(graph.alpha)
ggsave("./results/alpha/rf.png")


pp.settings <- readRDS("./pre-process/beta/3/pp.settings.RDS")
saveRDS(list(pp.settings, rf.beta), "./models/beta/rf.RDS")
saveRDS(tst.beta.df, "./results/beta/rf.RDS")
print(graph.beta)
ggsave("./results/beta/rf.png")

pp.settings <- readRDS("./pre-process/gamma/4/pp.settings.RDS")
saveRDS(list(pp.settings, rf.gamma), "./models/gamma/rf.RDS")
saveRDS(tst.gamma, "./results/gamma/rf.RDS")
print(graph.gamma)
ggsave("./results/gamma/rf.png")
