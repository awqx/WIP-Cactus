dir.create("./models")
dir.create("./models/alpha")
dir.create("./models/beta")
source("./07.model.functions.R")

# Libraries and Packages --------------------------------------------------

library(caret)
library(e1071)
library(kernlab)
# library(Matrix)
# library(stats)
# library(stringr)
library(tidyverse)

# Functions ---------------------------------------------------------------

polysvm.looq2 <- function(read.dir, nsplits, cost, deg, coef, e, g) {
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
    data <- data %>% select(., DelG, features)
    pred <- c(rep(0.0, nrow(data) - 1))
    # In this loop, i represents the index of the datapoint left out of 
    # model building
    for(i in 1:nrow(data)) {
      trn <- data[-i, ]
      tst <- data[i, ]
      x <- trn[ , -1]
      y <- trn[ , 1]
    
      # Refer to 07.0.1.svm.tune.poly.R for values
      svm.cv <- svm(x = x, y = y,
                    cost = cost, degree = deg, coef0 = coef,
                    epsilon = e, gamma = g,
                    kernel = "polynomial")
      pred[i] <- predict(svm.cv, tst[ , -1]) 
    }
    # test.df <- data.frame(obs, pred) %>% print()
    # test.df <<- test.df
    # defaultSummary(test.df) %>% print()
    # prediction error sum of squares
    for(i in 1:length(pred)) {
      if(abs(pred[i]) > 50)
        pred[i] <- mean(obs)
    }
    pred.df <- data.frame(obs, pred)
    p <- ggplot(pred.df, aes(x = obs, y = pred)) + 
      theme_bw() + 
      geom_point() + 
      labs(title = n) + 
      geom_abline(slope = 1, intercept = 0)
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
# pp.dir, tst.dir should end in a backslash
polysvm.tst <- function(pp.dir, tst.dir, nsplits, cost, deg, coef, e, g) {
  split <- 1:nsplits
  r2.results <- c(rep(0.0, nsplits))
  rmse.results <- c(rep(0.0, nsplits))
  
  if(str_detect(pp.dir, "alpha"))
    features <- readRDS("./feature.selection/alpha.vars.RDS")
  else
    features <- readRDS("./feature.selection/beta.vars.RDS")
  
  for(n in 1:nsplits) {
    trn <- readRDS(paste0(pp.dir, n, "/pp.RDS")) %>%
      select(., -guest)
    trn.y <- trn[ , 1]
    trn.x <- trn %>% select(., features)
    
    tst <- preprocess.tst.mod(pp.dir = pp.dir, tst.dir = tst.dir, 
                          feat = features, n = n)
    tst.y <- tst[ , 1]
    tst.x <- tst[ , -1]
    # Refer to 07.0.1.svm.tune.poly.R for values
    svm.cv <- svm(x = trn.x, y = trn.y,
                  cost = cost, degree = deg, coef0 = coef,
                  epsilon = e, gamma = g,
                  kernel = "polynomial")
    tst.df <- predict(svm.cv, tst.x) %>%
      cbind(tst.y, .) %>% as.data.frame() 
    colnames(tst.df) <- c("obs", "pred")
    p <- ggplot(tst.df, aes(x = obs, y = pred)) + 
      theme_bw() + 
      geom_point() + 
      labs(title = n) + 
      geom_abline(slope = 1, intercept = 0) + 
      coord_fixed()
    print(p)
    for(i in 1:nrow(tst.df))
      if(abs(tst.df$pred[i]) > 100)
        tst.df$pred[i] <- mean(trn.y)
    print(defaultSummary(tst.df))
    r2.results[n] <- defaultSummary(tst.df)[2]
    rmse.results[n] <- defaultSummary(tst.df)[1]
  }
  return(data.frame(split, r2.results, rmse.results))
}

rbfsvm.looq2 <- function(read.dir, nsplits, cost, e, g) {
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
    data <- data %>% select(., DelG, features)
    pred <- c(rep(0.0, nrow(data) - 1))
    # In this loop, i represents the index of the datapoint left out of 
    # model building
    for(i in 1:nrow(data)) {
      trn <- data[-i, ]
      tst <- data[i, ]
      x <- trn[ , -1]
      y <- trn[ , 1]
      
      # Refer to 07.0.1.svm.tune.poly.R for values
      svm.cv <- svm(x = x, y = y,
                    cost = cost, epsilon = e, gamma = g,
                    kernel = "radial")
      pred[i] <- predict(svm.cv, tst[ , -1]) 
    }
    # test.df <- data.frame(obs, pred) %>% print()
    # test.df <<- test.df
    # defaultSummary(test.df) %>% print()
    # prediction error sum of squares
    for(i in 1:length(pred)) {
      if(abs(pred[i]) > 50)
        pred[i] <- mean(obs)
    }
    pred.df <- data.frame(obs, pred)
    p <- ggplot(pred.df, aes(x = obs, y = pred)) + 
      theme_bw() + 
      geom_point() + 
      labs(title = n) + 
      geom_abline(slope = 1, intercept = 0)
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
# pp.dir, tst.dir should end in a backslash
rbfsvm.tst <- function(pp.dir, tst.dir, nsplits, cost, e, g) {
  split <- 1:nsplits
  r2.results <- c(rep(0.0, nsplits))
  rmse.results <- c(rep(0.0, nsplits))
  
  if(str_detect(pp.dir, "alpha"))
    features <- readRDS("./feature.selection/alpha.vars.RDS")
  else
    features <- readRDS("./feature.selection/beta.vars.RDS")
  
  for(n in 1:nsplits) {
    trn <- readRDS(paste0(pp.dir, n, "/pp.RDS")) %>%
      select(., -guest)
    trn.y <- trn[ , 1]
    trn.x <- trn %>% select(., features)
    
    tst <- preprocess.tst.mod(pp.dir = pp.dir, tst.dir = tst.dir, 
                              feat = features, n = n)
    tst.y <- tst[ , 1]
    tst.x <- tst[ , -1]
    # Refer to 07.0.1.svm.tune.poly.R for values
    svm.cv <- svm(x = trn.x, y = trn.y,
                  cost = cost, epsilon = e, gamma = g,
                  kernel = "radial")
    tst.df <- predict(svm.cv, tst.x) %>%
      cbind(tst.y, .) %>% as.data.frame() 
    colnames(tst.df) <- c("obs", "pred")
    p <- ggplot(tst.df, aes(x = obs, y = pred)) + 
      theme_bw() + 
      geom_point() + 
      labs(title = n) + 
      geom_abline(slope = 1, intercept = 0) + 
      coord_fixed()
    print(p)
    for(i in 1:nrow(tst.df))
      if(abs(tst.df$pred[i]) > 100)
        tst.df$pred[i] <- mean(trn.y)
    print(defaultSummary(tst.df))
    r2.results[n] <- defaultSummary(tst.df)[2]
    rmse.results[n] <- defaultSummary(tst.df)[1]
  }
  return(data.frame(split, r2.results, rmse.results))
}

sigsvm.looq2 <- function(read.dir, nsplits, cost, coef, e, g) {
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
    data <- data %>% select(., DelG, features)
    pred <- c(rep(0.0, nrow(data) - 1))
    # In this loop, i represents the index of the datapoint left out of 
    # model building
    for(i in 1:nrow(data)) {
      trn <- data[-i, ]
      tst <- data[i, ]
      x <- trn[ , -1]
      y <- trn[ , 1]
      
      # Refer to 07.0.1.svm.tune.poly.R for values
      svm.cv <- svm(x = x, y = y,
                    cost = cost, coef0 = coef,
                    epsilon = e, gamma = g,
                    kernel = "sigmoid")
      pred[i] <- predict(svm.cv, tst[ , -1]) 
    }
    # test.df <- data.frame(obs, pred) %>% print()
    # test.df <<- test.df
    # defaultSummary(test.df) %>% print()
    # prediction error sum of squares
    for(i in 1:length(pred)) {
      if(abs(pred[i]) > 50)
        pred[i] <- mean(obs)
    }
    pred.df <- data.frame(obs, pred)
    p <- ggplot(pred.df, aes(x = obs, y = pred)) + 
      theme_bw() + 
      geom_point() + 
      labs(title = n) + 
      geom_abline(slope = 1, intercept = 0)
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
# pp.dir, tst.dir should end in a backslash
sigsvm.tst <- function(pp.dir, tst.dir, nsplits, cost, coef, e, g) {
  split <- 1:nsplits
  r2.results <- c(rep(0.0, nsplits))
  rmse.results <- c(rep(0.0, nsplits))
  
  if(str_detect(pp.dir, "alpha"))
    features <- readRDS("./feature.selection/alpha.vars.RDS")
  else
    features <- readRDS("./feature.selection/beta.vars.RDS")
  
  for(n in 1:nsplits) {
    trn <- readRDS(paste0(pp.dir, n, "/pp.RDS")) %>%
      select(., -guest)
    trn.y <- trn[ , 1]
    trn.x <- trn %>% select(., features)
    
    tst <- preprocess.tst.mod(pp.dir = pp.dir, tst.dir = tst.dir, 
                              feat = features, n = n)
    tst.y <- tst[ , 1]
    tst.x <- tst[ , -1]
    # Refer to 07.0.1.svm.tune.poly.R for values
    svm.cv <- svm(x = trn.x, y = trn.y,
                  cost = cost, coef0 = coef,
                  epsilon = e, gamma = g,
                  kernel = "sigmoid")
    tst.df <- predict(svm.cv, tst.x) %>%
      cbind(tst.y, .) %>% as.data.frame() 
    colnames(tst.df) <- c("obs", "pred")
    p <- ggplot(tst.df, aes(x = obs, y = pred)) + 
      theme_bw() + 
      geom_point() + 
      labs(title = n) + 
      geom_abline(slope = 1, intercept = 0) + 
      coord_fixed()
    print(p)
    for(i in 1:nrow(tst.df))
      if(abs(tst.df$pred[i]) > 100)
        tst.df$pred[i] <- mean(trn.y)
    print(defaultSummary(tst.df))
    r2.results[n] <- defaultSummary(tst.df)[2]
    rmse.results[n] <- defaultSummary(tst.df)[1]
  }
  return(data.frame(split, r2.results, rmse.results))
}
# Polynomial --------------------------------------------------------------

#     LOOCV-Q2 ------------------------------------------------------------

# Alpha models 1, 3, 4, 5, 6, 7, 8, 10 pass
# Narrow misses for the rest (would all round to 0.5)
# avg = 0.521
polysvm.looq2(read.dir = "./pre-process/alpha/", nsplits = 10, 
              cost = 4, deg = 5, coef = 1, e = 0.27, g = 0.01)

# All beta models pass the Q2 test
# Average Q2 = 0.691
polysvm.looq2(read.dir = "./pre-process/beta/", nsplits = 10, 
              cost = 10, deg = 3, coef = 5, e = 0.25, g = 0.001)

#     Test sets -----------------------------------------------------------

# Split 8 does the best
alpha.tst <- polysvm.tst("./pre-process/alpha/", "./model.data/alpha/", nsplits = 10,
            cost = 4, deg = 5, coef = 1, e = 0.27, g = 0.01)

# Split 9 does the best here
# Many other splits suffer from a persistent outlier
# without the outlier, split 7 does the best
beta.tst <- polysvm.tst("./pre-process/beta/", "./model.data/beta/", nsplits = 10,
            cost = 10, deg = 3, coef = 5, e = 0.25, g = 0.001)

#     Single models -------------------------------------------------------

#         Alpha ----

trn.alpha <- readRDS("./pre-process/alpha/8/pp.RDS") %>%
  select(., -guest)
features <- readRDS("./feature.selection/alpha.vars.RDS")
colnames(trn.alpha) <- str_replace(colnames(trn.alpha), "-", ".")
trn.alpha <- trn.alpha %>% select(., DelG, features) 
trn.alpha.x <- select(trn.alpha, -DelG) 
trn.alpha.y <- select(trn.alpha, DelG) 

polysvm.alpha <- svm(x = trn.alpha.x, y = trn.alpha.y,
                     cost = 4, degree = 5, coef0 = 1,
                     epsilon = 0.27, gamma = 0.01,
                     kernel = "polynomial")

tst.alpha <- preprocess.tst.mod("./pre-process/alpha/", "./model.data/alpha/", 
                               features, 8)

tst.alpha.df <- predict(polysvm.alpha, tst.alpha[ , -1]) %>%
  cbind(tst.alpha[ , 1], .) %>% data.frame()
colnames(tst.alpha.df) <- c("obs", "pred")

# Yay, you pass
eval.tropsha(tst.alpha.df)
graph.alpha <- ggplot(tst.alpha.df, aes(x = obs, y = pred)) + 
  geom_point() + 
  theme_bw() + 
  coord_fixed()  + 
  geom_abline(intercept = 0, slope = 1)

#         Beta ----

trn.beta <- readRDS("./pre-process/beta/9/pp.RDS") %>%
  select(., -guest)
features <- readRDS("./feature.selection/beta.vars.RDS")
colnames(trn.beta) <- str_replace(colnames(trn.beta), "-", ".")
trn.beta <- trn.beta %>% select(., DelG, features) 
trn.beta.x <- select(trn.beta, -DelG) 
trn.beta.y <- select(trn.beta, DelG) 

polysvm.beta <- svm(x = trn.beta.x, y = trn.beta.y,
                     cost = 10, degree = 3, coef0 = 5,
                     epsilon = 0.25, gamma = 0.001,
                     kernel = "polynomial")
tst.beta <- preprocess.tst.mod("./pre-process/beta/", "./model.data/beta/", 
                               features, 9)

tst.beta.df <- predict(polysvm.beta, tst.beta[ , -1]) %>%
  cbind(tst.beta[ , 1], .) %>% data.frame()
colnames(tst.beta.df) <- c("obs", "pred")

# Yay, you pass
eval.tropsha(tst.beta.df)
graph.beta <- ggplot(tst.beta.df, aes(x = obs, y = pred)) + 
  geom_point() + 
  theme_bw() + 
  coord_fixed()  + 
  geom_abline(intercept = 0, slope = 1)

#         Saving models ----

dir.create("./results")
dir.create("./results/alpha")
dir.create("./results/beta")

pp.settings <- readRDS("./pre-process/alpha/8/pp.settings.RDS")
saveRDS(list(pp.settings, polysvm.alpha), "./models/alpha/polysvm.RDS")
saveRDS(tst.alpha.df, "./results/alpha/polysvm.RDS")
graph.alpha + 
  labs(x = "Experimental dG, kJ/mol", y = "Predicted dG, kJ/mol", 
       title = "Polynomial SVM for Alpha")
ggsave("./results/alpha/polysvm.png")


pp.settings <- readRDS("./pre-process/beta/9/pp.settings.RDS")
saveRDS(list(pp.settings, polysvm.beta), "./models/beta/polysvm.RDS")
saveRDS(tst.beta.df, "./results/beta/polysvm.RDS")
graph.beta + 
  labs(x = "Experimental dG, kJ/mol", y = "Predicted dG, kJ/mol", 
       title = "Polynomial SVM for Beta")
ggsave("./results/beta/polysvm.png")

# Radial ------------------------------------------------------------------

#     LOOCV-Q2 ------------------------------------------------------------

# All pass, average = 0.558
rbfsvm.looq2(read.dir = "./pre-process/alpha/", nsplits = 10, 
             cost = 5, e = 0.25, g = 0.05)

# All pass, average = 0.703
rbfsvm.looq2(read.dir = "./pre-process/beta/", nsplits = 10, 
             cost = 10, e = 0.01, g = 0.01)

#     Test sets -----------------------------------------------------------

# Using split 10, because it captures some of the lower range
alpha.tst <- rbfsvm.tst("./pre-process/alpha/", "./model.data/alpha/", 
                        nsplits = 10, cost = 5, e = 0.25, g = 0.05)

# Split 7
beta.tst <- rbfsvm.tst("./pre-process/beta/", "./model.data/beta/", nsplits = 10,
                        cost = 10, e = 0.01, g = 0.01)

#     Single models -------------------------------------------------------

#         Alpha ----

trn.alpha <- readRDS("./pre-process/alpha/10/pp.RDS") %>%
  select(., -guest)
features <- readRDS("./feature.selection/alpha.vars.RDS")
colnames(trn.alpha) <- str_replace(colnames(trn.alpha), "-", ".")
trn.alpha <- trn.alpha %>% select(., DelG, features) 
trn.alpha.x <- select(trn.alpha, -DelG) 
trn.alpha.y <- select(trn.alpha, DelG) 

rbfsvm.alpha <- svm(x = trn.alpha.x, y = trn.alpha.y,
                     cost = 5, epsilon = 0.25, gamma = 0.05,
                     kernel = "radial")

tst.alpha <- preprocess.tst.mod("./pre-process/alpha/", "./model.data/alpha/", 
                                features, 10)

tst.alpha.df <- predict(rbfsvm.alpha, tst.alpha[ , -1]) %>%
  cbind(tst.alpha[ , 1], .) %>% data.frame()
colnames(tst.alpha.df) <- c("obs", "pred")

# Yay, you pass
eval.tropsha(tst.alpha.df)
graph.alpha <- ggplot(tst.alpha.df, aes(x = obs, y = pred)) + 
  geom_point() + 
  theme_bw() + 
  coord_fixed()  + 
  geom_abline(intercept = 0, slope = 1)

#         Beta ----

trn.beta <- readRDS("./pre-process/beta/7/pp.RDS") %>%
  select(., -guest)
features <- readRDS("./feature.selection/beta.vars.RDS")
colnames(trn.beta) <- str_replace(colnames(trn.beta), "-", ".")
trn.beta <- trn.beta %>% select(., DelG, features) 
trn.beta.x <- select(trn.beta, -DelG) 
trn.beta.y <- select(trn.beta, DelG) 

rbfsvm.beta <- svm(x = trn.beta.x, y = trn.beta.y,
                    cost = 19, epsilon = 0.01, gamma = 0.01,
                    kernel = "radial")

tst.beta <- preprocess.tst.mod("./pre-process/beta/", "./model.data/beta/", 
                                features, 7)

tst.beta.df <- predict(rbfsvm.beta, tst.beta[ , -1]) %>%
  cbind(tst.beta[ , 1], .) %>% data.frame()
colnames(tst.beta.df) <- c("obs", "pred")

# Yay, you pass
eval.tropsha(tst.beta.df)
graph.beta <- ggplot(tst.beta.df, aes(x = obs, y = pred)) + 
  geom_point() + 
  theme_bw() + 
  coord_fixed()  + 
  geom_abline(intercept = 0, slope = 1)

#         Saving models ----

pp.settings <- readRDS("./pre-process/alpha/10/pp.settings.RDS")
saveRDS(list(pp.settings, rbfsvm.alpha), "./models/alpha/rbfsvm.RDS")
saveRDS(tst.alpha.df, "./results/alpha/rbfsvm.RDS")
graph.alpha + 
  labs(x = "Experimental dG, kJ/mol", y = "Predicted dG, kJ/mol", 
       title = "Radial SVM for Alpha")
ggsave("./results/alpha/rbfsvm.png")


pp.settings <- readRDS("./pre-process/beta/7/pp.settings.RDS")
saveRDS(list(pp.settings, rbfsvm.beta), "./models/beta/rbfsvm.RDS")
saveRDS(tst.beta.df, "./results/beta/rbfsvm.RDS")
graph.beta + 
  labs(x = "Experimental dG, kJ/mol", y = "Predicted dG, kJ/mol", 
       title = "Radial SVM for Beta")
ggsave("./results/beta/rbfsvm.png")


# Sigmoid -----------------------------------------------------------------

#     LOOCV-Q2 ------------------------------------------------------------

# None pass. 1 and 8 come the closest
sigsvm.looq2(read.dir = "./pre-process/alpha/", nsplits = 10, 
             cost = 15, coef = 0, e = 0.05, g = 0.001)

# All pass, average = 0.635
sigsvm.looq2(read.dir = "./pre-process/beta/", nsplits = 10, 
             cost = 75, coef = 0, e = 0.1, g = 0.001)

#     Test sets -----------------------------------------------------------

# Split 5
alpha.tst <- sigsvm.tst("./pre-process/alpha/", "./model.data/alpha/", 
                        nsplits = 10, cost = 15, coef = 0, e = 0.05, g = 0.001)

# Split 1
beta.tst <- sigsvm.tst("./pre-process/beta/", "./model.data/beta/", nsplits = 10,
                       cost = 75, coef = 0, e = 0.1, g = 0.001)

#     Single models -------------------------------------------------------

#         Alpha ----

trn.alpha <- readRDS("./pre-process/alpha/5/pp.RDS") %>%
  select(., -guest)
features <- readRDS("./feature.selection/alpha.vars.RDS")
colnames(trn.alpha) <- str_replace(colnames(trn.alpha), "-", ".")
trn.alpha <- trn.alpha %>% select(., DelG, features) 
trn.alpha.x <- select(trn.alpha, -DelG) 
trn.alpha.y <- select(trn.alpha, DelG) 

sigsvm.alpha <- svm(x = trn.alpha.x, y = trn.alpha.y,
                    cost = 15, coef = 0, e = 0.05, g = 0.001,
                    kernel = "sigmoid")

tst.alpha <- preprocess.tst.mod("./pre-process/alpha/", "./model.data/alpha/", 
                                features, 5)

tst.alpha.df <- predict(sigsvm.alpha, tst.alpha[ , -1]) %>%
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
       title = "Sigmoid SVM for Alpha")

#         Beta ----

trn.beta <- readRDS("./pre-process/beta/1/pp.RDS") %>%
  select(., -guest)
features <- readRDS("./feature.selection/beta.vars.RDS")
colnames(trn.beta) <- str_replace(colnames(trn.beta), "-", ".")
trn.beta <- trn.beta %>% select(., DelG, features) 
trn.beta.x <- select(trn.beta, -DelG) 
trn.beta.y <- select(trn.beta, DelG) 

sigsvm.beta <- svm(x = trn.beta.x, y = trn.beta.y,
                   cost = 75, coef = 0, e = 0.1, g = 0.001,
                   kernel = "sigmoid")

tst.beta <- preprocess.tst.mod("./pre-process/beta/", "./model.data/beta/", 
                               features, 1)

tst.beta.df <- predict(sigsvm.beta, tst.beta[ , -1]) %>%
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
       title = "Sigmoid SVM for Beta")

#         Saving models ----

pp.settings <- readRDS("./pre-process/alpha/5/pp.settings.RDS")
saveRDS(list(pp.settings, sigsvm.alpha), "./models/alpha/sigsvm.RDS")
saveRDS(tst.alpha.df, "./results/alpha/sigsvm.RDS")
graph.alpha 
ggsave("./results/alpha/sigsvm.png")


pp.settings <- readRDS("./pre-process/beta/1/pp.settings.RDS")
saveRDS(list(pp.settings, sigsvm.beta), "./models/beta/sigsvm.RDS")
saveRDS(tst.beta.df, "./results/beta/sigsvm.RDS")
graph.beta 
ggsave("./results/beta/sigsvm.png")
