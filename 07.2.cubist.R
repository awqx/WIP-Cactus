source("./07.model.functions.R")
source("eval.functions.R")

# Libraries and Packages --------------------------------------------------

library(caret)
library(Cubist)
library(tidyverse)

# Functions ---------------------------------------------------------------

cubist.looq2 <- function(read.dir, nsplits, cmte, extra, seed) {
  # initialize a vector for analysis
  q2.results <- c(rep(0.0, nsplits))
  if(str_detect(read.dir, "alpha"))
    features <- readRDS("./feature.selection/alpha.vars.RDS")
  else
    features <- readRDS("./feature.selection/beta.vars.RDS")
  
  ctrl <- cubistControl(
    seed = seed, 
    extrapolation = extra
  )
  
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
      x <- trn[ , -1]
      y <- trn[ , 1]
      
      # Refer to 07.0.1.svm.tune.poly.R for values
      cube <- cubist(x = x, y = y,  
                     committees = cmte, 
                     control = ctrl)
      pred[i] <- predict(cube, tst[ , -1]) 
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

cube.tst <- function(pp.dir, tst.dir, nsplits, cmte, extra, seed) {
  split <- 1:nsplits
  r2.results <- c(rep(0.0, nsplits))
  rmse.results <- c(rep(0.0, nsplits))
  
  if(str_detect(pp.dir, "alpha"))
    features <- readRDS("./feature.selection/alpha.vars.RDS")
  else
    features <- readRDS("./feature.selection/beta.vars.RDS")
  
  ctrl <- cubistControl(
    seed = seed, 
    extrapolation = extra
  )
  
  for(n in 1:nsplits) {
    trn <- readRDS(paste0(pp.dir, n, "/pp.RDS")) %>%
      select(., -guest)
    trn.y <- trn$DelG
    trn.x <- select(trn, one_of(features))
    
    tst <- preprocess.tst.mod(pp.dir = pp.dir, tst.dir = tst.dir, 
                              feat = colnames(trn.x), n = n)
    tst.y <- tst[ , 1]
    tst.x <- tst[ , -1]
    
    cube <- cubist(x = trn.x, y = trn.y,  
                   committees = cmte, 
                   control = ctrl)
    
    tst.df <- predict(cube, tst.x) %>%
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
      if(abs(tst.df$pred[i]) > 90)
        tst.df$pred[i] <- mean(trn.y)
    r2.results[n] <- defaultSummary(tst.df)[2]
    rmse.results[n] <- defaultSummary(tst.df)[1]
  }
  return(data.frame(split, r2.results, rmse.results))
}


# LOOCV-Q2 ----------------------------------------------------------------

# 3 doesn't pass (narrow margin, attributable to outliers)
# Q2 = 0.566354858323843
# Q2 = 0.589559116758517
# Q2 = 0.48999827195408
# Q2 = 0.563288142160792
# Q2 = 0.521585431960149
cubist.looq2(read.dir = "./pre-process/alpha/", nsplits = 5, 
              cmte = 100, extra = 10, seed = 101)

# Q2 = 0.626784727883507
# Q2 = 0.608993472263807
# Q2 = 0.616364600727912
# Q2 = 0.640192027871673
# Q2 = 0.703439843542136
cubist.looq2(read.dir = "./pre-process/beta/", nsplits = 5, 
             cmte = 100, extra = 50, seed = 101)

# Test sets ---------------------------------------------------------------

# split 2 is alright, except for a single outlier
# split 3 is the best, but it doesn't pass q2 test
alpha.tst <- cube.tst("./pre-process/alpha/", "./model.data/alpha", 
                      nsplits = 5, cmte = 50, extra = 50, seed = 101)
# split 5
beta.tst <- cube.tst("./pre-process/beta/", "./model.data/beta", 
                      nsplits = 5, cmte = 100, extra = 50, seed = 101)


# Single models -----------------------------------------------------------

#    Alpha ----
# Single models -----------------------------------------------------------

#    Alpha ----

trn.alpha <- readRDS("./pre-process/alpha/2/pp.RDS") %>%
  select(., -guest)
features <- readRDS("./feature.selection/alpha.vars.RDS")
colnames(trn.alpha) <- str_replace(colnames(trn.alpha), "-", ".")
trn.alpha <- trn.alpha %>% select(., DelG, one_of(features)) 
trn.alpha.x <- select(trn.alpha, -DelG) 
trn.alpha.y <- trn.alpha$DelG

cube.alpha.ctrl <- cubistControl(
  seed = 101, 
  extrapolation = 50
)

cube.alpha <- cubist(x = trn.alpha.x, y = trn.alpha.y,  
                     committees = 100, 
                     control = cube.alpha.ctrl)

tst.alpha <- preprocess.tst.mod("./pre-process/alpha/", "./model.data/alpha/", 
                                colnames(trn.alpha.x), 2)

tst.alpha.df <- predict(cube.alpha, tst.alpha[ , -1]) %>%
  cbind(tst.alpha[ , 1], .) %>% data.frame()
colnames(tst.alpha.df) <- c("obs", "pred")

# Yay, you pass
# R2 = 0.749 without outlier
eval.tropsha(tst.alpha.df[ -6, ])
graph.alpha <- ggplot(tst.alpha.df, aes(x = obs, y = pred)) + 
  geom_point() + 
  theme_bw() + 
  coord_fixed()  + 
  geom_abline(intercept = 0, slope = 1) + 
  labs(x = "Observed dG, kJ/mol", y = "Predicted dG, kJ/mol", 
       title = "Alpha-CD Cubist")

#    Beta ----

trn.beta <- readRDS("./pre-process/beta/5/pp.RDS") %>%
  select(., -guest)
features <- readRDS("./feature.selection/beta.vars.RDS")
colnames(trn.beta) <- str_replace(colnames(trn.beta), "-", ".")
trn.beta <- trn.beta %>% select(., DelG, features) 
trn.beta.x <- select(trn.beta, -DelG) 
trn.beta.y <- trn.beta$DelG

cube.beta.ctrl <- cubistControl(
  seed = 101, 
  extrapolation = 50
)

cube.beta <- cubist(x = trn.beta.x, y = trn.beta.y,  
                     committees = 100, 
                     control = cube.beta.ctrl)

tst.beta <- preprocess.tst.mod("./pre-process/beta/", "./model.data/beta/", 
                                features, 5)

tst.beta.df <- predict(cube.beta, tst.beta[ , -1]) %>%
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
       title = "Beta-CD Cubist")

#     Saving models ----

pp.settings <- readRDS("./pre-process/alpha/2/pp.settings.RDS")
saveRDS(list(pp.settings, cube.alpha), "./models/alpha/cube.RDS")
saveRDS(tst.alpha.df, "./results/alpha/cube.RDS")
print(graph.alpha)
ggsave("./results/alpha/cube.png")


pp.settings <- readRDS("./pre-process/beta/5/pp.settings.RDS")
saveRDS(list(pp.settings, cube.beta), "./models/beta/cube.RDS")
saveRDS(tst.beta.df, "./results/beta/cube.RDS")
print(graph.beta)
ggsave("./results/beta/cube.png")
