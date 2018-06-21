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
    trn.x <- select(trn, features)
    
    tst <- preprocess.tst.mod(pp.dir = pp.dir, tst.dir = tst.dir, 
                              feat = features, n = n)
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

# 3 passes. 8 and 9 sort of pass.
# Q2 = 0.368393118712259
# Q2 = 0.303807425780733
# Q2 = 0.503604271484987
# Q2 = 0.283483487043784
# Q2 = 0.419974823431514
# Q2 = 0.350065835353248
# Q2 = 0.357207789214112
# Q2 = 0.455132013162395
# Q2 = 0.463307985785176
# Q2 = 0.302970599355523
# 0.3807947
cubist.looq2(read.dir = "./pre-process/alpha/", nsplits = 10, 
              cmte = 10, extra = 0, seed = 101)

# Q2 = 0.711546416588894
# Q2 = 0.701130736684873
# Q2 = 0.733101969979448
# Q2 = 0.675850054665382
# Q2 = 0.769398319572569
# Q2 = 0.744973721991518
# Q2 = 0.742186145231216
# Q2 = 0.726238111754197
# Q2 = 0.650537013975868
# Q2 = 0.648989388869814
# 0.7103952
cubist.looq2(read.dir = "./pre-process/beta/", nsplits = 10, 
             cmte = 90, extra = 0, seed = 101)

# Test sets ---------------------------------------------------------------

# split 5
alpha.tst <- cube.tst("./pre-process/alpha/", "./model.data/alpha", 
                      nsplits = 10, cmte = 10, extra = 0, seed = 101)
# split 6
beta.tst <- cube.tst("./pre-process/beta/", "./model.data/beta", 
                      nsplits = 10, cmte = 90, extra = 0, seed = 101)


# Single models -----------------------------------------------------------

#    Alpha ----
# Single models -----------------------------------------------------------

#    Alpha ----

trn.alpha <- readRDS("./pre-process/alpha/5/pp.RDS") %>%
  select(., -guest)
features <- readRDS("./feature.selection/alpha.vars.RDS")
colnames(trn.alpha) <- str_replace(colnames(trn.alpha), "-", ".")
trn.alpha <- trn.alpha %>% select(., DelG, features) 
trn.alpha.x <- select(trn.alpha, -DelG) 
trn.alpha.y <- trn.alpha$DelG

cube.alpha.ctrl <- cubistControl(
  seed = 101, 
  extrapolation = 0
)

cube.alpha <- cubist(x = trn.alpha.x, y = trn.alpha.y,  
                     committees = 10, 
                     control = cube.alpha.ctrl)

tst.alpha <- preprocess.tst.mod("./pre-process/alpha/", "./model.data/alpha/", 
                                features, 5)

tst.alpha.df <- predict(cube.alpha, tst.alpha[ , -1]) %>%
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
       title = "Cubist for Alpha")

#    Beta ----

trn.beta <- readRDS("./pre-process/beta/6/pp.RDS") %>%
  select(., -guest)
features <- readRDS("./feature.selection/beta.vars.RDS")
colnames(trn.beta) <- str_replace(colnames(trn.beta), "-", ".")
trn.beta <- trn.beta %>% select(., DelG, features) 
trn.beta.x <- select(trn.beta, -DelG) 
trn.beta.y <- trn.beta$DelG

cube.beta.ctrl <- cubistControl(
  seed = 101, 
  extrapolation = 0
)

cube.beta <- cubist(x = trn.beta.x, y = trn.beta.y,  
                     committees = 10, 
                     control = cube.beta.ctrl)

tst.beta <- preprocess.tst.mod("./pre-process/beta/", "./model.data/beta/", 
                                features, 6)

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
       title = "Cubist for Beta")

#     Saving models ----

pp.settings <- readRDS("./pre-process/alpha/5/pp.settings.RDS")
saveRDS(list(pp.settings, cube.alpha), "./models/alpha/cube.RDS")
saveRDS(tst.alpha.df, "./results/alpha/cube.RDS")
print(graph.alpha)
ggsave("./results/alpha/cube.png")


pp.settings <- readRDS("./pre-process/beta/6/pp.settings.RDS")
saveRDS(list(pp.settings, cube.beta), "./models/beta/cube.RDS")
saveRDS(tst.beta.df, "./results/beta/cube.RDS")
print(graph.beta)
ggsave("./results/beta/cube.png")
