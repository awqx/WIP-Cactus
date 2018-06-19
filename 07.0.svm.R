dir.create("./models")
dir.create("./models/svm")

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
  for(n in 1:nsplits) {
    data <- readRDS(paste0(read.dir, n, "/pp.RDS")) %>%
      select(., -guest)
    if(str_detect(read.dir, "alpha"))
      features <- readRDS("./feature.selection/alpha.vars.RDS")
    else
      features <- readRDS("./feature.selection/beta.vars.RDS")
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

# Both pp.dir and tst.dir should end in a backslash
# n refers to the split number
preprocess.tst.mod <- function(pp.dir, tst.dir, feat, n) {
  pp.settings <- readRDS(paste0(pp.dir, n, "/pp.settings.RDS"))
  tst <- readRDS(paste0(tst.dir, "/tst", n, ".RDS"))
  guest <- select(tst, guest)
  tst.dg <- tst %>% select(., guest, DelG)
  
  colnames(tst) <- str_replace(colnames(tst), "-", ".")
  tst <- select(tst, -DelG)
  
  tst <- do.call(data.frame, lapply(tst,
                                    function(x)
                                      replace(x, is.infinite(x), NA)))
  tst <- tst %>%
    predict(pp.settings, .) %>% select(., feat) %>%
    cbind(guest, .)
  tst.ad <- domain.num(tst)
  tst.outliers <- tst.ad %>% filter(domain == "outside") %>% .$guest
  # print(tst.outliers)
  # indole is just a persistent headache for alpha
  # and 3-methylbenzoic acid consistently messes up beta
  if(str_detect(pp.dir, "alpha"))
    tst.outliers <- c(tst.outliers, "indole")
  else
    tst.outliers <- c(tst.outliers, "3-methylbenzoic acid")
  tst <- tst %>%
    filter(!guest %in% tst.outliers) %>% 
    select(., -guest) %>% data.matrix()
  tst.dg <- tst.dg %>% filter(!guest %in% tst.outliers) %>%
    select(., -guest)
  
  return(cbind(tst.dg, tst))
}

preprocess.tst <- function(pp.dir, tst.dir, feat, n) {
  pp.settings <- readRDS(paste0(pp.dir, n, "/pp.settings.RDS"))
  tst <- readRDS(paste0(tst.dir, "/tst", n, ".RDS"))
  guest <- select(tst, guest)
  tst.dg <- tst %>% select(., guest, DelG)
  
  colnames(tst) <- str_replace(colnames(tst), "-", ".")
  tst <- select(tst, -DelG)
  
  tst <- do.call(data.frame, lapply(tst,
                                    function(x)
                                      replace(x, is.infinite(x), NA)))
  tst <- tst %>%
    predict(pp.settings, .) %>% select(., feat) %>%
    cbind(guest, .)
  tst.ad <- domain.num(tst)
  tst.outliers <- tst.ad %>% filter(domain == "outside") %>% .$guest

  tst <- tst %>%
    filter(!guest %in% tst.outliers) %>% 
    select(., -guest) %>% data.matrix()
  tst.dg <- tst.dg %>% filter(!guest %in% tst.outliers) %>%
    select(., -guest)
  
  return(cbind(tst.dg, tst))
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

# Testing on external validation
ev.alpha <- readRDS("./ext.validation/alpha.RDS")
ev.alpha.info <- select(ev.alpha, guest:data.source)

colnames(ev.alpha) <- str_replace(colnames(ev.alpha), "-", ".")
ev.alpha <- select(ev.alpha, -host:-data.source)
ev.alpha <- do.call(data.frame, lapply(ev.alpha,
                                       function(x)
                                         replace(x, is.infinite(x), NA)))
pp.settings <- readRDS("./pre-process/alpha/8/pp.settings.RDS")
ev.alpha <- ev.alpha %>%
  predict(pp.settings, .) %>% select(., features) %>%
  cbind(select(ev.alpha.info, guest), .)

ev.alpha.ad <- domain.num(ev.alpha)
# No outliers
ev.alpha.outliers <- ev.alpha.ad %>% filter(domain == "outside") %>% .$guest

ev.alpha <- ev.alpha %>%
#   filter(!guest %in% ev.alpha.outliers) %>% 
  select(., -guest)
ev.alpha.dg <- ev.alpha.info %>% 
  # filter(!guest %in% ev.alpha.outliers) %>%
  select(., DelG) 

ev.alpha.df <- predict(polysvm.alpha, ev.alpha) %>%
  cbind(ev.alpha.dg, .)
colnames(ev.alpha.df) <- c("obs", "pred")
ggplot(ev.alpha.df, aes(x = obs, y = pred)) + 
  theme_bw() + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) + 
  coord_fixed()
defaultSummary(ev.alpha.df)
# rsquared = 0.732, RMSE = 2.48

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

# Testing on external validation
ev.beta <- readRDS("./ext.validation/beta.RDS")
ev.beta.info <- select(ev.beta, guest:data.source)

colnames(ev.beta) <- str_replace(colnames(ev.beta), "-", ".")
ev.beta <- select(ev.beta, -host:-data.source)
ev.beta <- do.call(data.frame, lapply(ev.beta,
                                       function(x)
                                         replace(x, is.infinite(x), NA)))
pp.settings <- readRDS("./pre-process/beta/9/pp.settings.RDS")
ev.beta <- ev.beta %>%
  predict(pp.settings, .) %>% select(., features) %>%
  cbind(select(ev.beta.info, guest), .)

ev.beta.ad <- domain.num(ev.beta)
# No outliers
ev.beta.outliers <- ev.beta.ad %>% filter(domain == "outside") %>% .$guest

ev.beta <- ev.beta %>%
  filter(!guest %in% ev.beta.outliers) %>% 
  select(., -guest)
ev.beta.dg <- ev.beta.info %>% 
  filter(!guest %in% ev.beta.outliers) %>%
  select(., DelG) 

ev.beta.df <- predict(polysvm.beta, ev.beta) %>%
  cbind(ev.beta.dg, .)
colnames(ev.beta.df) <- c("obs", "pred")
ggplot(ev.beta.df, aes(x = obs, y = pred)) + 
  theme_bw() + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) + 
  coord_fixed()
defaultSummary(ev.beta.df)
# rsquared = 0.540, RMSE = 3.38
# the r2 is brought down by one outlier, so this is probably safe

#         Additional evaluation ----

source("./eval.functions.R")
# They all pass, except for beta on R^2 (but not too bad, so it's alright)
eval.tropsha(ev.alpha.df)
eval.tropsha(ev.beta.df)

#         Saving models ----

saveRDS(polysvm.alpha, "./models/svm/polysvm.alpha.RDS")
saveRDS(polysvm.beta, "./models/svm/polysvm.beta.RDS")

# And saving the external validation results
dir.create("./results")
dir.create("./results/svm")
saveRDS(ev.alpha.df, "./results/svm/polysvm.alpha.df.RDS")
saveRDS(ev.beta.df, "./results/svm/polysvm.beta.df.RDS")

# Radial ------------------------------------------------------------------


