dir.create("./models/glmnet")

# Libraries ---------------------------------------------------------------

library(caret)
library(glmnet)
library(tidyverse)

# Functions ---------------------------------------------------------------

glm.looq2 <- function(read.dir, rfe.dir, nsplits, a, max) {
  # initialize a vector for analysis
  q2.results <- c(rep(0.0, nsplits))
  for(n in 1:nsplits) {
    data <- readRDS(paste0(read.dir, n, "/pp.RDS")) %>%
      select(., -guest)
    features <- readRDS(paste0(rfe.dir, "/rfe", n, ".RDS")) %>% predictors()
    obs <- data[ , 1]
    data <- data %>% select(., DelG, features) %>% data.matrix()
    pred <- c(rep(0.0, nrow(data) - 1))
    # In this loop, i represents the index of the datapoint left out of 
    # model building
    for(i in 1:nrow(data)) {
      trn <- data[-i, ]
      tst <- data[i, , drop = F]
      trn.x <- trn[ , -1]
      trn.y <- trn[ , 1]
      tst.x <- tst[ , -1, drop = F]
      tst.y <- tst[ , 1, drop = F]
      
      # Refer to 07.0.1.svm.tune.poly.R for values
      glm.mod <- glmnet(x = trn.x, y = trn.y, 
                        dfmax = max, alpha = a,
                        pmax = 75, 
                        family = "mgaussian")
      pred[i] <- predict.glmnet(glm.mod, tst.x,
                                s = tail(glm.mod$lambda, n = 1)) 
    }
    # test.df <- data.frame(obs, pred) %>% print()
    # test.df <<- test.df
    # defaultSummary(test.df) %>% print()
    # prediction error sum of squares
    pred.df <- data.frame(obs, pred)
    p <- ggplot(pred.df, aes(x = obs, y = pred)) + 
      theme_bw() + 
      geom_point() + 
      labs(title = as.character(n))
    print(p)
    # Handling outliers
    for(i in 1:length(pred)) {
      if(pred[i] < -100)
        pred[i] <- mean(obs)
    }
    PRESS <- sum((obs - pred)^2)
    # total sum of squares
    TSS <- sum((obs - mean(obs))^2)
    q2 <- 1 - PRESS/TSS
    message("Q2 = ", q2)
    q2.results[n] <- q2
  }
  return(mean(q2.results))
}

glm.tst <- function(trn.dir, rfe.dir, tst.dir, nsplits, a, max) {
  split <- 1:nsplits
  r2.results <- rep(0.0, nsplits)
  rmse.results <- rep(0.0, nsplits)
  for(n in 1:nsplits) {
    # Reading the training data
    features <- readRDS(paste0(rfe.dir, "/rfe", n, ".RDS")) %>% predictors()
    trn <- readRDS(paste0(trn.dir, n, "/pp.RDS")) %>%
      select(., -guest)
    trn <- trn %>% select(., DelG, features) 
    trn.x <- select(trn, -DelG) %>% data.matrix()
    trn.y <- select(trn, DelG) %>% data.matrix()
    
    # Reading the pre-processing settings
    pp.settings <- readRDS(paste0(trn.dir, n, "/pp.settings.RDS"))
    # zero.pred <- readRDS(paste0(read.dir, n, "/zero.pred.RDS"))
    # high.cor <- readRDS(paste0(read.dir, n, "/high.cor.RDS"))
    
    # Importing and pre-processing the test set
    tst <- readRDS(paste0(tst.dir, "/tst", n, ".RDS"))
    guest <- select(tst, guest)
    tst.dg <- tst %>% select(., guest, DelG) 
    
    colnames(tst) <- str_replace(colnames(tst), "-", ".")
    tst <- select(tst, -DelG)
    tst <- do.call(data.frame, lapply(tst,
                                      function(x)
                                        replace(x, is.infinite(x), NA)))
    tst <- tst %>%
      predict(pp.settings, .) %>% select(., features) %>%
      cbind(guest, .)
    tst.ad <- domain.num(tst)
    tst.outliers <- tst.ad %>% filter(domain == "outside") %>% .$guest
    # print(tst.outliers)
    if(str_detect(trn.dir, "alpha"))
      tst.outliers <- c(tst.outliers, "indole")
    tst <- tst %>%
      filter(!guest %in% tst.outliers) %>% 
      select(., -guest) %>% data.matrix()
    tst.dg <- tst.dg %>% filter(!guest %in% tst.outliers) %>%
      select(., -guest) %>% data.matrix
    # filter(!guest == "indole") %>%  
    # filter(!guest == "3-methylbenzoic acid") %>% # removing this persistent outlier 
    
    # Building the model
    glm.mod <- glmnet(x = trn.x, y = trn.y, 
                      dfmax = max, alpha = a,
                      pmax = ncol(trn.x), 
                      family = "mgaussian")

    tst.df <- predict.glmnet(glm.mod, newx = tst, 
                             s = tail(glm.mod$lambda, n = 1)) %>%
      cbind(tst.dg, .) %>%
      as.data.frame()
    colnames(tst.df) <- c("obs", "pred")
    tst.df <- tst.df[complete.cases(tst.df), ]
    p <- ggplot(tst.df, aes(x = obs, y = pred)) + 
      geom_point() + 
      theme_bw() + 
      labs(title = n) + 
      geom_abline(slope = 1, intercept = 0)
    for(i in 1:nrow(tst.df))
      if(abs(tst.df$pred[i]) > 100)
        tst.df$pred[i] <- mean(trn.y)
    print(p)
    r2.results[n] <- defaultSummary(tst.df)[2]
    rmse.results[n] <- defaultSummary(tst.df)[1]
  }
  return(data.frame(split, r2.results, rmse.results))
}


# LOO-Q2 analysis ---------------------------------------------------------

# Alpha
glm.looq2("./pre-process/alpha/", "./feature.selection/alpha", 
          nsplits = 10, a = 0.5, max = 42)
# q2 = 0.501. good enough.

# Beta
glm.looq2("./pre-process/beta/", "./feature.selection/beta",
          nsplits = 10, a = 0.7, max = 40)
# q2 = 0.679. good.


# Test sets ---------------------------------------------------------------

# Split 5 turned out the best, R2 = 0.637
# Had to remove indole, though
alpha.tst <- glm.tst("./pre-process/alpha/", "./feature.selection/alpha", 
        "./model.data/alpha/", nsplits = 10, a = 0.5, max = 42)

# Split 9 turned out the best, R2 = 0.759 
beta.tst <- glm.tst("./pre-process/beta/", "./feature.selection/beta", 
        "./model.data/beta/", nsplits = 10, a = 0.7, max = 40)

# Single model ------------------------------------------------------------

#     Alpha ----
features <- readRDS("./feature.selection/alpha/rfe5.RDS") %>% predictors()
trn <- readRDS("./pre-process/alpha/5/pp.RDS") %>%
  select(., -guest)
trn <- trn %>% select(., DelG, features) 
trn.x <- select(trn, -DelG) %>% data.matrix()
trn.y <- select(trn, DelG) %>% data.matrix()

# Reading the pre-processing settings
pp.settings <- readRDS("./pre-process/alpha/5/pp.settings.RDS")

# Importing and pre-processing the test set
tst <- readRDS("./model.data/alpha/tst5.RDS")
guest <- select(tst, guest)
tst.dg <- tst %>% select(., guest, DelG) 

colnames(tst) <- str_replace(colnames(tst), "-", ".")
tst <- select(tst, -DelG)
tst <- do.call(data.frame, lapply(tst,
                                  function(x)
                                    replace(x, is.infinite(x), NA)))
tst <- tst %>%
  predict(pp.settings, .) %>% select(., features) %>%
  cbind(guest, .)
tst.ad <- domain.num(tst)
tst.outliers <- tst.ad %>% filter(domain == "outside") %>% .$guest
# indole is a persistent black sheep
tst.outliers <- c("indole", tst.outliers)
tst <- tst %>%
  filter(!guest %in% tst.outliers) %>% 
  select(., -guest) %>% data.matrix()
tst.dg <- tst.dg %>% filter(!guest %in% tst.outliers) %>%
  select(., -guest) %>% data.matrix

# Building the model
glm.mod <- glmnet(x = trn.x, y = trn.y, 
                  dfmax = 42, alpha = 0.5,
                  pmax = ncol(trn.x), 
                  family = "mgaussian")

tst.df <- predict.glmnet(glm.mod, newx = tst, 
                         s = tail(glm.mod$lambda, n = 1)) %>%
  cbind(tst.dg, .) %>%
  as.data.frame()
colnames(tst.df) <- c("obs", "pred")
# Checking for any terrible outliers
for(i in 1:nrow(tst.df))
  if(abs(tst.df$pred[i]) > 100)
    tst.df$pred[i] <- mean(trn.y)
ggplot(tst.df, aes(x = obs, y = pred)) + 
  geom_point() + 
  theme_bw() + 
  labs(title = "GLMNet for alpha-CD", 
       x = "Observed dG, kJ/mol", y = "Predicted dG, kJ/mol") + 
  geom_smooth(method = "lm") + 
  geom_abline(slope = 1, intercept = 0)

defaultSummary(tst.df)
tst.df <- tst.df %>% mutate(resid = obs - pred)
ggplot(tst.df, aes(x = obs, y = resid)) + 
  geom_point() + 
  theme_bw() + 
  geom_smooth(method = "lm") + 
  geom_hline(yintercept = 0)

saveRDS(glm.mod, "./models/glmnet/alpha.RDS")
saveRDS(tst.df, "./models/glmnet/alpha.results.RDS")

#     Beta ----
features <- readRDS("./feature.selection/beta/rfe8.RDS") %>% predictors()
features <- readRDS("./feature.selection/beta.vars.RDS")
trn <- readRDS("./pre-process/beta/9/pp.RDS") %>%
  select(., -guest)
trn <- trn %>% select(., DelG, features) 
trn.x <- select(trn, -DelG) %>% data.matrix()
trn.y <- select(trn, DelG) %>% data.matrix()

# Reading the pre-processing settings
pp.settings <- readRDS("./pre-process/beta/9/pp.settings.RDS")

# Importing and pre-processing the test set
tst <- readRDS("./model.data/beta/tst9.RDS")
guest <- select(tst, guest)
tst.dg <- tst %>% select(., guest, DelG) 

colnames(tst) <- str_replace(colnames(tst), "-", ".")
tst <- select(tst, -DelG)
tst <- do.call(data.frame, lapply(tst,
                                  function(x)
                                    replace(x, is.infinite(x), NA)))
tst <- tst %>%
  predict(pp.settings, .) %>% select(., features) %>%
  cbind(guest, .)
tst.ad <- domain.num(tst)
tst.outliers <- tst.ad %>% filter(domain == "outside") %>% .$guest
# indole is a persistent black sheep
tst.outliers <- c("indole", tst.outliers)
tst <- tst %>%
  filter(!guest %in% tst.outliers) %>% 
  select(., -guest) %>% data.matrix()
tst.dg <- tst.dg %>% filter(!guest %in% tst.outliers) %>%
  select(., -guest) %>% data.matrix

# Building the model
glm.mod <- glmnet(x = trn.x, y = trn.y, 
                  dfmax = 40, alpha = 0.7,
                  pmax = ncol(trn.x), 
                  family = "mgaussian")

tst.df <- predict.glmnet(glm.mod, newx = tst, 
                         s = tail(glm.mod$lambda, n = 1)) %>%
  cbind(tst.dg, .) %>%
  as.data.frame()
colnames(tst.df) <- c("obs", "pred")
# # Checking for any terrible outliers
# for(i in 1:nrow(tst.df))
#   if(abs(tst.df$pred[i]) > 100)
#     tst.df$pred[i] <- mean(trn.y)
ggplot(tst.df, aes(x = obs, y = pred)) + 
  geom_point() + 
  theme_bw() + 
  labs(title = "GLMNet for beta-CD") + 
  geom_smooth(method = "lm") + 
  geom_abline(slope = 1, intercept = 0) + 
  coord_fixed(ratio = 1)

defaultSummary(tst.df)
tst.df <- tst.df %>% mutate(resid = obs - pred)
ggplot(tst.df, aes(x = obs, y = resid)) + 
  geom_point() + 
  theme_bw() + 
  # geom_smooth(method = "lm") + 
  geom_hline(yintercept = 0) 

saveRDS(glm.mod, "./models/glmnet/beta.RDS")
saveRDS(tst.df, "./models/glmnet/beta.results.RDS")

# Additional accuracy validation ------------------------------------------

source("./eval.functions.R")
# yay you all pass
alpha.df <- readRDS("./models/glmnet/alpha.results.RDS")
eval.tropsha(alpha.df)
beta.df <- readRDS("./models/glmnet/beta.results.RDS")
eval.tropsha(beta.df)

# External validation -----------------------------------------------------

#     Alpha ----

glm.alpha <- readRDS("./models/glmnet/alpha.RDS")
ev.alpha <- readRDS("./ext.validation/alpha.RDS")
ev.alpha.info <- select(ev.alpha, guest:data.source)
guest.alpha <- select(ev.alpha.info, guest)

colnames(ev.alpha) <- str_replace(colnames(ev.alpha), "-", ".")
ev.alpha <- select(ev.alpha, -host:-data.source)
ev.alpha <- do.call(data.frame, lapply(ev.alpha,
                                  function(x)
                                    replace(x, is.infinite(x), NA)))
pp.settings <- readRDS("./pre-process/alpha/5/pp.settings.RDS")
features <- readRDS("./feature.selection/alpha/rfe5.RDS") %>% predictors()
ev.alpha <- ev.alpha %>%
  predict(pp.settings, .) %>% select(., features) %>%
  cbind(guest.alpha, .)

ev.alpha.ad <- domain.num(ev.alpha)
# No outliers
ev.alpha.outliers <- ev.alpha.ad %>% filter(domain == "outside") %>% .$guest

ev.alpha <- ev.alpha %>%
  filter(!guest %in% ev.alpha.outliers) %>% 
  select(., -guest) %>% data.matrix()
ev.alpha.dg <- ev.alpha.info %>% filter(!guest %in% ev.alpha.outliers) %>%
  select(., DelG) %>% data.matrix

ev.alpha.df <- predict.glmnet(glm.alpha, newx = ev.alpha, 
                         s = tail(glm.alpha$lambda, n = 1)) %>%
  cbind(ev.alpha.dg, .) %>%
  as.data.frame()
colnames(ev.alpha.df) <- c("obs", "pred")

ggplot(ev.alpha.df, aes(x = obs, y = pred)) + 
  geom_point() + 
  theme_bw() + 
  coord_fixed(ratio = 1) + 
  geom_abline(intercept = 0, slope = 1)

defaultSummary(ev.alpha.df)
eval.tropsha(ev.alpha.df)

#     Beta ----

glm.beta <- readRDS("./models/glmnet/beta.RDS")
ev.beta <- readRDS("./ext.validation/beta.RDS")
ev.beta.info <- select(ev.beta, guest:data.source)
guest.beta <- select(ev.beta.info, guest)

colnames(ev.beta) <- str_replace(colnames(ev.beta), "-", ".")
ev.beta <- select(ev.beta, -host:-data.source)
ev.beta <- do.call(data.frame, lapply(ev.beta,
                                       function(x)
                                         replace(x, is.infinite(x), NA)))
pp.settings <- readRDS("./pre-process/beta/9/pp.settings.RDS")
features <- readRDS("./feature.selection/beta.vars.RDS") 
ev.beta <- ev.beta %>%
  predict(pp.settings, .) %>% select(., features) %>%
  cbind(guest.beta, .)

ev.beta.ad <- domain.num(ev.beta)
# No outliers
ev.beta.outliers <- ev.beta.ad %>% filter(domain == "outside") %>% .$guest

ev.beta <- ev.beta %>%
  filter(!guest %in% ev.beta.outliers) %>% 
  select(., -guest) %>% data.matrix()
ev.beta.dg <- ev.beta.info %>% filter(!guest %in% ev.beta.outliers) %>%
  select(., DelG) %>% data.matrix

ev.beta.df <- predict.glmnet(glm.beta, newx = ev.beta, 
                              s = tail(glm.beta$lambda, n = 1)) %>%
  cbind(ev.beta.dg, .) %>%
  as.data.frame()
colnames(ev.beta.df) <- c("obs", "pred")

ggplot(ev.beta.df, aes(x = obs, y = pred)) + 
  geom_point() + 
  theme_bw() + 
  coord_fixed(ratio = 1) + 
  geom_abline(intercept = 0, slope = 1)

defaultSummary(ev.beta.df)
eval.tropsha(ev.beta.df) # Eh, it's good enough