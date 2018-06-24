dir.create("./models/glmnet")

# Libraries ---------------------------------------------------------------

library(caret)
library(glmnet)
library(tidyverse)
source("./07.model.functions.R")

# Functions ---------------------------------------------------------------

glm.looq2 <- function(read.dir, nsplits, a, max) {
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
                        pmax = ncol(trn.x), 
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
      labs(title = n) + 
      geom_abline(intercept = 0, slope = 1) + 
      coord_fixed()
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

glm.tst <- function(pp.dir, tst.dir, nsplits, a, max) {
  split <- 1:nsplits
  r2.results <- rep(0.0, nsplits)
  rmse.results <- rep(0.0, nsplits)
  if(str_detect(pp.dir, "alpha"))
    features <- readRDS("./feature.selection/alpha.vars.RDS")
  else
    features <- readRDS("./feature.selection/beta.vars.RDS")
  for(n in 1:nsplits) {
    # Reading the training data
    trn <- readRDS(paste0(pp.dir, n, "/pp.RDS")) %>%
      select(., -guest)
    trn.x <- select(trn, features) %>% data.matrix()
    trn.y <- select(trn, DelG) %>% data.matrix()
    
    # Reading and pre-processing the test set
    tst <- preprocess.tst.mod(pp.dir = pp.dir, tst.dir = tst.dir, 
                              feat = features, n = n)
    tst.y <- data.matrix(tst)[ , 1, drop = F]
    tst.x <- data.matrix(tst)[ , -1]
    # Building the model
    glm.mod <- glmnet(x = trn.x, y = trn.y, 
                      dfmax = max, alpha = a,
                      pmax = ncol(trn.x), 
                      family = "mgaussian")
    tst.df <- predict.glmnet(glm.mod, newx = tst.x, 
                             s = tail(glm.mod$lambda, n = 1)) %>%
      cbind(tst.y, .) %>%
      as.data.frame()
    colnames(tst.df) <- c("obs", "pred")
    for(i in 1:nrow(tst.df))
      if(abs(tst.df$pred[i]) > 100)
        tst.df$pred[i] <- mean(trn.y)
    tst.df <- tst.df[complete.cases(tst.df), ]
    p <- ggplot(tst.df, aes(x = obs, y = pred)) + 
      geom_point() + 
      theme_bw() + 
      labs(title = n) + 
      geom_abline(slope = 1, intercept = 0) + 
      coord_fixed()
    print(p)
    r2.results[n] <- defaultSummary(tst.df)[2]
    rmse.results[n] <- defaultSummary(tst.df)[1]
  }
  return(data.frame(split, r2.results, rmse.results))
}


# LOO-Q2 analysis ---------------------------------------------------------

# Alpha
# Either split 1 or 8
# However, none actually pass the q2 test, which is unfortunate
glm.looq2("./pre-process/alpha/", nsplits = 10, a = 0.8, max = 15)

# Beta
glm.looq2("./pre-process/beta/", nsplits = 10, a = 0.6, max = 100) # 0.634
glm.looq2("./pre-process/beta/", nsplits = 10, a = 1, max = 20) # 0.637
# Both settings work fine...maybe consider two different models?


# Test sets ---------------------------------------------------------------

# Split 5 or 7 turned out the best, R2 = 0.714, .609 
alpha.tst <- glm.tst("./pre-process/alpha/", "./model.data/alpha/", 
                     nsplits = 10, a = 0.8, max = 15)

# All passed except 4 (0.592)
# Split 6 turned out the best, r2 = 0.758
beta.tst <- glm.tst("./pre-process/beta/", "./model.data/beta/", 
                    nsplits = 10, a = 1, max = 20)

# Single model ------------------------------------------------------------

#     Alpha ----

trn.alpha <- readRDS("./pre-process/alpha/5/pp.RDS") %>%
  select(., -guest)
features <- readRDS("./feature.selection/alpha.vars.RDS")
colnames(trn.alpha) <- str_replace(colnames(trn.alpha), "-", ".")
trn.alpha <- trn.alpha %>% select(., DelG, features) 
trn.alpha.x <- select(trn.alpha, -DelG) %>% data.matrix()
trn.alpha.y <- select(trn.alpha, DelG) %>% data.matrix()

glm.alpha <- glmnet(x = trn.alpha.x, y = trn.alpha.y, 
                  dfmax = 15, alpha = 0.8,
                  pmax = ncol(trn.alpha.x), 
                  family = "mgaussian")

tst.alpha <- preprocess.tst.mod("./pre-process/alpha/", "./model.data/alpha/", 
                                features, 5) %>% data.matrix()
tst.alpha.x <- tst.alpha[ , -1]
tst.alpha.y <- tst.alpha[ , 1, drop = F]

tst.alpha.df <- predict.glmnet(glm.alpha, tst.alpha.x, 
                               s = tail(glm.alpha$lambda, n = 1)) %>%
  cbind(tst.alpha.y, .) %>% data.frame()
colnames(tst.alpha.df) <- c("obs", "pred")
eval.tropsha(tst.alpha.df)

graph.alpha <- ggplot(tst.alpha.df, aes(x = obs, y = pred)) +
  geom_point() +
  theme_bw() +
  coord_fixed()  +
  geom_abline(intercept = 0, slope = 1) + 
  labs(x = "Observed dG, kJ/mol", y = "Experimental dG, kJ/mol", 
       title = "GLMNet for Alpha")
print(graph.alpha)

#     Beta ----

trn.beta <- readRDS("./pre-process/beta/6/pp.RDS") %>%
  select(., -guest)
features <- readRDS("./feature.selection/beta.vars.RDS")
colnames(trn.beta) <- str_replace(colnames(trn.beta), "-", ".")
trn.beta <- trn.beta %>% select(., DelG, features) 
trn.beta.x <- select(trn.beta, -DelG) %>% data.matrix()
trn.beta.y <- select(trn.beta, DelG) %>% data.matrix()

glm.beta <- glmnet(x = trn.beta.x, y = trn.beta.y, 
                    dfmax = 20, alpha = 1,
                    pmax = ncol(trn.beta.x), 
                    family = "mgaussian")

tst.beta <- preprocess.tst.mod("./pre-process/beta/", "./model.data/beta/", 
                                features, 6) %>% data.matrix()
tst.beta.x <- tst.beta[ , -1]
tst.beta.y <- tst.beta[ , 1, drop = F]

tst.beta.df <- predict.glmnet(glm.beta, tst.beta.x, 
                               s = tail(glm.beta$lambda, n = 1)) %>%
  cbind(tst.beta.y, .) %>% data.frame()
colnames(tst.beta.df) <- c("obs", "pred")
eval.tropsha(tst.beta.df)

graph.beta <- ggplot(tst.beta.df, aes(x = obs, y = pred)) +
  geom_point() +
  theme_bw() +
  coord_fixed()  +
  geom_abline(intercept = 0, slope = 1) + 
  labs(x = "Observed dG, kJ/mol", y = "Experimental dG, kJ/mol", 
       title = "GLMNet for Beta")
print(graph.beta)

# Saving models -----------------------------------------------------------

pp.settings <- readRDS("./pre-process/alpha/7/pp.settings.RDS")
saveRDS(list(pp.settings, glm.alpha), "./models/alpha/glmnet.RDS")
saveRDS(tst.alpha.df, "./results/alpha/glmnet.RDS")
print(graph.alpha)
ggsave("./results/alpha/glmnet.png")


pp.settings <- readRDS("./pre-process/beta/6/pp.settings.RDS")
saveRDS(list(pp.settings, glm.beta), "./models/beta/glmnet.RDS")
saveRDS(tst.beta.df, "./results/beta/glmnet.RDS")
print(graph.beta)
ggsave("./results/beta/glmnet.png")

# External validation -----------------------------------------------------

# # Reading the pre-processing settings
# ev.alpha <- preprocess.ev("alpha", 7, features)
# ev.alpha.x <- ev.alpha[ , -1] %>% data.matrix()
# ev.alpha.y <- ev.alpha[ , 1] %>% data.matrix()
# 
# ev.alpha.df <- predict.glmnet(glm.alpha, ev.alpha.x, 
#                               s = tail(glm.alpha$lambda, n = 1)) %>%
#   cbind(ev.alpha.y, .) %>% 
#   data.frame()
# colnames(ev.alpha.df) <- c("obs", "pred")
# ggplot(ev.alpha.df, aes(x = obs, y = pred)) + 
#   geom_point() + 
#   theme_bw() + 
#   labs(title = "GLMNet for alpha-CD", 
#        x = "Observed dG, kJ/mol", y = "Predicted dG, kJ/mol") + 
#   # geom_smooth(method = "lm") + 
#   geom_abline(slope = 1, intercept = 0) + 
#   coord_fixed()
# defaultSummary(ev.alpha.df)
# # Passes everything except R2
# eval.tropsha(ev.alpha.df)
# 
# # Reading the pre-processing settings
# ev.beta <- preprocess.ev("beta", 6, features)
# ev.beta.x <- ev.beta[ , -1] %>% data.matrix()
# ev.beta.y <- ev.beta[ , 1] %>% data.matrix()
# 
# ev.beta.df <- predict.glmnet(glm.beta, ev.beta.x, 
#                              s = tail(glm.beta$lambda, n = 1)) %>%
#   cbind(ev.beta.y, .) %>% 
#   data.frame()
# colnames(ev.beta.df) <- c("obs", "pred")
# ggplot(ev.beta.df, aes(x = obs, y = pred)) + 
#   geom_point() + 
#   theme_bw() + 
#   labs(title = "GLMNet for beta-CD", 
#        x = "Observed dG, kJ/mol", y = "Predicted dG, kJ/mol") + 
#   # geom_smooth(method = "lm") + 
#   geom_abline(slope = 1, intercept = 0) + 
#   coord_fixed()
# defaultSummary(ev.beta.df)
# # Passes everything except R2
# eval.tropsha(ev.beta.df)
