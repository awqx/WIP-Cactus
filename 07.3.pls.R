source("07.model.functions.R")

# Libraries and Packages --------------------------------------------------

library(pls)
library(tidyverse)

# Functions ---------------------------------------------------------------

pls.looq2 <- function(read.dir, rfe.dir, nsplits, method, ncomp) {
  trn.split <- 1:nsplits
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
    pred <- c()
    # In this loop, i represents the index of the datapoint left out of 
    # model building
    for(i in 1:nrow(data)) {
      trn <- data[-i, ]
      tst <- data[i, ]
      # trn.x <- trn[ , -1]
      # trn.y <- trn[ , 1]
      # tst.x <- tst[ , -1]
      # tst.y <- tst[ , 1]
      
      pls.mod <- plsr(DelG~., data = trn, 
                      ncomp = ncomp, method = method)
      pred[i] <- predict(pls.mod, tst[ , -1]) %>% .[2]
    }
    # Handling outliers
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

pls.tst <- function(pp.dir, tst.dir, nsplits, ncomp, method) {
  split <- 1:nsplits
  results.all <- data.frame()
  if(str_detect(pp.dir, "alpha"))
    features <- readRDS("./feature.selection/alpha.vars.RDS")
  else
    features <- readRDS("./feature.selection/beta.vars.RDS")
  pls.mod <- NULL
  for(i in 1:nsplits) {
    trn <- readRDS(paste0(pp.dir, i, "/pp.RDS")) %>%
      select(., DelG, one_of(features))
    pls.mod <- plsr(DelG~., data = trn, 
                    ncomp = ncomp, method = method)
    tst.all <- data.frame()
    
    for(j in 1:nsplits) {
      tst <- preprocess.tst.mod(pp.dir = pp.dir, tst.dir = tst.dir, 
                                feat = features, n = j)
      tst.y <- tst[ , 1]
      tst.x <- tst[ , -1]
      
      tst.df <- predict(pls.mod, tst[ , -1], ncomp = ncomp) %>%
        cbind(tst[ , 1], .) %>%
        as.data.frame()
      colnames(tst.df) <- c("obs", "pred")
      tst.df <- tst.df %>% mutate(split = j)
      # for(k in 1:nrow(tst.df))
      #   if(abs(tst.df$pred[k]) > 80)
      #     tst.df$pred[k] <- mean(trn.y)
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

# LOOCV-Q2 Evaluation -----------------------------------------------------

# 1, 4, 5, 6, 8, 9, 10
pls.looq2("./pre-process/alpha/", nsplits = 10, 
          ncomp = 8, method = "oscorespls")

# All pass. 0.526
pls.looq2("./pre-process/beta/", nsplits = 10, 
          ncomp = 16, method = "oscorespls")

# None
pls.looq2("./pre-process/gamma/", nsplits = 10, 
          ncomp = 2, method = "oscorespls")

# Test sets ---------------------------------------------------------------

# ALPHA - NA
alpha.tst <- pls.tst("./pre-process/alpha/", "./model.data/alpha/", 
                     nsplits = 10, ncomp = 8, method = "oscorespls")
alpha.1to1 <- alpha.tst %>% filter(trn.split == tst.split)
alpha.avg <- avg.tst(alpha.tst) %>% print()

# BETA - NA
beta.tst <- pls.tst("./pre-process/beta/", "./model.data/beta/", 
                    nsplits = 10, ncomp = 16, method = "oscorespls")
beta.1to1 <- beta.tst %>% filter(trn.split == tst.split)
beta.avg <- avg.tst(beta.tst) %>% print()

# GAMMA - NA
gamma.tst <- pls.tst("./pre-process/gamma/", "./model.data/gamma/", 
                    nsplits = 10, ncomp = 2, method = "oscorespls")
gamma.1to1 <- gamma.tst %>% filter(trn.split == tst.split)
gamma.avg <- avg.tst(gamma.tst) %>% print()

# Single models -----------------------------------------------------------

#    Alpha ----

trn.alpha <- readRDS("./pre-process/alpha/3/pp.RDS") %>%
  select(., -guest)
features <- readRDS("./feature.selection/alpha.vars.RDS")
colnames(trn.alpha) <- str_replace(colnames(trn.alpha), "-", ".")
trn.alpha <- trn.alpha %>% select(., DelG, features) 

pls.alpha <- plsr(DelG ~., data = trn.alpha, 
                  ncomp = 8, method = "oscorespls")

tst.alpha <- preprocess.tst.mod("./pre-process/alpha/", "./model.data/alpha/", 
                                features, 3)

tst.alpha.df <- predict(pls.alpha, tst.alpha[ , -1], ncomp = 8) %>%
  cbind(tst.alpha[ , 1], .) %>% data.frame()
colnames(tst.alpha.df) <- c("obs", "pred")

eval.tropsha(tst.alpha.df)
graph.alpha <- ggplot(tst.alpha.df, aes(x = obs, y = pred)) + 
  geom_point() + 
  theme_bw() + 
  coord_fixed()  + 
  geom_abline(intercept = 0, slope = 1) + 
  labs(x = "Observed dG, kJ/mol", y = "Predicted dG, kJ/mol", 
       title = "Alpha-CD PLS")

#    Beta ----

trn.beta <- readRDS("./pre-process/beta/3/pp.RDS") %>%
  select(., -guest)
features <- readRDS("./feature.selection/beta.vars.RDS")
colnames(trn.beta) <- str_replace(colnames(trn.beta), "-", ".")
trn.beta <- trn.beta %>% select(., DelG, features) 

pls.beta <- plsr(DelG ~., data = trn.beta, 
                 ncomp = 25, method = "oscorespls")

tst.beta <- preprocess.tst.mod("./pre-process/beta/", "./model.data/beta/", 
                               features, 3)

tst.beta.df <- predict(pls.beta, tst.beta[ , -1], ncomp = 25) %>%
  cbind(tst.beta[ , 1], .) %>% data.frame()
colnames(tst.beta.df) <- c("obs", "pred")

# Yay you pass
eval.tropsha(tst.beta.df)
graph.beta <- ggplot(tst.beta.df, aes(x = obs, y = pred)) + 
  geom_point() + 
  theme_bw() + 
  coord_fixed()  + 
  geom_abline(intercept = 0, slope = 1) + 
  labs(x = "Observed dG, kJ/mol", y = "Predicted dG, kJ/mol", 
       title = "Beta-CD PLS")

#     Saving models ----

pp.settings <- readRDS("./pre-process/alpha/3/pp.settings.RDS")
saveRDS(list(pp.settings, pls.alpha), "./models/alpha/pls.RDS")
saveRDS(tst.alpha.df, "./results/alpha/pls.RDS")
print(graph.alpha)
ggsave("./results/alpha/pls.png")


pp.settings <- readRDS("./pre-process/beta/3/pp.settings.RDS")
saveRDS(list(pp.settings, pls.beta), "./models/beta/pls.RDS")
saveRDS(tst.beta.df, "./results/beta/pls.RDS")
print(graph.beta)
ggsave("./results/beta/pls.png")
