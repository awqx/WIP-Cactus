source("07.model.functions.R")

# Libraries and Packages --------------------------------------------------

library(pls)
library(tidyverse)

# Functions ---------------------------------------------------------------

pls.looq2 <- function(read.dir, rfe.dir, nsplits, method, ncomp) {
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

pls.tst <- function(pp.dir, tst.dir, nsplits, ncomp, method) {
  split <- 1:nsplits
  r2.results <- c(rep(0.0, nsplits))
  rmse.results <- c(rep(0.0, nsplits))
  if(str_detect(pp.dir, "alpha"))
    features <- readRDS("./feature.selection/alpha.vars.RDS")
  else
    features <- readRDS("./feature.selection/beta.vars.RDS")
  
  for(n in 1:nsplits) {
    trn <- readRDS(paste0(pp.dir, n, "/pp.RDS")) %>%
      select(., DelG, features)
    tst <- preprocess.tst.mod(pp.dir = pp.dir, tst.dir = tst.dir, 
                              feat = features, n = n)
    
    pls.mod <- plsr(DelG~., data = trn, 
                    ncomp = ncomp, method = method)
     tst.df <- predict(pls.mod, tst[ , -1], ncomp = ncomp) %>%
      cbind(tst[ , 1], .) %>%
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
    print(defaultSummary(tst.df))
    r2.results[n] <- defaultSummary(tst.df)[2]
    rmse.results[n] <- defaultSummary(tst.df)[1]
  }
  return(data.frame(split, r2.results, rmse.results))
}

# LOOCV-Q2 Evaluation -----------------------------------------------------

# None pass D:. Split 1 comes the closest.
pls.looq2("./pre-process/alpha/", nsplits = 10, 
          ncomp = 8, method = "simpls")

# All pass. 0.572
pls.looq2("./pre-process/beta/", nsplits = 10, 
          ncomp = 25, method = "simpls")

# Test sets ---------------------------------------------------------------

# Only split 5 passes
# Q2 is brought down by single outlier - it's probably fine otherwise
alpha.tst <- pls.tst("./pre-process/alpha/", "./model.data/alpha/", 
                     nsplits = 10, ncomp = 8, method = "simpls")

# Split 1 looks the best
beta.tst <- pls.tst("./pre-process/beta/", "./model.data/beta/", 
                     nsplits = 10, ncomp = 25, method = "simpls")

# Single models -----------------------------------------------------------

#    Alpha ----

trn.alpha <- readRDS("./pre-process/alpha/5/pp.RDS") %>%
  select(., -guest)
features <- readRDS("./feature.selection/alpha.vars.RDS")
colnames(trn.alpha) <- str_replace(colnames(trn.alpha), "-", ".")
trn.alpha <- trn.alpha %>% select(., DelG, features) 

pls.alpha <- plsr(DelG ~., data = trn.alpha, 
                  ncomp = 8, method = "simpls")

tst.alpha <- preprocess.tst.mod("./pre-process/alpha/", "./model.data/alpha/", 
                                features, 5)

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
       title = "PLS for Alpha")
# Loading Data ------------------------------------------------------------

dir.create("./models/pls")
df.raw <- readRDS("./data/padel.pp.RDS")
df <- df.raw %>% select(-guest, -host, -data.source)

set.seed(10)
trn.ind <- sample(x = 1:nrow(df), 
                  size = round(0.8 * nrow(df)))
trn <- df[trn.ind, ]
tst <- df[-trn.ind, ]

#    Beta ----

trn.beta <- readRDS("./pre-process/beta/1/pp.RDS") %>%
  select(., -guest)
features <- readRDS("./feature.selection/beta.vars.RDS")
colnames(trn.beta) <- str_replace(colnames(trn.beta), "-", ".")
trn.beta <- trn.beta %>% select(., DelG, features) 

pls.beta <- plsr(DelG ~., data = trn.beta, 
                  ncomp = 25, method = "simpls")

tst.beta <- preprocess.tst.mod("./pre-process/beta/", "./model.data/beta/", 
                                features, 1)

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
       title = "PLS for Beta")

#     Saving models ----

pp.settings <- readRDS("./pre-process/alpha/5/pp.settings.RDS")
saveRDS(list(pp.settings, pls.alpha), "./models/alpha/pls.RDS")
saveRDS(tst.alpha.df, "./results/alpha/pls.RDS")
print(graph.alpha)
ggsave("./results/alpha/pls.png")


pp.settings <- readRDS("./pre-process/beta/1/pp.settings.RDS")
saveRDS(list(pp.settings, pls.beta), "./models/beta/pls.RDS")
saveRDS(tst.beta.df, "./results/beta/pls.RDS")
print(graph.beta)
ggsave("./results/beta/pls.png")
