# Libraries and Packages --------------------------------------------------

library(pls)
library(tidyverse)

# Loading Data ------------------------------------------------------------

setwd("~/SREP LAB/qsar")
dir.create("./tuning/pls")
source("./model.code/tuning.functions.R")

df.raw <- readRDS("./padel.pp.new.RDS")
df <- df.raw %>% select(-guest:-host) %>%
  select(-data.source)
mat <- sparse.model.matrix(~., df)

# Functions ---------------------------------------------------------------

tune.ncomp <- function(data, nfolds, comp, seed) {
  set.seed(seed)
  fold.list <- createFolds(y = data[ , 1], k = nfolds)
  results <- c(rep(0.0, nfolds))
  
  for(i in 1:nfolds) {
    fold <- fold.list[[i]]
    
    trn <- data[-fold, ]
    tst <- data[fold, ]
    
    pls <- plsr(DelG ~., data = trn, 
                ncomp = comp)
    pls.df <- predict(pls, tst[ , -1]) %>%
      cbind(tst[ , 1]) %>%
      data.frame() 
    
    colnames(pls.df)[1] <- "pred"
    colnames(pls.df)[2] <- "obs"
    
    results[i] <- defaultSummary(pls.df)[2]
  }
  
  return(data.frame( # Useful for records
    data = deparse(substitute(data)), # Turns the variable name into char
    nfolds = nfolds,
    ncomp = comp,
    rsquared = sum(results) / nfolds, 
    seed = seed))
}

# Method = "kernelpls", "widekernelpls", "simpls", or "oscorespls"
tune.method <- function(data, nfolds, met, seed) {
  set.seed(seed)
  fold.list <- createFolds(y = data[ , 1], k = nfolds)
  results <- c(rep(0.0, nfolds))
  
  for(i in 1:nfolds) {
    fold <- fold.list[[i]]
    
    trn <- data[-fold, ]
    tst <- data[fold, ]
    
    pls <- plsr(DelG ~., data = trn, 
                method = met #, ncomp = 11
                )
    pls.df <- predict(pls, tst[ , -1]) %>%
      cbind(tst[ , 1]) %>%
      data.frame() 
    
    colnames(pls.df)[1] <- "pred"
    colnames(pls.df)[2] <- "obs"
    
    results[i] <- defaultSummary(pls.df)[2]
    message("Fold ", i, " of ", met, " has been completed")
  }
  
  return(data.frame( # Useful for records
    data = deparse(substitute(data)), # Turns the variable name into char
    nfolds = nfolds,
    method = met,
    rsquared = sum(results) / nfolds, 
    seed = seed)) 
}


# Tuning ------------------------------------------------------------------

# "widekernelpls" seems to stall, so I'll remove it
method.options <- c("kernelpls", "simpls", "oscorespls")

#     Seed 1 --------------------------------------------------------------

ncomp <- do.call(rbind, lapply(
  seq(1, 29, 2), FUN = tune.ncomp,
  data = df, nfolds = 10, seed = 1
))

method.results <- do.call(rbind,
                          lapply(
                            method.options, FUN = tune.method,
                            data = df, nfolds = 10, seed = 1
                          ))

#    Seed 2 ---------------------------------------------------------------

ncomp2 <- do.call(rbind, lapply(
  seq(1, 29, 2), FUN = tune.ncomp,
  data = df, nfolds = 10, seed = 2
))

method.results2 <- do.call(rbind,
                          lapply(
                            method.options, FUN = tune.method,
                            data = df, nfolds = 10, seed = 2
                          ))

#     Seed 3 --------------------------------------------------------------

ncomp3 <- do.call(rbind, lapply(
  seq(1, 29, 2), FUN = tune.ncomp,
  data = df, nfolds = 10, seed = 3
))

method.results3 <- do.call(rbind,
                          lapply(
                            method.options, FUN = tune.method,
                            data = df, nfolds = 10, seed = 3
                          ))

#     Compilation ---------------------------------------------------------

ncomp.comp <- rbind(ncomp, ncomp2, ncomp3)
method.comp <- rbind(method.results, method.results2, method.results3)

saveRDS(ncomp.comp, "./tuning/pls/ncomp tuning.RDS")
saveRDS(method.comp, "./tuning/pls/method tuning.RDS")

# Graphs ------------------------------------------------------------------
dir.create("./tuning/pls")
temp <- ncomp.comp
temp$seed <- as.factor(temp$seed)
ggplot(temp, aes(x = ncomp, y = rsquared, group = seed, color = seed)) + 
  geom_line() + 
  theme_bw() + 
  labs(x = "Number of Components", y = "R-squared", 
       title = "PLS Tuning - Components", color = "Random Seed")
ggsave("./tuning/pls/2017-07-25 pls ncomp tune.png")

temp <- method.comp
temp$seed <- as.factor(temp$seed)
ggplot(temp, aes(x = method, y = rsquared, group = seed, fill = seed)) + 
  geom_col(position = "dodge") + 
  theme_bw() + 
  labs(x = "Method", y = "R-squared", 
       title = "PLS Tuning - Methods", color = "Random Seed")
ggsave("./tuning/pls/2017-07-25 pls method tune.png")
