# Libraries and Packages --------------------------------------------------

# install.packages("pls")
library(pls)
library(tidyverse)

# Functions ---------------------------------------------------------------

tune.pls.ncomp <- function(data, nfolds, comp, seed) {
  set.seed(seed)
  fold.list <- createFolds(y = data[ , 1], k = nfolds)
  r2.results <- c(rep(0.0, nfolds))
  rmse.results <- c(rep(0.0, nfolds))
  
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
    
    r2.results[i] <- defaultSummary(pls.df)[2]
    rmse.results[i] <- defaultSummary(pls.df)[1]
  }
  # message(comp)
  
  return(data.frame( # Useful for records
    data = deparse(substitute(data)), # Turns the variable name into char
    nfolds = nfolds,
    ncomp = comp,
    rsquared = sum(r2.results) / nfolds, 
    rmse = sum(rmse.results)/nfolds,
    seed = seed))
}

# Method = "kernelpls", "widekernelpls", "simpls", or "oscorespls"
tune.pls.method <- function(data, nfolds, met, seed) {
  set.seed(seed)
  fold.list <- createFolds(y = data[ , 1], k = nfolds)
  r2.results <- c(rep(0.0, nfolds))
  rmse.results <- c(rep(0.0), nfolds)
  
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
    
    rmse.results[i] <- defaultSummary(pls.df)[1]
    r2.results[i] <- defaultSummary(pls.df)[2]
  }
  # I have a message for this particular function in order to keep better
  # track, and to not stop the function too early
  message(met, " has been completed")
  return(data.frame( # Useful for records
    data = deparse(substitute(data)), # Turns the variable name into char
    nfolds = nfolds,
    method = met,
    rsquared = sum(r2.results)/nfolds, 
    rmse = sum(rmse.results)/nfolds,
    seed = seed)) 
}

tune.pls <- function(data, nfolds, comp, met) {
  fold.list <- createFolds(y = data[ , 1], k = nfolds)
  r2.results <- c(rep(0.0, nfolds))
  rmse.results <- c(rep(0.0, nfolds))
  
  for(i in 1:nfolds) {
    fold <- fold.list[[i]]
    
    trn <- data[-fold, ]
    tst <- data[fold, ]
    
    pls <- plsr(DelG ~., data = trn, 
                ncomp = comp, method = met)
    pls.df <- predict(pls, tst[ , -1]) %>%
      cbind(tst[ , 1]) %>%
      data.frame() 
    
    colnames(pls.df)[1] <- "pred"
    colnames(pls.df)[2] <- "obs"
    
    r2.results[i] <- defaultSummary(pls.df)[2]
    rmse.results[i] <- defaultSummary(pls.df)[1]
  }
  
  return(data.frame( 
    nfolds = nfolds,
    ncomp = comp,
    method = met,
    rsquared = sum(r2.results) / nfolds, 
    rmse = sum(rmse.results)/nfolds
    ))
}

# Alpha -------------------------------------------------------------------

#     Loading Data --------------------------------------------------------

dir.create("./tuning/pls")
dir.create("./tuning/pls/alpha")

# Reading data with all descriptors
trn.all <- readRDS("./pre-process/alpha/1/pp.RDS") 
colnames(trn.all) <- str_replace(colnames(trn.all), "-", ".")
trn.guest <- trn.all$guest
trn <- select(trn.all, -guest)

features <- readRDS("./feature.selection/alpha.vars.RDS")
trn <- trn[ , colnames(trn) %in% c("DelG", features)]

#     Method --------------------------------------------------------------

method.range <- c("kernelpls", "simpls", "oscorespls", "widekernelpls")
results1.methods <- do.call(rbind, lapply(method.range, FUN = tune.pls.method, 
                                          data = trn, nfolds = 10, seed = 101))
results2.methods <- do.call(rbind, lapply(method.range, FUN = tune.pls.method, 
                                          data = trn, nfolds = 10, seed = 102))
results.methods <- rbind(results1.methods, results2.methods)
saveRDS(results.methods, "./tuning/pls/alpha/methods.RDS")

#     Number of components ------------------------------------------------

ncomp.range <- c(1:5, 8, 13, 20)
results1.ncomp <- do.call(rbind, lapply(ncomp.range, FUN = tune.pls.ncomp, 
                                          data = trn, nfolds = 10, seed = 101))
results2.ncomp <- do.call(rbind, lapply(ncomp.range, FUN = tune.pls.ncomp, 
                                        data = trn, nfolds = 10, seed = 102))
results3.ncomp <- do.call(rbind, lapply(ncomp.range, FUN = tune.pls.ncomp, 
                                        data = trn, nfolds = 10, seed = 103))
results.ncomp <- rbind(results1.ncomp, results2.ncomp, results3.ncomp)
saveRDS(results.ncomp, "./tuning/pls/alpha/ncomp.RDS")  

#     Tuning --------------------------------------------------------------

pls.combos <- expand.grid(as.character(method.range), ncomp.range)
colnames(pls.combos) <- c("method", "ncomp")
method.combos <- pls.combos$method %>% as.character()
ncomp.combos <- pls.combos$ncomp

set.seed(1001)
system.time(
  results.combos <- do.call(
    rbind,
    mapply(
      FUN = tune.pls,
      comp = ncomp.combos,
      met = method.combos,
      MoreArgs = 
        list(nfolds = 10, data = trn), 
      SIMPLIFY = F
    )
  )
)

# system.time
# user  system elapsed 
# 4.77    0.02    4.84 

results.combos[order(results.combos$rsquared, decreasing = T), ] %>% head()
results.combos[order(results.combos$rmse), ] %>% head()

# r2 = 0.507, rmse = 3.46
# ncomp = 8, method = "simpls"

# r2 = 0.502, rmse = 3.42
# ncomp = 8. methpd = oscorelpls

saveRDS(results.combos, "./tuning/pls/alpha/tune.RDS")
results.combos <- results.combos %>%
  mutate(ncomp = as.factor(ncomp), method = as.factor(method))
ggplot(results.combos, aes(x = ncomp, y = method, fill = rsquared)) + 
  geom_raster() + 
  scale_fill_gradientn(colours = terrain.colors(20)) + 
  theme_bw() + 
  labs(title = "PLS Tuning for Alpha", x = "Number of Components", 
       y = "PLS Method", 
       fill = "R2")
ggsave("./tuning/pls/alpha/tune.png", dpi = 450)

# Beta -------------------------------------------------------------------

#     Loading Data --------------------------------------------------------

dir.create("./tuning/pls/beta")

# Reading data with all descriptors
trn.all <- readRDS("./pre-process/beta/1/pp.RDS") 
colnames(trn.all) <- str_replace(colnames(trn.all), "-", ".")
trn.guest <- trn.all$guest
trn <- select(trn.all, -guest)

features <- readRDS("./feature.selection/beta.vars.RDS")
trn <- trn[ , colnames(trn) %in% c("DelG", features)]

#     Method --------------------------------------------------------------

method.range <- c("kernelpls", "simpls", "oscorespls", "widekernelpls")
results1.methods <- do.call(rbind, lapply(method.range, FUN = tune.pls.method, 
                                          data = trn, nfolds = 10, seed = 101))
results2.methods <- do.call(rbind, lapply(method.range, FUN = tune.pls.method, 
                                          data = trn, nfolds = 10, seed = 102))
results.methods <- rbind(results1.methods, results2.methods)
saveRDS(results.methods, "./tuning/pls/beta/methods.RDS")

#     Number of components ------------------------------------------------

ncomp.range <- c(1:5, 10, 15, 25)
results1.ncomp <- do.call(rbind, lapply(ncomp.range, FUN = tune.pls.ncomp, 
                                        data = trn, nfolds = 10, seed = 101))
results2.ncomp <- do.call(rbind, lapply(ncomp.range, FUN = tune.pls.ncomp, 
                                        data = trn, nfolds = 10, seed = 102))
results3.ncomp <- do.call(rbind, lapply(ncomp.range, FUN = tune.pls.ncomp, 
                                        data = trn, nfolds = 10, seed = 103))
results.ncomp <- rbind(results1.ncomp, results2.ncomp, results3.ncomp)
saveRDS(results.ncomp, "./tuning/pls/beta/ncomp.RDS")  

#     Tuning --------------------------------------------------------------

pls.combos <- expand.grid(as.character(method.range), ncomp.range)
colnames(pls.combos) <- c("method", "ncomp")
method.combos <- pls.combos$method %>% as.character()
ncomp.combos <- pls.combos$ncomp

set.seed(1001)
system.time(
  results.combos <- do.call(
    rbind,
    mapply(
      FUN = tune.pls,
      comp = ncomp.combos,
      met = method.combos,
      MoreArgs = 
        list(nfolds = 10, data = trn), 
      SIMPLIFY = F
    )
  )
)

# system.time
#  user  system elapsed 
# 15.59    0.05   15.83 

results.combos[order(results.combos$rsquared, decreasing = T), ] %>% head()
results.combos[order(results.combos$rmse), ] %>% head()

# r2 = 0.644, rmse = 3.21
# ncomp = 25, method = simpls


saveRDS(results.combos, "./tuning/pls/beta/tune.RDS")
results.combos <- results.combos %>%
  mutate(ncomp = as.factor(ncomp), method = as.factor(method))
ggplot(results.combos, aes(x = ncomp, y = method, fill = rsquared)) + 
  geom_raster() + 
  scale_fill_gradientn(colours = terrain.colors(20)) + 
  theme_bw() + 
  labs(title = "GLMNet tuning for beta-CD", 
       x = "Number of components", y = "Maximum degrees of freedom", 
       fill = "R2")
ggsave("./tuning/pls/beta/tune.png", dpi = 450)

# General peak around 25 number of components. May be able to push higher. 