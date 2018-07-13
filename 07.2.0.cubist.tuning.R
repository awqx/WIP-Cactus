dir.create("./tuning/cubist")

# Libraries ---------------------------------------------------------------

# install.packages("Cubist")
library(caret)
library(Cubist)
library(tidyverse)

# Functions ---------------------------------------------------------------

# Cubist should be fed with matrix or data.frame
tune.cubist.cmte <- function(data, nfolds, cmte, seed) {
  set.seed(seed)
  fold.list <- createFolds(y = data[ , 1], k = nfolds)
  r2.results <- c(rep(0.0, nfolds))
  rmse.results <- c(rep(0.0, nfolds))
  
  ctrl <- cubistControl(
    seed = seed
  )
  
  for(i in 1:nfolds) {
    fold <- fold.list[[i]]
    
    trn.x <- data[-fold, -1]
    trn.y <- data[-fold, 1]
    tst.x <- data[fold, -1]
    tst.y <- data[fold, 1]
    
    cube <- cubist(x = trn.x, y = trn.y, control = ctrl,
                   committees = cmte)
    cube.df <- predict(cube, tst.x) %>%
      cbind(tst.y) %>%
      data.frame() 
    
    colnames(cube.df)[1] <- "pred"
    colnames(cube.df)[2] <- "obs"
    
    rmse.results[i] <- defaultSummary(cube.df)[1]
    r2.results[i] <- defaultSummary(cube.df)[2]
   
  }
  message(nfolds, "-fold cross-validation of ", cmte, " committees completed.")
  return(data.frame( 
    nfolds = nfolds,
    seed = seed,
    committees = cmte,
    rsquared = sum(r2.results)/nfolds, 
    rmse = sum(rmse.results)/nfolds))
}

tune.cubist.samp <- function(data, nfolds, samp, seed) {
  set.seed(seed)
  fold.list <- createFolds(y = data[ , 1], k = nfolds)
  r2.results <- c(rep(0.0, nfolds))
  rmse.results <- c(rep(0.0, nfolds))
  
  ctrl <- cubistControl(
    seed = seed, 
    sample = samp
  )
  
  for(i in 1:nfolds) {
    fold <- fold.list[[i]]
    
    trn.x <- data[-fold, -1]
    trn.y <- data[-fold, 1]
    tst.x <- data[fold, -1]
    tst.y <- data[fold, 1]
    
    cube <- cubist(x = trn.x, y = trn.y, control = ctrl)
    cube.df <- predict(cube, tst.x) %>%
      cbind(tst.y) %>%
      data.frame() %>%
      .[complete.cases(.), ]
    # for some reason, tuning samp causes NAs, so those cases 
    # must be removed before caret::defaultSummary
    
    colnames(cube.df)[1] <- "pred"
    colnames(cube.df)[2] <- "obs"
    
    rmse.results[i] <- defaultSummary(cube.df)[1]
    r2.results[i] <- defaultSummary(cube.df)[2]
  }
  
  message(nfolds, "-fold cross-validation of ", samp, "% sample size completed.")
  return(data.frame(
    nfolds = nfolds,
    seed = seed,
    samp = samp,
    rsquared = sum(r2.results)/nfolds, 
    rmse = sum(rmse.results)/nfolds))
}

tune.cubist.extra <- function(data, nfolds, extra, seed) {
  set.seed(seed)
  fold.list <- createFolds(y = data[ , 1], k = nfolds)
  r2.results <- c(rep(0.0, nfolds))
  rmse.results <- c(rep(0.0, nfolds))
  
  ctrl <- cubistControl(
    seed = seed, 
    extrapolation = extra
  )
  
  for(i in 1:nfolds) {
    fold <- fold.list[[i]]
    
    trn.x <- data[-fold, -1]
    trn.y <- data[-fold, 1]
    tst.x <- data[fold, -1]
    tst.y <- data[fold, 1]
    
    cube <- cubist(x = trn.x, y = trn.y, control = ctrl)
    cube.df <- predict(cube, tst.x) %>%
      cbind(tst.y) %>%
      data.frame() 
    
    colnames(cube.df)[1] <- "pred"
    colnames(cube.df)[2] <- "obs"
    
    rmse.results[i] <- defaultSummary(cube.df)[1]
    r2.results[i] <- defaultSummary(cube.df)[2]
  }
  message(nfolds, "-fold cross-validation of ", extra, " extrapolation completed.")
  return(data.frame(
    nfolds = nfolds,
    seed = seed,
    extrapolation = extra,
    rsquared = sum(r2.results)/nfolds, 
    rmse = sum(rmse.results)/nfolds))
}

# cubist is the only tuning function with a set.seed function
# because the entire thing needs a constant reminder to set.seed
tune.cubist <- function(data, nfolds, cmte, samp, extra, seed) {
  set.seed(seed)
  fold.list <- createFolds(y = data[ , 1], k = nfolds)
  r2.results <- c(rep(0.0, nfolds))
  rmse.results <- c(rep(0.0, nfolds))
  
  ctrl <- cubistControl(
    seed = seed, 
    extrapolation = extra,
    sample = samp
  )
  
  for(i in 1:nfolds) {
    fold <- fold.list[[i]]
    
    trn.x <- data[-fold, -1]
    trn.y <- data[-fold, 1]
    tst.x <- data[fold, -1]
    tst.y <- data[fold, 1]
    
    cube <- cubist(x = trn.x, y = trn.y,  
                   committees = cmte, 
                   control = ctrl)
    cube.df <- predict(cube, tst.x) %>%
      cbind(tst.y) %>%
      data.frame() 
    
    colnames(cube.df)[1] <- "pred"
    colnames(cube.df)[2] <- "obs"
    
    rmse.results[i] <- defaultSummary(cube.df)[1]
    r2.results[i] <- defaultSummary(cube.df)[2]
  }
  message(nfolds, "-fold cross-validation of ", cmte, " committees, ",
          samp, "% sample size, ",
          extra, " extrapolation completed.")
  return(data.frame(
    nfolds = nfolds,
    seed = seed,
    committees = cmte, 
    sample = samp, 
    extrapolation = extra,
    rsquared = sum(r2.results)/nfolds, 
    rmse = sum(rmse.results)/nfolds))
}

tune.cubist <- function(data, nfolds, cmte, extra, seed) {
  set.seed(seed)
  fold.list <- createFolds(y = data[ , 1], k = nfolds)
  r2.results <- c(rep(0.0, nfolds))
  rmse.results <- c(rep(0.0, nfolds))
  
  ctrl <- cubistControl(
    seed = seed, 
    extrapolation = extra
  )
  
  for(i in 1:nfolds) {
    fold <- fold.list[[i]]
    
    trn.x <- data[-fold, -1]
    trn.y <- data[-fold, 1]
    tst.x <- data[fold, -1]
    tst.y <- data[fold, 1]
    
    cube <- cubist(x = trn.x, y = trn.y,  
                   committees = cmte, 
                   control = ctrl)
    cube.df <- predict(cube, tst.x) %>%
      cbind(tst.y) %>%
      data.frame() 
    
    colnames(cube.df)[1] <- "pred"
    colnames(cube.df)[2] <- "obs"
    
    rmse.results[i] <- defaultSummary(cube.df)[1]
    r2.results[i] <- defaultSummary(cube.df)[2]
  }
  message(nfolds, "-fold cross-validation of ", 
          cmte, " committees, ",
          extra, " extrapolation completed.")
  return(data.frame(
    nfolds = nfolds,
    seed = seed,
    committees = cmte, 
    extrapolation = extra,
    rsquared = sum(r2.results)/nfolds, 
    rmse = sum(rmse.results)/nfolds))
}


# Alpha -------------------------------------------------------------------

dir.create("./tuning/cubist/alpha")

#     Data Organization ---------------------------------------------------

trn.all <- readRDS("./pre-process/alpha/1/pp.RDS") 
colnames(trn.all) <- str_replace(colnames(trn.all), "-", ".")
trn.guest <- trn.all$guest
trn <- select(trn.all, -guest)

features <- readRDS("./feature.selection/alpha.vars.RDS")
trn <- trn[ , colnames(trn) %in% c("DelG", features)]

#     Estimation ----------------------------------------------------------

#     Committee number ---
cmte.range <- c(10*seq(1, 9, 2), 100)
results1.cmte <- do.call(rbind, lapply(cmte.range, FUN = tune.cubist.cmte, 
                                       data = trn, nfolds = 10, seed = 101))
results2.cmte <- do.call(rbind, lapply(cmte.range, FUN = tune.cubist.cmte, 
                                       data = trn, nfolds = 10, seed = 102))
results3.cmte <- do.call(rbind, lapply(cmte.range, FUN = tune.cubist.cmte, 
                                       data = trn, nfolds = 10, seed = 103))
results.cmte <- rbind(results1.cmte, results2.cmte, results3.cmte) %>%
  mutate(seed = as.factor(seed)) 
ggplot(results.cmte, aes(x = committees, y = rsquared, color = seed)) + 
  geom_line() + 
  theme_bw()
saveRDS(results.cmte, "./tuning/cubist/alpha/cmte.RDS")

#     Sample percentage ---
# Keeps returning NAs. The only values that consistently work are 5, 10, 
# 85, 90, 95
# Will keep out of tuning for now
samp.range <- c(5, 10, 85, 90, 95)
 samp.range <- c(5, 10, 25, 50, 75, 95)
results1.samp <- do.call(rbind, lapply(samp.range, FUN = tune.cubist.samp, 
                                       data = trn, nfolds = 10, seed = 101)) %>% print()
results2.samp <- do.call(rbind, lapply(samp.range, FUN = tune.cubist.samp, 
                                       data = trn, nfolds = 10, seed = 102)) %>% print()
results3.samp <- do.call(rbind, lapply(samp.range, FUN = tune.cubist.samp, 
                                       data = trn, nfolds = 10, seed = 103)) %>% print()
results.samp <- rbind(results1.samp, results2.samp, results3.samp) %>%
  mutate(seed = as.factor(seed)) 
ggplot(results.samp, aes(x = samp, y = rsquared, color = seed)) + 
  geom_line() + 
  theme_bw()
saveRDS(results.samp, "./tuning/cubist/alpha/samp.RDS")

#     Extrapolation ---
extra.range <- 10 * c(0, 1, 4, 7, 10)
results1.extra <- do.call(rbind, lapply(extra.range, FUN = tune.cubist.extra, 
                                       data = trn, nfolds = 10, seed = 101))
results2.extra <- do.call(rbind, lapply(extra.range, FUN = tune.cubist.extra, 
                                       data = trn, nfolds = 10, seed = 102))
results3.extra <- do.call(rbind, lapply(extra.range, FUN = tune.cubist.extra, 
                                       data = trn, nfolds = 10, seed = 103))
results.extra <- rbind(results1.extra, results2.extra, results3.extra) %>%
  mutate(seed = as.factor(seed)) 
ggplot(results.extra, aes(x = extrapolation, y = rsquared, color = seed)) + 
  geom_line() + 
  theme_bw()

saveRDS(results.extra, "./tuning/cubist/alpha/extra.RDS")

#     Tuning --------------------------------------------------------------

# 6*5*6 = 150 combinations
cube.combos <- expand.grid(cmte.range, 
                           #samp.range, 
                           extra.range)
colnames(cube.combos) <- c("cmte",
                           # "samp", 
                           "extra") 
cmte.combos <- cube.combos$cmte
# samp.combos <- cube.combos$samp
extra.combos <- cube.combos$extra

set.seed(1001)
# I reduce the number of folds because Cubist takes a while
system.time(
  results.combos <- do.call(
    rbind,
    mapply(
      FUN = tune.cubist,
      cmte = cmte.combos,
      extra = extra.combos,
      MoreArgs = 
        list(nfolds = 5, data = trn, seed = 1001), 
      SIMPLIFY = F
    )
  )
)

# system.time output
# user  system elapsed 
# 41.85    0.05   42.13  
# 58.77    0.06   59.00

results.combos[order(results.combos$rsquared, decreasing = T), ] %>% head()
# nfolds seed committees extrapolation  rsquared     rmse
#     5 1001         10            10 0.6118653 3.003585
#     5 1001         10            70 0.6100001 3.058802
results.combos[order(results.combos$rmse), ] %>% head()
# nfolds seed committees extrapolation  rsquared     rmse
#      5 1001         10            10 0.6118653 3.003585
#      5 1001        100             0 0.6057734 3.019532

# r2 = 0.442, rmes = 3.83
# committees = 30, extra = 10

# r2 = 0.441, rmse = 3.81
# committees = 10, extrapolation = 0

saveRDS(results.combos, "./tuning/cubist/alpha/tune.RDS")

results.combos <- results.combos %>% 
  mutate(committees = as.factor(committees), 
         extrapolation = as.factor(extrapolation))
ggplot(results.combos, aes(x = committees, y = extrapolation, fill = rsquared)) + 
  geom_raster() + 
  scale_fill_gradientn(colors = terrain.colors(20)) + 
  theme_bw() + 
  labs(x = "Number of committees", y = "Degree of extrapolation", 
       title = "Cubist tuning for alpha-CD")
ggsave("./tuning/cubist/alpha/tune.png", dpi = 450)

# Beta -------------------------------------------------------------------

dir.create("./tuning/cubist/beta")

#     Data Organization ---------------------------------------------------

trn.all <- readRDS("./pre-process/beta/1/pp.RDS") 
colnames(trn.all) <- str_replace(colnames(trn.all), "-", ".")
trn.guest <- trn.all$guest
trn <- select(trn.all, -guest)

features <- readRDS("./feature.selection/beta.vars.RDS")
trn <- trn[ , colnames(trn) %in% c("DelG", features)]

#     Estimation ----------------------------------------------------------

#     Committee number ---
cmte.range <- c(10*seq(1, 9, 2), 100)
results1.cmte <- do.call(rbind, lapply(cmte.range, FUN = tune.cubist.cmte, 
                                       data = trn, nfolds = 10, seed = 101))
results2.cmte <- do.call(rbind, lapply(cmte.range, FUN = tune.cubist.cmte, 
                                       data = trn, nfolds = 10, seed = 102))
results3.cmte <- do.call(rbind, lapply(cmte.range, FUN = tune.cubist.cmte, 
                                       data = trn, nfolds = 10, seed = 103))
results.cmte <- rbind(results1.cmte, results2.cmte, results3.cmte) %>%
  mutate(seed = as.factor(seed)) 
ggplot(results.cmte, aes(x = committees, y = rsquared, color = seed)) + 
  geom_line() + 
  theme_bw()
saveRDS(results.cmte, "./tuning/cubist/beta/cmte.RDS")

# #     Sample percentage ---
# # Keeps returning NAs. The only values that consistently work are 5, 10, 
# # 85, 90, 95
# samp.range <- c(5, 10, 85, 90, 95)
# samp.range <- c(5, 10, 25, 50, 75, 95)
# results1.samp <- do.call(rbind, lapply(samp.range, FUN = tune.cubist.samp, 
#                                        data = trn, nfolds = 10, seed = 101)) %>% print()
# results2.samp <- do.call(rbind, lapply(samp.range, FUN = tune.cubist.samp, 
#                                        data = trn, nfolds = 10, seed = 102)) %>% print()
# results3.samp <- do.call(rbind, lapply(samp.range, FUN = tune.cubist.samp, 
#                                        data = trn, nfolds = 10, seed = 103)) %>% print()
# results.samp <- rbind(results1.samp, results2.samp, results3.samp) %>%
#   mutate(seed = as.factor(seed)) 
# ggplot(results.samp, aes(x = samp, y = rsquared, color = seed)) + 
#   geom_line() + 
#   theme_bw()
# saveRDS(results.samp, "./tuning/cubist/beta/samp.RDS")

#     Extrapolation ---
extra.range <- 10 * (0:10)
results1.extra <- do.call(rbind, lapply(extra.range, FUN = tune.cubist.extra, 
                                        data = trn, nfolds = 10, seed = 101))
results2.extra <- do.call(rbind, lapply(extra.range, FUN = tune.cubist.extra, 
                                        data = trn, nfolds = 10, seed = 102))
results3.extra <- do.call(rbind, lapply(extra.range, FUN = tune.cubist.extra, 
                                        data = trn, nfolds = 10, seed = 103))
results.extra <- rbind(results1.extra, results2.extra, results3.extra) %>%
  mutate(seed = as.factor(seed)) 
ggplot(results.extra, aes(x = extrapolation, y = rsquared, color = seed)) + 
  geom_line() + 
  theme_bw()

saveRDS(results.extra, "./tuning/cubist/beta/extra.RDS")

#     Tuning --------------------------------------------------------------

# 6*5*6 = 150 combinations
cube.combos <- expand.grid(cmte.range, 
                           #samp.range, 
                           extra.range)
colnames(cube.combos) <- c("cmte",
                           # "samp", 
                           "extra") 
cmte.combos <- cube.combos$cmte
# samp.combos <- cube.combos$samp
extra.combos <- cube.combos$extra

set.seed(1001)
# I reduce the number of folds because rforest takes a while
system.time(
  results.combos <- do.call(
    rbind,
    mapply(
      FUN = tune.cubist,
      cmte = cmte.combos,
      extra = extra.combos,
      MoreArgs = 
        list(nfolds = 5, data = trn, seed = 1001), 
      SIMPLIFY = F
    )
  )
)

# system.time output
# user  system elapsed 
# 189.24    0.36  199.34

results.combos[order(results.combos$rsquared, decreasing = T), ] %>% head()
# nfolds seed committees extrapolation  rsquared     rmse
# 29      5 1001         90            40 0.6485403 3.202934
# 30      5 1001        100            40 0.6483738 3.201895
results.combos[order(results.combos$rmse), ] %>% head()
# nfolds seed committees extrapolation  rsquared     rmse
# 36      5 1001        100            50 0.6479029 3.193700
# 35      5 1001         90            50 0.6469835 3.200482

saveRDS(results.combos, "./tuning/cubist/beta/tune.RDS")

results.combos <- results.combos %>% 
  mutate(committees = as.factor(committees), 
         extrapolation = as.factor(extrapolation))
ggplot(results.combos, aes(x = committees, y = extrapolation, fill = rsquared)) + 
  geom_raster() + 
  scale_fill_gradientn(colors = terrain.colors(20)) + 
  theme_bw() + 
  labs(x = "Number of committees", y = "Degree of extrapolation", 
       title = "Cubist tuning for beta-CD")
ggsave("./tuning/cubist/beta/tune.png", dpi = 450)

# gamma -------------------------------------------------------------------

dir.create("./tuning/cubist/gamma")

#     Data Organization ---------------------------------------------------

trn.all <- readRDS("./pre-process/gamma/1/pp.RDS") 
colnames(trn.all) <- str_replace(colnames(trn.all), "-", ".")
trn.guest <- trn.all$guest
trn <- select(trn.all, -guest)

features <- readRDS("./feature.selection/gamma.vars.RDS")
trn <- trn[ , colnames(trn) %in% c("DelG", features)]

#     Estimation ----------------------------------------------------------

#     Committee number ---
cmte.range <- c(10*seq(1, 9, 2), 100)
results1.cmte <- do.call(rbind, lapply(cmte.range, FUN = tune.cubist.cmte, 
                                       data = trn, nfolds = 10, seed = 101))
results2.cmte <- do.call(rbind, lapply(cmte.range, FUN = tune.cubist.cmte, 
                                       data = trn, nfolds = 10, seed = 102))
results3.cmte <- do.call(rbind, lapply(cmte.range, FUN = tune.cubist.cmte, 
                                       data = trn, nfolds = 10, seed = 103))
results.cmte <- rbind(results1.cmte, results2.cmte, results3.cmte) %>%
  mutate(seed = as.factor(seed)) 
ggplot(results.cmte, aes(x = committees, y = rsquared, color = seed)) + 
  geom_line() + 
  theme_bw()
saveRDS(results.cmte, "./tuning/cubist/gamma/cmte.RDS")

# #     Sample percentage ---
# # Keeps returning NAs. The only values that consistently work are 5, 10, 
# # 85, 90, 95
# samp.range <- c(5, 10, 85, 90, 95)
# samp.range <- c(5, 10, 25, 50, 75, 95)
# results1.samp <- do.call(rbind, lapply(samp.range, FUN = tune.cubist.samp, 
#                                        data = trn, nfolds = 10, seed = 101)) %>% print()
# results2.samp <- do.call(rbind, lapply(samp.range, FUN = tune.cubist.samp, 
#                                        data = trn, nfolds = 10, seed = 102)) %>% print()
# results3.samp <- do.call(rbind, lapply(samp.range, FUN = tune.cubist.samp, 
#                                        data = trn, nfolds = 10, seed = 103)) %>% print()
# results.samp <- rbind(results1.samp, results2.samp, results3.samp) %>%
#   mutate(seed = as.factor(seed)) 
# ggplot(results.samp, aes(x = samp, y = rsquared, color = seed)) + 
#   geom_line() + 
#   theme_bw()
# saveRDS(results.samp, "./tuning/cubist/gamma/samp.RDS")

#     Extrapolation ---
extra.range <- 10 * (0:10)
results1.extra <- do.call(rbind, lapply(extra.range, FUN = tune.cubist.extra, 
                                        data = trn, nfolds = 10, seed = 101))
results2.extra <- do.call(rbind, lapply(extra.range, FUN = tune.cubist.extra, 
                                        data = trn, nfolds = 10, seed = 102))
results3.extra <- do.call(rbind, lapply(extra.range, FUN = tune.cubist.extra, 
                                        data = trn, nfolds = 10, seed = 103))
results.extra <- rbind(results1.extra, results2.extra, results3.extra) %>%
  mutate(seed = as.factor(seed)) 
ggplot(results.extra, aes(x = extrapolation, y = rsquared, color = seed)) + 
  geom_line() + 
  theme_bw()

saveRDS(results.extra, "./tuning/cubist/gamma/extra.RDS")

#     Tuning --------------------------------------------------------------

# 6*5*6 = 150 combinations
cube.combos <- expand.grid(cmte.range, 
                           #samp.range, 
                           extra.range)
colnames(cube.combos) <- c("cmte",
                           # "samp", 
                           "extra") 
cmte.combos <- cube.combos$cmte
# samp.combos <- cube.combos$samp
extra.combos <- cube.combos$extra

set.seed(1001)
# I reduce the number of folds because cubist takes a while
system.time(
  results.combos <- do.call(
    rbind,
    mapply(
      FUN = tune.cubist,
      cmte = cmte.combos,
      extra = extra.combos,
      MoreArgs = 
        list(nfolds = 5, data = trn, seed = 1001), 
      SIMPLIFY = F
    )
  )
)

# system.time output
# user  system elapsed 
# 189.24    0.36  199.34

results.combos[order(results.combos$rsquared, decreasing = T), ] %>% head()
# nfolds seed committees extrapolation  rsquared     rmse
# 37      5 1001         10            60 0.1911192 1.713634
# 40      5 1001         70            60 0.1854841 1.757859
results.combos[order(results.combos$rmse), ] %>% head()
# nfolds seed committees extrapolation  rsquared     rmse
# 13      5 1001         10            20 0.1129040 1.677326
# 19      5 1001         10            30 0.1566335 1.684599

saveRDS(results.combos, "./tuning/cubist/gamma/tune.RDS")

results.combos <- results.combos %>% 
  mutate(committees = as.factor(committees), 
         extrapolation = as.factor(extrapolation))
ggplot(results.combos, aes(x = committees, y = extrapolation, fill = rsquared)) + 
  geom_raster() + 
  scale_fill_gradientn(colors = terrain.colors(20)) + 
  theme_bw() + 
  labs(x = "Number of committees", y = "Degree of extrapolation", 
       title = "Cubist tuning for gamma-CD")
ggsave("./tuning/cubist/gamma/tune.png", dpi = 450)
