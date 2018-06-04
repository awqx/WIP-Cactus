
# Libraries ---------------------------------------------------------------

library(Cubist)
library(tidyverse)

# Data Organization -------------------------------------------------------

rpt <- readRDS("./data/padel.pp.RDS")
mat.dg <- rpt %>% 
  dplyr::select(., -guest:-host, -data.source)

set.seed(25)
trn.ind <- sample(x = 1:nrow(mat.dg), size = round(0.7 * nrow(mat.dg)))
trn <- mat.dg[trn.ind, ]
tst <- mat.dg[-trn.ind, ]

dir.create("./tuning/cubist")

# Functions ---------------------------------------------------------------

# Cubist should be fed with matrix or data.frame
tune.cubist.cmte <- function(data, nfolds, cmte, seed) {
  set.seed(seed)
  fold.list <- createFolds(y = data[ , 1], k = nfolds)
  results <- c(rep(0.0, nfolds))
  
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
    
    results[i] <- defaultSummary(cube.df)[2]
    message("Fold ", i, " of ", cmte, " committees completed.")
  }
  
  return(data.frame( # Useful for records
    data = deparse(substitute(data)), # Turns the variable name into char
    nfolds = nfolds,
    seed = seed,
    committees = cmte,
    rsquared = sum(results) / nfolds))
}

tune.cubist.samp <- function(data, nfolds, samp, seed) {
  set.seed(seed)
  fold.list <- createFolds(y = data[ , 1], k = nfolds)
  results <- c(rep(0.0, nfolds))
  
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
      data.frame() 
    
    colnames(cube.df)[1] <- "pred"
    colnames(cube.df)[2] <- "obs"
    
    results[i] <- defaultSummary(cube.df)[2]
    message("Fold ", i, " of ", samp, " % sample size completed.")
  }
  
  return(data.frame( # Useful for records
    data = deparse(substitute(data)), # Turns the variable name into char
    nfolds = nfolds,
    seed = seed,
    sample.perc = samp,
    rsquared = sum(results) / nfolds))
}

tune.cubist.extra <- function(data, nfolds, extra, seed) {
  set.seed(seed)
  fold.list <- createFolds(y = data[ , 1], k = nfolds)
  results <- c(rep(0.0, nfolds))
  
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
    
    results[i] <- defaultSummary(cube.df)[2]
    message("Fold ", i, " of ", extra, " extrapolation completed.")
  }
  
  return(data.frame( # Useful for records
    data = deparse(substitute(data)), # Turns the variable name into char
    nfolds = nfolds,
    seed = seed,
    extrapolation = extra, 
    rsquared = sum(results) / nfolds))
}


# Committees --------------------------------------------------------------

#     Seed 1 --------------------------------------------------------------

cube.cmte <- do.call(
  rbind,
  lapply(
    10 * c(seq(1, 9, 2), 10),
    FUN = tune.cubist.cmte,
    data = trn,
    nfolds = 10, 
    seed = 1
  )
)

cube.samp <- do.call(
  rbind,
  lapply(
    10 * c(seq(1, 9, 2), 9.9),
    FUN = tune.cubist.samp,
    data = trn,
    nfolds = 10, 
    seed = 1
  )
)

cube.extra <- do.call(
  rbind,
  lapply(
    10 * c(0, seq(1, 9, 2), 10),
    FUN = tune.cubist.extra,
    data = trn,
    nfolds = 10, 
    seed = 1
  )
)

#     Seed 2 --------------------------------------------------------------

cube.cmte2 <- do.call(
  rbind,
  lapply(
    10 * c(seq(1, 9, 2), 10),
    FUN = tune.cubist.cmte,
    data = trn,
    nfolds = 10, 
    seed = 2
  )
)

cube.samp2 <- do.call(
  rbind,
  lapply(
    10 * c(seq(1, 9, 2), 9.9),
    FUN = tune.cubist.samp,
    data = trn,
    nfolds = 10, 
    seed = 2
  )
)

cube.extra2 <- do.call(
  rbind,
  lapply(
    10 * c(0, seq(1, 9, 2), 10),
    FUN = tune.cubist.extra,
    data = trn,
    nfolds = 10, 
    seed = 2
  )
)

#     Seed 3 --------------------------------------------------------------

cube.cmte3 <- do.call(
  rbind,
  lapply(
    10 * c(seq(1, 9, 2), 10),
    FUN = tune.cubist.cmte,
    data = trn,
    nfolds = 10, 
    seed = 3
  )
)

cube.samp3 <- do.call(
  rbind,
  lapply(
    10 * c(seq(1, 9, 2), 9.9),
    FUN = tune.cubist.samp,
    data = trn,
    nfolds = 10, 
    seed = 3
  )
)

cube.extra3 <- do.call(
  rbind,
  lapply(
    10 * c(0, seq(1, 9, 2), 10),
    FUN = tune.cubist.extra,
    data = trn,
    nfolds = 10, 
    seed = 3 
  )
)

#     Seed 4 --------------------------------------------------------------

cube.cmte4 <- do.call(
  rbind,
  lapply(
    10 * c(seq(1, 9, 2), 10),
    FUN = tune.cubist.cmte,
    data = trn,
    nfolds = 10, 
    seed = 4
  )
)

cube.samp4 <- do.call(
  rbind,
  lapply(
    10 * c(seq(1, 9, 2), 9.9),
    FUN = tune.cubist.samp,
    data = trn,
    nfolds = 10, 
    seed = 4
  )
)

cube.extra4 <- do.call(
  rbind,
  lapply(
    10 * c(0, seq(1, 9, 2), 10),
    FUN = tune.cubist.extra,
    data = trn,
    nfolds = 10, 
    seed = 4
  )
)

#     Seed 5 --------------------------------------------------------------

cube.cmte5 <- do.call(
  rbind,
  lapply(
    10 * c(seq(1, 9, 2), 10),
    FUN = tune.cubist.cmte,
    data = trn,
    nfolds = 10, 
    seed = 5
  )
)

cube.samp5 <- do.call(
  rbind,
  lapply(
    10 * c(seq(1, 9, 2), 9.9),
    FUN = tune.cubist.samp,
    data = trn,
    nfolds = 10, 
    seed = 5
  )
)

cube.extra5 <- do.call(
  rbind,
  lapply(
    10 * c(0, seq(1, 9, 2), 10),
    FUN = tune.cubist.extra,
    data = trn,
    nfolds = 10, 
    seed = 5
  )
)

#     Seed 6 --------------------------------------------------------------

cube.cmte6 <- do.call(
  rbind,
  lapply(
    10 * c(seq(1, 9, 2), 10),
    FUN = tune.cubist.cmte,
    data = trn,
    nfolds = 10, 
    seed = 6
  )
)

cube.samp6 <- do.call(
  rbind,
  lapply(
    10 * c(seq(1, 9, 2), 9.9),
    FUN = tune.cubist.samp,
    data = trn,
    nfolds = 10, 
    seed = 6
  )
)

cube.extra6 <- do.call(
  rbind,
  lapply(
    10 * c(0, seq(1, 9, 2), 10),
    FUN = tune.cubist.extra,
    data = trn,
    nfolds = 10, 
    seed = 6
  )
)

#     Compilation ---------------------------------------------------------

cmte.comp <- rbind(cube.cmte, cube.cmte2, cube.cmte3, 
                   cube.cmte4, cube.cmte5, cube.cmte6)
samp.comp <- rbind(cube.samp, cube.samp2, cube.samp3, 
                   cube.samp4, cube.samp5, cube.samp6)
extra.comp <- rbind(cube.extra, cube.extra2, cube.extra3, 
                    cube.extra4, cube.extra5, cube.extra6)

saveRDS(cmte.comp, "./tuning/cubist/cmte.results.RDS")
saveRDS(samp.comp, "./tuning/cubist/samp.results.RDS")
saveRDS(extra.comp, "./tuning/cubist/extra.results.RDS")

# Graphs ------------------------------------------------------------------

# General trend is the committees should be > 25
temp <- cmte.comp
temp$seed <- as.factor(temp$seed)
ggplot(temp, aes(x = committees, y = rsquared, 
                      group = seed, color = seed)) + 
  geom_line() + 
  theme_bw() + 
  labs(x = "Committees", y = "R-squared", color = "Random Seed", 
      title = "Cubist - Tuning Number of Committees")
# ggsave("./tuning/cubist/2017-07-27 cubist commitees.png")

temp <- samp.comp
temp$seed <- as.factor(temp$seed)
ggplot(temp, aes(x = sample.perc, y = rsquared, 
                 group = seed, color = seed)) + 
  geom_line() + 
  theme_bw() + 
  labs(x = "Sample Percentage", y = "R-squared", color = "Random Seed", 
       title = "Cubist - Tuning Sample Percentage")
# ggsave("./tuning/cubist/2017-07-27 cubist sample percentage.png")

temp <- extra.comp
temp$seed <- as.factor(temp$seed)
ggplot(temp, aes(x = extrapolation, y = rsquared, 
                 group = seed, color = seed)) + 
  geom_line() + 
  theme_bw() + 
  labs(x = "Extrapolation", y = "R-squared", color = "Random Seed", 
       title = "Cubist - Tuning Extrapolation")
# ggsave("./tuning/cubist/2017-07-27 cubist extrapolation.png")
