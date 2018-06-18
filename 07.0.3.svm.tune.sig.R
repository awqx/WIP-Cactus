# Runs tuning for alpha- and beta-CD models with a sigmoid kernel
# SVM code sourced from package e1071
source("./07.0.0.svm.functions.R")

# Alpha-CD ----------------------------------------------------------------

#     Loading Data --------------------------------------------------------

# Reading data with all descriptors
trn.all <- readRDS("./pre-process/alpha/1/pp.RDS") 
colnames(trn.all) <- str_replace(colnames(trn.all), "-", ".")
trn.guest <- trn.all$guest
trn <- select(trn.all, -guest)

rfe1 <- readRDS("./feature.selection/alpha/rfe1.RDS")
trn.pred <- c("DelG", predictors(rfe1))

trn <- trn[ , colnames(trn) %in% trn.pred]

#     Estimation ----------------------------------------------------------

# Working out where an acceptable range of values exists
# Trying to keep ranges to 7 values (arbitrary)

#     Cost ---

cost.range <- c(1, 2, 5, 10, 25, 50, 100)
results1.cost <- do.call(rbind, lapply(cost.range, FUN = tune.svm.cost, 
                                       data = trn, kerneltype = "sigmoid",
                                       nfolds = 10, seed = 101)) 
results2.cost <- do.call(rbind, lapply(cost.range, FUN = tune.svm.cost,
                                       data = trn, kerneltype = "sigmoid",
                                       nfolds = 10, seed = 102)) 
results3.cost <- do.call(rbind, lapply(cost.range, FUN = tune.svm.cost,
                                       data = trn, kerneltype = "sigmoid", 
                                       nfolds = 10, seed = 103)) 
results.cost <- rbind(results1.cost, results2.cost, results3.cost) %>%
  mutate(seed = as.factor(seed))
ggplot(results.cost, aes(x = cost, color = seed, group = seed)) + 
  geom_line(aes(y = rsquared)) + 
  theme_bw()

#     Gamma ---

gamma.range <- c(0.001, 0.01, 0.05, 0.1, 0.25, 0.5, 1)
results1.gamma <- do.call(rbind, lapply(gamma.range, FUN = tune.svm.gamma,
                                        data = trn, kerneltype = "sigmoid", 
                                        nfolds = 10, seed = 101))
results2.gamma <- do.call(rbind, lapply(gamma.range, FUN = tune.svm.gamma,
                                        data = trn, kerneltype = "sigmoid", 
                                        nfolds = 10, seed = 102))
results3.gamma <- do.call(rbind, lapply(gamma.range, FUN = tune.svm.gamma,
                                        data = trn, kerneltype = "sigmoid", 
                                        nfolds = 10, seed = 103))
results.gamma <- rbind(results1.gamma, results2.gamma, results3.gamma) %>%
  mutate(seed = as.factor(seed))
ggplot(results.gamma, aes(x = gamma, color = seed, group = seed)) + 
  geom_line(aes(y = rsquared)) + 
  theme_bw()

#     Epsilon ---

epsilon.range <- c(0, 0.01, 0.05, 0.1, 0.15, 0.25, 0.5, 1)
results1.epsilon <- do.call(rbind, lapply(epsilon.range, FUN = tune.svm.epsilon,
                                          data = trn, kerneltype = "sigmoid",
                                          nfolds = 10, seed = 101)) 
results2.epsilon <- do.call(rbind, lapply(epsilon.range, FUN = tune.svm.epsilon,
                                          data = trn, kerneltype = "sigmoid",
                                          nfolds = 10, seed = 102))
results3.epsilon <- do.call(rbind, lapply(epsilon.range, FUN = tune.svm.epsilon,
                                          data = trn, kerneltype = "sigmoid",
                                          nfolds = 10, seed = 103))
results4.epsilon <- do.call(rbind, lapply(epsilon.range, FUN = tune.svm.epsilon,
                                          data = trn, kerneltype = "sigmoid",
                                          nfolds = 10, seed = 104))
results.epsilon <- rbind(results1.epsilon, results2.epsilon, 
                         results3.epsilon, results4.epsilon) %>%
  mutate(seed = as.factor(seed))
ggplot(results.epsilon, aes(x = epsilon, color = seed, group = seed)) + 
  geom_line(aes(y = rsquared)) + 
  # scale_x_continuous(trans = "log2") + 
  theme_bw()

#     Coef ---

coef.range <- c(0, 1, 2, 3, 4, 7, 10)
results1.coef <- do.call(rbind, lapply(coef.range, FUN = tune.svm.coef,
                                       data = trn, kerneltype = "sigmoid", 
                                       nfolds = 10, seed = 101)) 
results2.coef <- do.call(rbind, lapply(coef.range, FUN = tune.svm.coef,
                                       data = trn, kerneltype = "sigmoid", 
                                       nfolds = 10, seed = 102))
results3.coef <- do.call(rbind, lapply(coef.range, FUN = tune.svm.coef,
                                       data = trn, kerneltype = "sigmoid", 
                                       nfolds = 10, seed = 103))
results4.coef <- do.call(rbind, lapply(coef.range, FUN = tune.svm.coef,
                                       data = trn, kerneltype = "sigmoid", 
                                       nfolds = 10, seed = 104))
results.coef <- rbind(results1.coef, results2.coef, 
                      results3.coef, results4.coef) %>%
  mutate(seed = as.factor(seed))
ggplot(results.coef, aes(x = coef, color = seed, group = seed)) + 
  geom_line(aes(y = rsquared)) + 
  theme_bw()

saveRDS(results.cost, "./tuning/svm/alpha/sig.cost.RDS")
saveRDS(results.gamma, "./tuning/svm/alpha/sig.gamma.RDS")
saveRDS(results.epsilon, "./tuning/svm/alpha/sig.epsilon.RDS")
saveRDS(results.coef, "./tuning/svm/alpha/sig.coef.RDS")

#     Tuning --------------------------------------------------------------

# 2744 combos
svm.combos <- expand.grid(cost.range, gamma.range, 
                          epsilon.range, coef.range)
colnames(svm.combos) <- c("cost", "gamma", "epsilon")
cost.combos <- svm.combos$cost
gamma.combos <- svm.combos$gamma
eps.combos <- svm.combos$epsilon

set.seed(1001)
system.time(
  results.combos <- do.call(
    rbind,
    mapply(
      FUN = tune.svm.sig,
      cost = cost.combos, 
      e = eps.combos, 
      g = gamma.combos,
      MoreArgs = 
        list(nfolds = 5, data = trn, coef = 0), 
      SIMPLIFY = F
    )
  )
)
# system.time output
#  user  system elapsed 
# 88.03    0.05   89.92 

saveRDS(results.combos, "./tuning/svm/alpha/sig.tuning.RDS")

results.combos[order(results.combos$rsquared, decreasing = T), ] %>% head()
results.combos[order(results.combos$rmse), ] %>% head()

# Best rsquared = 0.679, rmse = 5.74
# cost = 2, eps = 0.25, gamma = 0.001, coef0 = 0
# Best rmse = 4.28, r2 = 0.430
# cost = 1, eps = 1, gamma = 0.05, coef0 = 0 

# ========================================================================
# Beta-CD ----------------------------------------------------------------

#     Loading Data --------------------------------------------------------

# Reading data with all descriptors
trn.all <- readRDS("./pre-process/beta/1/pp.RDS") 
colnames(trn.all) <- str_replace(colnames(trn.all), "-", ".")
trn.guest <- trn.all$guest
trn <- select(trn.all, -guest)

rfe1 <- readRDS("./feature.selection/beta/rfe1.RDS")
trn.pred <- c("DelG", predictors(rfe1))

trn <- trn[ , colnames(trn) %in% trn.pred]

#     Estimation ----------------------------------------------------------

# Working out where an acceptable range of values exists
# Trying to keep ranges to 7 values (arbitrary)

#     Cost ---

cost.range <- c(1, 2, 5, 10, 25, 75, 150, 250)
results1.cost <- do.call(rbind, lapply(cost.range, FUN = tune.svm.cost, 
                                       data = trn, kerneltype = "sigmoid",
                                       nfolds = 10, seed = 101)) 
results2.cost <- do.call(rbind, lapply(cost.range, FUN = tune.svm.cost,
                                       data = trn, kerneltype = "sigmoid",
                                       nfolds = 10, seed = 102)) 
results3.cost <- do.call(rbind, lapply(cost.range, FUN = tune.svm.cost,
                                       data = trn, kerneltype = "sigmoid", 
                                       nfolds = 10, seed = 103)) 
results.cost <- rbind(results1.cost, results2.cost, results3.cost) %>%
  mutate(seed = as.factor(seed))
ggplot(results.cost, aes(x = cost, color = seed, group = seed)) + 
  geom_line(aes(y = rsquared)) + 
  theme_bw()

#     Gamma ---

gamma.range <- c(0.001, 0.01, 0.1, 0.25, 0.5, 1)
results1.gamma <- do.call(rbind, lapply(gamma.range, FUN = tune.svm.gamma,
                                        data = trn, kerneltype = "sigmoid", 
                                        nfolds = 10, seed = 101))
results2.gamma <- do.call(rbind, lapply(gamma.range, FUN = tune.svm.gamma,
                                        data = trn, kerneltype = "sigmoid", 
                                        nfolds = 10, seed = 102))
results3.gamma <- do.call(rbind, lapply(gamma.range, FUN = tune.svm.gamma,
                                        data = trn, kerneltype = "sigmoid", 
                                        nfolds = 10, seed = 103))
results.gamma <- rbind(results1.gamma, results2.gamma, results3.gamma) %>%
  mutate(seed = as.factor(seed))
ggplot(results.gamma, aes(x = gamma, color = seed, group = seed)) + 
  geom_line(aes(y = rsquared)) + 
  theme_bw() 

#     Epsilon ---

epsilon.range <- c(0, 0.01, 0.1, 0.25, 0.5, 0.75, 1)
results1.epsilon <- do.call(rbind, lapply(epsilon.range, FUN = tune.svm.epsilon,
                                          data = trn, kerneltype = "sigmoid",
                                          nfolds = 10, seed = 101)) 
results2.epsilon <- do.call(rbind, lapply(epsilon.range, FUN = tune.svm.epsilon,
                                          data = trn, kerneltype = "sigmoid",
                                          nfolds = 10, seed = 102))
results3.epsilon <- do.call(rbind, lapply(epsilon.range, FUN = tune.svm.epsilon,
                                          data = trn, kerneltype = "sigmoid",
                                          nfolds = 10, seed = 103))
results4.epsilon <- do.call(rbind, lapply(epsilon.range, FUN = tune.svm.epsilon,
                                          data = trn, kerneltype = "sigmoid",
                                          nfolds = 10, seed = 104))
results.epsilon <- rbind(results1.epsilon, results2.epsilon, 
                         results3.epsilon, results4.epsilon) %>%
  mutate(seed = as.factor(seed))
ggplot(results.epsilon, aes(x = epsilon, color = seed, group = seed)) + 
  geom_line(aes(y = rsquared)) + 
  theme_bw()

#     Coef ---

coef.range <- c(0, 1, 2, 3, 7, 10)
results1.coef <- do.call(rbind, lapply(coef.range, FUN = tune.svm.coef,
                                       data = trn, kerneltype = "sigmoid", 
                                       nfolds = 10, seed = 101)) 
results2.coef <- do.call(rbind, lapply(coef.range, FUN = tune.svm.coef,
                                       data = trn, kerneltype = "sigmoid", 
                                       nfolds = 10, seed = 102))
results3.coef <- do.call(rbind, lapply(coef.range, FUN = tune.svm.coef,
                                       data = trn, kerneltype = "sigmoid", 
                                       nfolds = 10, seed = 103))
results4.coef <- do.call(rbind, lapply(coef.range, FUN = tune.svm.coef,
                                       data = trn, kerneltype = "sigmoid", 
                                       nfolds = 10, seed = 104))
results.coef <- rbind(results1.coef, results2.coef, 
                      results3.coef, results4.coef) %>%
  mutate(seed = as.factor(seed))
ggplot(results.coef, aes(x = coef, color = seed, group = seed)) + 
  geom_line(aes(y = rsquared)) + 
  theme_bw()

saveRDS(results.cost, "./tuning/svm/beta/sig.cost.RDS")
saveRDS(results.gamma, "./tuning/svm/beta/sig.gamma.RDS")
saveRDS(results.epsilon, "./tuning/svm/beta/sig.epsilon.RDS")
saveRDS(results.coef, "./tuning/svm/beta/sig.coef.RDS")

#     Tuning --------------------------------------------------------------

# 2016 tuning combinations
svm.combos <- expand.grid(cost.range, gamma.range, 
                          epsilon.range, coef.range)
colnames(svm.combos) <- c("cost", "gamma", "epsilon", "coef")
cost.combos <- svm.combos$cost
gamma.combos <- svm.combos$gamma
eps.combos <- svm.combos$epsilon
coef.combos <- svm.combos$coef


set.seed(1001)
system.time(
  results.combos <- do.call(
    rbind,
    mapply(
      FUN = tune.svm.sig,
      cost = cost.combos, 
      e = eps.combos, 
      g = gamma.combos,
      coef = coef.combos,
      MoreArgs = 
        list(nfolds = 5, data = trn), 
      SIMPLIFY = F
    )
  )
)
# system.time output
# user  system elapsed 
# 83.59    0.03   84.06

saveRDS(results.combos, "./tuning/svm/beta/sig.tuning.RDS")

results.combos[order(results.combos$rsquared, decreasing = T), ] %>% head()
results.combos[order(results.combos$rmse), ] %>% head()

# Best rsquared = 0.682 (rmse = 2.94)
# cost = 250, epsilon = 0.25, gamma = 0.001, coef0 = 1
# Best rmse = 2.84 (r2 = 0.634)
# cost = 5, epsilon = 0.25, gamma = 0.01, coef0 = 0
