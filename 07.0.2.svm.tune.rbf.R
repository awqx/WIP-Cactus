# Runs tuning for alpha- and beta-CD models with radial kernel
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

cost.range <- 2 ^ (0:7)
results1.cost <- do.call(rbind, lapply(cost.range, FUN = tune.svm.cost, 
                                       data = trn, kerneltype = "radial",
                                       nfolds = 10, seed = 101)) 
results2.cost <- do.call(rbind, lapply(cost.range, FUN = tune.svm.cost,
                                       data = trn, kerneltype = "radial",
                                       nfolds = 10, seed = 102)) 
results3.cost <- do.call(rbind, lapply(cost.range, FUN = tune.svm.cost,
                                       data = trn, kerneltype = "radial", 
                                       nfolds = 10, seed = 103)) 
results.cost <- rbind(results1.cost, results2.cost, results3.cost) %>%
  mutate(seed = as.factor(seed))
ggplot(results.cost, aes(x = cost, color = seed, group = seed)) + 
  geom_line(aes(y = rsquared)) + 
  theme_bw() + 
  scale_x_continuous(trans = "log2")

#     Gamma ---

gamma.range <- 2^(-8:-1)
results1.gamma <- do.call(rbind, lapply(gamma.range, FUN = tune.svm.gamma,
                                        data = trn, kerneltype = "radial", 
                                        nfolds = 10, seed = 101)) 
results2.gamma <- do.call(rbind, lapply(gamma.range, FUN = tune.svm.gamma,
                                        data = trn, kerneltype = "radial", 
                                        nfolds = 10, seed = 102))
results3.gamma <- do.call(rbind, lapply(gamma.range, FUN = tune.svm.gamma,
                                        data = trn, kerneltype = "radial", 
                                        nfolds = 10, seed = 103))
results.gamma <- rbind(results1.gamma, results2.gamma, results3.gamma) %>%
  mutate(seed = as.factor(seed))
ggplot(results.gamma, aes(x = gamma, color = seed, group = seed)) + 
  geom_line(aes(y = rsquared)) + 
  theme_bw() + 
  scale_x_continuous(trans = "log2")

#     Epsilon ---

epsilon.range <- 2^(-6:0)
results1.epsilon <- do.call(rbind, lapply(epsilon.range, FUN = tune.svm.epsilon,
                                          data = trn, kerneltype = "radial",
                                          nfolds = 10, seed = 101))
results2.epsilon <- do.call(rbind, lapply(epsilon.range, FUN = tune.svm.epsilon,
                                          data = trn, kerneltype = "radial",
                                          nfolds = 10, seed = 102))
results3.epsilon <- do.call(rbind, lapply(epsilon.range, FUN = tune.svm.epsilon,
                                          data = trn, kerneltype = "radial",
                                          nfolds = 10, seed = 103))
results4.epsilon <- do.call(rbind, lapply(epsilon.range, FUN = tune.svm.epsilon,
                                          data = trn, kerneltype = "radial",
                                          nfolds = 10, seed = 104))
results.epsilon <- rbind(results1.epsilon, results2.epsilon, 
                         results3.epsilon, results4.epsilon) %>%
  mutate(seed = as.factor(seed))
ggplot(results.epsilon, aes(x = epsilon, color = seed, group = seed)) + 
  geom_line(aes(y = rsquared)) + 
  scale_x_continuous(trans = "log2") + 
  theme_bw()

saveRDS(results.cost, "./tuning/svm/alpha/rbf.cost.RDS")
saveRDS(results.gamma, "./tuning/svm/alpha/rbf.gamma.RDS")
saveRDS(results.epsilon, "./tuning/svm/alpha/rbf.epsilon.RDS")

#     Tuning --------------------------------------------------------------

# 8*8*7 = 448 combinations

svm.combos <- expand.grid(cost.range, gamma.range, 
                          epsilon.range)
colnames(svm.combos) <- c("cost", "gamma", "epsilon")
cost.combos <- svm.combos$cost
gamma.combos <- svm.combos$gamma
eps.combos <- svm.combos$epsilon

set.seed(1001)
system.time(
  results.combos <- do.call(
    rbind,
    mapply(
      FUN = tune.svm.rbf,
      cost = cost.combos, 
      e = eps.combos, 
      g = gamma.combos,
      MoreArgs = 
        list(nfolds = 5, data = trn), 
      SIMPLIFY = F
    )
  )
)

# system.time output
# user  system elapsed 
# 52.35    0.14   54.61 

saveRDS(results.combos, "./tuning/svm/alpha/rbf.tuning.RDS")

# results.combos[order(results.combos$rsquared, decreasing = T), ] %>% head()
# results.combos[order(results.combos$rmse), ] %>% head()

# Best rsquared = 0.755 (rmse = 3.64)
# cost = 2, eps = 0.5, gamma = 0.00390625 
# Best rmse = 3.19 (rsquared = 0.629)
# cost = 8, eps = 0.0625, gamma = 0.0078125

# middle ground: r2 = 0.663, rmse = 3.22
# cost = 2, eps = 0.03125, gamma = 0.00390625

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

cost.range <- 2^(1:7)
results1.cost <- do.call(rbind, lapply(cost.range, FUN = tune.svm.cost, 
                                       data = trn, kerneltype = "radial",
                                       nfolds = 10, seed = 101)) 
results2.cost <- do.call(rbind, lapply(cost.range, FUN = tune.svm.cost,
                                       data = trn, kerneltype = "radial",
                                       nfolds = 10, seed = 102)) 
results3.cost <- do.call(rbind, lapply(cost.range, FUN = tune.svm.cost,
                                       data = trn, kerneltype = "radial", 
                                       nfolds = 10, seed = 103)) 
results.cost <- rbind(results1.cost, results2.cost, results3.cost)
ggplot(results.cost, aes(x = cost, color = seed, group = seed)) + 
  geom_line(aes(y = rsquared)) + 
  theme_bw()

#     Gamma ---

gamma.range <- 2^(-8:-2)
results1.gamma <- do.call(rbind, lapply(gamma.range, FUN = tune.svm.gamma,
                                        data = trn, kerneltype = "radial", 
                                        nfolds = 10, seed = 101))
results2.gamma <- do.call(rbind, lapply(gamma.range, FUN = tune.svm.gamma,
                                        data = trn, kerneltype = "radial", 
                                        nfolds = 10, seed = 102))
results3.gamma <- do.call(rbind, lapply(gamma.range, FUN = tune.svm.gamma,
                                        data = trn, kerneltype = "radial", 
                                        nfolds = 10, seed = 103))
results.gamma <- rbind(results1.gamma, results2.gamma, results3.gamma) %>%
  mutate(seed = as.factor(seed))
ggplot(results.gamma, aes(x = gamma, color = seed, group = seed)) + 
  geom_line(aes(y = rsquared)) + 
  theme_bw() + 
  scale_x_continuous(trans = "log2")

#     Epsilon ---

epsilon.range <- 2^(-5:0)
results1.epsilon <- do.call(rbind, lapply(epsilon.range, FUN = tune.svm.epsilon,
                                          data = trn, kerneltype = "radial",
                                          nfolds = 10, seed = 101))
results2.epsilon <- do.call(rbind, lapply(epsilon.range, FUN = tune.svm.epsilon,
                                          data = trn, kerneltype = "radial",
                                          nfolds = 10, seed = 102))
results3.epsilon <- do.call(rbind, lapply(epsilon.range, FUN = tune.svm.epsilon,
                                          data = trn, kerneltype = "radial",
                                          nfolds = 10, seed = 103))
results4.epsilon <- do.call(rbind, lapply(epsilon.range, FUN = tune.svm.epsilon,
                                          data = trn, kerneltype = "radial",
                                          nfolds = 10, seed = 104))
results.epsilon <- rbind(results1.epsilon, results2.epsilon, 
                         results3.epsilon, results4.epsilon) %>%
  mutate(seed = as.factor(seed))
ggplot(results.epsilon, aes(x = epsilon, color = seed, group = seed)) + 
  geom_line(aes(y = rsquared)) + 
  scale_x_continuous(trans = "log2") + 
  theme_bw()

dir.create("./tuning/svm/beta")
saveRDS(results.cost, "./tuning/svm/beta/rbf.cost.RDS")
saveRDS(results.gamma, "./tuning/svm/beta/rbf.gamma.RDS")
saveRDS(results.epsilon, "./tuning/svm/beta/rbf.epsilon.RDS")

#     Tuning --------------------------------------------------------------

# 7^3 = 343 tuning combinations
svm.combos <- expand.grid(cost.range, gamma.range, 
                          epsilon.range)
colnames(svm.combos) <- c("cost", "gamma", "epsilon")
cost.combos <- svm.combos$cost
gamma.combos <- svm.combos$gamma
eps.combos <- svm.combos$epsilon


set.seed(1001)
system.time(
  results.combos <- do.call(
    rbind,
    mapply(
      FUN = tune.svm.rbf,
      cost = cost.combos, 
      e = eps.combos, 
      g = gamma.combos,
      MoreArgs = 
        list(nfolds = 5, data = trn), 
      SIMPLIFY = F
    )
  )
)

# system.time output
#   user  system elapsed 
# 14.37    0.05   15.49 

saveRDS(results.combos, "./tuning/svm/beta/rbf.tuning.RDS")

# results.combos[order(results.combos$rsquared, decreasing = T), ] %>% head()
# results.combos[order(results.combos$rmse), ] %>% head()

# Best rsquared = 0.709, rmse = 3.67
# cost = 2, eps = 0.0625, gamma = 0.00390625 (rmse = 3.96)
# Best rmse = 3.64, r2 = 0.639
# cost = 32, eps = 0.25, gamma = 0.0078125 (r2 = 0.420)