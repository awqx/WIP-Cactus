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

features <- readRDS("./feature.selection/alpha.vars.RDS")
trn <- trn[ , colnames(trn) %in% c("DelG", features)]

#     Estimation ----------------------------------------------------------

# Working out where an acceptable range of values exists
# Trying to keep ranges to 7 values (arbitrary)

#     Cost ---

cost.range <- c(1:5, 10, 25)
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
  theme_bw()

#     Gamma ---

gamma.range <- c(0.01, 0.05, 0.1, 0.25)
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
  theme_bw() 

#     Epsilon ---

epsilon.range <- c(0.01, 0.05, 0.1, 0.25, 0.5)
results1.epsilon <- do.call(rbind, lapply(epsilon.range, FUN = tune.svm.epsilon,
                                          data = trn, kerneltype = "radial",
                                          nfolds = 10, seed = 101))
results2.epsilon <- do.call(rbind, lapply(epsilon.range, FUN = tune.svm.epsilon,
                                          data = trn, kerneltype = "radial",
                                          nfolds = 10, seed = 102))
results3.epsilon <- do.call(rbind, lapply(epsilon.range, FUN = tune.svm.epsilon,
                                          data = trn, kerneltype = "radial",
                                          nfolds = 10, seed = 103))
results.epsilon <- rbind(results1.epsilon, results2.epsilon, 
                         results3.epsilon) %>%
  mutate(seed = as.factor(seed))
ggplot(results.epsilon, aes(x = epsilon, color = seed, group = seed)) + 
  geom_line(aes(y = rsquared)) + 
  theme_bw()

saveRDS(results.cost, "./tuning/svm/alpha/rbf.cost.RDS")
saveRDS(results.gamma, "./tuning/svm/alpha/rbf.gamma.RDS")
saveRDS(results.epsilon, "./tuning/svm/alpha/rbf.epsilon.RDS")

#     Tuning --------------------------------------------------------------

# 7*8*7 = 392 combinations

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
        list(nfolds = 10, data = trn), 
      SIMPLIFY = F
    )
  )
)

# system.time output
# user  system elapsed 
# 22.18    0.06   22.49 
# 19.84    0.09   20.19

saveRDS(results.combos, "./tuning/svm/alpha/rbf.tuning.RDS")

results.combos[order(results.combos$rsquared, decreasing = T), ] %>% head()
# nfolds kernel cost epsilon gamma  rsquared     rmse
#    10    rbf    2    0.10  0.01 0.6834522 2.694921
#    10    rbf    1    0.25  0.01 0.6745141 2.808424
results.combos[order(results.combos$rmse), ] %>% head()
# nfolds kernel cost epsilon gamma  rsquared     rmse
#    10    rbf    2    0.10  0.01 0.6834522 2.694921
#    10    rbf    2    0.05  0.01 0.6689463 2.695471

# ========================================================================
# Beta-CD ----------------------------------------------------------------

#     Loading Data --------------------------------------------------------

trn.all <- readRDS("./pre-process/beta/1/pp.RDS") 
colnames(trn.all) <- str_replace(colnames(trn.all), "-", ".")
trn.guest <- trn.all$guest
trn <- select(trn.all, -guest)

features <- readRDS("./feature.selection/beta.vars.RDS")
trn <- trn[ , colnames(trn) %in% c("DelG", features)]

#     Estimation ----------------------------------------------------------

# Working out where an acceptable range of values exists
# Trying to keep ranges to 7 values (arbitrary)

#     Cost ---

cost.range <- c(1:5, 10, 25)
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
  theme_bw()

#     Gamma ---

gamma.range <- c(0.01, 0.05, 0.05, 0.1, 0.25)
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
  theme_bw() 

#     Epsilon ---

epsilon.range <- c(0.01, 0.05, 0.1, 0.25, 0.5)
results1.epsilon <- do.call(rbind, lapply(epsilon.range, FUN = tune.svm.epsilon,
                                          data = trn, kerneltype = "radial",
                                          nfolds = 10, seed = 101))
results2.epsilon <- do.call(rbind, lapply(epsilon.range, FUN = tune.svm.epsilon,
                                          data = trn, kerneltype = "radial",
                                          nfolds = 10, seed = 102))
results3.epsilon <- do.call(rbind, lapply(epsilon.range, FUN = tune.svm.epsilon,
                                          data = trn, kerneltype = "radial",
                                          nfolds = 10, seed = 103))
results.epsilon <- rbind(results1.epsilon, results2.epsilon, 
                         results3.epsilon) %>%
  mutate(seed = as.factor(seed))
ggplot(results.epsilon, aes(x = epsilon, color = seed, group = seed)) + 
  geom_line(aes(y = rsquared)) + 
  theme_bw()

dir.create("./tuning/svm/beta")
saveRDS(results.cost, "./tuning/svm/beta/rbf.cost.RDS")
saveRDS(results.gamma, "./tuning/svm/beta/rbf.gamma.RDS")
saveRDS(results.epsilon, "./tuning/svm/beta/rbf.epsilon.RDS")

#     Tuning --------------------------------------------------------------

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
        list(nfolds = 10, data = trn), 
      SIMPLIFY = F
    )
  )
)

# system.time output
# user  system elapsed 
#  42.37    0.21   45.09 
# 39.52    0.21   41.27 
saveRDS(results.combos, "./tuning/svm/beta/rbf.tuning.RDS")

results.combos[order(results.combos$rsquared, decreasing = T), ] %>% head()
# nfolds kernel cost epsilon gamma  rsquared     rmse
#    10    rbf    3    0.50  0.01 0.7231572 3.099700
#    10    rbf    5    0.01  0.01 0.7168283 3.103063
results.combos[order(results.combos$rmse), ] %>% head()
# nfolds kernel cost epsilon gamma  rsquared     rmse
#     10    rbf   10    0.05  0.01 0.7075812 2.992280
#    10    rbf    5    0.25  0.01 0.6980419 3.012280