# Runs tuning for alpha- and beta-CD models with polynomial kernel
# SVM code sourced from package e1071
source("./07.0.0.svm.functions.R")

# Alpha-CD ----------------------------------------------------------------

#     Loading Data --------------------------------------------------------

dir.create("./tuning")
dir.create("./tuning/svm")
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

cost.range <- c(1, 5, 10, 20, 40, 75, 150)
results1.cost <- do.call(rbind, lapply(cost.range, FUN = tune.svm.cost, 
                         data = trn, kerneltype = "polynomial",
                         nfolds = 10, seed = 101)) 
results2.cost <- do.call(rbind, lapply(cost.range, FUN = tune.svm.cost,
                         data = trn, kerneltype = "polynomial",
                         nfolds = 10, seed = 102)) 
results3.cost <- do.call(rbind, lapply(cost.range, FUN = tune.svm.cost,
                         data = trn, kerneltype = "polynomial", 
                         nfolds = 10, seed = 103)) 
results.cost <- rbind(results1.cost, results2.cost, results3.cost) %>%
  mutate(seed = as.factor(seed))
ggplot(results.cost, aes(x = cost, color = seed, group = seed)) + 
  geom_line(aes(y = rsquared)) + 
 # scale_x_continuous(tran = "log2") + 
  theme_bw()

#     Gamma ---

gamma.range <- c(0.001, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1)
results1.gamma <- do.call(rbind, lapply(gamma.range, FUN = tune.svm.gamma,
                          data = trn, kerneltype = "polynomial", 
                          nfolds = 10, seed = 101))
results2.gamma <- do.call(rbind, lapply(gamma.range, FUN = tune.svm.gamma,
                                        data = trn, kerneltype = "polynomial", 
                                        nfolds = 10, seed = 102))
results3.gamma <- do.call(rbind, lapply(gamma.range, FUN = tune.svm.gamma,
                                        data = trn, kerneltype = "polynomial", 
                                        nfolds = 10, seed = 103))
results.gamma <- rbind(results1.gamma, results2.gamma, results3.gamma) %>%
  mutate(seed = as.factor(seed))
ggplot(results.gamma, aes(x = gamma, color = seed, group = seed)) + 
  geom_line(aes(y = rsquared)) + 
  theme_bw() 

#     Epsilon ---

epsilon.range <- 2^(-6:0)
results1.epsilon <- do.call(rbind, lapply(epsilon.range, FUN = tune.svm.epsilon,
                            data = trn, kerneltype = "polynomial",
                            nfolds = 10, seed = 101))
results2.epsilon <- do.call(rbind, lapply(epsilon.range, FUN = tune.svm.epsilon,
                                          data = trn, kerneltype = "polynomial",
                                          nfolds = 10, seed = 102))
results3.epsilon <- do.call(rbind, lapply(epsilon.range, FUN = tune.svm.epsilon,
                                          data = trn, kerneltype = "polynomial",
                                          nfolds = 10, seed = 103))
results4.epsilon <- do.call(rbind, lapply(epsilon.range, FUN = tune.svm.epsilon,
                                          data = trn, kerneltype = "polynomial",
                                          nfolds = 10, seed = 104))
results.epsilon <- rbind(results1.epsilon, results2.epsilon, 
                         results3.epsilon, results4.epsilon) %>%
  mutate(seed = as.factor(seed))
ggplot(results.epsilon, aes(x = epsilon, color = seed, group = seed)) + 
  geom_line(aes(y = rsquared)) + 
  # scale_x_continuous(trans = "log2") + 
  theme_bw()

#     Coef ---

coef.range <- c(0, 1, 5, 10, 20, 40, 75)
results1.coef <- do.call(rbind, lapply(coef.range, FUN = tune.svm.coef,
                         data = trn, kerneltype = "polynomial", 
                         nfolds = 10, seed = 101))
results2.coef <- do.call(rbind, lapply(coef.range, FUN = tune.svm.coef,
                                       data = trn, kerneltype = "polynomial", 
                                       nfolds = 10, seed = 102))
results3.coef <- do.call(rbind, lapply(coef.range, FUN = tune.svm.coef,
                                       data = trn, kerneltype = "polynomial", 
                                       nfolds = 10, seed = 103))
results4.coef <- do.call(rbind, lapply(coef.range, FUN = tune.svm.coef,
                                       data = trn, kerneltype = "polynomial", 
                                       nfolds = 10, seed = 104))
results.coef <- rbind(results1.coef, results2.coef, 
                      results3.coef, results4.coef) %>%
  mutate(seed = as.factor(seed))
ggplot(results.coef, aes(x = coef, color = seed, group = seed)) + 
  geom_line(aes(y = rsquared)) + 
  theme_bw()

#     Degree ---

deg.range <- 1:5
results1.deg <- do.call(rbind, lapply(deg.range, FUN = tune.svm.degree,
                                      data = trn, nfolds = 10, seed = 101)) 
results2.deg <- do.call(rbind, lapply(deg.range, FUN = tune.svm.degree,
                                   data = trn,  nfolds = 10, seed = 102))
results3.deg <- do.call(rbind, lapply(deg.range, FUN = tune.svm.degree,
                                      data = trn, nfolds = 10, seed = 103))
results.deg <- rbind(results1.deg, results2.deg, results3.deg) %>%
  mutate(seed = as.factor(seed))
ggplot(results.deg, aes(x = degree, color = seed, group = seed)) + 
  geom_line(aes(y = rsquared)) + 
  scale_x_continuous(trans = "log2") + 
  theme_bw()

dir.create("./tuning/svm/alpha")
saveRDS(results.cost, "./tuning/svm/alpha/poly.cost.RDS")
saveRDS(results.gamma, "./tuning/svm/alpha/poly.gamma.RDS")
saveRDS(results.epsilon, "./tuning/svm/alpha/poly.epsilon.RDS")
saveRDS(results.coef, "./tuning/svm/alpha/poly.coef.RDS")
saveRDS(results.coef, "./tuning/svm/alpha/poly.deg.RDS")

#     Tuning --------------------------------------------------------------

coef.range <- c(0, 1, 5, 10, 20, 40, 75)
cost.range <- c(1, 5, 10, 20, 40, 75)
deg.range <- 1:4
epsilon.range <- 2^(-6:0)
gamma.range <- c(0.001, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75)
# 7*5*7*6*4 = 9604 tuning combinations
svm.combos <- expand.grid(cost.range, gamma.range, 
                          epsilon.range, coef.range, deg.range)
colnames(svm.combos) <- c("cost", "gamma", "epsilon", "coef", "degree")
cost.combos <- svm.combos$cost
gamma.combos <- svm.combos$gamma
eps.combos <- svm.combos$epsilon
coef.combos <- svm.combos$coef
deg.combos <- svm.combos$degree

set.seed(1001)
system.time(
  results.combos <- do.call(
    rbind,
    mapply(
      FUN = tune.svm.poly,
      cost = cost.combos, 
      e = eps.combos, 
      g = gamma.combos,
      deg = deg.combos, 
      coef = coef.combos,
      MoreArgs = 
        list(nfolds = 5, data = trn), 
      SIMPLIFY = F
    )
  )
)
# system.time output
# user  system elapsed 
# 7908.05    7.72 8145.98 

saveRDS(results.combos, "./tuning/svm/alpha/poly.tuning.RDS")

results.combos[order(results.combos$rsquared, decreasing = T), ] %>% head()
results.combos[order(results.combos$rmse), ] %>% head()

# Best rsquared (0.775) (rmse = 604.5, no that isn't a typo)
# degree = 4, cost = 20, eps = 0.25, gamma = 0.05, coef0 = 1 
# Best rmse = 3.55 (r2 = 0.665)
# deg = 4, cost = 1, epsilon = 0.5, gamma = 0.1, coef0 = 10
# middle ground: r2 = 0.684, rmse = 4.02
# degree = 4, cost = 10, epsilon = 0.125, gamma = 1, coef0 = 20

# r2 = 0.604, rmse = 3.08
# degree = 4, cost = 5, eps = 0.25, gammma = 0.01, coef0 = 1

# ========================================================================
# Beta-CD ----------------------------------------------------------------

#     Loading Data --------------------------------------------------------

dir.create("./tuning/svm/beta")
# Reading data with all descriptors
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

cost.range <- 2^(1:6)
results1.cost <- do.call(rbind, lapply(cost.range, FUN = tune.svm.cost, 
                                       data = trn, kerneltype = "polynomial",
                                       nfolds = 10, seed = 101)) 
results2.cost <- do.call(rbind, lapply(cost.range, FUN = tune.svm.cost,
                                       data = trn, kerneltype = "polynomial",
                                       nfolds = 10, seed = 102)) 
results3.cost <- do.call(rbind, lapply(cost.range, FUN = tune.svm.cost,
                                       data = trn, kerneltype = "polynomial", 
                                       nfolds = 10, seed = 103)) 
results.cost <- rbind(results1.cost, results2.cost, results3.cost) %>%
  mutate(seed = as.factor(seed))
ggplot(results.cost, aes(x = cost, color = seed)) + 
  geom_line(aes(y = rsquared)) + 
  scale_x_continuous(tran = "log2") + 
  theme_bw()

#     Gamma ---

gamma.range <- c(0.001, 0.01, 0.05, 0.1, 0.25, 0.5, 1)
results1.gamma <- do.call(rbind, lapply(gamma.range, FUN = tune.svm.gamma,
                                        data = trn, kerneltype = "polynomial", 
                                        nfolds = 10, seed = 101))
results2.gamma <- do.call(rbind, lapply(gamma.range, FUN = tune.svm.gamma,
                                        data = trn, kerneltype = "polynomial", 
                                        nfolds = 10, seed = 102))
results3.gamma <- do.call(rbind, lapply(gamma.range, FUN = tune.svm.gamma,
                                        data = trn, kerneltype = "polynomial", 
                                        nfolds = 10, seed = 103))
results.gamma <- rbind(results1.gamma, results2.gamma, results3.gamma) %>%
  mutate(seed = as.factor(seed))
ggplot(results.gamma, aes(x = gamma, color = seed, group = seed)) + 
  geom_line(aes(y = rsquared)) + 
  theme_bw() 

#     Epsilon ---

epsilon.range <- 2^(-4:0)
results1.epsilon <- do.call(rbind, lapply(epsilon.range, FUN = tune.svm.epsilon,
                                          data = trn, kerneltype = "polynomial",
                                          nfolds = 10, seed = 101))
results2.epsilon <- do.call(rbind, lapply(epsilon.range, FUN = tune.svm.epsilon,
                                          data = trn, kerneltype = "polynomial",
                                          nfolds = 10, seed = 102))
results3.epsilon <- do.call(rbind, lapply(epsilon.range, FUN = tune.svm.epsilon,
                                          data = trn, kerneltype = "polynomial",
                                          nfolds = 10, seed = 103))
results4.epsilon <- do.call(rbind, lapply(epsilon.range, FUN = tune.svm.epsilon,
                                          data = trn, kerneltype = "polynomial",
                                          nfolds = 10, seed = 104))
results.epsilon <- rbind(results1.epsilon, results2.epsilon, 
                         results3.epsilon, results4.epsilon) %>%
  mutate(seed = as.factor(seed))
ggplot(results.epsilon, aes(x = epsilon, color = seed, group = seed)) + 
  geom_line(aes(y = rsquared)) + 
  scale_x_continuous(trans = "log2") + 
  theme_bw()

#     Coef ---

coef.range <- c(0, 1, 5, 10, 25, 50, 75)
results1.coef <- do.call(rbind, lapply(coef.range, FUN = tune.svm.coef,
                                       data = trn, kerneltype = "polynomial", 
                                       nfolds = 10, seed = 101)) 
results2.coef <- do.call(rbind, lapply(coef.range, FUN = tune.svm.coef,
                                       data = trn, kerneltype = "polynomial", 
                                       nfolds = 10, seed = 102))
results3.coef <- do.call(rbind, lapply(coef.range, FUN = tune.svm.coef,
                                       data = trn, kerneltype = "polynomial", 
                                       nfolds = 10, seed = 103))
results4.coef <- do.call(rbind, lapply(coef.range, FUN = tune.svm.coef,
                                       data = trn, kerneltype = "polynomial", 
                                       nfolds = 10, seed = 104))
results.coef <- rbind(results1.coef, results2.coef, 
                      results3.coef, results4.coef) %>%
  mutate(seed = as.factor(seed))
ggplot(results.coef, aes(x = coef, color = seed, group = seed)) + 
  geom_line(aes(y = rsquared)) + 
  theme_bw()

#     Degree ---

deg.range <- 1:5
results1.deg <- do.call(rbind, lapply(deg.range, FUN = tune.svm.degree,
                                       data = trn, nfolds = 10, seed = 101)) 
results2.deg <- do.call(rbind, lapply(deg.range, FUN = tune.svm.degree,
                                       data = trn,  nfolds = 10, seed = 102))
results3.deg <- do.call(rbind, lapply(deg.range, FUN = tune.svm.degree,
                                       data = trn, nfolds = 10, seed = 103))
results.deg <- rbind(results1.deg, results2.deg, results3.deg) %>%
  mutate(seed = as.factor(seed))
ggplot(results.deg, aes(x = degree, color = seed, group = seed)) + 
  geom_line(aes(y = rsquared)) + 
  scale_x_continuous(trans = "log2") + 
  theme_bw()

saveRDS(results.cost, "./tuning/svm/beta/poly.cost.RDS")
saveRDS(results.gamma, "./tuning/svm/beta/poly.gamma.RDS")
saveRDS(results.epsilon, "./tuning/svm/beta/poly.epsilon.RDS")
saveRDS(results.coef, "./tuning/svm/beta/poly.coef.RDS")
saveRDS(results.coef, "./tuning/svm/beta/poly.deg.RDS")

#     Tuning --------------------------------------------------------------

# 5040 combinations
cost.range <- c(1, 5, 10, 20, 50)
deg.range <- 1:4
coef.range <- c(0, 1, 2, 5, 10, 25)
gamma.range <- c(0.001, 0.01, 0.05, 0.1, 0.25, 0.5)
epsilon.range <- 2^(-4:0)
# 3600 combos
svm.combos <- expand.grid(cost.range, gamma.range, 
                          epsilon.range, coef.range, deg.range)
colnames(svm.combos) <- c("cost", "gamma", "epsilon", "coef", "degree")
cost.combos <- svm.combos$cost
gamma.combos <- svm.combos$gamma
eps.combos <- svm.combos$epsilon
coef.combos <- svm.combos$coef
deg.combos <- svm.combos$degree

set.seed(1001)
system.time(
  results.combos <- do.call(
    rbind,
    mapply(
      FUN = tune.svm.poly,
      cost = cost.combos, 
      e = eps.combos, 
      g = gamma.combos,
      deg = deg.combos, 
      coef = coef.combos,
      MoreArgs = 
        list(nfolds = 5, data = trn), 
      SIMPLIFY = F
    )
  )
)

# system.time output
# user  system elapsed 
# 1633.64    1.28 1682.16 

saveRDS(results.combos, "./tuning/svm/beta/poly.tuning.RDS")

results.combos[order(results.combos$rsquared, decreasing = T), ] %>% head()
results.combos[order(results.combos$rmse), ] %>% head()

# Best rsquared = 0.761, rmse = 2.70
# degree = 3, cost = 10, eps = 0.25, gamma = 0.001, coef0 = 25
# Best rmse = 2.65 (r2 = 0.757)
# deg = 4, cost = 10, eps = 0.25, gamma = 0.001, coef0 = 5