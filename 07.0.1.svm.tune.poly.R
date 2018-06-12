# Runs tuning for alpha- and beta-CD models with polynomial kernel
# SVM code sourced from package e1071
source("./07.0.0.svm.functions.R")

# Alpha-CD ----------------------------------------------------------------

#     Loading Data --------------------------------------------------------

dir.create("./tuning")
dir.create("./tuning/svm")
# Reading data with all descriptors
trn.all <- readRDS("./model.data/alpha/trn1.RDS") 
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

cost.range <- 2*(1:7)
results1.cost <- do.call(rbind, lapply(cost.range, FUN = tune.svm.cost, 
                         data = trn, kerneltype = "polynomial",
                         nfolds = 10, seed = 101)) 
results2.cost <- do.call(rbind, lapply(cost.range, FUN = tune.svm.cost,
                         data = trn, kerneltype = "polynomial",
                         nfolds = 10, seed = 102)) 
results3.cost <- do.call(rbind, lapply(cost.range, FUN = tune.svm.cost,
                         data = trn, kerneltype = "polynomial", 
                         nfolds = 10, seed = 103)) 
results.cost <- rbind(results1.cost, results2.cost, results3.cost)
ggplot(results.cost, aes(x = cost, color = seed, group = seed)) + 
  geom_line(aes(y = rsquared)) + 
  theme_bw()

#     Gamma ---

gamma.range <- 2^(-8:-2)
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
  theme_bw() + 
  scale_x_continuous(trans = "log2")

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

coef.range <- 2^(0:6)
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
  scale_x_continuous(trans = "log2") + 
  theme_bw()

#     Degree ---

deg.range <- 1:4
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

# 7^4 = 2401 tuning combinations

cost.range <- 2*(1:7)
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
# user        system      elapsed 
# 1322.01     1.55        1370.91 

saveRDS(results.combos, "./tuning/svm/alpha/poly.tuning.RDS")

# results.combos[order(results.combos$rsquared, decreasing = T), ] %>% head()
# results.combos[order(results.combos$rmse), ] %>% head()

# Best rsquared (0.504)
# degree = 3, cost = 6, eps = 0.125, gamma = 0.-3125, coef0 = 2 (rmse = 3.78)
# Best rmse (3.68)
# deg = 1, cost = 8, eps = 0.25, gamma = 0.0078125, coef0 = 64 (r2 = 0.475) 
# or, for a middle ground
# rmse = 3.72, r2 = 0.498
# deg = 3, cost = 4, eps = 0.125, gamma = 0.00390625, coef0 = 1

# ========================================================================
# Beta-CD ----------------------------------------------------------------

#     Loading Data --------------------------------------------------------

dir.create("./tuning/svm/beta")
# Reading data with all descriptors
trn.all <- readRDS("./model.data/beta/trn8.RDS") 
colnames(trn.all) <- str_replace(colnames(trn.all), "-", ".")
trn.guest <- trn.all$guest
trn <- select(trn.all, -guest)

rfe1 <- readRDS("./feature.selection/beta/rfe8.RDS")
trn.pred <- c("DelG", predictors(rfe1))

trn <- trn[ , colnames(trn) %in% trn.pred]

#     Estimation ----------------------------------------------------------

# Working out where an acceptable range of values exists
# Trying to keep ranges to 7 values (arbitrary)

#     Cost ---

cost.range <- 2^(1:7)
results1.cost <- do.call(rbind, lapply(cost.range, FUN = tune.svm.cost, 
                                       data = trn, kerneltype = "polynomial",
                                       nfolds = 10, seed = 101)) 
results2.cost <- do.call(rbind, lapply(cost.range, FUN = tune.svm.cost,
                                       data = trn, kerneltype = "polynomial",
                                       nfolds = 10, seed = 102)) 
results3.cost <- do.call(rbind, lapply(cost.range, FUN = tune.svm.cost,
                                       data = trn, kerneltype = "polynomial", 
                                       nfolds = 10, seed = 103)) 
results.cost <- rbind(results1.cost, results2.cost, results3.cost)
ggplot(results.cost, aes(x = cost, color = seed, group = seed)) + 
  geom_line(aes(y = rsquared)) + 
  theme_bw()

#     Gamma ---

gamma.range <- 2^(-8:-2)
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
  theme_bw() + 
  scale_x_continuous(trans = "log2")

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
  scale_x_continuous(trans = "log2") + 
  theme_bw()

#     Coef ---

coef.range <- c(2^(-1:5))
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
  scale_x_continuous(trans = "log2") + 
  theme_bw()

#     Degree ---

deg.range <- 1:4
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

# 7^4 = 2401 tuning combinations
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
# user        system      elapsed 
# 4403.22     1.95        4499.24

saveRDS(results.combos, "./tuning/svm/beta/poly.tuning.RDS")

# results.combos[order(results.combos$rsquared, decreasing = T), ] %>% head()
# results.combos[order(results.combos$rmse), ] %>% head()

# Best rsquared (0.486)
# degree = 2, cost = 2, eps = 0.03125, gamma = 0.0078125, coef0 = 8 (rmse = 3.880)
# Best rmse = same parameters as above