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

cost.range <- c(1:5, 10, 20, 40)
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

gamma.range <- c(0.001, 0.01, 0.05, 0.1, 0.25, 0.5)
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

epsilon.range <- c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75)
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

coef.range <- c(0:5, 10, 20)
results1.coef <- do.call(rbind, lapply(coef.range, FUN = tune.svm.coef,
                         data = trn, kerneltype = "polynomial", 
                         nfolds = 10, seed = 101))
results2.coef <- do.call(rbind, lapply(coef.range, FUN = tune.svm.coef,
                                       data = trn, kerneltype = "polynomial", 
                                       nfolds = 10, seed = 102))
results3.coef <- do.call(rbind, lapply(coef.range, FUN = tune.svm.coef,
                                       data = trn, kerneltype = "polynomial", 
                                       nfolds = 10, seed = 103))
results.coef <- rbind(results1.coef, results2.coef, 
                      results3.coef) %>%
  mutate(seed = as.factor(seed))
ggplot(results.coef, aes(x = coef, color = seed, group = seed)) + 
  geom_line(aes(y = rsquared)) + 
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
  theme_bw()

dir.create("./tuning/svm/alpha")
saveRDS(results.cost, "./tuning/svm/alpha/poly.cost.RDS")
saveRDS(results.gamma, "./tuning/svm/alpha/poly.gamma.RDS")
saveRDS(results.epsilon, "./tuning/svm/alpha/poly.epsilon.RDS")
saveRDS(results.coef, "./tuning/svm/alpha/poly.coef.RDS")
saveRDS(results.coef, "./tuning/svm/alpha/poly.deg.RDS")

#     Tuning --------------------------------------------------------------

coef.range <- c(0:5, 10, 20)
cost.range <- c(1:5, 10, 20)
deg.range <- 1:4
epsilon.range <- c(0.01, 0.05, 0.1, 0.25, 0.5)
gamma.range <- c(0.001, 0.01, 0.05, 0.1, 0.25, 0.5)
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
# 1055.41    2.88 1082.83 

saveRDS(results.combos, "./tuning/svm/alpha/poly.tuning.RDS")

results.combos[order(results.combos$rsquared, decreasing = T), ] %>% head()
# nfolds     kernel degree cost epsilon gamma coef0  rsquared     rmse
#     5 polynomial      2    4    0.05  0.05     5 0.6946747 2.684865
#     5 polynomial      2    5    0.05  0.05     2 0.6945467 2.626579
results.combos[order(results.combos$rmse), ] %>% head()
# nfolds     kernel degree cost epsilon gamma coef0  rsquared     rmse
#     5 polynomial      2   20    0.25  0.01     2 0.6462556 2.584697
#     5 polynomial      3    5    0.25  0.01     1 0.6481750 2.604591

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

cost.range <- c(1:5, 10, 25)
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
#   scale_x_continuous(tran = "log2") + 
  theme_bw()

#     Gamma ---

gamma.range <- c(0.01, 0.05, 0.1, 0.25, 0.5)
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

epsilon.range <- c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75)
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
  theme_bw()

#     Coef ---

coef.range <- c(0:5, 10)
results1.coef <- do.call(rbind, lapply(coef.range, FUN = tune.svm.coef,
                                       data = trn, kerneltype = "polynomial", 
                                       nfolds = 10, seed = 101)) 
results2.coef <- do.call(rbind, lapply(coef.range, FUN = tune.svm.coef,
                                       data = trn, kerneltype = "polynomial", 
                                       nfolds = 10, seed = 102))
results3.coef <- do.call(rbind, lapply(coef.range, FUN = tune.svm.coef,
                                       data = trn, kerneltype = "polynomial", 
                                       nfolds = 10, seed = 103))
results.coef <- rbind(results1.coef, results2.coef, 
                      results3.coef) %>%
  mutate(seed = as.factor(seed))
ggplot(results.coef, aes(x = coef, color = seed, group = seed)) + 
  geom_line(aes(y = rsquared)) + 
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
  theme_bw()

saveRDS(results.cost, "./tuning/svm/beta/poly.cost.RDS")
saveRDS(results.gamma, "./tuning/svm/beta/poly.gamma.RDS")
saveRDS(results.epsilon, "./tuning/svm/beta/poly.epsilon.RDS")
saveRDS(results.coef, "./tuning/svm/beta/poly.coef.RDS")
saveRDS(results.coef, "./tuning/svm/beta/poly.deg.RDS")

#     Tuning --------------------------------------------------------------

# # 5040 combinations
# cost.range <- c(1:5, 10)
# deg.range <- 1:4
# coef.range <- c(1:5, 10)
# gamma.range <- c(0.01, 0.05, 0.1, 0.25, 0.5)
# epsilon.range <- c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75)

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
# 1165.10    3.78 1230.36 

saveRDS(results.combos, "./tuning/svm/beta/poly.tuning.RDS")

results.combos[order(results.combos$rsquared, decreasing = T), ] %>% head()
# nfolds     kernel degree cost epsilon gamma coef0  rsquared     rmse
#       5 polynomial      4   25    0.75  0.50    10 0.6931928 2.528484
#       5 polynomial      4   25    0.75  0.25     2 0.6914132 2.667870
results.combos[order(results.combos$rmse), ] %>% head()
# nfolds     kernel degree cost epsilon gamma coef0  rsquared     rmse
#       5 polynomial      4   25    0.75  0.50    10 0.6931928 2.528484
#       5 polynomial      4    2    0.75  0.50    10 0.6733786 2.570205

# ========================================================================
# Gamma-CD ----------------------------------------------------------------

#     Loading Data --------------------------------------------------------

dir.create("./tuning/svm/gamma")
# Reading data with all descriptors
trn.all <- readRDS("./pre-process/gamma/1/pp.RDS") 
colnames(trn.all) <- str_replace(colnames(trn.all), "-", ".")
trn.guest <- trn.all$guest
trn <- select(trn.all, -guest)

features <- readRDS("./feature.selection/gamma.vars.RDS")
trn <- trn[ , colnames(trn) %in% c("DelG", features)]

#     Estimation ----------------------------------------------------------

# Working out where an acceptable range of values exists
# Trying to keep ranges to 7 values (arbitrary)

#     Cost ---

cost.range <- c(1, 3, 5, 10, 25, 50, 100)
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

gamma.range <- c(0.01, 0.05, 0.1, 0.25, 0.5)
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

epsilon.range <- c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1)
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
  theme_bw()

#     Coef ---

coef.range <- c(0, 1, 2, 5, 10, 25)
results1.coef <- do.call(rbind, lapply(coef.range, FUN = tune.svm.coef,
                                       data = trn, kerneltype = "polynomial", 
                                       nfolds = 10, seed = 101)) 
results2.coef <- do.call(rbind, lapply(coef.range, FUN = tune.svm.coef,
                                       data = trn, kerneltype = "polynomial", 
                                       nfolds = 10, seed = 102))
results3.coef <- do.call(rbind, lapply(coef.range, FUN = tune.svm.coef,
                                       data = trn, kerneltype = "polynomial", 
                                       nfolds = 10, seed = 103))
results.coef <- rbind(results1.coef, results2.coef, 
                      results3.coef) %>%
  mutate(seed = as.factor(seed))
ggplot(results.coef, aes(x = coef, color = seed, group = seed)) + 
  geom_line(aes(y = rsquared)) + 
  theme_bw()

#     Degree ---

deg.range <- 2:5
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
  theme_bw()

saveRDS(results.cost, "./tuning/svm/gamma/poly.cost.RDS")
saveRDS(results.gamma, "./tuning/svm/gamma/poly.gamma.RDS")
saveRDS(results.epsilon, "./tuning/svm/gamma/poly.epsilon.RDS")
saveRDS(results.coef, "./tuning/svm/gamma/poly.coef.RDS")
saveRDS(results.coef, "./tuning/svm/gamma/poly.deg.RDS")

#     Tuning --------------------------------------------------------------

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
# 1114.51    4.29 1181.49

saveRDS(results.combos, "./tuning/svm/gamma/poly.tuning.RDS")

results.combos[order(results.combos$rsquared, decreasing = T), ] %>% head()
# Clearly outliers, but whatever
# nfolds     kernel degree cost epsilon gamma coef0  rsquared     rmse
#     5 polynomial      3    3    0.01  0.01     2 0.3887948 1.980814
#     5 polynomial      2   25    1.00  0.50     5 0.3655031 1.821173
results.combos[order(results.combos$rmse), ] %>% head()
# nfolds     kernel degree cost epsilon gamma coef0   rsquared     rmse
#     5 polynomial      3    1    0.50  0.01     2 0.3168359 1.529669
#     5 polynomial      2    3    0.10  0.01     0 0.2212526 1.551559