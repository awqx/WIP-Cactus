#     Radial Kernel -------------------------------------------------------
#         Seed 1 ----------------------------------------------------------
# Peaks at 8
results.rbf.cost <-
  do.call(
    rbind,
    lapply(
      2 ^ (0:12),
      FUN = tune.svm.cost,
      data = sprse.trn,
      kerneltype = "radial",
      nfolds = 10,
      seed = 1
    )
  ) 

# Peaks around 0.0005
results.rbf.gamma <-
  do.call(
    rbind,
    lapply(
      2 ^ (-13:2),
      FUN = tune.svm.gamma,
      data = sprse.trn,
      kerneltype = "radial",
      nfolds = 10,
      seed = 1
    )
  )

# Peaks around 1
results.rbf.epsilon <-
  do.call(
    rbind,
    lapply(
      2 ^ (-6:3),
      FUN = tune.svm.epsilon,
      data = sprse.trn,
      kerneltype = "radial",
      nfolds = 10,
      seed = 1
    )
  )

#         Seed 2 ----------------------------------------------------------
# Peaks at 8
results2.rbf.cost <-
  do.call(
    rbind,
    lapply(
      2 ^ (0:12),
      FUN = tune.svm.cost,
      data = sprse.trn,
      kerneltype = "radial",
      nfolds = 10,
      seed = 2
    )
  ) 

# Peaks around 0.0005
results2.rbf.gamma <-
  do.call(
    rbind,
    lapply(
      2 ^ (-13:2),
      FUN = tune.svm.gamma,
      data = sprse.trn,
      kerneltype = "radial",
      nfolds = 10,
      seed = 2
    )
  )

# Peaks around 1
results2.rbf.epsilon <-
  do.call(
    rbind,
    lapply(
      2 ^ (-6:3),
      FUN = tune.svm.epsilon,
      data = sprse.trn,
      kerneltype = "radial",
      nfolds = 10,
      seed = 2
    )
  )

#         Seed 3 ----------------------------------------------------------
# Peaks at 8
results3.rbf.cost <-
  do.call(
    rbind,
    lapply(
      2 ^ (0:12),
      FUN = tune.svm.cost,
      data = sprse.trn,
      kerneltype = "radial",
      nfolds = 10,
      seed = 3
    )
  ) 

# Peaks around 0.0005
results3.rbf.gamma <-
  do.call(
    rbind,
    lapply(
      2 ^ (-13:2),
      FUN = tune.svm.gamma,
      data = sprse.trn,
      kerneltype = "radial",
      nfolds = 10,
      seed = 3
    )
  )

# Peaks around 1
results3.rbf.epsilon <-
  do.call(
    rbind,
    lapply(
      2 ^ (-6:3),
      FUN = tune.svm.epsilon,
      data = sprse.trn,
      kerneltype = "radial",
      nfolds = 10,
      seed = 3
    )
  )

#         Seed 4 ----------------------------------------------------------
# Peaks at 8
results4.rbf.cost <-
  do.call(
    rbind,
    lapply(
      2 ^ (0:12),
      FUN = tune.svm.cost,
      data = sprse.trn,
      kerneltype = "radial",
      nfolds = 10,
      seed = 4
    )
  ) 

# Peaks around 0.0005
results4.rbf.gamma <-
  do.call(
    rbind,
    lapply(
      2 ^ (-13:2),
      FUN = tune.svm.gamma,
      data = sprse.trn,
      kerneltype = "radial",
      nfolds = 10,
      seed = 4
    )
  )

# Peaks around 1
results4.rbf.epsilon <-
  do.call(
    rbind,
    lapply(
      2 ^ (-6:3),
      FUN = tune.svm.epsilon,
      data = sprse.trn,
      kerneltype = "radial",
      nfolds = 10,
      seed = 4
    )
  )

#         Seed 5 ----------------------------------------------------------
# Peaks at 8
results5.rbf.cost <-
  do.call(
    rbind,
    lapply(
      2 ^ (0:12),
      FUN = tune.svm.cost,
      data = sprse.trn,
      kerneltype = "radial",
      nfolds = 10,
      seed = 5
    )
  ) 

# Peaks around 0.0005
results5.rbf.gamma <-
  do.call(
    rbind,
    lapply(
      2 ^ (-13:2),
      FUN = tune.svm.gamma,
      data = sprse.trn,
      kerneltype = "radial",
      nfolds = 10,
      seed = 5
    )
  )

# Peaks around 1
results5.rbf.epsilon <-
  do.call(
    rbind,
    lapply(
      2 ^ (-6:3),
      FUN = tune.svm.epsilon,
      data = sprse.trn,
      kerneltype = "radial",
      nfolds = 10,
      seed = 5
    )
  )
#         Seed 6 ----------------------------------------------------------

# Peaks at 8
results6.rbf.cost <-
  do.call(
    rbind,
    lapply(
      2 ^ (0:12),
      FUN = tune.svm.cost,
      data = sprse.trn,
      kerneltype = "radial",
      nfolds = 10,
      seed = 6
    )
  ) 

# Peaks around 0.0005
results6.rbf.gamma <-
  do.call(
    rbind,
    lapply(
      2 ^ (-13:2),
      FUN = tune.svm.gamma,
      data = sprse.trn,
      kerneltype = "radial",
      nfolds = 10,
      seed = 6
    )
  )

# Peaks around 1
results6.rbf.epsilon <-
  do.call(
    rbind,
    lapply(
      2 ^ (-6:3),
      FUN = tune.svm.epsilon,
      data = sprse.trn,
      kerneltype = "radial",
      nfolds = 10,
      seed = 6
    )
  )

#         Compilation -----------------------------------------------------

rbf.cost.comp <- rbind(results.rbf.cost, results2.rbf.cost, 
                       results3.rbf.cost, results4.rbf.cost, 
                       results5.rbf.cost, results6.rbf.cost)
rbf.gamma.comp <- rbind(results.rbf.gamma, results2.rbf.gamma, 
                        results3.rbf.gamma, results4.rbf.gamma, 
                        results5.rbf.gamma, results6.rbf.gamma)
rbf.epsilon.comp <- rbind(results.rbf.epsilon, results2.rbf.epsilon, 
                          results3.rbf.epsilon, results4.rbf.epsilon, 
                          results5.rbf.epsilon, results6.rbf.epsilon)

saveRDS(rbf.cost.comp, "./tuning/svm/rbf.tune.cost.RDS")
saveRDS(rbf.gamma.comp, "./tuning/svm/rbf.tune.gamma.RDS")
saveRDS(rbf.epsilon.comp, "./tuning/svm/rbf.tune.epsilon.RDS")

#     RBF graphs ----------------------------------------------------------

temp.data <- rbf.cost.comp
temp.data$seed <- as.factor(temp.data$seed)
ggplot(data = temp.data, aes(x = cost, y = rsquared, 
                             group = seed, color = seed)) + 
  geom_line(size = 1) + 
  scale_x_continuous(trans = 'log2') + 
  labs(x = "Cost", y = "R-squared", 
       title = "Radial Basis Kernel - Cost", 
       color = "Random Seed") + 
  theme_bw()
ggsave(filename = "./tuning/svm/2017-07-20 rbf cost.png")

temp.data <- rbf.gamma.comp
temp.data$seed <- as.factor(temp.data$seed)
ggplot(data = temp.data, aes(x = gamma, y = rsquared,
                             group = seed, color = seed)) + 
  geom_line(size = 1) + 
  scale_x_continuous(trans = 'log2') + 
  labs(x = "Gamma", y = "R-squared", 
       title = "Radial Basis Kernel - Gamma", 
       color = "Random Seed") + 
  theme_bw()
ggsave(filename = "./tuning/svm/2017-07-20 rbf gamma.png")

temp.data <- rbf.epsilon.comp
temp.data$seed <- as.factor(temp.data$seed)
ggplot(data = temp.data, aes(x = epsilon, y = rsquared, 
                             group = seed, color = seed)) + 
  geom_line(size = 1) + 
  scale_x_continuous(trans = 'log2') + 
  labs(x = "Epsilon", y = "R-squared", 
       title = "Radial Basis Kernel - Epsilon", 
       color = "Random Seed") + 
  theme_bw()
ggsave("./tuning/svm/2017-07-20 rbf epsilon.png")
# Final parameters:
# for polynomial:coef0 = 2, gamma = 0.0625, cost = 8, epsilon = 1
# for radial: cost = 8, gamma = 0.0005, epsilon = 1
