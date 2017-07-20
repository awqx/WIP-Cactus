#     Seed 1 ----------------------------------------------------------

sig.cost <-
  do.call(
    rbind,
    lapply(
      2 ^ (-2:6),
      FUN = tune.svm.cost,
      data = sprse.trn,
      kerneltype = "sigmoid",
      nfolds = 10,
      seed = 1
    )
  )

sig.gamma <-
  do.call(
    rbind,
    lapply(
      2 ^ (-12:1),
      FUN = tune.svm.gamma,
      data = sprse.trn,
      kerneltype = "sigmoid",
      nfolds = 10,
      seed = 1
    )
  )

sig.epsilon <-
  do.call(
    rbind,
    lapply(
      2 ^ (-8:3),
      FUN = tune.svm.epsilon,
      data = sprse.trn,
      kerneltype = "sigmoid",
      nfolds = 10,
      seed = 1
    )
  )

#     Seed 2 --------------------------------------------------------------

sig2.cost <-
  do.call(
    rbind,
    lapply(
      2 ^ (-2:6),
      FUN = tune.svm.cost,
      data = sprse.trn,
      kerneltype = "sigmoid",
      nfolds = 10,
      seed = 2
    )
  )

sig2.gamma <-
  do.call(
    rbind,
    lapply(
      2 ^ (-12:1),
      FUN = tune.svm.gamma,
      data = sprse.trn,
      kerneltype = "sigmoid",
      nfolds = 10,
      seed = 2
    )
  )

sig2.epsilon <-
  do.call(
    rbind,
    lapply(
      2 ^ (-8:3),
      FUN = tune.svm.epsilon,
      data = sprse.trn,
      kerneltype = "sigmoid",
      nfolds = 10,
      seed = 2
    )
  )

#     Seed 3 --------------------------------------------------------------

sig3.cost <-
  do.call(
    rbind,
    lapply(
      2 ^ (-2:6),
      FUN = tune.svm.cost,
      data = sprse.trn,
      kerneltype = "sigmoid",
      nfolds = 10,
      seed = 3
    )
  )

sig3.gamma <-
  do.call(
    rbind,
    lapply(
      2 ^ (-12:1),
      FUN = tune.svm.gamma,
      data = sprse.trn,
      kerneltype = "sigmoid",
      nfolds = 10,
      seed = 3
    )
  )

sig3.epsilon <-
  do.call(
    rbind,
    lapply(
      2 ^ (-8:3),
      FUN = tune.svm.epsilon,
      data = sprse.trn,
      kerneltype = "sigmoid",
      nfolds = 10,
      seed = 3
    )
  )

#     Seed 4 --------------------------------------------------------------

sig4.cost <-
  do.call(
    rbind,
    lapply(
      2 ^ (-2:6),
      FUN = tune.svm.cost,
      data = sprse.trn,
      kerneltype = "sigmoid",
      nfolds = 10,
      seed = 4
    )
  )

sig4.gamma <-
  do.call(
    rbind,
    lapply(
      2 ^ (-12:1),
      FUN = tune.svm.gamma,
      data = sprse.trn,
      kerneltype = "sigmoid",
      nfolds = 10,
      seed = 4
    )
  )

sig4.epsilon <-
  do.call(
    rbind,
    lapply(
      2 ^ (-8:3),
      FUN = tune.svm.epsilon,
      data = sprse.trn,
      kerneltype = "sigmoid",
      nfolds = 10,
      seed = 4
    )
  )

#     Seed 5 --------------------------------------------------------------

sig5.cost <-
  do.call(
    rbind,
    lapply(
      2 ^ (-2:6),
      FUN = tune.svm.cost,
      data = sprse.trn,
      kerneltype = "sigmoid",
      nfolds = 10,
      seed = 5
    )
  )

sig5.gamma <-
  do.call(
    rbind,
    lapply(
      2 ^ (-12:1),
      FUN = tune.svm.gamma,
      data = sprse.trn,
      kerneltype = "sigmoid",
      nfolds = 10,
      seed = 5
    )
  )

sig5.epsilon <-
  do.call(
    rbind,
    lapply(
      2 ^ (-8:3),
      FUN = tune.svm.epsilon,
      data = sprse.trn,
      kerneltype = "sigmoid",
      nfolds = 10,
      seed = 5
    )
  )

#     Seed 6 --------------------------------------------------------------

sig6.cost <-
  do.call(
    rbind,
    lapply(
      2 ^ (-2:6),
      FUN = tune.svm.cost,
      data = sprse.trn,
      kerneltype = "sigmoid",
      nfolds = 10,
      seed = 6
    )
  )

sig6.gamma <-
  do.call(
    rbind,
    lapply(
      2 ^ (-12:1),
      FUN = tune.svm.gamma,
      data = sprse.trn,
      kerneltype = "sigmoid",
      nfolds = 10,
      seed = 6
    )
  )

sig6.epsilon <-
  do.call(
    rbind,
    lapply(
      2 ^ (-8:3),
      FUN = tune.svm.epsilon,
      data = sprse.trn,
      kerneltype = "sigmoid",
      nfolds = 10,
      seed = 6
    )
  )

#     Compilation ---------------------------------------------------------

sig.cost.comp <- rbind(sig.cost, sig2.cost, 
                       sig3.cost, sig4.cost, 
                       sig5.cost, sig6.cost)
sig.gamma.comp <- rbind(sig.gamma, sig2.gamma, 
                       sig3.gamma, sig4.gamma, 
                       sig5.gamma, sig6.gamma)
sig.epsilon.comp <- rbind(sig.epsilon, sig2.epsilon, 
                       sig3.epsilon, sig4.epsilon, 
                       sig5.epsilon, sig6.epsilon)

#     Sigmoid graphs ----------------------------------------------------------
temp.data <- sig.cost.comp
temp.data$seed <- temp.data$seed %>% as.factor()
ggplot(data = temp.data, aes(x = cost, y = rsquared, 
                                 group = seed, color = seed)) + 
  geom_line(size = 1) + 
  scale_x_continuous(trans = 'log2') + 
  theme_bw() + 
  labs(x = "Cost", y = "R-squared", 
       title = "Sigmoid Kernel - Cost", 
       color = "Random Seed")

temp.data <- sig.gamma.comp
temp.data$seed <- temp.data$seed %>% as.factor()
ggplot(data = temp.data, aes(x = gamma, y = rsquared, 
                             group = seed, color = seed)) + 
  geom_line(size = 1) + 
  scale_x_continuous(trans = 'log2') + 
  theme_bw() + 
  labs(x = "Gamma", y = "R-squared", 
       title = "Sigmoid Kernel - Gamma", 
       color = "Random Seed")

temp.data <- sig.epsilon.comp
temp.data$seed <- temp.data$seed %>% as.factor()
ggplot(data = temp.data, aes(x = epsilon, y = rsquared, 
                             group = seed, color = seed)) + 
  geom_line(size = 1) + 
  scale_x_continuous(trans = 'log2') + 
  theme_bw() + 
  labs(x = "Epsilon", y = "R-squared", 
       title = "Sigmoid Kernel - Epsilon", 
       color = "Random Seed")


