#     Sigmoid kernel ----------------------------------------------------------
# Cost = 2
results.sig.cost <-
  do.call(
    rbind,
    lapply(
      2 ^ (-6:6),
      FUN = tune.svm.cost,
      data = sprse.trn,
      kerneltype = "sigmoid",
      nfolds = 10,
      seed = 5
    )
  ) 

# Peaks around 0.0005
results.sig.gamma <-
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

# Epsilon = 0.25
results.sig.epsilon <-
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

# Peaks around 0.004
results.sig.ga.cost <-
  do.call(
    rbind,
    lapply(
      2 ^ (-14:-4),
      FUN = tune.svm.gamma,
      data = sprse.ga.trn,
      kerneltype = "sigmoid",
      nfolds = 10,
      seed = 5
    )
  )

#     Sigmoid graphs ----------------------------------------------------------
ggplot(data = results.sig.cost, aes(x = cost, y = rsquared)) + 
  geom_line() + 
  scale_x_continuous(trans = 'log2') + 
  theme_bw()

ggplot(data = results.sig.gamma, aes(x = gamma, y = rsquared)) + 
  geom_line() + 
  scale_x_continuous(trans = 'log2') + 
  theme_bw()

ggplot(data = results.sig.epsilon, aes(x = epsilon, y = rsquared)) + 
  geom_line() + 
  scale_x_continuous(trans = 'log2') + 
  theme_bw()

ggplot(data = results.sig.ga.cost, aes(x = gamma, y = rsquared)) + 
  geom_line() + 
  scale_x_continuous(trans = 'log2') + 
  theme_bw()


