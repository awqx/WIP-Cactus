source('tuning.R')

# Y-randomization ---------------------------------------------------------

#     Analysis ------------------------------------------------------------

cube.q2.yr <- function(vars, trial.path, nsplits, cmte, extra, seed) {
  trn.split <- 1:nsplits # Used for first col of saved table
  q2.results <- c(rep(0.0, nsplits))
  
  ctrl <- cubistControl(
    seed = seed, 
    extrapolation = extra
  )
  
  for(n in 1:nsplits) {
    data <- readRDS(paste0(trial.path, 'pp/', n, '/pp.RDS'))
    data <- data %>% select(., info, one_of(vars)) 
    obs <- data$info
    pred <- c() # initializing
    # In this loop, i represents the index of the datapoint left out of 
    # model building
    for(i in 1:nrow(data)) {
      trn <- data[-i, ]
      tst <- data[i, ]
      x <- trn[ , -1]
      y <- trn[ , 1]
      
      cube <- cubist(x = x, y = y,  
                     committees = cmte, 
                     control = ctrl)
      pred[i] <- predict(cube, tst[ , -1]) 
    }
    pred.df <- data.frame(obs, pred)
    colnames(pred.df) <- c("obs", "pred")
    PRESS <- sum((obs - pred)^2)
    TSS <- sum((obs - mean(obs))^2)
    q2 <- 1 - PRESS/TSS
    q2.results[n] <- q2
  }
  results <- data.frame(trn.split, q2.results)
  print(results)
  return(results)
}

rf.q2.yr <- function(vars, trial.path, nsplits, ntree, node, m) {
  trn.split <- 1:nsplits # Used for first col of saved table
  q2.results <- c(rep(0.0, nsplits))
  
  for(n in 1:nsplits) {
    data <- readRDS(paste0(trial.path, 'pp/', n, '/pp.RDS'))
    data <- data %>% select(., info, one_of(vars)) 
    obs <- data$info
    pred <- c() # initializing
    # In this loop, i represents the index of the datapoint left out of 
    # model building
    for(i in 1:nrow(data)) {
      trn <- data[-i, ]
      tst <- data[i, ]
      x <- trn[ , -1]
      y <- trn[ , 1]
      
      rf <- svm(x = x, y = y,
                    ntree = ntree,
                    nodesize = node,
                    mtry = m,  
                    importance = T)
      pred[i] <- predict(rf, tst[ , -1]) 
    }
    pred.df <- data.frame(obs, pred)
    colnames(pred.df) <- c("obs", "pred")
    PRESS <- sum((obs - pred)^2)
    TSS <- sum((obs - mean(obs))^2)
    q2 <- 1 - PRESS/TSS
    q2.results[n] <- q2
  }
  results <- data.frame(trn.split, q2.results)
  print(results)
  return(results)
}

polysvm.q2.yr <- function(vars, trial.path, nsplits, d, cost, e, g, c0) {
  trn.split <- 1:nsplits # Used for first col of saved table
  q2.results <- c(rep(0.0, nsplits))
  
  for(n in 1:nsplits) {
    data <- readRDS(paste0(trial.path, 'pp/', n, '/pp.RDS'))
    data <- data %>% select(., info, one_of(vars)) 
    obs <- data$info
    pred <- c() # initializing
    # In this loop, i represents the index of the datapoint left out of 
    # model building
    for(i in 1:nrow(data)) {
      trn <- data[-i, ]
      tst <- data[i, ]
      x <- trn[ , -1]
      y <- trn[ , 1]
      
      svm.cv <- svm(x = x, y = y,
                    cost = cost, degree = d, coef0 = c0,
                    epsilon = e, gamma = g,
                    kernel = 'polynomial')
      pred[i] <- predict(svm.cv, tst[ , -1]) 
    }
    # Handling outliers
    # for(i in 1:length(pred)) {
    #   if(abs(pred[i]) > 85)
    #     pred[i] <- mean(obs)
    # }
    pred.df <- data.frame(obs, pred)
    colnames(pred.df) <- c("obs", "pred")
    # p <- ggplot(pred.df, aes(x = obs, y = pred)) + 
    #   theme_bw() + 
    #   geom_point() + 
    #   labs(title = n) + 
    #   geom_abline(slope = 1, intercept = 0) + 
    #   coord_fixed()
    # print(p)
    PRESS <- sum((obs - pred)^2)
    TSS <- sum((obs - mean(obs))^2)
    q2 <- 1 - PRESS/TSS
    # message("Q2 = ", q2)
    q2.results[n] <- q2
  }
  results <- data.frame(trn.split, q2.results)
  print(results)
  return(results)
}

rbfsvm.q2.yr <- function(vars, trial.path, nsplits, cost, e, g) {
  trn.split <- 1:nsplits # Used for first col of saved table
  q2.results <- c(rep(0.0, nsplits))
  
  for(n in 1:nsplits) {
    data <- readRDS(paste0(trial.path, 'pp/', n, '/pp.RDS'))
    data <- data %>% select(., info, one_of(vars)) 
    obs <- data$info
    pred <- c() # initializing
    # In this loop, i represents the index of the datapoint left out of 
    # model building
    for(i in 1:nrow(data)) {
      trn <- data[-i, ]
      tst <- data[i, ]
      x <- trn[ , -1]
      y <- trn[ , 1]
      
      svm.cv <- svm(x = x, y = y,
                    cost = cost, epsilon = e, gamma = g,
                    kernel = 'radial')
      pred[i] <- predict(svm.cv, tst[ , -1]) 
    }
    # Handling outliers
    # for(i in 1:length(pred)) {
    #   if(abs(pred[i]) > 85)
    #     pred[i] <- mean(obs)
    # }
    pred.df <- data.frame(obs, pred)
    colnames(pred.df) <- c("obs", "pred")
    # p <- ggplot(pred.df, aes(x = obs, y = pred)) + 
    #   theme_bw() + 
    #   geom_point() + 
    #   labs(title = n) + 
    #   geom_abline(slope = 1, intercept = 0) + 
    #   coord_fixed()
    # print(p)
    PRESS <- sum((obs - pred)^2)
    TSS <- sum((obs - mean(obs))^2)
    q2 <- 1 - PRESS/TSS
    # message("Q2 = ", q2)
    q2.results[n] <- q2
  }
  results <- data.frame(trn.split, q2.results)
  print(results)
  return(results)
}

sigsvm.q2.yr <- function(vars, trial.path, nsplits, cost, e, g, c0) {
  trn.split <- 1:nsplits # Used for first col of saved table
  q2.results <- c(rep(0.0, nsplits))
  
  for(n in 1:nsplits) {
    data <- readRDS(paste0(trial.path, 'pp/', n, '/pp.RDS'))
    data <- data %>% select(., info, one_of(vars)) 
    obs <- data$info
    pred <- c() # initializing
    # In this loop, i represents the index of the datapoint left out of 
    # model building
    for(i in 1:nrow(data)) {
      trn <- data[-i, ]
      tst <- data[i, ]
      x <- trn[ , -1]
      y <- trn[ , 1]
      
      svm.cv <- svm(x = x, y = y,
                    cost = cost, coef0 = c0,
                    epsilon = e, gamma = g,
                    kernel = 'sigmoid')
      pred[i] <- predict(svm.cv, tst[ , -1]) 
    }
    # Handling outliers
    # for(i in 1:length(pred)) {
    #   if(abs(pred[i]) > 85)
    #     pred[i] <- mean(obs)
    # }
    pred.df <- data.frame(obs, pred)
    colnames(pred.df) <- c("obs", "pred")
    # p <- ggplot(pred.df, aes(x = obs, y = pred)) + 
    #   theme_bw() + 
    #   geom_point() + 
    #   labs(title = n) + 
    #   geom_abline(slope = 1, intercept = 0) + 
    #   coord_fixed()
    # print(p)
    PRESS <- sum((obs - pred)^2)
    TSS <- sum((obs - mean(obs))^2)
    q2 <- 1 - PRESS/TSS
    # message("Q2 = ", q2)
    q2.results[n] <- q2
  }
  results <- data.frame(trn.split, q2.results)
  print(results)
  return(results)
}

pls.q2.yr <- function(vars, trial.path, nsplits, method, ncomp) {
  trn.split <- 1:nsplits # Used for first col of saved table
  q2.results <- c(rep(0.0, nsplits))
  
  for(n in 1:nsplits) {
    data <- readRDS(paste0(trial.path, 'pp/', n, '/pp.RDS'))
    data <- data %>% select(., info, one_of(vars)) 
    obs <- data$info
    pred <- c() # initializing
    # In this loop, i represents the index of the datapoint left out of 
    # model building
    for(i in 1:nrow(data)) {
      trn <- data[-i, ]
      tst <- data[i, ]
      
      pls.mod <- plsr(info ~ ., data = trn, 
                      ncomp = ncomp, method = method)
      pred[i] <- predict(pls.mod, tst[ , -1]) %>% .[2]
    }
    # Handling outliers
    # for(i in 1:length(pred)) {
    #   if(abs(pred[i]) > 85)
    #     pred[i] <- mean(obs)
    # }
    pred.df <- data.frame(obs, pred)
    colnames(pred.df) <- c("obs", "pred")
    # p <- ggplot(pred.df, aes(x = obs, y = pred)) + 
    #   theme_bw() + 
    #   geom_point() + 
    #   labs(title = n) + 
    #   geom_abline(slope = 1, intercept = 0) + 
    #   coord_fixed()
    # print(p)
    PRESS <- sum((obs - pred)^2)
    TSS <- sum((obs - mean(obs))^2)
    q2 <- 1 - PRESS/TSS
    # message("Q2 = ", q2)
    q2.results[n] <- q2
  }
  results <- data.frame(trn.split, q2.results)
  print(results)
  return(results)
}

#     Model building ------------------------------------------------------

# Cubist === 
build.cube <- function(host, ntrial, nsplit, seed) {
  trial.path <- paste0('yrand/', host, '/', ntrial, '/')
  # Tuning
  trn <- readRDS(paste0(trial.path, 'pp/1/pp.RDS'))
  colnames(trn) <- str_replace(colnames(trn), "-", ".")
  features <- readRDS(paste0(trial.path, 'vars.RDS'))
  trn <- trn %>% select(., info, one_of(features)) 
    # Ranges
  cmte.range <- c(10, 30, 60, 80, 100)
  extra.range <- c(0, 10, 30, 60, 100)
  
  cube.combos <- expand.grid(cmte.range, extra.range)
  cmte.combos <- cube.combos$Var1
  extra.combos <- cube.combos$Var2
  
  set.seed(seed)
  system.time(
    results.combos <- do.call(
      rbind,
      mapply(
        FUN = tune.cubist,
        cmte= cmte.combos, 
        extra = extra.combos, 
        MoreArgs = 
          list(nfolds = 5, data = trn, seed = seed), 
        SIMPLIFY = F
      )
    )
  )
  temp <<- results.combos
  message('Tuning of trial ', ntrial, ' completed.')
  thecmte <- results.combos[order(results.combos$rsquared, decreasing = T), ] %>% 
    .[1, 'committees']
  theextra <- results.combos[order(results.combos$rsquared, decreasing = T), ] %>% 
    .[1, 'extrapolation']
  message('cmte = ', thecmte, ' || extra = ', theextra)
  
  q2 <- cube.q2.yr(vars = features, trial.path = trial.path, nsplits = 5, 
                   cmte = thecmte, 
                   extra = theextra, seed = seed)
  best.split <- which.max(q2$q2.results)
  saveRDS(q2, paste0(trial.path, 'results/cubist.q2.RDS'))
  
  # Building a Cubist model
  if(is.null(best.split) || length(best.split) == 0) # Just in case
    best.split <- 1
  trn <- readRDS(paste0(trial.path, 'pp/', best.split, '/pp.RDS'))
  colnames(trn) <- str_replace(colnames(trn), "-", ".")
  features <- readRDS(paste0(trial.path, 'vars.RDS'))
  trn <- trn[ , colnames(trn) %in% c("info", features)]
  
  ctrl <- cubistControl(
    seed = seed, 
    extrapolation = theextra
  )
  
  cube <- cubist(x = trn[ , -1], y = trn[ , 1],  
                 committees = thecmte, 
                 control = ctrl)
  
  saveRDS(cube, paste0(trial.path, 'models/cubist.RDS'))
  message('Cubist model of yrand trial ', ntrial, ' completed.')
}

# Random forest ===
build.rf <- function(host, ntrial, nsplit, seed) {
  trial.path <- paste0('yrand/', host, '/', ntrial, '/')
  # Tuning
  # Import the data
  # Loading split 1 is pretty arbitrary
  trn <- readRDS(paste0(trial.path, 'pp/1/pp.RDS'))
  colnames(trn) <- str_replace(colnames(trn), "-", ".")
  features <- readRDS(paste0(trial.path, 'vars.RDS'))
  trn <- trn %>% select(., info, one_of(features)) 
  # Establishing ranges for tuning
  ntree.range <- c(25, 50, 75, 150, 250)
  node.range <- c(1, 3, 5, 8, 15)
  mtry.range <- c(1, 2, 3, 5, 8) %>%
    .[. <= length(features)]
  
  rf.combos <- expand.grid(ntree.range, node.range, mtry.range)
  ntree.combos <- rf.combos$Var1
  node.combos <- rf.combos$Var2
  mtry.combos <- rf.combos$Var3
  
  set.seed(seed)
  system.time(
    results.combos <- do.call(
      rbind,
      mapply(
        FUN = tune.rf,
        ntree = ntree.combos, 
        node = node.combos, 
        m = mtry.combos,
        MoreArgs = 
          list(nfolds = 5, data = trn), 
        SIMPLIFY = F
      )
    )
  )
  
  message('Tuning of trial ', ntrial, ' completed.')
  thentree <- results.combos[order(results.combos$rsquared, decreasing = T), ] %>% 
    .[1, 'ntree'] 
  thenodesize <- results.combos[order(results.combos$rsquared, decreasing = T), ] %>% 
    .[1, 'nodesize']
  themtry <- results.combos[order(results.combos$rsquared, decreasing = T), ] %>% 
    .[1, 'mtry'] 
  message('ntree = ', thentree, ' || node = ', thenodesize, ' || mtry = ', themtry)
  
  q2 <- rf.q2.yr(vars = features, trial.path = trial.path, nsplits = 5, 
                 ntree = thentree, 
                 node = thenodesize, 
                 m = themtry)
  best.split <- which.max(q2$q2.results)
  saveRDS(q2, paste0(trial.path, 'results/rf.q2.RDS'))
  
  # Building and saving a SVM model
  # Importing training data
  if(is.null(best.split) || length(best.split) == 0) # Just in case
    best.split <- 1
  trn <- readRDS(paste0(trial.path, 'pp/', best.split, '/pp.RDS'))
  colnames(trn) <- str_replace(colnames(trn), "-", ".")
  features <- readRDS(paste0(trial.path, 'vars.RDS'))
  trn <- trn[ , colnames(trn) %in% c("info", features)]
  
  rf.mod <- randomForest(x = trn[ , -1], y = trn[ , 1],
                         ntree = thentree, 
                         nodesize = thenodesize, 
                         mtry = themtry) 
  
  saveRDS(rf.mod, paste0(trial.path, 'models/rf.RDS'))
  message('RF model of yrand trial ', ntrial, ' completed.')
}

# SVM ===
build.polysvm <- function(host, ntrial, nsplit, seed) {
  trial.path <- paste0('yrand/', host, '/', ntrial, '/')
  # Tuning
    # Import the data
    # Loading split 1 is pretty arbitrary
  trn <- readRDS(paste0(trial.path, 'pp/1/pp.RDS'))
  colnames(trn) <- str_replace(colnames(trn), "-", ".")
  features <- readRDS(paste0(trial.path, 'vars.RDS'))
  trn <- trn %>% select(., info, one_of(features)) 
    # Establishing ranges for tuning
  coef.range <- 0:4
  cost.range <- c(1, 2, 5, 10)
  deg.range <- 2:4
  epsilon.range <- c(0.1, 0.4, 0.8)
  gamma.range <- c(0.01, 0.1, 0.5)
  
  svm.combos <- expand.grid(cost.range, gamma.range, 
                            epsilon.range, coef.range, deg.range)
  cost.combos <- svm.combos$Var1
  gamma.combos <- svm.combos$Var2
  eps.combos <- svm.combos$Var3
  coef.combos <- svm.combos$Var4
  deg.combos <- svm.combos$Var5
  
  set.seed(seed)
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
  
  message('Tuning of trial ', ntrial, ' completed.')
  thedegree <- results.combos[order(results.combos$rsquared, decreasing = T), ] %>% 
    .[1, 'degree'] %>% as.character()
  thecost <- results.combos[order(results.combos$rsquared, decreasing = T), ] %>% 
    .[1, 'cost'] %>% as.character()
  theepsilon <- results.combos[order(results.combos$rsquared, decreasing = T), ] %>% 
    .[1, 'epsilon'] %>% as.character()
  thegamma <- results.combos[order(results.combos$rsquared, decreasing = T), ] %>% 
    .[1, 'gamma'] %>% as.character()
  thecoef0 <- results.combos[order(results.combos$rsquared, decreasing = T), ] %>% 
    .[1, 'coef0'] %>% as.character()
  message('degree = ', thedegree, ' || cost = ', thecost, 
          ' || eps = ', theepsilon, ' || gamma = ', thegamma, 
          ' || coef0 = ', thecoef0)
  
  q2 <- polysvm.q2.yr(vars = features, trial.path = trial.path, nsplits = 5, 
                      d = thedegree, cost = thecost, e = theepsilon, 
                      g = thegamma, c0 = thecoef0)
  best.split <- which.max(q2$q2.results)
  saveRDS(q2, paste0(trial.path, 'results/svm-poly.q2.RDS'))
  
  # Building and saving a SVM model
  # Importing training data
  if(is.null(best.split) || length(best.split) == 0) # Just in case
    best.split <- 1
  trn <- readRDS(paste0(trial.path, 'pp/', best.split, '/pp.RDS'))
  colnames(trn) <- str_replace(colnames(trn), "-", ".")
  features <- readRDS(paste0(trial.path, 'vars.RDS'))
  trn <- trn[ , colnames(trn) %in% c("info", features)]
  
  svm.mod <- svm(x = trn[ , -1], y = trn[ , 1],
                  degree = thedegree, 
                  cost = thecost, 
                  epsilon = theepsilon, 
                  gamma = thegamma, 
                  coef0 = thecoef0,
                  kernel = 'polynomial') 
  
  saveRDS(svm.mod, paste0(trial.path, 'models/svm-poly.RDS'))
  message('SVM-Poly model of yrand trial ', ntrial, ' completed.')
}
build.rbfsvm <- function(host, ntrial, nsplit, seed) {
  trial.path <- paste0('yrand/', host, '/', ntrial, '/')
  # Tuning
  # Import the data
  # Loading split 1 is pretty arbitrary
  trn <- readRDS(paste0(trial.path, 'pp/1/pp.RDS'))
  colnames(trn) <- str_replace(colnames(trn), "-", ".")
  features <- readRDS(paste0(trial.path, 'vars.RDS'))
  trn <- trn %>% select(., info, one_of(features)) 
  # Establishing ranges for tuning
  cost.range <- c(1, 2, 5, 10, 20)
  epsilon.range <- c(0.01, 0.1, 0.25, 0.5, 0.75)
  gamma.range <- c(0.01, 0.1, 0.25, 0.5)
  
  svm.combos <- expand.grid(cost.range, gamma.range, epsilon.range)
  cost.combos <- svm.combos$Var1
  gamma.combos <- svm.combos$Var2
  eps.combos <- svm.combos$Var3
  
  set.seed(seed)
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
  
  message('Tuning of trial ', ntrial, ' completed.')
  thecost <- results.combos[order(results.combos$rsquared, decreasing = T), ] %>% 
    .[1, 'cost'] %>% as.character()
  theepsilon <- results.combos[order(results.combos$rsquared, decreasing = T), ] %>% 
    .[1, 'epsilon'] %>% as.character()
  thegamma <- results.combos[order(results.combos$rsquared, decreasing = T), ] %>% 
    .[1, 'gamma'] %>% as.character()
  message('cost = ', thecost, ' || eps = ', theepsilon, ' || gamma = ', thegamma)
  
  q2 <- rbfsvm.q2.yr(vars = features, trial.path = trial.path, nsplits = 5, 
                  cost = thecost, e = theepsilon, g = thegamma)
  best.split <- which.max(q2$q2.results)
  saveRDS(q2, paste0(trial.path, 'results/svm-rbf.q2.RDS'))
  
  # Building and saving a SVM model
  # Importing training data
  if(is.null(best.split) || length(best.split) == 0) # Just in case
    best.split <- 1
  trn <- readRDS(paste0(trial.path, 'pp/', best.split, '/pp.RDS'))
  colnames(trn) <- str_replace(colnames(trn), "-", ".")
  features <- readRDS(paste0(trial.path, 'vars.RDS'))
  trn <- trn[ , colnames(trn) %in% c("info", features)]
  
  svm.mod <- svm(x = trn[ , -1], y = trn[ , 1],
                 cost = thecost, 
                 epsilon = theepsilon, 
                 gamma = thegamma, 
                 kernel = 'radial') 
  
  saveRDS(svm.mod, paste0(trial.path, 'models/svm-rbf.RDS'))
  message('SVM-RBF model of yrand trial ', ntrial, ' completed.')
}
build.sigsvm <- function(host, ntrial, nsplit, seed) {
  trial.path <- paste0('yrand/', host, '/', ntrial, '/')
  # Tuning
  # Import the data
  # Loading split 1 is pretty arbitrary
  trn <- readRDS(paste0(trial.path, 'pp/1/pp.RDS'))
  colnames(trn) <- str_replace(colnames(trn), "-", ".")
  features <- readRDS(paste0(trial.path, 'vars.RDS'))
  trn <- trn %>% select(., info, one_of(features)) 
  # Establishing ranges for tuning
  coef.range <- 0:4
  cost.range <- c(1, 2, 5, 10)
  epsilon.range <- c(0.1, 0.4, 0.8)
  gamma.range <- c(0.01, 0.1, 0.5)
  
  svm.combos <- expand.grid(cost.range, gamma.range, 
                            epsilon.range, coef.range)
  cost.combos <- svm.combos$Var1
  gamma.combos <- svm.combos$Var2
  eps.combos <- svm.combos$Var3
  coef.combos <- svm.combos$Var4
  
  set.seed(seed)
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
  
  message('Tuning of trial ', ntrial, ' completed.')
  thecost <- results.combos[order(results.combos$rsquared, decreasing = T), ] %>% 
    .[1, 'cost'] %>% as.character()
  theepsilon <- results.combos[order(results.combos$rsquared, decreasing = T), ] %>% 
    .[1, 'epsilon'] %>% as.character()
  thegamma <- results.combos[order(results.combos$rsquared, decreasing = T), ] %>% 
    .[1, 'gamma'] %>% as.character()
  thecoef0 <- results.combos[order(results.combos$rsquared, decreasing = T), ] %>% 
    .[1, 'coef0'] %>% as.character()
  message('cost = ', thecost, ' || coef0 = ', thecoef0,
          ' || eps = ', theepsilon, ' || gamma = ', thegamma)
  
  q2 <- sigsvm.q2.yr(vars = features, trial.path = trial.path, nsplits = 5, 
                     cost = thecost, e = theepsilon, 
                     g = thegamma, c0 = thecoef0)
  best.split <- which.max(q2$q2.results)
  saveRDS(q2, paste0(trial.path, 'results/svm-sig.q2.RDS'))
  
  # Building and saving a SVM model
  # Importing training data
  if(is.null(best.split) || length(best.split) == 0) # Just in case
    best.split <- 1
  trn <- readRDS(paste0(trial.path, 'pp/', best.split, '/pp.RDS'))
  colnames(trn) <- str_replace(colnames(trn), "-", ".")
  features <- readRDS(paste0(trial.path, 'vars.RDS'))
  trn <- trn[ , colnames(trn) %in% c("info", features)]
  
  svm.mod <- svm(x = trn[ , -1], y = trn[ , 1],
                 cost = thecost, 
                 epsilon = theepsilon, 
                 gamma = thegamma, 
                 coef0 = thecoef0,
                 kernel = 'sigmoid') 
  
  saveRDS(svm.mod, paste0(trial.path, 'models/svm-sig.RDS'))
  message('SVM-sig model of yrand trial ', ntrial, ' completed.')
}
# PLS ===
# Modling functions only used for yrand
# Tunes and builds a MARS model
# ntrial: the number of the trial (not total)
# nsplit: total number of splits
build.pls <- function(host, ntrial, nsplit, seed) {
  trial.path <- paste0('yrand/', host, '/', ntrial, '/')
  # Tuning
  # Import the data
      # Loading split 1 is pretty arbitrary
  trn <- readRDS(paste0(trial.path, 'pp/1/pp.RDS'))
  colnames(trn) <- str_replace(colnames(trn), "-", ".")
  
  features <- readRDS(paste0(trial.path, 'vars.RDS'))
  trn <- trn %>% select(., info, one_of(features)) 
  
  # ncomp.range <- c(1, round(ncol(trn) * 0.1), round(ncol(trn) * 0.25), 
  #                  round(ncol(trn) * 0.5), round(ncol(trn) * 0.85), 
  #                  ncol(trn)) %>%
  #   unique() # in case variables are few
  # ncomp.range <- unique(ncomp.range) 
  ncomp.range <- c(2:5, 8, 15)
  if(length(features) > 1)
    ncomp.range <- ncomp.range[ncomp.range < length(features)] %>%
      .[. > 0] # Filtering for invalid ncomp values
  else
    ncomp.range <- c(1)
  method.range <- c("kernelpls", "simpls", "oscorespls", "widekernelpls")
  
  combos <- expand.grid(method.range, 
                        ncomp.range)
  method.combos <- combos$Var1 %>% as.character()
  ncomp.combos <- combos$Var2
  
  # Creating a table of results of tuning
  set.seed(seed)
  results.combos <- do.call(
    rbind,
    mapply(
      FUN = tune.pls,
      comp = ncomp.combos,
      met = method.combos,
      MoreArgs =
        list(nfolds = 10, data = trn),
      SIMPLIFY = F
    )
  )
  message('Tuning of trial ', ntrial, ' completed.')
  # Finding the tuning parameters that correspond to highest R2
  if(is.null(results.combos)) { # in case there is only one feature
    themethod <- 'oscorespls'
    thencomp <- 1
  } else {
    themethod <- results.combos[order(results.combos$rsquared, decreasing = T), ] %>% 
      .[1, 'method'] %>% as.character()
    thencomp <- results.combos[order(results.combos$rsquared, decreasing = T), ] %>% 
      .[1, 'ncomp'] 
  }
  print(themethod)
  print(thencomp)
  
  # best.split <- 0 # initializing this variable, in case q2.results are NA
  q2 <- pls.q2.yr(vars = features, trial.path = trial.path, 
                          nsplits = 5, ncomp = thencomp, method = themethod)
  best.split <- which.max(q2$q2.results)
  saveRDS(q2, paste0(trial.path, 'results/pls.q2.RDS'))
  
  # Building and saving a PLS model
  # Importing training data
  if(is.null(best.split) || length(best.split) == 0)
    best.split <- 1
  trn <- readRDS(paste0(trial.path, 'pp/', best.split, '/pp.RDS'))
  colnames(trn) <- str_replace(colnames(trn), "-", ".")
  
  features <- readRDS(paste0(trial.path, 'vars.RDS'))
  trn <- trn[ , colnames(trn) %in% c("info", features)]
  
  pls.mod <- plsr(info~., data = trn, 
                  ncomp = thencomp, method = themethod) 
  
  saveRDS(pls.mod, paste0(trial.path, 'models/pls.RDS'))
  message('PLS model of yrand trial ', ntrial, ' completed.')
}
