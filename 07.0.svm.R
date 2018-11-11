dir.create("models")
dir.create("models/alpha")
dir.create("models/beta")
dir.create("models/gamma")
dir.create("results")
dir.create("results/alpha")
dir.create("results/beta")
dir.create("results/gamma")

# library(Matrix)
# library(stats)
# library(stringr)
source("07.model.functions.R")
p_load(e1071, kernlab)

# Functions ---------------------------------------------------------------

# For building the finished model
# Returns the results of the model on the test splits as a data.frame
svm.tst.splits <- function(pp.dir, tst.dir, feat, nsplits, model) {
  tst.all <- data.frame()
  for (j in 1:nsplits) {
    tst <- preprocess.tst.mod(pp.dir = pp.dir, tst.dir = tst.dir, 
                              feat = feat, n = j)
    tst.y <- tst[ , 1]
    tst.x <- tst[ , -1]
    
    tst.df <- predict(model, tst.x) %>%
      cbind(tst.y, .) %>% as.data.frame() 
    colnames(tst.df) <- c("obs", "pred")
    tst.df <- tst.df %>% mutate(split = j)
    
    for(k in 1:nrow(tst.df))
      if(abs(tst.df$pred[k]) > 100)
        tst.df$pred[k] <- mean(trn.y)
    
    tst.all <- rbind(tst.all, tst.df)
  }
  tst.all$split <- as.factor(tst.all$split)
  return(tst.all)
}
  
  # Polynomial -----
polysvm.looq2 <- function(read.dir, nsplits, cost, deg, coef, e, g) {
  # initialize a vector for analysis
  trn.split <- 1:nsplits
  q2.results <- c(rep(0.0, nsplits))
  if(str_detect(read.dir, "alpha"))
    features <- readRDS("./feature.selection/alpha.vars.RDS")
  else
    features <- readRDS("./feature.selection/beta.vars.RDS")
  for(n in 1:nsplits) {
    data <- readRDS(paste0(read.dir, n, "/pp.RDS")) %>%
      select(., -guest)
    obs <- data[ , 1]
    # Using dplyr::one_of for redundancy (in case pre-processing 
    # accidentally dropped a couple variables)
    data <- data %>% select(., DelG, one_of(features))
    pred <- c(rep(0.0, nrow(data) - 1))
    # In this loop, i represents the index of the datapoint left out of 
    # model building
    for(i in 1:nrow(data)) {
      trn <- data[-i, ]
      tst <- data[i, ]
      x <- trn[ , -1]
      y <- trn[ , 1]
    
      # Refer to 07.0.1.svm.tune.poly.R for values
      svm.cv <- svm(x = x, y = y,
                    cost = cost, degree = deg, coef0 = coef,
                    epsilon = e, gamma = g,
                    kernel = "polynomial")
      pred[i] <- predict(svm.cv, tst[ , -1]) 
    }
    # test.df <- data.frame(obs, pred) %>% print()
    # test.df <<- test.df
    # defaultSummary(test.df) %>% print()
    # prediction error sum of squares
    for(i in 1:length(pred)) {
      if(abs(pred[i]) > 50)
        pred[i] <- mean(obs)
    }
    pred.df <- data.frame(obs, pred)
    p <- ggplot(pred.df, aes(x = obs, y = pred)) + 
      theme_bw() + 
      geom_point() + 
      labs(title = n) + 
      geom_abline(slope = 1, intercept = 0)
    print(p)
    PRESS <- sum((obs - pred)^2)
    # total sum of squares
    TSS <- sum((obs - mean(obs))^2)
    q2 <- 1 - PRESS/TSS
    message("Q2 = ", q2)
    q2.results[n] <- q2
  }
  return(data.frame(trn.split, q2.results))
}

# pp.dir, tst.dir should end in a backslash
polysvm.tst <- function(pp.dir, tst.dir, nsplits, cost, deg, coef, e, g) {
  split <- 1:nsplits
  # r2.results <- c(rep(0.0, nsplits))
  # rmse.results <- c(rep(0.0, nsplits))
  
  if(str_detect(pp.dir, "alpha"))
    features <- readRDS("./feature.selection/alpha.vars.RDS")
  else
    features <- readRDS("./feature.selection/beta.vars.RDS")
  
  results.all <- data.frame() 
  
  for(i in 1:nsplits) {
    trn <- readRDS(paste0(pp.dir, i, "/pp.RDS")) %>%
      select(., -guest)
    trn.y <- trn[ , 1]
    trn.x <- trn %>% select(., one_of(features))
    
    # Refer to 07.0.1.svm.tune.poly.R for values
    svm.cv <- svm(x = trn.x, y = trn.y,
                  cost = cost, degree = deg, coef0 = coef,
                  epsilon = e, gamma = g,
                  kernel = "polynomial")
    
    tst.all <- data.frame()
    
    # Testing on each test split
    for(j in 1:nsplits) {
      tst <- preprocess.tst.mod(pp.dir = pp.dir, tst.dir = tst.dir, 
                                feat = colnames(trn.x), n = j)
      tst.y <- tst[ , 1]
      tst.x <- tst[ , -1]
      
      tst.df <- predict(svm.cv, tst.x) %>%
        cbind(tst.y, .) %>% as.data.frame() 
      colnames(tst.df) <- c("obs", "pred")
      tst.df <- tst.df %>% mutate(split = j)
      # Controlling for extreme outliers
      for(k in 1:nrow(tst.df))
        if(abs(tst.df$pred[k]) > 100)
          tst.df$pred[k] <- mean(trn.y)
      
      results.tst <- data.frame(trn.split = i, tst.split = j, 
                                r2 = defaultSummary(tst.df)[2], 
                                rmse = defaultSummary(tst.df)[1])
      results.all <- rbind(results.all, results.tst)
      tst.all <- rbind(tst.all, tst.df)
    }
    
    tst.all$split <- as.factor(tst.all$split)
    p <- ggplot(tst.all, aes(x = obs, y = pred, color = split)) + 
      theme_bw() + 
      geom_point() + 
      labs(title = i) + 
      geom_abline(slope = 1, intercept = 0) + 
      coord_fixed()
    print(p)
  }
  row.names(results.all) <- NULL
  return(results.all)
}

  # RBF -----
rbfsvm.looq2 <- function(read.dir, nsplits, cost, e, g) {
  # initialize a vector for analysis
  trn.split <- 1:nsplits
  q2.results <- c(rep(0.0, nsplits))
  if(str_detect(read.dir, "alpha"))
    features <- readRDS("./feature.selection/alpha.vars.RDS")
  else
    features <- readRDS("./feature.selection/beta.vars.RDS")
  for(n in 1:nsplits) {
    data <- readRDS(paste0(read.dir, n, "/pp.RDS")) %>%
      select(., -guest)
    obs <- data[ , 1]
    # Using dplyr::one_of for redundancy (in case pre-processing 
    # accidentally dropped a couple variables)
    data <- data %>% select(., DelG, one_of(features))
    pred <- c(rep(0.0, nrow(data) - 1))
    # In this loop, i represents the index of the datapoint left out of 
    # model building
    for(i in 1:nrow(data)) {
      trn <- data[-i, ]
      tst <- data[i, ]
      x <- trn[ , -1]
      y <- trn[ , 1]
      
      # Refer to 07.0.1.svm.tune.poly.R for values
      svm.cv <- svm(x = x, y = y,
                    cost = cost, epsilon = e, gamma = g,
                    kernel = "radial")
      pred[i] <- predict(svm.cv, tst[ , -1]) 
    }
    # test.df <- data.frame(obs, pred) %>% print()
    # test.df <<- test.df
    # defaultSummary(test.df) %>% print()
    # prediction error sum of squares
    for(i in 1:length(pred)) {
      if(abs(pred[i]) > 50)
        pred[i] <- mean(obs)
    }
    pred.df <- data.frame(obs, pred)
    p <- ggplot(pred.df, aes(x = obs, y = pred)) + 
      theme_bw() + 
      geom_point() + 
      labs(title = n) + 
      geom_abline(slope = 1, intercept = 0)
    print(p)
    PRESS <- sum((obs - pred)^2)
    # total sum of squares
    TSS <- sum((obs - mean(obs))^2)
    q2 <- 1 - PRESS/TSS
    message("Q2 = ", q2)
    q2.results[n] <- q2
  }
  return(data.frame(trn.split, q2.results))
}

# pp.dir, tst.dir should end in a backslash
rbfsvm.tst <- function(pp.dir, tst.dir, nsplits, cost, e, g) {
  split <- 1:nsplits
  # r2.results <- c(rep(0.0, nsplits))
  # rmse.results <- c(rep(0.0, nsplits))
  
  if(str_detect(pp.dir, "alpha"))
    features <- readRDS("./feature.selection/alpha.vars.RDS")
  else
    features <- readRDS("./feature.selection/beta.vars.RDS")
  
  results.all <- data.frame() 
  
  for(i in 1:nsplits) {
    trn <- readRDS(paste0(pp.dir, i, "/pp.RDS")) %>%
      select(., -guest)
    trn.y <- trn[ , 1]
    trn.x <- trn %>% select(., one_of(features))
    
    # Refer to 07.0.1.svm.tune.poly.R for values
    svm.cv <- svm(x = trn.x, y = trn.y,
                  cost = cost, epsilon = e, gamma = g,
                  kernel = "radial")
    
    tst.all <- data.frame()
    
    # Testing on each test split
    for(j in 1:nsplits) {
      tst <- preprocess.tst.mod(pp.dir = pp.dir, tst.dir = tst.dir, 
                                feat = colnames(trn.x), n = j)
      tst.y <- tst[ , 1]
      tst.x <- tst[ , -1]
      
      tst.df <- predict(svm.cv, tst.x) %>%
        cbind(tst.y, .) %>% as.data.frame() 
      colnames(tst.df) <- c("obs", "pred")
      tst.df <- tst.df %>% mutate(split = j)
      # Controlling for extreme outliers
      for(k in 1:nrow(tst.df))
        if(abs(tst.df$pred[k]) > 100)
          tst.df$pred[k] <- mean(trn.y)
      
      results.tst <- data.frame(trn.split = i, tst.split = j, 
                                r2 = defaultSummary(tst.df)[2], 
                                rmse = defaultSummary(tst.df)[1])
      results.all <- rbind(results.all, results.tst)
      tst.all <- rbind(tst.all, tst.df)
    }
    
    tst.all$split <- as.factor(tst.all$split)
    p <- ggplot(tst.all, aes(x = obs, y = pred, color = split)) + 
      theme_bw() + 
      geom_point() + 
      labs(title = i) + 
      geom_abline(slope = 1, intercept = 0) + 
      coord_fixed()
    print(p)
  }
  row.names(results.all) <- NULL
  return(results.all)
}

  # Sigmoid ----
sigsvm.looq2 <- function(read.dir, nsplits, cost, coef, e, g) {
  trn.split <- 1:nsplits
  q2.results <- c(rep(0.0, nsplits))
  if(str_detect(read.dir, "alpha"))
    features <- readRDS("./feature.selection/alpha.vars.RDS")
  else
    features <- readRDS("./feature.selection/beta.vars.RDS")
  for(n in 1:nsplits) {
    data <- readRDS(paste0(read.dir, n, "/pp.RDS")) %>%
      select(., -guest)
    obs <- data[ , 1]
    data <- data %>% select(., DelG, one_of(features))
    pred <- c(rep(0.0, nrow(data) - 1))
    for(i in 1:nrow(data)) {
      trn <- data[-i, ]
      tst <- data[i, ]
      x <- trn[ , -1]
      y <- trn[ , 1]

      svm.cv <- svm(x = x, y = y,
                    cost = cost, coef0 = coef,
                    epsilon = e, gamma = g,
                    kernel = "sigmoid")
      pred[i] <- predict(svm.cv, tst[ , -1]) 
    }
    for(i in 1:length(pred)) {
      if(abs(pred[i]) > 50)
        pred[i] <- mean(obs)
    }
    pred.df <- data.frame(obs, pred)
    p <- ggplot(pred.df, aes(x = obs, y = pred)) + 
      theme_bw() + 
      geom_point() + 
      labs(title = n) + 
      geom_abline(slope = 1, intercept = 0)
    print(p)
    PRESS <- sum((obs - pred)^2)
    TSS <- sum((obs - mean(obs))^2)
    q2 <- 1 - PRESS/TSS
    message("Q2 = ", q2)
    q2.results[n] <- q2
  }
  return(data.frame(trn.split, q2.results))
}
# pp.dir, tst.dir should end in a backslash
sigsvm.tst <- function(pp.dir, tst.dir, nsplits, cost, coef, e, g) {
  split <- 1:nsplits
  # r2.results <- c(rep(0.0, nsplits))
  # rmse.results <- c(rep(0.0, nsplits))
  
  if(str_detect(pp.dir, "alpha"))
    features <- readRDS("./feature.selection/alpha.vars.RDS")
  else
    features <- readRDS("./feature.selection/beta.vars.RDS")
  
  results.all <- data.frame() 
  
  for(i in 1:nsplits) {
    trn <- readRDS(paste0(pp.dir, i, "/pp.RDS")) %>%
      select(., -guest)
    trn.y <- trn[ , 1]
    trn.x <- trn %>% select(., one_of(features))
    
    # Refer to 07.0.1.svm.tune.poly.R for values
    svm.cv <- svm(x = trn.x, y = trn.y,
                  cost = cost, coef0 = coef,
                  epsilon = e, gamma = g,
                  kernel = "sigmoid")
    
    tst.all <- data.frame()
    
    # Testing on each test split
    for(j in 1:nsplits) {
      tst <- preprocess.tst.mod(pp.dir = pp.dir, tst.dir = tst.dir, 
                                feat = colnames(trn.x), n = j)
      tst.y <- tst[ , 1]
      tst.x <- tst[ , -1]
      
      tst.df <- predict(svm.cv, tst.x) %>%
        cbind(tst.y, .) %>% as.data.frame() 
      colnames(tst.df) <- c("obs", "pred")
      tst.df <- tst.df %>% mutate(split = j)
      # Controlling for extreme outliers
      for(k in 1:nrow(tst.df))
        if(abs(tst.df$pred[k]) > 100)
          tst.df$pred[k] <- mean(trn.y)
      
      results.tst <- data.frame(trn.split = i, tst.split = j, 
                                r2 = defaultSummary(tst.df)[2], 
                                rmse = defaultSummary(tst.df)[1])
      results.all <- rbind(results.all, results.tst)
      tst.all <- rbind(tst.all, tst.df)
    }
    
    tst.all$split <- as.factor(tst.all$split)
    p <- ggplot(tst.all, aes(x = obs, y = pred, color = split)) + 
      theme_bw() + 
      geom_point() + 
      labs(title = i) + 
      geom_abline(slope = 1, intercept = 0) + 
      coord_fixed()
    print(p)
  }
  row.names(results.all) <- NULL
  return(results.all)
}

#####
# Polynomial --------------------------------------------------------------
#   Alpha -----------------------------------------------------------------

#     LOOCV-Q2 -----

# 1, 6, 8, and 9
alpha.q2 <- polysvm.looq2(read.dir = "./pre-process/alpha/", nsplits = 10, 
                          deg = 2, cost = 5, e = 0.25, 
                          g = 0.05, coef = 5) %>% print()
saveRDS(alpha.q2, 'results/alpha/polysvm.q2.RDS')

#     Test sets -----

# Decided on split 6
alpha.tst <- polysvm.tst("./pre-process/alpha/", "./model.data/alpha/", 
                         nsplits = 10,
                         deg = 2, cost = 5, e = 0.25, 
                         g = 0.05, coef = 5)
# This is the data that would have been created with only a single
# train-test split
# You can see, no. 2 is very optimistic
alpha.1to1 <- alpha.tst %>% filter(trn.split == tst.split) %>% print()
# Taking the averages
# (Here, No 2, while still good, is much more realistic)
alpha.avg <- data.table(alpha.tst, key = 'trn.split')
alpha.avg <- alpha.avg[ , list(r2 = mean(r2), 
                               rmse = mean(rmse)), 
                        by = trn.split] %>% print()

#     Single model -----

trn.alpha <- readRDS("./pre-process/alpha/6/pp.RDS") %>%
  select(., -guest)
features <- readRDS("./feature.selection/alpha.vars.RDS")
colnames(trn.alpha) <- str_replace(colnames(trn.alpha), "-", ".")
trn.alpha <- trn.alpha %>% select(., DelG, features) 
trn.alpha.x <- select(trn.alpha, -DelG) 
trn.alpha.y <- select(trn.alpha, DelG) 

polysvm.alpha <- svm(x = trn.alpha.x, y = trn.alpha.y,
                     degree = 2, cost = 5, epsilon = 0.25, 
                     gamma = 0.05, coef0 = 5,
                     kernel = "polynomial")
tst.alpha.df <- svm.tst.splits("pre-process/alpha/", "model.data/alpha/", 
                               features, 10, polysvm.alpha)

# Passes all criteria
eval.tropsha(tst.alpha.df)

# Graph of results on all test splits
graph.alpha <- ggplot(tst.alpha.df, aes(x = obs, y = pred, color = split)) + 
  geom_point() + 
  theme_bw() + 
  coord_fixed()  + 
  geom_abline(intercept = 0, slope = 1) + 
  labs(x = "Experimental dG, kJ/mol", 
       y = "Predicted dG, kJ/mol", 
       color = "Test split",
       title = "Alpha-CD Polynomial SVM") 
print(graph.alpha)

#     Saving results -----
# Apprpriate pre-processing settings
pp.settings <- readRDS("./pre-process/beta/6/pp.settings.RDS")
saveRDS(list(pp.settings, polysvm.alpha), "./models/alpha/polysvm.RDS")
# Results from all the test splits
saveRDS(tst.alpha.df, "./results/alpha/polysvm.all.RDS")
# Results from the relevant split
saveRDS(tst.alpha.df %>% filter(split == "6"),
        "./results/alpha/polysvm.RDS")
graph.alpha
ggsave("./results/alpha/polysvm.png")

#   Beta ------------------------------------------------------------------

#     LOO-CV -----
# 1, 2, 5, 9, 10
beta.q2 <- polysvm.looq2(read.dir = "./pre-process/beta/", nsplits = 10, 
                         deg = 4, cost = 25, e = 0.75, 
                         g = 0.5, coef = 10) %>% print()
saveRDS(beta.q2, 'results/beta/polysvm.q2.RDS')

#     Test sets -----
# SPLIT 2
beta.tst <- polysvm.tst("./pre-process/beta/", "./model.data/beta/", 
                        nsplits = 10,
                        deg = 4, cost = 25, e = 0.75, 
                        g = 0.5, coef = 10) 
# Again, 1-to-1, split 2 does very well, possibly as a result of
# a fortuitous data split
beta.1to1 <- beta.tst %>% filter(trn.split == tst.split)
beta.avg <- data.table(beta.tst, key = 'trn.split')
beta.avg <- beta.avg[ , list(r2 = mean(r2), 
                             rmse = mean(rmse)), 
                      by = trn.split] %>% print()

#     Model -----

trn.beta <- readRDS("./pre-process/beta/2/pp.RDS") %>%
  select(., -guest)
features <- readRDS("./feature.selection/beta.vars.RDS")
colnames(trn.beta) <- str_replace(colnames(trn.beta), "-", ".")
trn.beta <- trn.beta %>% select(., DelG, features) 
trn.beta.x <- select(trn.beta, -DelG) 
trn.beta.y <- select(trn.beta, DelG) 

polysvm.beta <- svm(x = trn.beta.x, y = trn.beta.y,
                    deg = 4, cost = 25, epsilon = 0.75, 
                    gamma = 0.5, coef0 = 10,
                    kernel = "polynomial")

tst.beta.df <- svm.tst.splits("pre-process/beta/", "model.data/beta/", 
                              features, 10, polysvm.beta)

# Yay, you pass
eval.tropsha(tst.beta.df)
graph.beta <- ggplot(tst.beta.df, aes(x = obs, y = pred, color = split)) + 
  geom_point() + 
  theme_bw() + 
  coord_fixed()  + 
  geom_abline(intercept = 0, slope = 1) + 
  labs(x = "Experimental dG, kJ/mol", 
       y = "Predicted dG, kJ/mol", 
       color = "Test split", 
       title = "Beta-CD Polynomial SVM")
print(graph.beta)

#     Save -----

pp.settings <- readRDS("./pre-process/beta/2/pp.settings.RDS")
saveRDS(list(pp.settings, polysvm.beta), "./models/beta/polysvm.RDS")
saveRDS(tst.beta.df, "./results/beta/polysvm.all.RDS")
saveRDS(tst.beta.df %>% filter(split == "2"), 
        "./results/beta/polysvm.RDS")
graph.beta 
ggsave("./results/beta/polysvm.png")

#   Gamma -----------------------------------------------------------------

#    LOO-CV -----
# As expected, none pass
gamma.q2 <- polysvm.looq2(read.dir = "./pre-process/gamma/", nsplits = 10, 
                          deg = 3, cost =3, e = 0.01, 
                          g = 0.01, coef = 2) %>% print()
saveRDS(gamma.q2, 'results/gamma/polysvm.q2.RDS')

#     Test -----

# SPLIT 6
gamma.tst <- polysvm.tst("./pre-process/gamma/", "./model.data/gamma/", 
                        nsplits = 10,
                        deg = 4, cost = 25, e = 0.75, 
                        g = 0.5, coef = 10) 
gamma.1to1 <- gamma.tst %>% filter(trn.split == tst.split)
gamma.avg <- data.table(gamma.tst, key = 'trn.split')
gamma.avg <- gamma.avg[ , list(r2 = mean(r2), 
                             rmse = mean(rmse)), 
                      by = trn.split] %>% print()

#     Model -----

trn.gamma <- readRDS("./pre-process/gamma/6/pp.RDS") %>%
  select(., -guest)
features <- readRDS("./feature.selection/gamma.vars.RDS")
colnames(trn.gamma) <- str_replace(colnames(trn.gamma), "-", ".")
trn.gamma <- trn.gamma %>% select(., DelG, features) 
trn.gamma.x <- select(trn.gamma, -DelG) 
trn.gamma.y <- select(trn.gamma, DelG) 

polysvm.gamma <- svm(x = trn.gamma.x, y = trn.gamma.y,
                     degree = 3, cost =3, epsilon = 0.01, 
                     gamma = 0.01, coef0 = 2,
                     kernel = "polynomial")

tst.gamma.df <- svm.tst.splits("pre-process/gamma/", "model.data/gamma/", 
                              features, 10, polysvm.gamma)

# Does not pass criteria (already failed Q2 test anyway)
eval.tropsha(tst.gamma.df)
graph.gamma <- ggplot(tst.gamma.df, aes(x = obs, y = pred, color = split)) + 
  geom_point() + 
  theme_bw() + 
  coord_fixed()  + 
  geom_abline(intercept = 0, slope = 1) + 
  labs(x = "Experimental dG, kJ/mol", 
       y = "Predicted dG, kJ/mol", 
       color = "Test split", 
       title = "Gamma-CD Polynomial SVM")
print(graph.gamma)

#     Save -----

pp.settings <- readRDS("./pre-process/gamma/6/pp.settings.RDS")
saveRDS(list(pp.settings, polysvm.gamma), "./models/gamma/polysvm.RDS")
saveRDS(tst.gamma.df, "./results/gamma/polysvm.all.RDS")
saveRDS(tst.gamma.df %>% filter(split == "6"), 
        "./results/gamma/polysvm.RDS")
graph.gamma 
ggsave("./results/gamma/polysvm.png")

#####
# Radial ------------------------------------------------------------------
#   Alpha-CD -----  
#     LOOCV-Q2 -----
# All pass
alpha.rbf.q2 <- rbfsvm.looq2(read.dir = "./pre-process/alpha/", 
                             nsplits = 10, cost = 3, 
                             e = .1, g = .05) %>% print()
saveRDS(alpha.rbf.q2, 'results/alpha/rbfsvm.q2.RDS')

#     Test ------
# Split 10
alpha.tst <- rbfsvm.tst("./pre-process/alpha/", "./model.data/alpha/", 
                        nsplits = 10, cost = 3, 
                        e = 0.1, g = 0.05) 
alpha.1to1 <- alpha.tst %>% filter(trn.split == tst.split)
alpha.avg <- data.table(alpha.tst, key = 'trn.split')
alpha.avg <- alpha.avg[ , list(r2 = mean(r2), 
                               rmse = mean(rmse)), 
                        by = trn.split] %>% print()

#     Model -----

trn.alpha <- readRDS("./pre-process/alpha/10/pp.RDS") %>%
  select(., -guest)
features <- readRDS("./feature.selection/alpha.vars.RDS")
colnames(trn.alpha) <- str_replace(colnames(trn.alpha), "-", ".")
trn.alpha <- trn.alpha %>% select(., DelG, features) 
trn.alpha.x <- select(trn.alpha, -DelG) 
trn.alpha.y <- select(trn.alpha, DelG) 

rbfsvm.alpha <- svm(x = trn.alpha.x, y = trn.alpha.y,
                    cost = 3, epsilon = 0.1, gamma = 0.05,
                    kernel = "radial")
tst.alpha.df <- svm.tst.splits("pre-process/alpha/", "model.data/alpha/", 
                               features, 10, rbfsvm.alpha)

# Passes evaluation
eval.tropsha(tst.alpha.df)
graph.alpha <- ggplot(tst.alpha.df, aes(x = obs, y = pred, color = split)) + 
  geom_point() + 
  theme_bw() + 
  coord_fixed()  + 
  geom_abline(intercept = 0, slope = 1) + 
  labs(x = "Experimental dG, kJ/mol", y = "Predicted dG, kJ/mol", 
       title = "Alpha-CD Radial SVM", color = "Test split")
graph.alpha

#     Saving ----

pp.settings <- readRDS("./pre-process/alpha/10/pp.settings.RDS")
saveRDS(list(pp.settings, rbfsvm.alpha), "./models/alpha/rbfsvm.RDS")
saveRDS(tst.alpha.df, "./results/alpha/rbfsvm.all.RDS")
saveRDS(tst.alpha.df %>% filter(split == "10"), 
        "./results/alpha/rbfsvm.RDS")
graph.alpha 
ggsave("./results/alpha/rbfsvm.png")

#   Beta-CD -----
#     LOO-CV -----

# All pass
beta.rbf.q2 <- rbfsvm.looq2(read.dir = "./pre-process/beta/", nsplits = 10, 
                            cost = 10, e = .1, g = 0.05)
saveRDS(beta.rbf.q2, 'results/beta/rbfsvm.q2.RDS')

#     Test -----
# Split 7
beta.tst <- rbfsvm.tst("./pre-process/beta/", "./model.data/beta/", 
                       nsplits = 10, cost = 10, e = 0.1, g = 0.05)
beta.1to1 <- beta.tst %>% filter(trn.split == tst.split)
beta.avg <- data.table(beta.tst, key = 'trn.split')
beta.avg <- beta.avg[ , list(r2 = mean(r2), 
                             rmse = mean(rmse)), 
                      by = trn.split] %>% print()

#     Model ----

trn.beta <- readRDS("./pre-process/beta/7/pp.RDS") %>%
  select(., -guest)
features <- readRDS("./feature.selection/beta.vars.RDS")
colnames(trn.beta) <- str_replace(colnames(trn.beta), "-", ".")
trn.beta <- trn.beta %>% select(., DelG, features) 
trn.beta.x <- select(trn.beta, -DelG) 
trn.beta.y <- select(trn.beta, DelG) 

rbfsvm.beta <- svm(x = trn.beta.x, y = trn.beta.y,
                   cost = 10, epsilon = 0.1, gamma = 0.05,
                   kernel = "radial")
tst.beta.df <- svm.tst.splits("pre-process/beta/", "model.data/beta/", 
                              features, 10, rbfsvm.beta)

# All pass
eval.tropsha(tst.beta.df)
graph.beta <- ggplot(tst.beta.df, aes(x = obs, y = pred, color = split)) + 
  geom_point() + 
  theme_bw() + 
  coord_fixed()  + 
  geom_abline(intercept = 0, slope = 1) + 
  labs(x = "Experimental dG, kJ/mol", y = "Predicted dG, kJ/mol", 
       title = "Beta-CD Polynomial SVM", color = "Test split")
graph.beta

#     Save ----

pp.settings <- readRDS("./pre-process/beta/7/pp.settings.RDS")
saveRDS(list(pp.settings, rbfsvm.beta), "./models/beta/rbfsvm.RDS")
saveRDS(tst.beta.df, "./results/beta/rbfsvm.all.RDS")
saveRDS(tst.beta.df %>% filter(split == "7"), 
        "./results/beta/rbfsvm.RDS")
graph.beta 
ggsave("./results/beta/rbfsvm.png")

#   Gamma-CD -----
#     LOO-CV -----
gamma.rbf.q2 <- rbfsvm.looq2(read.dir = "./pre-process/gamma/", nsplits = 10, 
                            cost = 3, e = .15, g = 0.01)
saveRDS(gamma.rbf.q2, "./results/gamma/rbfsvm.q2.RDS")

#     Test ----
# None pass the 0.6 line
# split 1 is best
gamma.tst <- rbfsvm.tst("./pre-process/gamma/", "./model.data/gamma/", 
                       nsplits = 10, cost = 3, e = 0.15, g = 0.01)
gamma.1to1 <- gamma.tst %>% filter(trn.split == tst.split)
gamma.avg <- data.table(gamma.tst, key = 'trn.split')
gamma.avg <- gamma.avg[ , list(r2 = mean(r2), 
                             rmse = mean(rmse)), 
                      by = trn.split] %>% print()

#     Model ----

# Again, none pass, but I thought that it might be interesting to keep 
# these models around
trn.gamma <- readRDS("./pre-process/gamma/1/pp.RDS") %>%
  select(., -guest)
features <- readRDS("./feature.selection/gamma.vars.RDS")
colnames(trn.gamma) <- str_replace(colnames(trn.gamma), "-", ".")
trn.gamma <- trn.gamma %>% select(., DelG, features) 
trn.gamma.x <- select(trn.gamma, -DelG) 
trn.gamma.y <- select(trn.gamma, DelG) 

rbfsvm.gamma <- svm(x = trn.gamma.x, y = trn.gamma.y,
                   cost = 3, epsilon = 0.15, gamma = 0.01,
                   kernel = "radial")
tst.gamma.df <- svm.tst.splits("pre-process/gamma/", "model.data/gamma/", 
                              features, 10, rbfsvm.gamma)

# None pass
eval.tropsha(tst.gamma.df)
graph.gamma <- ggplot(tst.gamma.df, aes(x = obs, y = pred, color = split)) + 
  geom_point() + 
  theme_bw() + 
  coord_fixed()  + 
  geom_abline(intercept = 0, slope = 1) + 
  labs(x = "Experimental dG, kJ/mol", y = "Predicted dG, kJ/mol", 
       title = "Gamma-CD Polynomial SVM", color = "Test split")
graph.gamma

#     Save ----

pp.settings <- readRDS("./pre-process/gamma/10/pp.settings.RDS")
saveRDS(list(pp.settings, rbfsvm.gamma), "./models/gamma/rbfsvm.RDS")
saveRDS(tst.gamma.df, "./results/gamma/rbfsvm.all.RDS")
saveRDS(tst.gamma.df %>% filter(split == "1"), 
        "./results/gamma/rbfsvm.RDS")
graph.gamma 
ggsave("./results/gamma/rbfsvm.png")

# Sigmoid -----------------------------------------------------------------
#   Alpha ----
#     LOOCV-Q2 ----

# Alpha - 5, 6, 9 pass
sig.alpha.q2 <- sigsvm.looq2(read.dir = "./pre-process/alpha/", nsplits = 10, 
             cost = 3, coef = 0, e = 0.15, g = 0.01)
saveRDS(sig.alpha.q2, "./results/alpha/sigsvm.q2.RDS")

#     Test ----

# Highest is 5, at 0.5178
# None pass lowest threshold
alpha.tst <- sigsvm.tst("./pre-process/alpha/", "./model.data/alpha/", 
                        nsplits = 10, cost = 3, 
                        coef = 0, e = 0.15, g = 0.01)
alpha.1to1 <- alpha.tst %>% filter(trn.split == tst.split)
alpha.avg <- data.table(alpha.tst, key = 'trn.split')
alpha.avg <- alpha.avg[ , list(r2 = mean(r2), 
                               rmse = mean(rmse)), 
                        by = trn.split] %>% print()
#     Model ----

trn.alpha <- readRDS("./pre-process/alpha/5/pp.RDS") %>%
  select(., -guest)
features <- readRDS("./feature.selection/alpha.vars.RDS")
colnames(trn.alpha) <- str_replace(colnames(trn.alpha), "-", ".")
trn.alpha <- trn.alpha %>% select(., DelG, features)
trn.alpha.x <- select(trn.alpha, -DelG)
trn.alpha.y <- select(trn.alpha, DelG)

sigsvm.alpha <- svm(x = trn.alpha.x, y = trn.alpha.y,
                    cost = 15, coef = 0, e = 0.05, g = 0.001,
                    kernel = "sigmoid")

tst.alpha <- preprocess.tst.mod("./pre-process/alpha/", "./model.data/alpha/",
                                features, 5)

tst.alpha.df <- svm.tst.splits("pre-process/alpha/", "model.data/alpha/", 
                               features, 10, sigsvm.alpha)

eval.tropsha(tst.alpha.df)
graph.alpha <- ggplot(tst.alpha.df, aes(x = obs, y = pred, color = split)) +
  geom_point() +
  theme_bw() +
  coord_fixed()  +
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Observed dG, kJ/mol", y = "Predicted dG, kJ/mol",
       title = "Sigmoid SVM for Alpha", color = "Test split")
graph.alpha

#     Save ----
# No model saved because none passed evaluation
ggsave("./results/alpha/sigsvm.png")
saveRDS(tst.alpha.df, "./results/alpha/sigsvm.all.RDS")
saveRDS(tst.alpha.df %>% filter(split == "5"), 
        "./results/alpha/sigsvm.RDS")

#   Beta ----
#     LOO-CV ----
# Beta - 1, 4, 5, 6, 7, 9, 10
sig.beta.q2 <- sigsvm.looq2(read.dir = "./pre-process/beta/", nsplits = 10, 
                            cost = 150, coef = 0, e = 0.5, g = 0.001)
saveRDS(sig.beta.q2, "./results/beta/sigsvm.q2.RDS")

#     Test ----

# None pass
# Highest is 4, at 0.5306
beta.tst <- sigsvm.tst("./pre-process/beta/", "./model.data/beta/", 
                       nsplits = 10, cost = 150, coef = 0, 
                       e = 0.5, g = 0.001)
beta.1to1 <- beta.tst %>% filter(trn.split == tst.split)
beta.avg <- data.table(beta.tst, key = 'trn.split')
beta.avg <- beta.avg[ , list(r2 = mean(r2), 
                             rmse = mean(rmse)), 
                      by = trn.split] %>% print()

#     Model ----

trn.beta <- readRDS("./pre-process/beta/1/pp.RDS") %>%
  select(., -guest)
features <- readRDS("./feature.selection/beta.vars.RDS")
colnames(trn.beta) <- str_replace(colnames(trn.beta), "-", ".")
trn.beta <- trn.beta %>% select(., DelG, features)
trn.beta.x <- select(trn.beta, -DelG)
trn.beta.y <- select(trn.beta, DelG)

sigsvm.beta <- svm(x = trn.beta.x, y = trn.beta.y,
                   cost = 75, coef = 0, e = 0.1, g = 0.001,
                   kernel = "sigmoid")

tst.beta.df <- svm.tst.splits("pre-process/beta/", "model.data/beta/", 
                               features, 10, sigsvm.beta)

eval.tropsha(tst.beta.df)
graph.beta <- ggplot(tst.beta.df, aes(x = obs, y = pred, color = split)) +
  geom_point() +
  theme_bw() +
  coord_fixed()  +
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Observed dG, kJ/mol", y = "Predicted dG, kJ/mol",
       title = "Sigmoid SVM for Beta", color = "Test split")
graph.beta

#     Save ----

ggsave("./results/beta/sigsvm.png")
saveRDS(tst.beta.df, "./results/beta/sigsvm.all.RDS")
saveRDS(tst.beta.df %>% filter(split == "1"), 
        "./results/beta/sigsvm.RDS")

#   Gamma ----
#      LOO-CV ----
sig.gamma.q2 <- sigsvm.looq2(read.dir = "./pre-process/gamma/", nsplits = 10, 
                             cost = 4, coef = 0, e = 0.05, g = 0.01)
saveRDS(sig.gamma.q2, "./results/gamma/sigsvm.q2.RDS")

#      Test ----

# Highest is 10, at 0.1695
gamma.tst <- sigsvm.tst("./pre-process/gamma/", "./model.data/gamma/", 
                       nsplits = 10, cost = 150, coef = 0, 
                       e = 0.5, g = 0.001)
gamma.1to1 <- gamma.tst %>% filter(trn.split == tst.split)
gamma.avg <- data.table(gamma.tst, key = 'trn.split')
gamma.avg <- gamma.avg[ , list(r2 = mean(r2), 
                             rmse = mean(rmse)), 
                      by = trn.split] %>% print()

#     Model ----

trn.gamma <- readRDS("./pre-process/gamma/10/pp.RDS") %>%
  select(., -guest)
features <- readRDS("./feature.selection/gamma.vars.RDS")
colnames(trn.gamma) <- str_replace(colnames(trn.gamma), "-", ".")
trn.gamma <- trn.gamma %>% select(., DelG, features)
trn.gamma.x <- select(trn.gamma, -DelG)
trn.gamma.y <- select(trn.gamma, DelG)

sigsvm.gamma <- svm(x = trn.gamma.x, y = trn.gamma.y,
                   cost = 150, coef = 0, e = 0.5, g = 0.001,
                   kernel = "sigmoid")

tst.gamma.df <- svm.tst.splits("pre-process/gamma/", "model.data/gamma/", 
                              features, 10, sigsvm.gamma)

eval.tropsha(tst.gamma.df)
graph.gamma <- ggplot(tst.gamma.df, aes(x = obs, y = pred, color = split)) +
  geom_point() +
  theme_bw() +
  coord_fixed()  +
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Observed dG, kJ/mol", y = "Predicted dG, kJ/mol",
       title = "Sigmoid SVM for Gamma", color = "Test split")
graph.gamma

#     Save ----

ggsave("./results/gamma/sigsvm.png")
saveRDS(tst.gamma.df, "./results/gamma/sigsvm.all.RDS")
saveRDS(tst.gamma.df %>% filter(split == "10"), 
        "./results/gamma/sigsvm.RDS")
