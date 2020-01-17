source("./07.model.functions.R")
p_load(Cubist)

# Functions ---------------------------------------------------------------

cubist.looq2 <- function(read.dir, nsplits, cmte, extra, seed) {
  trn.split <- 1:nsplits
  q2.results <- c(rep(0.0, nsplits))
  if(str_detect(read.dir, "alpha"))
    features <- readRDS("./feature.selection/alpha.vars.RDS")
  else
    features <- readRDS("./feature.selection/beta.vars.RDS")
  
  ctrl <- cubistControl(
    seed = seed, 
    extrapolation = extra
  )
  
  for(n in 1:nsplits) {
    data <- readRDS(paste0(read.dir, n, "/pp.RDS")) %>%
      select(., -guest)
    obs <- data[ , 1]
    data <- data %>% select(., DelG, one_of(features))
    pred <- c(rep(0.0, nrow(data) - 1))
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

cube.tst <- function(pp.dir, tst.dir, nsplits, cmte, extra, seed) {
  if(str_detect(pp.dir, "alpha"))
    features <- readRDS("./feature.selection/alpha.vars.RDS")
  else
    features <- readRDS("./feature.selection/beta.vars.RDS")
  
  ctrl <- cubistControl(
    seed = seed, 
    extrapolation = extra
  )
  results.all <- data.frame()
  
  for(n in 1:nsplits) {
    trn <- readRDS(paste0(pp.dir, n, "/pp.RDS")) %>%
      select(., -guest)
    trn.y <- trn$DelG
    trn.x <- select(trn, one_of(features))
    cube <- cubist(x = trn.x, y = trn.y,  
                   committees = cmte, 
                   control = ctrl)
    tst.all <- data.frame()
    
    for(j in 1:nsplits) {
      tst <- preprocess.tst.mod(pp.dir = pp.dir, tst.dir = tst.dir, 
                                feat = colnames(trn.x), n = j)
      tst.y <- tst[ , 1]
      tst.x <- tst[ , -1]
      
      tst.df <- predict(cube, tst.x) %>%
        cbind(tst.y, .) %>%
        as.data.frame()
      colnames(tst.df) <- c("obs", "pred")
      tst.df <- tst.df %>% mutate(split = j)
      
      for(k in 1:nrow(tst.df))
        if(abs(tst.df$pred[k]) > 100)
          tst.df$pred[k] <- mean(trn.y)
      
      results.tst <- data.frame(trn.split = n, tst.split = j, 
                                r2 = defaultSummary(tst.df)[2], 
                                rmse = defaultSummary(tst.df)[1])
      results.all <- rbind(results.all, results.tst)
      tst.all <- rbind(tst.all, tst.df)
    }
    tst.all <- tst.all %>% mutate(split = as.factor(split))
    p <- ggplot(tst.all, aes(x = obs, y = pred, color = split)) + 
      theme_bw() + 
      geom_point() + 
      labs(title = n) + 
      geom_abline(slope = 1, intercept = 0) + 
      coord_fixed()
    print(p)
  }
  row.names(results.all) <- NULL
  return(results.all)
}

cube.tst.splits <- function(pp.dir, tst.dir, feat, nsplits, model) {
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

# Alpha ----
#   LOO-CV ----
#  1, 4, 5, 8, 10
cube.alpha.q2 <- cubist.looq2(read.dir = "./pre-process/alpha/", 
                              nsplits = 10, seed = 101,
                              cmte = 90, extra = 10)

#   Test ----

# 1 is the best
alpha.tst <- cube.tst("pre-process/alpha/", "model.data/alpha", 
                      nsplits = 10, seed = 101,
                      cmte = 90, extra = 10)
avg.tst(alpha.tst)

#   Model ----

trn.alpha <- readRDS("./pre-process/alpha/1/pp.RDS") %>%
  select(., -guest)
features <- readRDS("./feature.selection/alpha.vars.RDS")
colnames(trn.alpha) <- str_replace(colnames(trn.alpha), "-", ".")
trn.alpha <- trn.alpha %>% select(., DelG, one_of(features)) 
trn.alpha.x <- select(trn.alpha, -DelG) 
trn.alpha.y <- trn.alpha$DelG

cube.alpha.ctrl <- cubistControl(
  seed = 101, 
  extrapolation = 10
)
cube.alpha <- cubist(x = trn.alpha.x, y = trn.alpha.y,  
                     committees = 90, 
                     control = cube.alpha.ctrl)

tst.alpha.df <- cube.tst.splits('pre-process/alpha/', 'model.data/alpha/', 
                                features, 10, cube.alpha)
tst.alpha.df2 <- tst.alpha.df %>% filter(split == "1")

# All pass
eval.tropsha(tst.alpha.df)
eval.tropsha(tst.alpha.df2)

graph.alpha <- ggplot(tst.alpha.df, aes(x = obs, y = pred, color = split)) + 
  geom_point() + 
  theme_bw() + 
  coord_fixed()  + 
  geom_abline(intercept = 0, slope = 1) + 
  labs(x = "Observed dG, kJ/mol", y = "Predicted dG, kJ/mol", 
       title = "Alpha-CD Cubist", color = "Test split")
graph.alpha

#   Save ----

pp.settings <- readRDS("./pre-process/alpha/1/pp.settings.RDS")
saveRDS(list(pp.settings, cube.alpha), "./models/alpha/cube.RDS")
saveRDS(tst.alpha.df, "./results/alpha/cube.all.RDS")
saveRDS(tst.alpha.df2, "./results/alpha/cube.RDS")

# Beta ----
#   LOO-CV ----
# All pass
cube.beta.q2 <- cubist.looq2(read.dir = "./pre-process/beta/", 
                             nsplits = 10, seed = 101, 
                             cmte = 50, extra = 40)

#   Test ----

beta.tst <- cube.tst("./pre-process/beta/", "./model.data/beta", 
                     nsplits = 10, seed = 101, 
                     cmte = 50, extra = 40)
avg.tst(beta.tst)

#   Model ----

trn.beta <- readRDS("./pre-process/beta/9/pp.RDS") %>%
  select(., -guest)
features <- readRDS("./feature.selection/beta.vars.RDS")
colnames(trn.beta) <- str_replace(colnames(trn.beta), "-", ".")
trn.beta <- trn.beta %>% select(., DelG, features) 
trn.beta.x <- select(trn.beta, -DelG) 
trn.beta.y <- trn.beta$DelG

cube.beta.ctrl <- cubistControl(
  seed = 101, 
  extrapolation = 40
)

cube.beta <- cubist(x = trn.beta.x, y = trn.beta.y,  
                    committees = 50, 
                    control = cube.beta.ctrl)

tst.beta.df <- cube.tst.splits('pre-process/beta/', 'model.data/beta/', 
                               features, 10, cube.beta)
tst.beta.df2 <- tst.beta.df %>% filter(split == "9")

# All pass
eval.tropsha(tst.beta.df)
eval.tropsha(tst.beta.df2)

graph.beta <- ggplot(tst.beta.df, aes(x = obs, y = pred, color = split)) + 
  geom_point() + 
  theme_bw() + 
  coord_fixed()  + 
  geom_abline(intercept = 0, slope = 1) + 
  labs(x = "Observed dG, kJ/mol", y = "Predicted dG, kJ/mol", 
       title = "Beta-CD Cubist", color = "Test split")
graph.beta

#   Save ----

pp.settings <- readRDS("./pre-process/beta/9/pp.settings.RDS")
saveRDS(list(pp.settings, cube.beta), "./models/beta/cube.RDS")
saveRDS(tst.beta.df, "./results/beta/cube.all.RDS")
saveRDS(tst.beta.df2, "./results/beta/cube.RDS")

# Gamma ----

#   LOO-CV ----
# Gamma - none pass
cube.gamma.q2 <- cubist.looq2(read.dir = "./pre-process/gamma/", 
                             nsplits = 10, seed = 101, 
                             cmte = 100, extra = 20)

#   Test ----

# GAMMA - 6
gamma.tst <- cube.tst("./pre-process/gamma/", "./model.data/gamma", 
                     nsplits = 10, seed = 101, 
                     cmte = 100, extra = 20)
avg.tst(gamma.tst)

#   Model ----

trn.gamma <- readRDS("./pre-process/gamma/6/pp.RDS") %>%
  select(., -guest)
features <- readRDS("./feature.selection/gamma.vars.RDS")
colnames(trn.gamma) <- str_replace(colnames(trn.gamma), "-", ".")
trn.gamma <- trn.gamma %>% select(., DelG, features) 
trn.gamma.x <- select(trn.gamma, -DelG) 
trn.gamma.y <- trn.gamma$DelG

cube.gamma.ctrl <- cubistControl(
  seed = 101, 
  extrapolation = 20
)

cube.gamma <- cubist(x = trn.gamma.x, y = trn.gamma.y,  
                    committees = 100, 
                    control = cube.gamma.ctrl)

tst.gamma.df <- cube.tst.splits('pre-process/gamma/', 'model.data/gamma/', 
                               features, 10, cube.gamma)
tst.gamma.df2 <- tst.gamma.df %>% filter(split == "6")

eval.tropsha(tst.gamma.df)
eval.tropsha(tst.gamma.df2)
graph.gamma <- ggplot(tst.gamma.df, aes(x = obs, y = pred, color = split)) + 
  geom_point() + 
  theme_bw() + 
#   coord_fixed()  + 
  geom_abline(intercept = 0, slope = 1) + 
  labs(x = "Observed dG, kJ/mol", y = "Predicted dG, kJ/mol", 
       title = "Gamma-CD Cubist", color = "Test split")
graph.gamma

#     Save ----

pp.settings <- readRDS("./pre-process/gamma/6/pp.settings.RDS")
saveRDS(list(pp.settings, cube.gamma), "./models/gamma/cube.RDS")
saveRDS(tst.gamma.df, "./results/gamma/cube.all.RDS")
saveRDS(tst.gamma.df2, "./results/gamma/cube.RDS")
