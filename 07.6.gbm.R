source("07.model.functions.R")
p_load(data.table, gbm)

# Functions ---------------------------------------------------------------

gbm.looq2 <- function(read.dir, nsplits, num, d, s, n) {
  trn.split <- 1:nsplits
  q2.results <- c(rep(0.0, nsplits))
  
  features <- find.features(read.dir)
  
  for(n in 1:nsplits) {
    data <- readRDS(paste0(read.dir, n, "/pp.RDS")) %>%
      select(., -guest)
    obs <- data[ , 1]
    data <- data %>% select(., DelG, one_of(features)) 
    pred <- c(rep(0.0, nrow(data) - 1))
    
    for(i in 1:nrow(data)) {
      trn <- data[-i, ]
      tst <- data[i, ]
      trn.x <- trn[ , -1]
      trn.y <- trn[ , 1]
      tst.x <- tst[ , -1]
      tst.y <- tst[ , 1]
      
      gbm <- gbm.fit(x = trn.x, 
                     y = trn.y,
                     n.trees = num,
                     interaction.depth = d, 
                     shrinkage = s, 
                     n.minobsinnode = n,
                     verbose = F, 
                     distribution = 'gaussian')
      pred[i] <- predict(gbm, tst.x, 
                         n.trees = num) 
    }
    # # Handling outliers
    # for(i in 1:length(pred)) {
    #   if(abs(pred[i]) > 85)
    #     pred[i] <- mean(obs)
    # }
    pred.df <- data.frame(obs, pred)
    colnames(pred.df) <- c("obs", "pred")
    p <- ggplot(pred.df, aes(x = obs, y = pred)) + 
      theme_bw() + 
      geom_point() + 
      labs(title = n) + 
      geom_abline(slope = 1, intercept = 0) + 
      coord_fixed()
    print(p)
    
    PRESS <- sum((obs - pred)^2)
    TSS <- sum((obs - mean(obs))^2)
    q2 <- 1 - PRESS/TSS
    message("Q2 = ", q2)
    q2.results[n] <- q2
  }
  return(data.frame(trn.split, q2.results))
}

gbm.tst <- function(pp.dir, tst.dir, nsplits,
                    num, d, s, n) {
  split <- 1:nsplits
  features <- find.features(pp.dir)
  
  results.all <- data.frame() 
  
  for(i in 1:nsplits) {
    trn <- readRDS(paste0(pp.dir, i, "/pp.RDS")) %>%
      select(., -guest)
    trn.y <- trn[ , 1]
    trn.x <- trn %>% select(., one_of(features))
    
    # Refer to 07.0.1.svm.tune.poly.R for values
    gbm <- gbm.fit(x = trn.x, 
                   y = trn.y,
                   n.trees = num,
                   interaction.depth = d, 
                   shrinkage = s, 
                   n.minobsinnode = n,
                   verbose = F, 
                   distribution = 'gaussian')
    
    tst.all <- data.frame()
    for(j in 1:nsplits) {
      tst <- preprocess.tst.mod(pp.dir = pp.dir, tst.dir = tst.dir, 
                                feat = colnames(trn.x), n = j)
      tst.y <- tst[ , 1]
      tst.x <- tst[ , -1]
      
      tst.df <- predict(gbm, tst.x, n.trees = num) %>%
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

gbm.tst.splits <- function(pp.dir, tst.dir, feat, nsplits, num, model) {
  tst.all <- data.frame()
  for (j in 1:nsplits) {
    tst <- preprocess.tst.mod(pp.dir = pp.dir, tst.dir = tst.dir, 
                              feat = feat, n = j)
    tst.y <- tst[ , 1]
    tst.x <- tst[ , -1]
    
    tst.df <- predict(model, tst.x,
                      n.trees = num) %>%
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
# Alpha -------------------------------------------------------------------

# All pass 
alpha.q2 <- gbm.looq2("./pre-process/alpha/", 10, 500, 4, .01, 4)
alpha.tst <- gbm.tst("./pre-process/alpha/", "./model.data/alpha/", 10,
                     500, 4, .01, 4)
alpha.1to1 <- alpha.tst %>% filter(trn.split == tst.split)
alpha.avg <- data.table(alpha.tst, key = 'trn.split')
# 2
alpha.avg <- alpha.avg[ , list(r2 = mean(r2), 
                               rmse = mean(rmse)), 
                        by = trn.split] %>% print()

trn.alpha <- readRDS("./pre-process/alpha/2/pp.RDS") %>%
  select(., -guest)
features <- readRDS("./feature.selection/alpha.vars.RDS")
colnames(trn.alpha) <- str_replace(colnames(trn.alpha), "-", ".")
trn.alpha <- trn.alpha %>% select(., DelG, features) 
trn.alpha.x <- select(trn.alpha, -DelG) 
trn.alpha.y <- trn.alpha$DelG

gbm.alpha <- gbm.fit(trn.alpha.x, 
                     trn.alpha.y, 
                     distribution = 'gaussian', 
                     verbose = F, 
                     n.trees = 500, 
                     interaction.depth = 4, 
                     shrinkage = 0.01, 
                     n.minobsinnode = 4)
alpha.single <- gbm.tst.splits("pre-process/alpha/", "model.data/alpha/", 
                                   features, 10, 500, gbm.alpha)

eval.tropsha(alpha.single)
graph.alpha <- ggplot(alpha.single, aes(x = obs, y = pred, color = split)) + 
  geom_point() + 
  theme_bw() + 
  coord_fixed()  + 
  geom_abline(intercept = 0, slope = 1) + 
  labs(x = "Observed dG, kJ/mol", y = "Predicted dG, kJ/mol", 
       title = "Alpha-CD GBM", color = "Test split")
print(graph.alpha)

pp.settings <- readRDS("./pre-process/alpha/2/pp.settings.RDS")
saveRDS(list(pp.settings, gbm.alpha), "./models/alpha/gbm.RDS")
saveRDS(alpha.single, "./results/alpha/gbm.RDS")
print(graph.alpha)
ggsave("./results/alpha/gbm.png")

# Beta --------------------------------------------------------------------

beta.q2 <- gbm.looq2("./pre-process/beta/", 10,
                     500, 5, .01, 1)

beta.tst <- gbm.tst("./pre-process/beta/", "./model.data/beta/", 10,
                     500, 5, .01, 1)
beta.1to1 <- beta.tst %>% filter(trn.split == tst.split)
beta.avg <- data.table(beta.tst, key = 'trn.split')
# 10
beta.avg <- beta.avg[ , list(r2 = mean(r2), 
                               rmse = mean(rmse)), 
                        by = trn.split] %>% print()


trn.beta <- readRDS("./pre-process/beta/10/pp.RDS") %>%
  select(., -guest)
features <- readRDS("./feature.selection/beta.vars.RDS")
colnames(trn.beta) <- str_replace(colnames(trn.beta), "-", ".")
trn.beta <- trn.beta %>% select(., DelG, features) 
trn.beta.x <- select(trn.beta, -DelG) 
trn.beta.y <- trn.beta$DelG

gbm.beta <- gbm.fit(trn.beta.x, 
                     trn.beta.y, 
                     distribution = 'gaussian', 
                     verbose = F, 
                     n.trees = 500, 
                     interaction.depth = 5, 
                     shrinkage = 0.01, 
                     n.minobsinnode = 1)
beta.single <- gbm.tst.splits("pre-process/beta/", "model.data/beta/", 
                                   features, 10, 500, gbm.beta) 
eval.tropsha(beta.single)
graph.beta <- ggplot(beta.single, aes(x = obs, y = pred, color = split)) + 
  geom_point() + 
  theme_bw() + 
  coord_fixed()  + 
  geom_abline(intercept = 0, slope = 1) + 
  labs(x = "Observed dG, kJ/mol", y = "Predicted dG, kJ/mol", 
       title = "Beta-CD GBM", color = "Test split")
print(graph.beta)

pp.settings <- readRDS("./pre-process/beta/10/pp.settings.RDS")
saveRDS(list(pp.settings, gbm.beta), "./models/beta/gbm.RDS")
saveRDS(beta.single, "./results/beta/gbm.RDS")
print(graph.beta)
ggsave("./results/beta/gbm.png")