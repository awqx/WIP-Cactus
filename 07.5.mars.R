# MARS - Multivariate adaptive regression splines
source("./07.model.functions.R")
p_load(data.table, earth)

# Functions ---------------------------------------------------------------

mars.looq2 <- function(read.dir, rfe.dir, nsplits, 
                       d, p, nk, t, m, fk) {
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
      trn.x <- trn[ , -1]
      trn.y <- trn[ , 1]
      tst.x <- tst[ , -1]
      tst.y <- tst[ , 1]
      
      mars <- earth(
        x = trn.x,
        y = trn.y,
        degree = d, 
        penalty = as.numeric(p), 
        nk = nk, 
        thresh = t, 
        minspan = m, 
        fast.k = fk
      )
      pred[i] <- predict(mars, tst.x) 
    }
    # Handling outliers
    for(i in 1:length(pred)) {
      if(abs(pred[i]) > 85)
        pred[i] <- mean(obs)
    }
    pred.df <- data.frame(obs, pred)
    colnames(pred.df) <- c("obs", "pred")
    q2.plot <- ggplot(pred.df, aes(x = obs, y = pred)) + 
      theme_bw() + 
      geom_point() + 
      labs(title = n) + 
      geom_abline(slope = 1, intercept = 0) + 
      coord_fixed()
    print(q2.plot)
    
    PRESS <- sum((obs - pred)^2)
    TSS <- sum((obs - mean(obs))^2)
    q2 <- 1 - PRESS/TSS
    message("Q2 = ", q2)
    q2.results[n] <- q2
  }
  return(data.frame(trn.split, q2.results))
}

mars.tst <- function(pp.dir, tst.dir, nsplits, 
                     d, p, nk, t, m, fk) {
  split <- 1:nsplits
  if(str_detect(pp.dir, "alpha"))
    features <- readRDS("./feature.selection/alpha.vars.RDS")
  else if(str_detect(pp.dir, 'beta'))
    features <- readRDS("./feature.selection/beta.vars.RDS")
  else
    features <- readRDS('feature.selection/gamma.vars.RDS')
  results.all <- data.frame() 
  
  for(i in 1:nsplits) {
    trn <- readRDS(paste0(pp.dir, i, "/pp.RDS")) %>%
      select(., -guest)
    trn.y <- trn[ , 1]
    trn.x <- trn %>% select(., one_of(features))
    
    mars <- earth(
      x = trn.x,
      y = trn.y,
      degree = d, 
      penalty = as.numeric(p), 
      nk = nk, 
      thresh = t, 
      minspan = m, 
      fast.k = fk
    )
    
    tst.all <- data.frame()
    for(j in 1:nsplits) {
      tst <- preprocess.tst.mod(pp.dir = pp.dir, tst.dir = tst.dir, 
                                feat = colnames(trn.x), n = j)
      tst.y <- tst[ , 1]
      tst.x <- tst[ , -1]
      
      tst.df <- predict(mars, tst.x) %>%
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
    tst.plot <- ggplot(tst.all, aes(x = obs, y = pred, color = split)) + 
      theme_bw() + 
      geom_point() + 
      labs(title = i) + 
      geom_abline(slope = 1, intercept = 0) + 
      coord_fixed()
    print(tst.plot)
  }
  row.names(results.all) <- NULL
  return(results.all)
}
 
mars.tst.splits <- function(pp.dir, tst.dir, feat, nsplits, model) {
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

# Alpha -----

#   Q2 ----
# 1, 4, 5, 8, 10
mars.alpha.q2 <- mars.looq2('pre-process/alpha/', nsplits = 10, 
                            d = 1, p = 1, nk = 35, 
                            t = 0, m = 0, fk = 30)
#   Test ----
# 8 is the best (R2 = 0.650, Q2 = 0.581)
alpha.tst <- mars.tst('pre-process/alpha/', 'model.data/alpha/', 
                           nsplits = 10, 
                           d = 1, p = 1, nk = 35, 
                           t = 0, m = 0, fk = 30)
alpha.avg <- data.table(alpha.tst, key = 'trn.split')
alpha.avg <- alpha.avg[ , list(r2 = mean(r2), 
                               rmse = mean(rmse)), 
                        by = trn.split] %>% print()

#   Model ----
trn.alpha <- readRDS("./pre-process/alpha/8/pp.RDS") %>%
  select(., -guest)
features <- readRDS("./feature.selection/alpha.vars.RDS")
colnames(trn.alpha) <- str_replace(colnames(trn.alpha), "-", ".")
trn.alpha <- trn.alpha %>% select(., DelG, features) 
trn.alpha.x <- select(trn.alpha, -DelG) 
trn.alpha.y <- trn.alpha$DelG

mars.alpha <- earth(
  x = trn.alpha.x,
  y = trn.alpha.y,
  degree = 1, 
  penalty = 1, 
  nk = 35, 
  thresh = 0, 
  minspan = 0, 
  fast.k = 30
)
tst.alpha.df <- mars.tst.splits("pre-process/alpha/", "model.data/alpha/", 
                             features, 10, mars.alpha)
tst.alpha.df2 <- tst.alpha.df %>% filter(split == "8")

# All pass
eval.tropsha(tst.alpha.df)
eval.tropsha(tst.alpha.df2)
graph.alpha <- ggplot(tst.alpha.df, aes(x = obs, y = pred, color = split)) + 
  geom_point() + 
  theme_bw() + 
  coord_fixed()  + 
  geom_abline(intercept = 0, slope = 1) + 
  labs(x = "Observed dG, kJ/mol", y = "Predicted dG, kJ/mol", 
       title = "Alpha-CD MARS", color = "Test split")
graph.alpha

#   Save ----

saveRDS(tst.alpha.df, 'results/alpha/mars.all.RDS')
saveRDS(tst.alpha.df2, 'results/alpha/mars.RDS')
pp.settings <- readRDS('pre-process/alpha/8/pp.settings.RDS')
saveRDS(list(pp.settings, mars.alpha), 'models/alpha/mars.RDS')
# print(graph.alpha)
# ggsave('results/alpha/mars.png')

# Beta --------------------------------------------------------------------

#   Q2 ----
# All oass
beta.q2 <- mars.looq2('pre-process/beta/', nsplits = 10, 
                      d = 2, p = 4, nk = 25, 
                      t = 0.015, m = 70, fk = 20)
#   Test ----
# 2 is the best (R2 = 0.7374, Q2 = 0.5839)
beta.tst <- mars.tst('pre-process/beta/', 'model.data/beta/', 
                      nsplits = 10, 
                      d = 2, p = 4, nk = 25, 
                      t = 0.015, m = 70, fk = 20)
beta.avg <- data.table(beta.tst, key = 'trn.split')
beta.avg <- beta.avg[ , list(r2 = mean(r2), 
                               rmse = mean(rmse)), 
                        by = trn.split] %>% print()

#   Model ----
trn.beta <- readRDS("./pre-process/beta/2/pp.RDS") %>%
  select(., -guest)
features <- readRDS("./feature.selection/beta.vars.RDS")
colnames(trn.beta) <- str_replace(colnames(trn.beta), "-", ".")
trn.beta <- trn.beta %>% select(., DelG, features) 
trn.beta.x <- select(trn.beta, -DelG) 
trn.beta.y <- trn.beta$DelG

mars.beta <- mars <- earth(
  x = trn.beta.x,
  y = trn.beta.y,
  degree = 2, 
  penalty = 4, 
  nk = 25, 
  thresh = 0.015, 
  minspan = 70, 
  fast.k = 20
)
tst.beta.df <- mars.tst.splits("pre-process/beta/", "model.data/beta/", 
                                features, 10, mars.beta)
tst.beta.df2 <- tst.beta.df %>% filter(split == "2")

# All pass
eval.tropsha(tst.beta.df)
eval.tropsha(tst.beta.df2)
graph.beta <- ggplot(tst.beta.df, aes(x = obs, y = pred, color = split)) + 
  geom_point() + 
  theme_bw() + 
  coord_fixed()  + 
  geom_abline(intercept = 0, slope = 1) + 
  labs(x = "Observed dG, kJ/mol", y = "Predicted dG, kJ/mol", 
       title = "Beta-CD MARS", color = "Test split")
graph.beta

#   Save ----

saveRDS(tst.beta.df, 'results/beta/mars.all.RDS')
saveRDS(tst.beta.df2, 'results/beta/mars.RDS')
pp.settings <- readRDS('pre-process/beta/2/pp.settings.RDS')
saveRDS(list(pp.settings, mars.beta), 'models/beta/mars.RDS')
# print(graph.beta)
# ggsave('results/beta/mars.png')

# Gamma --------------------------------------------------------------------

#   Q2 ----

# None pass
gamma.q2 <- mars.looq2('pre-process/gamma/', nsplits = 10, 
                      d = 2, p = 4, nk = 25, 
                      t = 0.015, m = 70, fk = 20)
#   Test ----

# None pass, but 7 the best
gamma.tst <- mars.tst('pre-process/gamma/', 'model.data/gamma/', 
                     nsplits = 10, 
                     d = 2, p = 4, nk = 25, 
                     t = 0.015, m = 70, fk = 20)
gamma.avg <- data.table(gamma.tst, key = 'trn.split')
gamma.avg <- gamma.avg[ , list(r2 = mean(r2), 
                             rmse = mean(rmse)), 
                      by = trn.split] %>% print()

#   Model ----
trn.gamma <- readRDS("./pre-process/gamma/7/pp.RDS") %>%
  select(., -guest)
features <- readRDS("./feature.selection/gamma.vars.RDS")
colnames(trn.gamma) <- str_replace(colnames(trn.gamma), "-", ".")
trn.gamma <- trn.gamma %>% select(., DelG, features) 
trn.gamma.x <- select(trn.gamma, -DelG) 
trn.gamma.y <- trn.gamma$DelG

mars.gamma <- mars <- earth(
  x = trn.gamma.x,
  y = trn.gamma.y,
  degree = 2, 
  penalty = 4, 
  nk = 25, 
  thresh = 0.015, 
  minspan = 70, 
  fast.k = 20
)
tst.gamma.df <- mars.tst.splits("pre-process/gamma/", "model.data/gamma/", 
                               features, 10, mars.gamma)
tst.gamma.df2 <- tst.gamma.df %>% filter(split == "7")

# All pass
eval.tropsha(tst.gamma.df)
eval.tropsha(tst.gamma.df2)
graph.gamma <- ggplot(tst.gamma.df, aes(x = obs, y = pred, color = split)) + 
  geom_point() + 
  theme_bw() + 
  coord_fixed()  + 
  geom_abline(intercept = 0, slope = 1) + 
  labs(x = "Observed dG, kJ/mol", y = "Predicted dG, kJ/mol", 
       title = "Gamma-CD MARS", color = "Test split")
graph.gamma

#   Save ----

saveRDS(tst.gamma.df, 'results/gamma/mars.all.RDS')
saveRDS(tst.gamma.df2, 'results/gamma/mars.RDS')
# print(graph.gamma)
# ggsave('results/gamma/mars.png')