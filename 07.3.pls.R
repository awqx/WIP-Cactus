source("07.model.functions.R")
p_load(data.table, pls)

# Functions ---------------------------------------------------------------

pls.looq2 <- function(read.dir, rfe.dir, nsplits, method, ncomp) {
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
    pred <- c()
    # In this loop, i represents the index of the datapoint left out of 
    # model building
    for(i in 1:nrow(data)) {
      trn <- data[-i, ]
      tst <- data[i, ]
      # trn.x <- trn[ , -1]
      # trn.y <- trn[ , 1]
      # tst.x <- tst[ , -1]
      # tst.y <- tst[ , 1]
      
      pls.mod <- plsr(DelG~., data = trn, 
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

pls.tst <- function(pp.dir, tst.dir, nsplits, ncomp, method) {
  split <- 1:nsplits
  results.all <- data.frame()
  if(str_detect(pp.dir, "alpha"))
    features <- readRDS("./feature.selection/alpha.vars.RDS")
  else
    features <- readRDS("./feature.selection/beta.vars.RDS")
  pls.mod <- NULL
  for(i in 1:nsplits) {
    trn <- readRDS(paste0(pp.dir, i, "/pp.RDS")) %>%
      select(., DelG, one_of(features))
    pls.mod <- plsr(DelG~., data = trn, 
                    ncomp = ncomp, method = method)
    tst.all <- data.frame()
    
    for(j in 1:nsplits) {
      tst <- preprocess.tst.mod(pp.dir = pp.dir, tst.dir = tst.dir, 
                                feat = features, n = j)
      tst.y <- tst[ , 1]
      tst.x <- tst[ , -1]
      
      tst.df <- predict(pls.mod, tst[ , -1], ncomp = ncomp) %>%
        cbind(tst[ , 1], .) %>%
        as.data.frame()
      colnames(tst.df) <- c("obs", "pred")
      tst.df <- tst.df %>% mutate(split = j)
      # for(k in 1:nrow(tst.df))
      #   if(abs(tst.df$pred[k]) > 80)
      #     tst.df$pred[k] <- mean(trn.y)
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
# Alpha ----
#   LOO-CV ----

# 1, 4, 5, 6, 8, 9, 10
pls.looq2("./pre-process/alpha/", nsplits = 10, 
          ncomp = 8, method = "oscorespls")

#   Test ----

# None pass, 3 does the best
alpha.tst <- pls.tst("./pre-process/alpha/", "./model.data/alpha/", 
                     nsplits = 10, ncomp = 8, method = "oscorespls")
alpha.1to1 <- alpha.tst %>% filter(trn.split == tst.split)
alpha.avg <- avg.tst(alpha.tst) %>% print()

#   Model ----

trn.alpha <- readRDS("./pre-process/alpha/3/pp.RDS") %>%
  select(., -guest)
features <- readRDS("./feature.selection/alpha.vars.RDS")
colnames(trn.alpha) <- str_replace(colnames(trn.alpha), "-", ".")
trn.alpha <- trn.alpha %>% select(., DelG, features) 
# trn.alpha.x <- select(trn.alpha, -DelG) 
# trn.alpha.y <- trn.alpha$DelG

pls.alpha <- plsr(DelG~., data = trn.alpha, 
                  ncomp = 8, method = "oscorespls")
tst.alpha.df <- tst.splits("pre-process/alpha/", "model.data/alpha/", 
                           features, 10, pls.alpha) %>%
  data.table(., key = c("obs", "split"))
tst.alpha.df <- tst.alpha.df[ , list(pred = mean(pred)), 
                              by = list(obs, split)] %>%
  as.data.frame()

tst.alpha.df2 <- tst.alpha.df %>% filter(split == "3")

eval.tropsha(tst.alpha.df)
eval.tropsha(tst.alpha.df2)

#   Save ----
saveRDS(tst.alpha.df, "./results/alpha/pls.all.RDS")
saveRDS(tst.alpha.df2, "./results/alpha/pls.RDS")

# Beta ----
#   LOO-CV ----

# All pass. 0.526
pls.looq2("./pre-process/beta/", nsplits = 10, 
          ncomp = 16, method = "oscorespls")

#   Test ----

# None pass
# 10 does the best
beta.tst <- pls.tst("./pre-process/beta/", "./model.data/beta/", 
                    nsplits = 10, ncomp = 16, method = "oscorespls")
beta.1to1 <- beta.tst %>% filter(trn.split == tst.split)
beta.avg <- avg.tst(beta.tst) %>% print()

#   Model ----

trn.beta <- readRDS("./pre-process/beta/10/pp.RDS") %>%
  select(., -guest)
features <- readRDS("./feature.selection/beta.vars.RDS")
colnames(trn.beta) <- str_replace(colnames(trn.beta), "-", ".")
trn.beta <- trn.beta %>% select(., DelG, features) 
# trn.beta.x <- select(trn.beta, -DelG) 
# trn.beta.y <- trn.beta$DelG

pls.beta <- plsr(DelG~., data = trn.beta, 
                  ncomp = 16, method = "oscorespls")
tst.beta.df <- tst.splits("pre-process/beta/", "model.data/beta/", 
                           features, 10, pls.beta) %>%
  data.table(., key = c("obs", "split"))
tst.beta.df <- tst.beta.df[ , list(pred = mean(pred)), 
                              by = list(obs, split)] %>%
  as.data.frame()

tst.beta.df2 <- tst.beta.df %>% filter(split == "3")

eval.tropsha(tst.beta.df)
eval.tropsha(tst.beta.df2)

#   Save ----
saveRDS(tst.beta.df, "./results/beta/pls.all.RDS")
saveRDS(tst.beta.df2, "./results/beta/pls.RDS")

#   Save ----

saveRDS(tst.beta.df, "./results/beta/pls.all.RDS")
saveRDS(tst.beta.df2, "./results/beta/pls.RDS")

# Gamma ----
#   LOOCV ----

# None
pls.looq2("./pre-process/gamma/", nsplits = 10, 
          ncomp = 2, method = "oscorespls")

#   Test ----

# None pass
# 6 does the best
gamma.tst <- pls.tst("./pre-process/gamma/", "./model.data/gamma/", 
                    nsplits = 10, ncomp = 2, method = "oscorespls")
gamma.1to1 <- gamma.tst %>% filter(trn.split == tst.split)
gamma.avg <- avg.tst(gamma.tst) %>% print()

#   Model ----

trn.gamma <- readRDS("./pre-process/gamma/6/pp.RDS") %>%
  select(., -guest)
features <- readRDS("./feature.selection/gamma.vars.RDS")
colnames(trn.gamma) <- str_replace(colnames(trn.gamma), "-", ".")
trn.gamma <- trn.gamma %>% select(., DelG, features) 
# trn.gamma.x <- select(trn.gamma, -DelG) 
# trn.gamma.y <- trn.gamma$DelG

pls.gamma <- plsr(DelG~., data = trn.gamma, 
                 ncomp = 2, method = "oscorespls")
tst.gamma.df <- tst.splits("pre-process/gamma/", "model.data/gamma/", 
                           features, 10, pls.gamma) %>%
  data.table(., key = c("obs", "split"))
tst.gamma.df <- tst.gamma.df[ , list(pred = mean(pred)), 
                              by = list(obs, split)] %>%
  as.data.frame()

tst.gamma.df2 <- tst.gamma.df %>% filter(split == "3")

eval.tropsha(tst.gamma.df)
eval.tropsha(tst.gamma.df2)

#   Save ----
saveRDS(tst.gamma.df, "./results/gamma/pls.all.RDS")
saveRDS(tst.gamma.df2, "./results/gamma/pls.RDS")
