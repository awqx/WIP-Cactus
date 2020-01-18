# makes use of the `randomForest` package
source("helpers/training.R")

# Tuning ------------------------------------------------------------------

  # Alpha ----

trn <- readRDS("modeling-data/alpha/split.RDS")[[10]][["trn"]]
fold_list <- split_train(trn, nfolds = 10, seed = 20200117)

ntree <- c(25, 50, 75, 100, 250, 500)
nodesize <- c(1:5, 10, 25)
mtry <- c(1:5, 7, 10)
combo_grid <- expand.grid(ntree, nodesize, mtry)
names(combo_grid) <- c("ntree", "nodesize", "mtry")
combo_list <- split(combo_grid, seq(nrow(combo_grid)))

alpha_tune <- do.call(
  rbind, 
  lapply(
    combo_list, 
    train_rf, 
    folds = fold_list
  )
)
  # Beta ----

trn <- readRDS("modeling-data/beta/split.RDS")[[10]][["trn"]]
fold_list <- split_train(trn, nfolds = 10, seed = 20200117)

ntree <- c(25, 50, 75, 100, 250, 500)
nodesize <- c(1:5, 10, 25)
mtry <- c(1:5, 8, 12, 17)
combo_grid <- expand.grid(ntree, nodesize, mtry)
names(combo_grid) <- c("ntree", "nodesize", "mtry")
combo_list <- split(combo_grid, seq(nrow(combo_grid)))

beta_tune <- do.call(
  rbind, 
  lapply(
    combo_list, 
    train_rf, 
    folds = fold_list
  )
)

  # Gamma ----

trn <- readRDS("modeling-data/gamma/split.RDS")[[10]][["trn"]]
fold_list <- split_train(trn, nfolds = 10, seed = 20200117)

ntree <- c(25, 50, 75, 100, 250, 500)
nodesize <- c(1:5, 10, 25)
mtry <- c(1:5, 8, 20, 40)
combo_grid <- expand.grid(ntree, nodesize, mtry)
names(combo_grid) <- c("ntree", "nodesize", "mtry")
combo_list <- split(combo_grid, seq(nrow(combo_grid)))

gamma_tune <- do.call(
  rbind, 
  lapply(
    combo_list, 
    train_rf, 
    folds = fold_list
  )
)

  # Save ----

saveRDS(alpha_tune, "training/alpha/rf.RDS")
saveRDS(beta_tune, "training/beta/rf.RDS")
saveRDS(gamma_tune, "training/gamma/rf.RDS")

# Model -------------------------------------------------------------------

  # Alpha ----

# refer to `05-svm.R` for comments and further explanation
alpha_split <- readRDS("modeling-data/alpha/split.RDS")
names(alpha_split) <- c(1:10)
list(
  "ntree" = 25, 
  "nodesize" = 2, 
  "mtry" = 1
) %>%
  train_rf(
    folds = alpha_split, 
    param = ., 
    fold_avg = F
) # best is 4

# saving the best model
trn <- readRDS("modeling-data/alpha/split.RDS")[[4]]$trn
saveRDS(
  randomForest(
    x = trn[, -1:-2], 
    y = trn$dG, 
    ntree = 25, 
    nodesize = 2, 
    mtry = 1
  ),
  "model/alpha/rf.RDS"
)

  # Beta ----

beta_split <- readRDS("modeling-data/beta/split.RDS")
names(beta_split) <- c(1:10)
list(
  "ntree" = 25, 
  "nodesize" = 5, 
  "mtry" = 3
) %>%
  train_rf(
    folds = beta_split, 
    param = ., 
    fold_avg = F
) # best is 2
trn <- readRDS("modeling-data/beta/split.RDS")[[2]]$trn
saveRDS(
  randomForest(
    x = trn[, -1:-2], 
    y = trn$dG, 
    ntree = 25, 
    nodesize = 2, 
    mtry = 1
  ),
  "model/beta/rf.RDS"
)

  # Gamma ----

gamma_split <- readRDS("modeling-data/gamma/split.RDS")
names(gamma_split) <- c(1:10)
list(
  "ntree" = 25, 
  "nodesize" = 10, 
  "mtry" = 8
) %>% 
  train_rf(
    folds = gamma_split, 
    param = ., 
    fold_avg = F
  )
# the best is split 8
trn <- readRDS("modeling-data/gamma/split.RDS")[[8]]$trn
saveRDS(
  randomForest(
    x = trn[, -1:-2], 
    y = trn$dG, 
    ntree = 25, 
    nodesize = 10, 
    mtry = 8
  ),
  "model/gamma/rf.RDS"
)