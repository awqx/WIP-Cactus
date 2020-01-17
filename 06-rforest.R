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
