source("helpers/training.R")
if (!dir.exists("training")) {
  dir.create("training")
  dir.create("training/alpha")
  dir.create("training/beta")
  dir.create("training/gamma")
}
if (!dir.exists("model")) {
  dir.create("model")
  dir.create("model/alpha")
  dir.create("model/beta")
  dir.create("model/gamma")
}

# I'm tuning kernels separately to more easily run the code. Each kernel
# requires a different number of parameters

# Tuning ------------------------------------------------------------------

  # Alpha ----
  
    # Importing data ----

# tuning will be performed on one of the splits. Though this is less 
# optimal than tuning on all 10 splits, the latter is too resource-
# intensive. It may be possible for simple models with few parameters, 
# such as GLM.
trn <- readRDS("modeling-data/alpha/split.RDS")[[10]][["trn"]]
fold_list <- split_train(trn, nfolds = 10, seed = 20200116)

    # Polynomial ----

# the relevant parameters are coef0, cost, degree, epsilon, gamma.
# generate all possible combinations of a reasonable range for each
# parameters to test
coef0 <- c(-10, -5, -1, 0, 1, 5, 10)
cost <- c(1:5, 10, 20)
degree <- c(1:4)
epsilon <- c(0.01, 0.05, 0.1, 0.25, 0.5)
gamma <- c(0.001, 0.01, 0.05, 0.1, 0.25, 0.5)
combo_grid <- expand.grid(coef0, cost, degree, epsilon, gamma)
names(combo_grid) <- c("coef0", "cost", "degree", "epsilon", "gamma")
combo_list <- split(
  combo_grid, 
  seq(nrow(combo_grid))
)
alpha_poly_tune <- do.call(
  rbind, 
  lapply(
    combo_list, 
    train_svm,
    kernel = "polynomial", 
    folds = fold_list
  )
)

    # Radial ----

# parameters: cost, epsilon, gamma
cost <- c(1:5, 10, 20, 40)
epsilon <- c(0.001, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75)
gamma <- c(0.001, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9)
combo_grid <- expand.grid(cost, epsilon, gamma)
names(combo_grid) <- c("cost", "epsilon", "gamma")
combo_list <- split(
  combo_grid, 
  seq(nrow(combo_grid))
)
alpha_rbf_tune <- do.call(
  rbind, 
  lapply(
    combo_list, 
    train_svm,
    kernel = "radial", 
    folds = fold_list
  )
)

    # Sigmoid ----

# parameters: coef0, cost, epsilon, gamma
coef0 <- c(0, 1:5, 10, 20)
cost <- c(1:5, 10, 20)
epsilon <- c(0.01, 0.05, 0.1, 0.25, 0.5)
gamma <- c(0.001, 0.01, 0.05, 0.1, 0.25, 0.5)
combo_grid <- expand.grid(coef0, cost, epsilon, gamma)
names(combo_grid) <- c("coef0", "cost", "epsilon", "gamma")
combo_list <- split(
  combo_grid, 
  seq(nrow(combo_grid))
)
alpha_sig_tune <- do.call(
  rbind, 
  lapply(
    combo_list, 
    train_svm,
    kernel = "sigmoid", 
    folds = fold_list
  )
)

    # Linear ----

# parameters: cost, epsilon
cost <- c(1:5, 10, 20, 50, 100)
epsilon <- c(0.01, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9)
combo_grid <- expand.grid(cost, epsilon)
names(combo_grid) <- c("cost", "epsilon")
combo_list <- split(
  combo_grid, 
  seq(nrow(combo_grid))
)
alpha_lin_tune <- do.call(
  rbind, 
  lapply(
    combo_list, 
    train_svm,
    kernel = "linear", 
    folds = fold_list
  )
)

  # Beta ----

    # Importing data ----

trn <- readRDS("modeling-data/beta/split.RDS")[[10]][["trn"]]
fold_list <- split_train(trn, nfolds = 10, seed = 20200116)

    # Polynomial ----

coef0 <- c(-10, -5, -1, 0, 1, 5, 10)
cost <- c(1:5, 10, 20)
degree <- c(1:4)
epsilon <- c(0.01, 0.05, 0.1, 0.25, 0.5)
gamma <- c(0.001, 0.01, 0.05, 0.1, 0.25, 0.5)
combo_grid <- expand.grid(coef0, cost, degree, epsilon, gamma)
names(combo_grid) <- c("coef0", "cost", "degree", "epsilon", "gamma")
combo_list <- split(
  combo_grid, 
  seq(nrow(combo_grid))
)
beta_poly_tune <- do.call(
  rbind, 
  lapply(
    combo_list, 
    train_svm,
    kernel = "polynomial", 
    folds = fold_list
  )
)

    # Radial ----

# parameters: cost, epsilon, gamma
cost <- c(1:5, 10, 20, 40)
epsilon <- c(0.001, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75)
gamma <- c(0.001, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9)
combo_grid <- expand.grid(cost, epsilon, gamma)
names(combo_grid) <- c("cost", "epsilon", "gamma")
combo_list <- split(
  combo_grid, 
  seq(nrow(combo_grid))
)
beta_rbf_tune <- do.call(
  rbind, 
  lapply(
    combo_list, 
    train_svm,
    kernel = "radial", 
    folds = fold_list
  )
)

    # Sigmoid ----

# parameters: coef0, cost, epsilon, gamma
coef0 <- c(-10, -5, -1, 0, 1, 5, 10)
cost <- c(1:5, 10, 20)
epsilon <- c(0.01, 0.05, 0.1, 0.25, 0.5)
gamma <- c(0.001, 0.01, 0.05, 0.1, 0.25, 0.5)
combo_grid <- expand.grid(coef0, cost, epsilon, gamma)
names(combo_grid) <- c("coef0", "cost", "epsilon", "gamma")
combo_list <- split(
  combo_grid,
  seq(nrow(combo_grid))
)
beta_sig_tune <- do.call(
  rbind,
  lapply(
    combo_list,
    train_svm,
    kernel = "sigmoid",
    folds = fold_list
  )
)

    # Linear ----

# parameters: cost, epsilon
cost <- c(1:5, 10, 20, 50, 100)
epsilon <- c(0.01, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9)
combo_grid <- expand.grid(cost, epsilon)
names(combo_grid) <- c("cost", "epsilon")
combo_list <- split(
  combo_grid, 
  seq(nrow(combo_grid))
)
beta_lin_tune <- do.call(
  rbind, 
  lapply(
    combo_list, 
    train_svm,
    kernel = "linear", 
    folds = fold_list
  )
)
  # Gamma ----

    # Importing data ----

trn <- readRDS("modeling-data/gamma/split.RDS")[[10]][["trn"]]
fold_list <- split_train(trn, nfolds = 10, seed = 20200116)

    # Polynomial ----

coef0 <- c(-10, -5, -1, 0, 1, 5, 10)
cost <- c(1:5, 10, 20)
degree <- c(1:4)
epsilon <- c(0.01, 0.05, 0.1, 0.25, 0.5)
gamma <- c(0.001, 0.01, 0.05, 0.1, 0.25, 0.5)
combo_grid <- expand.grid(coef0, cost, degree, epsilon, gamma)
names(combo_grid) <- c("coef0", "cost", "degree", "epsilon", "gamma")
combo_list <- split(
  combo_grid, 
  seq(nrow(combo_grid))
)
gamma_poly_tune <- do.call(
  rbind, 
  lapply(
    combo_list, 
    train_svm,
    kernel = "polynomial", 
    folds = fold_list
  )
)

    # Radial ----

# parameters: cost, epsilon, gamma
cost <- c(1:5, 10, 20, 40)
epsilon <- c(0.001, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75)
gamma <- c(0.001, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9)
combo_grid <- expand.grid(cost, epsilon, gamma)
names(combo_grid) <- c("cost", "epsilon", "gamma")
combo_list <- split(
  combo_grid, 
  seq(nrow(combo_grid))
)
gamma_rbf_tune <- do.call(
  rbind, 
  lapply(
    combo_list, 
    train_svm,
    kernel = "radial", 
    folds = fold_list
  )
)

    # Sigmoid ----

# parameters: coef0, cost, epsilon, gamma
coef0 <- c(-10, -5, -1, 0, 1, 5, 10)
cost <- c(1:5, 10, 20)
epsilon <- c(0.01, 0.05, 0.1, 0.25, 0.5)
gamma <- c(0.001, 0.01, 0.05, 0.1, 0.25, 0.5)
combo_grid <- expand.grid(coef0, cost, epsilon, gamma)
names(combo_grid) <- c("coef0", "cost", "epsilon", "gamma")
combo_list <- split(
  combo_grid,
  seq(nrow(combo_grid))
)
gamma_sig_tune <- do.call(
  rbind,
  lapply(
    combo_list,
    train_svm,
    kernel = "sigmoid",
    folds = fold_list
  )
)

    # Linear ----

# parameters: cost, epsilon
cost <- c(1:5, 10, 20, 50, 100)
epsilon <- c(0.01, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9)
combo_grid <- expand.grid(cost, epsilon)
names(combo_grid) <- c("cost", "epsilon")
combo_list <- split(
  combo_grid, 
  seq(nrow(combo_grid))
)
gamma_lin_tune <- do.call(
  rbind, 
  lapply(
    combo_list, 
    train_svm,
    kernel = "linear", 
    folds = fold_list
  )
)

  # Saving ----

saveRDS(alpha_lin_tune, "training/alpha/svm-lin.RDS")
saveRDS(alpha_poly_tune, "training/alpha/svm-poly.RDS")
saveRDS(alpha_rbf_tune, "training/alpha/svm-rbf.RDS")
saveRDS(alpha_sig_tune, "training/alpha/svm-sig.RDS")
saveRDS(beta_lin_tune, "training/beta/svm-lin.RDS")
saveRDS(beta_poly_tune, "training/beta/svm-poly.RDS")
saveRDS(beta_rbf_tune, "training/beta/svm-rbf.RDS")
saveRDS(beta_sig_tune, "training/beta/svm-sig.RDS")
saveRDS(gamma_lin_tune, "training/gamma/svm-lin.RDS")
saveRDS(gamma_poly_tune, "training/gamma/svm-poly.RDS")
saveRDS(gamma_rbf_tune, "training/gamma/svm-rbf.RDS")
saveRDS(gamma_sig_tune, "training/gamma/svm-sig.RDS")

# Model -------------------------------------------------------------------

  # Alpha ----

# determing split with best performance --
# we use the results of tuning to inform model building on each fold
# of the training set. we then select the best performing model as the
# final model to go toward the ensemble.
alpha_split <- readRDS("modeling-data/alpha/split.RDS")
names(alpha_split) <- c(1:10)
list("cost" = 3,
     "epsilon" = 0.1,
     "gamma" = 0.25) %>%
  train_svm(
    folds = alpha_split,
    kernel = "radial",
    param = .,
    fold_avg = F
  ) # best is 3
trn <- readRDS("modeling-data/alpha/split.RDS")[[3]]$trn
saveRDS(
  svm(
  x = trn[, -1:-2], 
  y = trn$dG, 
  cost = 3, 
  epsilon = 0.1, 
  gamma = 0.25
  ),
  "model/alpha/svm.RDS"
)

  # Beta ----

beta_split <- readRDS("modeling-data/beta/split.RDS")
names(beta_split) <- c(1:10)
# because the highest r2 occurred between cost = 20 and cost = 40, 
# i tested values until settling on 30
list(
  "cost" = 30, 
  "epsilon" = 0.01, 
  "gamma" = 0.01
) %>%
  train_svm(
    folds = beta_split, 
    kernel = "radial", 
    param = ., 
    fold_avg = F
  ) # split 9
trn <- readRDS("modeling-data/beta/split.RDS")[[9]]$trn
saveRDS(
  svm(
    x = trn[, -1:-2], 
    y = trn$dG, 
    cost = 30, 
    epsilon = 0.01, 
    gamma = 0.01
  ),
  "model/beta/svm.RDS"
)

  # Gamma ----

gamma_split <- readRDS("modeling-data/gamma/split.RDS")
names(gamma_split) <- c(1:10)
list(
  "cost" = 50, 
  "epsilon" = 0.01, 
  "gamma" = 0.01
) %>%
  train_svm(
    folds = gamma_split, 
    kernel = "radial", 
    param = ., 
    fold_avg = F
  )
# none of the splits reached r2 > 0.6, so nothing will be saved
# split 7 is the best but only has r2 of 0.45
# trn <- readRDS("modeling-data/gamma/split.RDS")[[7]]$trn
# saveRDS(
#   svm(
#     x = trn[, -1:-2], 
#     y = trn$dG, 
#     cost = 1, 
#     epsilon = 0.1, 
#     gamma = 0.01
#   ),
#   "model/gamma/svm.RDS"
# )