source("helpers/tuning.R")
if(!dir.exists("tuning")) {
  dir.create("tuning")
  dir.create("tuning/alpha")
  dir.create("tuning/beta")
  dir.create("tuning/gamma")
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
x <- trn[, -1:-2]
y <- trn[, 2]
# create folds in training set
# due to the number of parameters, instead of taking 10 folds, 
# calculate 6 folds
nfold <- 6
set.seed(2020)
fold_list <- createFolds(y, k = nfold)

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
    tune_svm,
    x = x, 
    y = y, 
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
    tune_svm,
    x = x, 
    y = y, 
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
    tune_svm,
    x = x, 
    y = y, 
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
    tune_svm,
    x = x, 
    y = y, 
    kernel = "linear", 
    folds = fold_list
  )
)

  # Beta ----

    # Importing data ----

trn <- readRDS("modeling-data/beta/split.RDS")[[10]][["trn"]]
x <- trn[, -1:-2]
y <- trn[, 2]

nfold <- 6
set.seed(2020)
fold_list <- createFolds(y, k = nfold)

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
    tune_svm,
    x = x, 
    y = y, 
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
    tune_svm,
    x = x, 
    y = y, 
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
    tune_svm,
    x = x,
    y = y,
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
    tune_svm,
    x = x, 
    y = y, 
    kernel = "linear", 
    folds = fold_list
  )
)
  # Gamma ----

    # Importing data ----

trn <- readRDS("modeling-data/gamma/split.RDS")[[10]][["trn"]]
x <- trn[, -1:-2]
y <- trn[, 2]

nfold <- 6
set.seed(2020)
fold_list <- createFolds(y, k = nfold)

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
    tune_svm,
    x = x, 
    y = y, 
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
    tune_svm,
    x = x, 
    y = y, 
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
    tune_svm,
    x = x,
    y = y,
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
    tune_svm,
    x = x, 
    y = y, 
    kernel = "linear", 
    folds = fold_list
  )
)

  # Saving ----

saveRDS(alpha_lin_tune, "tuning/alpha/svm-lin.RDS")
saveRDS(alpha_poly_tune, "tuning/alpha/svm-poly.RDS")
saveRDS(alpha_rbf_tune, "tuning/alpha/svm-rbf.RDS")
saveRDS(alpha_sig_tune, "tuning/alpha/svm-sig.RDS")
saveRDS(beta_lin_tune, "tuning/beta/svm-lin.RDS")
saveRDS(beta_poly_tune, "tuning/beta/svm-poly.RDS")
saveRDS(beta_rbf_tune, "tuning/beta/svm-rbf.RDS")
saveRDS(beta_sig_tune, "tuning/beta/svm-sig.RDS")
saveRDS(gamma_lin_tune, "tuning/gamma/svm-lin.RDS")
saveRDS(gamma_poly_tune, "tuning/gamma/svm-poly.RDS")
saveRDS(gamma_rbf_tune, "tuning/gamma/svm-rbf.RDS")
saveRDS(gamma_sig_tune, "tuning/gamma/svm-sig.RDS")

# Model -------------------------------------------------------------------


