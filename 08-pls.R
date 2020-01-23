source("helpers/training.R")

trn <- readRDS("modeling-data/alpha/split.RDS")[[10]][["trn"]]
tst <- readRDS("modeling-data/alpha/split.RDS")[[10]][["tst"]]
tst <- tst[, -1:-2]
trn <- trn[, -1]
pls_model <- plsr(
  dG~., 
  data = trn,
  method = "kernelpls", 
  ncomp = 1
)
pls_model <- do.call(
  plsr, 
  args = 
    list(
      "formula" = dG~., 
      "data" = trn, 
      "method" = "kernelpls", 
      "ncomp" = 1
    )
)
predict(pls_model, tst)

# Tuning ------------------------------------------------------------------

  # Alpha ----

combo_grid <- expand.grid(
  c("kernelpls", "simpls", "oscorespls", "widekernelpls"), # method
  c(1:8), # ncomp
  stringsAsFactors = F
)
names(combo_grid) <- c("method", "ncomp")
combo_list <- split(combo_grid, seq(nrow(combo_grid)))

trn <- readRDS("modeling-data/alpha/split.RDS")[[10]][["trn"]]
fold_list <- split_train(trn, nfolds = 10, seed = 20200123)
alpha_tune <- do.call(
  rbind, 
  lapply(
    combo_list, 
    train_pls, 
    folds = fold_list
  )
)

  # Beta ----

combo_grid <- expand.grid(
  c("kernelpls", "simpls", "oscorespls", "widekernelpls"), # method
  c(1:17), # ncomp
  stringsAsFactors = F
)
names(combo_grid) <- c("method", "ncomp")
combo_list <- split(combo_grid, seq(nrow(combo_grid)))

trn <- readRDS("modeling-data/beta/split.RDS")[[10]][["trn"]]
fold_list <- split_train(trn, nfolds = 10, seed = 20200123)
beta_tune <- do.call(
  rbind, 
  lapply(
    combo_list, 
    train_pls, 
    folds = fold_list
  )
)

  # Gamma ----

combo_grid <- expand.grid(
  c("kernelpls", "simpls", "oscorespls", "widekernelpls"), # method
  c(1:8, 11, 15, 20), # ncomp
  stringsAsFactors = F
)
names(combo_grid) <- c("method", "ncomp")
combo_list <- split(combo_grid, seq(nrow(combo_grid)))

trn <- readRDS("modeling-data/gamma/split.RDS")[[10]][["trn"]]
fold_list <- split_train(trn, nfolds = 10, seed = 20200123)
gamma_tune <- do.call(
  rbind, 
  lapply(
    combo_list, 
    train_pls, 
    folds = fold_list
  )
)

  # Save ----

saveRDS(alpha_tune, "training/alpha/pls.RDS")
saveRDS(beta_tune, "training/beta/pls.RDS")
saveRDS(gamma_tune, "training/gamma/pls.RDS")

# Models ------------------------------------------------------------------

  # Alpha ----

alpha_split <- readRDS("modeling-data/alpha/split.RDS")
names(alpha_split) <- c(1:10)
list(
  "ncomp" = 2, 
  "kernel" = "simpls"
) %>%
  train_pls(
    folds = alpha_split, 
    param = ., 
    fold_avg = F
  ) # best is 10
# saving the best model
trn <- readRDS("modeling-data/alpha/split.RDS")[[10]]$trn
saveRDS(
  plsr(
    dG~.,
    data = trn[, -1],
    ncomp = 2, 
    kernel = "simpls"
  ),
  "model/alpha/pls.RDS"
)

  # Beta ----

beta_split <- readRDS("modeling-data/beta/split.RDS")
names(beta_split) <- c(1:10)
list(
  "ncomp" = 2, 
  "kernel" = "simpls"
) %>%
  train_pls(
    folds = beta_split, 
    param = ., 
    fold_avg = F
  ) # best is 6
# saving the best model
trn <- readRDS("modeling-data/beta/split.RDS")[[6]]$trn
saveRDS(
  plsr(
    dG~.,
    data = trn[, -1],
    ncomp = 2, 
    kernel = "simpls"
  ),
  "model/beta/pls.RDS"
)

  # Gamma ----

gamma_split <- readRDS("modeling-data/gamma/split.RDS")
names(gamma_split) <- c(1:10)
list(
  "ncomp" = 1, 
  "kernel" = "simpls"
) %>%
  train_pls(
    folds = gamma_split, 
    param = ., 
    fold_avg = F
  ) # best is 2
# saving the best model
trn <- readRDS("modeling-data/gamma/split.RDS")[[2]]$trn
saveRDS(
  plsr(
    dG~.,
    data = trn[, -1],
    ncomp = 2, 
    kernel = "simpls"
  ),
  "model/gamma/pls.RDS"
)