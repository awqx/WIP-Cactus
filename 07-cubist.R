source("helpers/training.R")

# Tuning ------------------------------------------------------------------

combo_grid <- expand.grid(
  c(5, 10, 20, 50, 75, 100), # committees 
  c(0, 5, 25, 50, 75, 100), # extrapolation
  c(5, 10, 25, 50, 80), # rules
  c(0, 5, 20, 40, 75, 95) # sample
)
names(combo_grid) <-
  c("committees", "extrapolation", "rules", "sample")
combo_list <- split(combo_grid, seq(nrow(combo_grid)))

  # Alpha ----

trn <- readRDS("modeling-data/alpha/split.RDS")[[10]][["trn"]]
fold_list <- split_train(trn, nfolds = 10, seed = 20200118)
alpha_tune <- do.call(
  rbind, 
  lapply(
    combo_list, 
    train_cubist, 
    folds = fold_list, 
    seed = 20200118
  )
)

test_param <- list(
  "committees" = 10, 
  "extrapolation" = 50, 
  "rules" = 10, 
  "sample" = 50
)
train_cubist(fold_list, test_param, 20200118)

  # Beta ----

trn <- readRDS("modeling-data/beta/split.RDS")[[10]][["trn"]]
fold_list <- split_train(trn, nfolds = 10, seed = 20200118)
beta_tune <- do.call(
  rbind, 
  lapply(
    combo_list, 
    train_cubist, 
    folds = fold_list, 
    seed = 20200118
  )
)

  # Gamma ----

# adjusting tuning parameters due to stalling on sample
combo_grid <- expand.grid(
  c(5, 10, 20, 50, 75, 100), # committees 
  c(0, 5, 25, 50, 75, 100), # extrapolation
  c(5, 10, 25, 50, 80), # rules
  c(0) # sample default
)
names(combo_grid) <-
  c("committees", "extrapolation", "rules", "sample")
combo_list <- split(combo_grid, seq(nrow(combo_grid)))
trn <- readRDS("modeling-data/gamma/split.RDS")[[10]][["trn"]]
# 10 folds did not work with gamma, so trying 5 instead
fold_list <- split_train(trn, nfolds = 10, seed = 20200118)
gamma_tune <- do.call(
  rbind, 
  lapply(
    combo_list, 
    train_cubist, 
    folds = fold_list, 
    seed = 20200118
  )
)

  # Save ----

saveRDS(alpha_tune, "training/alpha/cubist.RDS")
saveRDS(beta_tune, "training/beta/cubist.RDS")
saveRDS(gamma_tune, "training/gamma/cubist.RDS")

# Model -------------------------------------------------------------------

  # Alpha ----

# refer to `05-svm.R` for comments and further explanation
alpha_split <- readRDS("modeling-data/alpha/split.RDS")
names(alpha_split) <- c(1:10)
list(
  "committees" = 100,
  "extrapolation" = 80,
  "rules" = 80, 
  "sample" = 40
) %>%
  train_cubist(
    folds = alpha_split, 
    param = ., 
    fold_avg = F, 
    seed = 20200118
  ) # best is 4

# saving the best model
trn <- readRDS("modeling-data/alpha/split.RDS")[[4]]$trn
alpha_ctrl <- cubistControl(
  seed = 20200118, 
  extrapolation = 100,
  sample = 40
)
saveRDS(
  cubist(
    x = trn[, -1:-2], 
    y = trn$dG, 
    committees = 100,
    control = alpha_ctrl
  ),
  "model/alpha/cubist.RDS"
)

# Beta ----

beta_split <- readRDS("modeling-data/beta/split.RDS")
names(beta_split) <- c(1:10)
list(
  "committees" = 100,
  "extrapolation" = 5,
  "rules" = 25, 
  "sample" = 0
) %>%
  train_cubist(
    folds = beta_split, 
    param = ., 
    fold_avg = F, 
    seed = 20200118
  ) # best is 4

# saving the best model
trn <- readRDS("modeling-data/beta/split.RDS")[[2]]$trn
beta_ctrl <- cubistControl(
  seed = 20200118, 
  extrapolation = 100,
  sample = 0
)
saveRDS(
  cubist(
    x = trn[, -1:-2], 
    y = trn$dG, 
    committees = 100,
    control = beta_ctrl
  ),
  "model/beta/cubist.RDS"
)

# Gamma ----

gamma_split <- readRDS("modeling-data/gamma/split.RDS")
names(gamma_split) <- c(1:10)
list(
  "committees" = 20,
  "extrapolation" = 25,
  "rules" = 80, 
  "sample" = 0
) %>%
  train_cubist(
    folds = gamma_split, 
    param = ., 
    fold_avg = F, 
    seed = 20200118
  ) # best is 4

# saving the best model
trn <- readRDS("modeling-data/gamma/split.RDS")[[8]]$trn
gamma_ctrl <- cubistControl(
  seed = 20200118, 
  extrapolation = 25,
  sample = 0
)
saveRDS(
  cubist(
    x = trn[, -1:-2], 
    y = trn$dG, 
    committees = 20,
    control = gamma_ctrl
  ),
  "model/gamma/cubist.RDS"
)