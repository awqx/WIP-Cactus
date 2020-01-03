# Basic preprocessing of descriptors, including centering, scaling, and
# removal of highly correlated variables
source("helpers/preprocess.R")
source("helpers/applicability.R")
# if(!dir.exists("preprocess")) {
#   dir.create("preprocess")
#   dir.create("preprocess/outliers")
# }
if(!dir.exists("modeling-data")) {
  dir.create("modeling-data")
  dir.create("modeling-data/alpha")
  dir.create("modeling-data/beta")
  dir.create("modeling-data/gamma")
  dir.create("modeling-data/unprocessed")
}
if(!dir.exists("ext-val")) dir.create("ext-val")

# Activity outliers ------------------------------------------------------

alpha <- readRDS("descriptors/alpha-dg.RDS")
alpha <- alpha[-find_outlier(alpha, 2.5), ] # -2 outliers
beta <- readRDS("descriptors/beta-dg.RDS")
beta <- beta[-find_outlier(beta, 2.5), ] # - 10 outliers
gamma <- readRDS("descriptors/gamma-dg.RDS")
gamma <- gamma[-find_outlier(gamma, 2.5), ] # -1 outlier

# External validation -----------------------------------------------------

# Citing Tropsha, Best practices for QSAR (DOI:10.1002/minf.201000061)
# An external validation set (constituting 10-20% of the original data) 
# should be created to be tested after model building. \

set.seed(101)
alpha_ev_index <- sample(nrow(alpha), round(nrow(alpha) * 0.15))
beta_ev_index <- sample(nrow(beta), round(nrow(beta) * 0.15))
gamma_ev_index <- sample(nrow(gamma), round(nrow(gamma) * 0.15))
saveRDS(alpha[alpha_ev_index, ], "ext-val/alpha-ev.RDS")
saveRDS(beta[beta_ev_index, ], "ext-val/beta-ev.RDS")
saveRDS(gamma[gamma_ev_index, ], "ext-val/gamma-ev.RDS")
alpha <- alpha[-alpha_ev_index, ]
beta <- beta[-beta_ev_index, ]
gamma <- gamma[-gamma_ev_index, ]
saveRDS(alpha, "modeling-data/alpha/unprocessed.RDS")
saveRDS(beta, "modeling-data/beta/unprocessed.RDS")
saveRDS(gamma, "modeling-data/gamma/unprocessed.RDS")

# Structural outliers -----------------------------------------------------
    # Alpha ----
# requires preprocessing of the values
# separating guest and dG b/c preprocessing functions require all numeric
alpha_guest_dg <- select(alpha, guest, dG)
alpha <- select(alpha, -guest, -dG)
# centering and scaling
alpha_pp <- center_scale(alpha)[[1]]
# remove variables with near zero variance
alpha_zero_var <- nearZeroVar(alpha_pp)
alpha_pp <- alpha_pp[, -alpha_zero_var]
# removing highly correlated variables
alpha_high_cor <- findCorrelation(
  cor(alpha_pp, use = "pairwise.complete.obs"), 
  cutoff = 0.95
)
alpha_pp <- alpha_pp[, -alpha_high_cor]
# replacing the names of the guests
alpha_pp_guest <- data.frame(alpha_guest_dg, alpha_pp)
# testing for outliers
alpha_outliers_guests <- alpha_pp_guest %>%
  select(-dG) %>% # dG is not a predictor
  domain_num() %>%
  filter(domain == "outside") %>%
  .$guest # 5 outliers
# table of outliers, to potentially experiment with later
alpha_outliers <- filter(
  alpha_pp_guest, 
  guest %in% alpha_outliers_guests
  )
alpha_pp <- filter(
  alpha_pp_guest, 
  !guest %in% alpha_outliers
  )

    # Beta ----

beta_guest_dg <- select(beta, guest, dG)
beta <- select(beta, -guest, -dG)
beta_pp <- center_scale(beta)[[1]]
beta_pp <- beta_pp[, -nearZeroVar(beta_pp)]
beta_high_cor <- findCorrelation(
  cor(beta_pp, use = "pairwise.complete.obs"), 
  cutoff = 0.95
)
beta_pp <- beta_pp[, -beta_high_cor]
beta_pp_guest <- data.frame(beta_guest_dg, beta_pp)

beta_outliers_guests <- beta_pp_guest %>%
  select(-dG) %>% 
  domain_num() %>%
  filter(domain == "outside") %>%
  .$guest # 5 outliers
beta_outliers <- filter(
  beta_pp_guest, 
  guest %in% beta_outliers_guests
)
beta_pp <- filter(
  beta_pp_guest, 
  !guest %in% beta_outliers
)
    
    # Gamma ----

gamma_guest_dg <- select(gamma, guest, dG)
gamma <- select(gamma, -guest, -dG)
gamma_pp <- center_scale(gamma)[[1]]
gamma_pp <- gamma_pp[, -nearZeroVar(gamma_pp)]
gamma_high_cor <- findCorrelation(
  cor(gamma_pp, use = "pairwise.complete.obs"), 
  cutoff = 0.95
)
gamma_pp <- gamma_pp[, -gamma_high_cor]
gamma_pp_guest <- data.frame(gamma_guest_dg, gamma_pp)

gamma_outliers_guests <- gamma_pp_guest %>%
  select(-dG) %>% 
  domain_num() %>%
  filter(domain == "outside") %>%
  .$guest # 5 outliers
gamma_outliers <- filter(
  gamma_pp_guest, 
  guest %in% gamma_outliers_guests
)
gamma_pp <- filter(
  gamma_pp_guest, 
  !guest %in% gamma_outliers
)

# Feature selection -------------------------------------------------------

set.seed(10104)
nsplits <- 10
alpha_splits <- split_train_test(nsplits, alpha_pp)
beta_splits <- split_train_test(nsplits, beta_pp)
gamma_splits <- split_train_test(nsplits, gamma_pp)

    # Alpha ----

# running recursive feature elimination on the alpha data
alpha_rfe <- lapply(
  alpha_splits, 
  use_rfe
)
# retrieiving a list of all predictors that appeared (w/ repeats)
alpha_all_predictors <- unlist(
  lapply(
    alpha_rfe, 
    predictors
  ) 
)
# sorting based on number of occurences
alpha_patterns <- paste0("^", unique(alpha_all_predictors), "$")
alpha_predictors <- lapply(
  alpha_patterns, 
  str_count, 
  alpha_all_predictors
) %>%
  lapply(sum) %>%
  unlist()
# selecting the features that appear in all runs. 
# RFE replaces dashes in predictor names with dots, so str_replace 
# reverts this change
alpha_features <- 
  unique(alpha_all_predictors)[alpha_predictors == nsplits] %>%
  str_replace("\\.", "-")
saveRDS(alpha_features, "modeling-data/alpha/predictors.RDS")

    # Beta ----

beta_rfe <- lapply(
  beta_splits, 
  use_rfe
)
beta_all_predictors <- unlist(
  lapply(
    beta_rfe, 
    predictors
  ) 
)
beta_patterns <- paste0("^", unique(beta_all_predictors), "$")
beta_predictors <- lapply(
  beta_patterns, 
  str_count, 
  beta_all_predictors
) %>%
  lapply(sum) %>%
  unlist()

beta_features <- 
  unique(beta_all_predictors)[beta_predictors == nsplits] %>%
  str_replace("\\.", "-")
saveRDS(beta_features, "modeling-data/beta/predictors.RDS")

    # Gamma ----

gamma_rfe <- lapply(
  gamma_splits, 
  use_rfe
)
gamma_all_predictors <- unlist(
  lapply(
    gamma_rfe, 
    predictors
  ) 
)
gamma_patterns <- paste0("^", unique(gamma_all_predictors), "$")
gamma_predictors <- lapply(
  gamma_patterns, 
  str_count, 
  gamma_all_predictors
) %>%
  lapply(sum) %>%
  unlist()

gamma_features <- 
  unique(gamma_all_predictors)[gamma_predictors >= (nsplits - 1)] %>%
  str_replace("\\.", "-")
saveRDS(gamma_features, "modeling-data/gamma/predictors.RDS")

# Pre-processing ----------------------------------------------------------

# redoing centering and scaling with features from previous step.
# removing zero variance and highly correlated descriptors not necessary
# as those were removed for feature selection to work

# this step may seem redundant, but it will make predicting on new 
# molecules much easier, as the descriptors will be quickly narrowed down
# to the relevant features and will then only require a center-scale.

# reading features saved from previous step.
# this becomes useful as RFE may span multiple work sessions.
alpha <- readRDS("modeling-data/alpha/unprocessed.RDS") %>%
  select(readRDS("modeling-data/alpha/predictors.RDS")) %>%
  center_scale()
alpha_pp <- data.frame(alpha_guest_dg, alpha[[1]])
saveRDS(alpha_pp, "modeling-data/alpha/preprocessed.RDS")
saveRDS(alpha[[2]], "modeling-data/alpha/preprocess-settings.RDS")

beta <- readRDS("modeling-data/beta/unprocessed.RDS")
beta <- beta %>%
  select(readRDS("modeling-data/beta/predictors.RDS")) %>%
  center_scale()
beta_pp <- data.frame(beta_guest_dg, beta[[1]])
saveRDS(beta_pp, "modeling-data/beta/preprocessed.RDS")
saveRDS(beta[[2]], "modeling-data/beta/preprocess-settings.RDS")

gamma <- readRDS("modeling-data/gamma/unprocessed.RDS") %>%
  select(readRDS("modeling-data/gamma/predictors.RDS")) %>%
  center_scale()
gamma_pp <- data.frame(gamma_guest_dg, gamma[[1]])
saveRDS(gamma_pp, "modeling-data/gamma/preprocessed.RDS")
saveRDS(gamma[[2]], "modeling-data/gamma/preprocess-settings.RDS")

# Splitting Data ----------------------------------------------------------

# Multiple combinations of test and train sets should be created in order
# to fully validate the models. No specification was made in the paper
# as to how many different splits should be created, exactly, so I
# decided 10 sets would be created

save_splits(10, alpha_pp, "modeling-data/alpha/split")
save_splits(10, beta_pp, "modeling-data/beta/split")
save_splits(10, gamma_pp, "modeling-data/gamma/split")