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
alpha_guest_dg <- select(alpha, guest, dG)
alpha <- select(alpha, -guest, -dG)
alpha_pp_settings <- preProcess(
  alpha, 
  na.remove = T, 
  method = c("center", "scale"), 
  verbose = F
)
alpha_pp <- predict(alpha_pp_settings, alpha)
alpha_zero_var <- nearZeroVar(alpha_pp)
alpha_pp <- alpha_pp[, -alpha_zero_var]
alpha_high_cor <- findCorrelation(
  cor(alpha_pp, use = "pairwise.complete.obs"), 
  cutoff = 0.95
)
alpha_pp <- alpha_pp[, -alpha_high_cor]
alpha_pp_guest <- data.frame(alpha_guest_dg, alpha_pp)
alpha_outliers_guests <- alpha_pp_guest %>%
  domain_num() %>%
  filter(domain == "outside") %>%
  .$guest # 6 outliers
alpha_outliers <- filter(
  alpha_pp_guest, 
  guest %in% alpha_outliers_guests)
alpha_modeling_data <- filter(
  alpha_pp_guest, 
  !guest %in% alpha_outliers)

saveRDS(alpha_pp_settings, "modeling-data/alpha/pp-settings.RDS")
saveRDS(alpha_zero_var, "modeling-data/alpha/zero-var.RDS")
saveRDS(alpha_high_cor, "modeling-data/alpha/high-cor.RDS")
saveRDS(alpha_outliers, "modeling-data/alpha/outliers.RDS")
saveRDS(alpha_modeling_data, "modeling-data/alpha/all.RDS")

    # Beta ----
beta_guest_dg <- select(beta, guest, dG)
beta <- select(beta, -guest, -dG)
beta_pp_settings <- preProcess(
  beta, 
  na.remove = T, 
  method = c("center", "scale"), 
  verbose = F
)
beta_pp <- predict(beta_pp_settings, beta)
beta_zero_var <- nearZeroVar(beta_pp)
beta_pp <- beta_pp[, -beta_zero_var]
beta_high_cor <- findCorrelation(
  cor(beta_pp, use = "pairwise.complete.obs"), 
  cutoff = 0.95
)
beta_pp <- beta_pp[, -beta_high_cor]
beta_pp_guest <- data.frame(beta_guest_dg, beta_pp)
beta_outliers_guests <- beta_pp_guest %>%
  domain_num() %>%
  filter(domain == "outside") %>%
  .$guest # 6 outliers
beta_outliers <- filter(
  beta_pp_guest, 
  guest %in% beta_outliers_guests)
beta_modeling_data <- filter(
  beta_pp_guest, 
  !guest %in% beta_outliers_guests)

saveRDS(beta_pp_settings, "modeling-data/beta/pp-settings.RDS")
saveRDS(beta_zero_var, "modeling-data/beta/zero-var.RDS")
saveRDS(beta_high_cor, "modeling-data/beta/high-cor.RDS")
saveRDS(beta_outliers, "modeling-data/beta/outliers.RDS")
saveRDS(beta_modeling_data, "modeling-data/beta/all.RDS")
    
    # Gamma ----
gamma_guest_dg <- select(gamma, guest, dG)
gamma <- select(gamma, -guest, -dG)
gamma_pp_settings <- preProcess(
  gamma, 
  na.remove = T, 
  method = c("center", "scale"), 
  verbose = F
)
gamma_pp <- predict(gamma_pp_settings, gamma)
gamma_zero_var <- nearZeroVar(gamma_pp)
gamma_pp <- gamma_pp[, -gamma_zero_var]
gamma_high_cor <- findCorrelation(
  cor(gamma_pp, use = "pairwise.complete.obs"), 
  cutoff = 0.95
)
gamma_pp <- gamma_pp[, -gamma_high_cor]
gamma_pp_guest <- data.frame(gamma_guest_dg, gamma_pp)
gamma_outliers_guests <- gamma_pp_guest %>%
  domain_num() %>%
  filter(domain == "outside") %>%
  .$guest # 6 outliers
gamma_outliers <- filter(gamma_pp_guest, guest %in% gamma_outliers_guests)
gamma_modeling_data <- filter(
  gamma_pp_guest, 
  !guest %in% gamma_outliers_guests)

saveRDS(gamma_pp_settings, "modeling-data/gamma/pp-settings.RDS")
saveRDS(gamma_zero_var, "modeling-data/gamma/zero-var.RDS")
saveRDS(gamma_high_cor, "modeling-data/gamma/high-cor.RDS")
saveRDS(gamma_outliers, "modeling-data/gamma/outliers.RDS")
saveRDS(gamma_modeling_data, "modeling-data/gamma/all.RDS")

# Splitting Data ----------------------------------------------------------

# Multiple combinations of test and train sets should be created in order
# to fully validate the models. No specification was made in the paper
# as to how many different splits should be created, exactly, so I
# decided 10 sets would be created

split_train_test(10, alpha_modeling_data, "modeling-data/alpha/split")
split_train_test(10, beta_modeling_data, "modeling-data/beta/split")
split_train_test(10, gamma_modeling_data, "modeling-data/gamma/split")