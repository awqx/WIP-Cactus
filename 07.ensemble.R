source("07.model.functions.R")

# Libraries and Packages -------------------------------------------------

if(!require("pacman")) { 
  install.packages("pacman")
  library(pacman)
} else
  library(pacman)

p_load(caret,
       Cubist,
       e1071, 
       earth,
       gbm,
       glmnet,
       kernlab, 
       # pls, 
       randomForest, 
       tidyverse
       )

# Functions ---------------------------------------------------------------

# df: data.frame containing names of guests as well as 1376 PaDEL-desc
predict.alpha <- function(df) {
  outliers <- c()
  features <- readRDS("./feature.selection/alpha.vars.RDS")
  name.models <- list.files("./models/alpha") %>% str_remove(".RDS$")
  path.models <- paste0("./models/alpha/", list.files("models/alpha"))
  results <- data.frame()
  
  for(n in 1:length(name.models)) {
    pp.settings <- readRDS(path.models[n])[[1]]
    pp.results <- preprocess.desc(df, features, pp.settings)
    outliers <- c(outliers, pp.results[[2]])
    desc.pp <- pp.results[[1]]
    guests <- desc.pp$guests
    desc.pp <- desc.pp[ , -1]
    qsar <- readRDS(path.models[n])[[2]]
    if(str_detect(name.models[n], "gbm")) 
      pred <- predict(qsar, desc.pp, n.trees = 500)
    # else if(str_detect(name.models[n], "cubist"))
    #   pred <- Cubist::predict(qsar, desc.pp)
    else
      pred <- predict(qsar, desc.pp)
    results.qsar <- data.frame(guests, pred) 
    colnames(results.qsar)[2] <- name.models[n]
    if(n == 1) # initialize the data.frame if it's the first
      results <- results.qsar
    else
      results <- inner_join(results, results.qsar, by = "guests")
  }
  outliers <- unique(outliers)
  full.results <- results %>% 
    mutate(ensemble = rowMeans(results[ , -1]))
  predictions <- full.results %>% select(., guests, ensemble) %>%
    rename(dG = ensemble)
  
  return(list(predictions, full.results, outliers))
}

# Predicting using a list of models and features
# the list of models should be named
predict.ensemble <- function(df, feat, models) {
  
  # Initializing some variables
  outliers <- c()
  results <- data.frame()
  
  for(n in 1:length(models)) {
    message(names(models)[n], " starting")
    pp.settings <- models[[n]][[1]]
    pp.results <- preprocess.desc(df, feat, pp.settings)
    outliers <- c(outliers, pp.results[[2]])
    desc.pp <- pp.results[[1]]
    guests <- desc.pp$guests
    desc.pp <- desc.pp[ , -1]
    
    qsar <- models[[n]][[2]]
    # if statement acceptable here b/c there is only one case to evaluate
    # would use switch otherwise
    if(str_detect(names(models)[n], "gbm")) 
      pred <- predict(qsar, desc.pp, n.trees = 500)
    else
      pred <- predict(qsar, desc.pp)
    results.qsar <- data.frame(guests, pred) 
    colnames(results.qsar)[2] <- names(models)[n]
    if(n == 1) # initialize the data.frame if it's the first
      results <- results.qsar
    else
      results <- inner_join(results, results.qsar, by = "guests")
  }
  outliers <- unique(outliers)
  full.results <- results %>% 
    mutate(ensemble = rowMeans(results[ , -1]))
  predictions <- full.results %>% select(., guests, ensemble) %>%
    rename(dG = ensemble)
  
  return(list(predictions, full.results, outliers))
}

# predict.ensemble that allows outliers
predict.ensemble.all <- function(df, feat, models) {
  
  # Initializing some variables
  outliers <- c()
  results <- data.frame()
  
  for(n in 1:length(models)) {
    message(names(models)[n], " starting")
    pp.settings <- models[[n]][[1]]
    pp.results <- preprocess.desc.all(df, feat, pp.settings)
    outliers <- c(outliers, pp.results[[2]])
    desc.pp <- pp.results[[1]]
    guests <- desc.pp$guests
    desc.pp <- desc.pp[ , -1]
    
    qsar <- models[[n]][[2]]
    # if statement acceptable here b/c there is only one case to evaluate
    # would use switch otherwise
    if(str_detect(names(models)[n], "gbm")) 
      pred <- predict(qsar, desc.pp, n.trees = 500)
    else
      pred <- predict(qsar, desc.pp)
    results.qsar <- data.frame(guests, pred) 
    colnames(results.qsar)[2] <- names(models)[n]
    if(n == 1) # initialize the data.frame if it's the first
      results <- results.qsar
    else
      results <- inner_join(results, results.qsar, by = "guests")
  }
  outliers <- unique(outliers)
  full.results <- results %>% 
    mutate(ensemble = rowMeans(results[ , -1]))
  predictions <- full.results %>% select(., guests, ensemble) %>%
    rename(dG = ensemble)
  
  return(list(predictions, full.results, outliers))
}

predict.beta <- function(df) {
  outliers <- c()
  features <- readRDS("./feature.selection/beta.vars.RDS")
  name.models <- list.files("./models/beta") %>% str_remove(".RDS$")
  path.models <- paste0("./models/beta/", list.files("models/beta"))
  results <- data.frame()
  
  for(n in 1:length(name.models)) {
    pp.settings <- readRDS(path.models[n])[[1]]
    pp.results <- preprocess.desc(df, features, pp.settings)
    outliers <- c(outliers, pp.results[[2]])
    desc.pp <- pp.results[[1]]
    guests <- desc.pp$guests
    desc.pp <- desc.pp[ , -1]
    qsar <- readRDS(path.models[n])[[2]]
    if(str_detect(name.models[n], "gbm")) 
      pred <- predict(qsar, desc.pp, n.trees = 500)
    else
      pred <- predict(qsar, desc.pp)
    results.qsar <- data.frame(guests, pred) 
    colnames(results.qsar)[2] <- name.models[n]
    if(n == 1) # initialize the data.frame if it's the first
      results <- results.qsar
    else
      results <- inner_join(results, results.qsar, by = "guests")
  }
  outliers <- unique(outliers)
  full.results <- results %>% 
    mutate(ensemble = rowMeans(results[ , -1]))
  predictions <- full.results %>% select(., guests, ensemble) %>%
    rename(dG = ensemble)
  
  return(list(predictions, full.results, outliers))
}


# Alpha -------------------------------------------------------------------

#   Fill-in ----

# Taking the average of each descriptor from modeling data
# For missing values, the mean will be inputted
alpha.md <- readRDS("./model.data/alpha.md.RDS") 
colnames(alpha.md) <- str_replace(colnames(alpha.md), "-", "\\.")
alpha.data <- alpha.md %>% select(-guest:-data.source)
alpha.fill <- apply(alpha.data, 2, mean, na.rm = T) 

# alpha.fill should be 1376 descriptors long

# Saving
saveRDS(alpha.fill, "./models/alpha.fill.RDS")

#     Loading models and data ----
# Creating a list of all models
alpha.models <- lapply(list.files("./models/alpha", full.names = T), 
                       readRDS)
# Naming the list appropriately
alpha.names <- list.files("./models/alpha") %>% str_remove("\\.RDS")
names(alpha.models) <- alpha.names

# Reading the variables for alpha 
alpha.feat <- readRDS("./feature.selection/alpha.vars.RDS")

# Read ext val data 
ev.alpha <- readRDS("./ext.validation/alpha.RDS") %>%
  select(., -host)
# Save the actual dG values for later evaluation
ev.alpha.dg <- ev.alpha %>% select(., guest, DelG) %>%
  rename(obs = DelG, guests = guest)
ev.alpha <- ev.alpha %>% select(., -DelG)

#    Ensemble ----
ev.alpha.pred <- predict.ensemble(ev.alpha, alpha.feat, alpha.models)
# Outliers: iodobenzene, chlorcyclizine
alpha.outliers <- ev.alpha.pred[[3]]

ev.alpha.df <- inner_join(ev.alpha.pred[[1]], ev.alpha.dg, by = "guests") %>%
  rename(pred = dG)
# R2 = 0.69
defaultSummary(ev.alpha.df)

ggplot(ev.alpha.df, aes(x = obs, y = pred)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) + 
  theme_bw() + 
  coord_fixed() 

# Everything passes
eval.tropsha(ev.alpha.df)

# Beta --------------------------------------------------------------------

#   Fill-in ----

# Taking the average of each descriptor from modeling data
# For missing values, the mean will be inputted
beta.md <- readRDS("./model.data/beta.md.RDS") 
colnames(beta.md) <- str_replace(colnames(beta.md), "-", "\\.")
beta.data <- beta.md %>% select(-guest:-data.source)
beta.fill <- apply(beta.data, 2, mean, na.rm = T) 

# beta.fill should be 1376 descriptors long

# Saving
saveRDS(beta.fill, "./models/beta.fill.RDS")

#     Loading models and data ----

# Creating a list of all models
beta.models <- lapply(list.files("./models/beta", full.names = T), 
                       readRDS)
# Naming the list appropriately
beta.names <- list.files("./models/beta") %>% str_remove("\\.RDS")
names(beta.models) <- beta.names

# Reading the variables for beta 
beta.feat <- readRDS("./feature.selection/beta.vars.RDS")

# External validation data
ev.beta <- readRDS("./ext.validation/beta.RDS") %>%
  select(., -host)
# Save the actual dG values for later evaluation
ev.beta.dg <- ev.beta %>% select(., guest, DelG) %>%
  rename(obs = DelG, guests = guest)
ev.beta <- ev.beta %>% select(., -DelG)

# 1, 4-diiodobenzene
beta.outliers <- ev.beta.pred[[3]]
# For some reason, flufenamic acid decided to duplicate itself
ev.beta.temp <- ev.beta.pred[[1]][!duplicated(ev.beta.pred[[1]]), ]
ev.beta.df <- inner_join(ev.beta.temp, ev.beta.dg, by = "guests") %>%
  rename(pred = dG)
# R2 = 0.735
defaultSummary(ev.beta.df)
ggplot(ev.beta.df, aes(x = obs, y = pred)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) + 
  theme_bw() + 
  coord_fixed() 
# Everything passes
eval.tropsha(ev.beta.df)

combined.df <- rbind(
  ev.alpha.df %>% mutate(host = "alpha"), 
  ev.beta.df %>% mutate(host = "beta")
)

ggplot(combined.df, aes(x = obs, y = pred, color = host)) + 
  geom_point() + 
  theme_bw() + 
  coord_fixed() + 
  geom_abline(intercept = 0, slope = 1) + 
  labs(x = "Observed dG, kJ/mol", y = "Predicted dG, kJ/mol", 
       title = "Ensemble QSAR results", 
       color = "CD type")
ggsave("./results/2018-08 ensemble.png", dpi = 450)
saveRDS(combined.df, "./results/ensemble.RDS")


# Temp --------------------------------------------------------------------

cisplatin <- read.csv("~/SREP LAB/qsar/cisplatin_DB00515.csv", 
                      na.strings = NA) %>%
  rename(guest = Name)
# replacing infinite values w/ NA
is.na(cisplatin) <- sapply(cisplatin, is.infinite)
colnames(cisplatin) <- str_replace(colnames(cisplatin), 
                                   "-", "\\.")
cisplatin <- desc.fill(cisplatin)

# cisplatin <- rbind(cisplatin, cisplatin, cisplatin) %>%
#   mutate(Name = as.character(Name))
# doubling the rows b/c of a quirk of an applicability domain function
# will debug later
cisplatin <- rbind(cisplatin, cisplatin)

# ends up as 64 observations...I don't know why
# using 3 rows has 729 observations
# row ^ 6 = # of observations? Why?
cisplatin.alpha <- predict.ensemble(cisplatin, alpha.feat, alpha.models)
# -9.256442 kJ/mol

# beta.models1 <- beta.models
# beta.models1[[1]] <- NULL # removing cubist
temp <- ev.beta[ , -2]
colnames(temp) <- colnames(cisplatin)
temp <- rbind(temp, cisplatin)
cisplatin.beta <- predict.ensemble.all(temp, beta.feat, beta.models)
cisplatin.beta <- predict.beta(cisplatin)
# -12.92012 kJ/mol