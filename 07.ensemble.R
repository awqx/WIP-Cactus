source("07.model.functions.R")

# Libraries and Packages -------------------------------------------------

library(caret)
library(Cubist)
library(e1071)
library(glmnet)
library(kernlab)
library(pls)
library(randomForest)
library(tidyverse)

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
    if(str_detect(name.models[n], "glm")) {
      desc.pp <- data.matrix(desc.pp)
      pred <- predict.glmnet(qsar, desc.pp, 
                             s = tail(qsar$lambda, n = 1))
    } else if (str_detect(name.models[n], "pls"))
      pred <- predict(qsar, desc.pp, ncomp = 8)
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
    if(str_detect(name.models[n], "glm")) {
      desc.pp <- data.matrix(desc.pp)
      pred <- predict.glmnet(qsar, desc.pp, 
                             s = tail(qsar$lambda, n = 1))
    } else if (str_detect(name.models[n], "pls"))
      pred <- predict(qsar, desc.pp, ncomp = 25)
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

# Read the file
ev.alpha <- readRDS("./ext.validation/alpha.RDS") %>%
  select(., -host)
# Save the actual dG values for later evaluation
ev.alpha.dg <- ev.alpha %>% select(., guest, DelG) %>%
  rename(obs = DelG, guests = guest)
ev.alpha <- ev.alpha %>% select(., -DelG)

ev.alpha.pred <- predict.alpha(ev.alpha)
# No outliers
alpha.outliers <- ev.alpha.pred[[3]]
  
ev.alpha.df <- inner_join(ev.alpha.pred, ev.alpha.dg, by = "guests") %>%
  rename(obs = ev.alpha.dg, pred = dG)
defaultSummary(ev.alpha.df)
ggplot(ev.alpha.df, aes(x = obs, y = pred)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) + 
  theme_bw() + 
  coord_fixed() 
# Everything passes
# R2 = 0.66
eval.tropsha(ev.alpha.df)

# Beta --------------------------------------------------------------------

ev.beta <- readRDS("./ext.validation/beta.RDS") %>%
  select(., -host)
# Save the actual dG values for later evaluation
ev.beta.dg <- ev.beta %>% select(., guest, DelG) %>%
  rename(obs = DelG, guests = guest)
ev.beta <- ev.beta %>% select(., -DelG)

ev.beta.pred <- predict.beta(ev.beta)
# No outliers
beta.outliers <- ev.beta.pred[[3]]

ev.beta.df <- inner_join(ev.beta.pred[[1]], ev.beta.dg, by = "guests") %>%
  rename(pred = dG)
defaultSummary(ev.beta.df)
ggplot(ev.beta.df, aes(x = obs, y = pred)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) + 
  theme_bw() + 
  coord_fixed() 
# Everything passes
# R2 = 0.62
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

saveRDS(combined.df, "./results/ensemble.RDS")
 