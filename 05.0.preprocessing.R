# Libraries and Packages --------------------------------------------------

library(caret)
library(data.table)
library(Matrix)
library(stringr)
library(tidyverse)

# Pre-processing and cleaning ---------------------------------------------

setwd("~/SREP LAB/qsar")
padel <- readRDS("./molecules/descriptors/all.padel.RDS") # prev: 04.all.padel
# Removing predictors with near zero variance - 1302 predictors
zero.pred <- nearZeroVar(padel)
zero.pred.names <- colnames(padel)[zero.pred]
padel.no.zero <- padel[ , -zero.pred]

# Binning alpha, beta, and gamma
padel.cd <- padel.no.zero %>% 
  mutate(alpha = ifelse(str_detect(host, "alpha"), 1, 0)) %>%
  mutate(beta = ifelse(str_detect(host, "beta"), 1, 0)) %>%
  mutate(gamma = ifelse(str_detect(host, "gamma"), 1, 0))

# Separating non-predictors from data
padel.split1 <- padel.cd %>% dplyr::select(., guest:data.source) # %>%
  # filter(data.source == "rekharsky.inoue")
padel.split2 <- padel.cd %>% # filter(data.source == "rekharsky.inoue") %>% 
  dplyr::select(., -guest:-data.source)

padel.split3 <- padel.cd %>% dplyr::select(., guest:data.source) %>%
  filter(data.source == "suzuki")
padel.split4 <- padel.cd %>% filter(data.source == "suzuki") %>% 
  dplyr::select(., -guest:-data.source)
  

# padel.pp = padel.preprocess, just shortened
padel.temp <- do.call(data.frame,lapply(padel.split2, 
                                        function(x) replace(x, is.infinite(x),NA)))
pp.settings <- preProcess(padel.temp, na.remove = T, 
                          method = c("knnImpute", "center", "scale"))
saveRDS(pp.settings, "preprocess.settings.RDS")
padel.pp <- predict(pp.settings, padel.temp)

# padel.pp1 <- preProcess(padel.pp, method = "BoxCox") %>%
#   predict(., padel.pp) # Not sure if any effect from BoxCox

# More nearZeroVar analysis after preprocessing, to be safe
zero.pred2 <- nearZeroVar(padel.pp)
zero.pred2.names <- colnames(padel.pp)[zero.pred2]  # Originally gmin and gamma, 
                                                    # but I wanted to keep gamma
zero.pred2 <- zero.pred2[1]
padel.pp <- padel.pp[ , -zero.pred2]

# Removing highly correlated predictors
too.high <- findCorrelation(cor(padel.pp), 0.95) # 0.95 mostly arbitrary
corr <- names(padel.pp)[too.high]
padel.pp <- padel.pp[ , -too.high]

padel.pp <- cbind(padel.split1, padel.pp)

colnames(padel.pp) <- str_replace(colnames(padel.pp), "-", ".")

#     Creating External Validation Set ----------------------------------------

set.seed(4)
ext.val.ind <- sample(x = 1:nrow(padel.pp), 
                      size = round(0.15 * nrow(padel.pp)))
ext.val <- padel.pp[ext.val.ind, ]
padel.pp.all <- padel.pp
padel.pp <- padel.pp[-ext.val.ind, ]
#     Data Organization and Saving --------------------------------------------

sprse.padel <- sparse.model.matrix(~., padel.pp)
mat.padel <- as.matrix(padel.pp)

saveRDS(padel.pp, "./padel.pp.new.RDS")
saveRDS(sprse.padel, "./sprse.padel.RDS")
saveRDS(mat.padel, "./mat.padel.RDS")

saveRDS(ext.val, "./external validation set new.RDS")

# Processing GAFS Data ----------------------------------------------------

ga.padel <- readRDS("./feature.select/GAFS df.RDS")
# Removing predictors with near zero variance - 1302 predictors
# ga.zero.pred <- nearZeroVar(ga.padel) # Good news: no zero variances

# Putting "gamma" back into the data
ga.cd <- cbind(ga.padel, padel.pp[ , "gamma"]) %>%
  rename(gamma = `padel.pp[, \"gamma\"]`)

# Removing highly correlated predictors
ga.too.high <- findCorrelation(cor(ga.cd[ , -1:-4]), 0.95) # No high corr

ga.pp <- ga.cd

# Modifying external validation set to match genetic alg
ga.extval <- readRDS("./external validation set.RDS")
ga.extval <- ga.extval[ , colnames(ga.extval) %in% colnames(ga.pp)]

#     Data Organization and Saving ----------------------------------------

sprse.ga <- sparse.model.matrix(~., ga.pp)
mat.ga <- as.matrix(ga.pp)

saveRDS(ga.pp, "./feature.select/ga.pp.RDS")
saveRDS(sprse.ga, "./feature.select/sprse.ga.RDS")
saveRDS(mat.ga, "./feature.select/mat.ga.RDS")

saveRDS(ga.extval, "./feature.select/gafs.extval.RDS")

# Processing 2D Only Data -------------------------------------------------

twod <- readRDS("./molecules/descriptors/04.all.2d.RDS")
# Removing predictors with near zero variance - 1302 predictors
zero.pred <- nearZeroVar(twod)
zero.pred.names <- colnames(twod)[zero.pred]
twod.no.zero <- twod[ , -zero.pred]

# Binning alpha, beta, and gamma
twod.cd <- twod.no.zero %>% 
  mutate(alpha = ifelse(str_detect(host, "alpha"), 1, 0)) %>%
  mutate(beta = ifelse(str_detect(host, "beta"), 1, 0)) %>%
  mutate(gamma = ifelse(str_detect(host, "gamma"), 1, 0))

# Separating non-predictors from data
twod.split1 <- twod.cd %>% dplyr::select(., guest:data.source)
twod.split2 <- twod.cd %>% dplyr::select(., -guest:-data.source)

# twod.pp = twod.preprocess, just shortened
twod.pp <- preProcess(twod.split2, na.remove = T, 
                       method = c("medianImpute", "center", "scale")) %>%
  predict(., twod.split2)
twod.pp <- lapply(twod.pp, as.numeric) %>% data.frame()
twod.pp <- preProcess(twod.pp, na.remove = T, 
                      method = c("medianImpute")) %>%
  predict(., twod.pp)

# Removing highly correlated predictors
too.high <- findCorrelation(cor(twod.pp), 0.95) # 0.95 mostly arbitrary
corr <- names(twod.pp)[too.high]
twod.pp <- twod.pp[ , -too.high]
twod.pp <- cbind(twod.split1, twod.pp)

colnames(twod.pp) <- str_replace(colnames(twod.pp), "-", ".")

#     External Validation ----
set.seed(4)
ext.val.ind <- sample(x = 1:nrow(twod.pp), 
                      size = round(0.15 * nrow(twod.pp)))
ext.val <- twod.pp[ext.val.ind, ]

#     Data Organization and Saving --------------------------------------------

saveRDS(twod.pp, "./2d.pp.RDS")
saveRDS(ext.val, "./2d.extval.RDS")

# Suzuki Only -------------------------------------------------------------

suz <- readRDS("./molecules/descriptors/suz.all.padel.RDS")
suz.nz <- suz[ , -zero.pred]
suz.cd <- suz.nz  %>% 
  mutate(alpha = ifelse(str_detect(host, "alpha"), 1, 0)) %>%
  mutate(beta = ifelse(str_detect(host, "beta"), 1, 0)) %>%
  mutate(gamma = ifelse(str_detect(host, "gamma"), 1, 0))

suz.split <- suz.cd %>% dplyr::select(., -guest:-data.source)
suz.temp <- do.call(data.frame,
                    lapply(suz.split, function(x) replace(x, is.infinite(x),NA)))
suz.pp <- predict(pp.settings, suz.temp)
suz.pp <- suz.pp[ , -zero.pred2]
suz.pp <- suz.pp[ , -too.high]
suz.pp <- cbind(suz.nz[ , 1:4], suz.pp)
colnames(suz.pp) <- str_replace(colnames(suz.pp), "-", ".")

#     External Validation -------------------------------------------------

set.seed(4)
ext.val.ind <- sample(x = 1:nrow(suz.pp), 
                      size = round(0.15 * nrow(suz.pp)))
ext.val <- suz.pp[ext.val.ind, ]
suz.pp.all <- suz.pp
suz.pp <- suz.pp[-ext.val.ind, ]

#     Saving --------------------------------------------------------------

saveRDS(suz.pp, "./suz.pp.RDS")
saveRDS(ext.val, "./suz.extval.RDS")
