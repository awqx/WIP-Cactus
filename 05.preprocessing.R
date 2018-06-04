# Libraries and Packages --------------------------------------------------

library(caret)
library(data.table)
library(Matrix)
library(stringr)
library(tidyverse)

# Pre-processing and cleaning ---------------------------------------------

padel <- readRDS("./descriptors/all.padel.RDS") 

# 1. Removing predictors with near zero variance - 1302 predictors
zero.pred <- nearZeroVar(padel)
zero.pred.names <- colnames(padel)[zero.pred]
padel.no.zero <- padel[ , -zero.pred]

# 2. Binning alpha, beta, and gamma
padel.cd <- padel.no.zero %>% 
  mutate(alpha = ifelse(str_detect(host, "alpha"), 1, 0)) %>%
  mutate(beta = ifelse(str_detect(host, "beta"), 1, 0)) %>%
  mutate(gamma = ifelse(str_detect(host, "gamma"), 1, 0))

# 3. Separating non-predictors from data
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
padel.pp <- predict(pp.settings, padel.temp)

# padel.pp1 <- preProcess(padel.pp, method = "BoxCox") %>%
#   predict(., padel.pp) # Not sure if any effect from BoxCox

# 4. More nearZeroVar analysis after preprocessing, to be safe
zero.pred2 <- nearZeroVar(padel.pp)
zero.pred2.names <- colnames(padel.pp)[zero.pred2]  # Originally gmin and gamma, 
                                                    # but I wanted to keep gamma
zero.pred2 <- zero.pred2[1]
padel.pp <- padel.pp[ , -zero.pred2]

# 5. Removing highly correlated predictors
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

dir.create("./pre-process")
saveRDS(pp.settings, "./pre-process/pp.settings.RDS")

dir.create("./data")
saveRDS(padel.pp, "./data/padel.pp.RDS")
saveRDS(sprse.padel, "./data/sprse.padel.RDS")
saveRDS(mat.padel, "./data/mat.padel.RDS")

saveRDS(ext.val, "./data/ext.val.RDS")

# Suzuki Only -------------------------------------------------------------

suz <- readRDS("./descriptors/suz.all.padel.RDS")
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

saveRDS(suz.pp, "./data/suz.pp.RDS")
saveRDS(ext.val, "./data/suz.extval.RDS")
