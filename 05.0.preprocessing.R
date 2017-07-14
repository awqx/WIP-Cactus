library(caret)
library(data.table)
library(Matrix)
library(stringr)
library(tidyverse)

# Pre-processing and cleaning ---------------------------------------------

setwd("~/SREP LAB/qsar")
padel <- readRDS("./molecules/descriptors/04.all.padel.RDS")
# Removing predictors with near zero variance - 1302 predictors
zero.pred <- nearZeroVar(padel)
zero.pred.names <- colnames(padel)[zero.pred]
padel.no.zero <- padel[ , -zero.pred]

# Binning alpha, beta, and gamma
padel.cd <- padel.no.zero %>% 
  mutate(alpha = ifelse(str_detect(host, "alpha"), 1, 0)) %>%
  mutate(beta = ifelse(str_detect(host, "beta"), 1, 0)) %>%
  mutate(gamma = ifelse(str_detect(host, "gamma"), 1, 0))

# Separating nnon-predictors from the data
padel.split1 <- padel.cd %>% dplyr::select(., guest:data.source)
padel.split2 <- padel.cd %>% dplyr::select(., -guest:-data.source)

# padel.pp = padel.preprocess, just shortened
padel.pp <- preProcess(padel.split2, na.remove = T, 
                  method = c("knnImpute", "center", "scale")) %>%
  predict(., padel.split2)
# padel.pp1 <- preProcess(padel.pp, method = "BoxCox") %>%
#   predict(., padel.pp) # Not sure if any effect from BoxCox

# More nearZeroVar analysis after preprocessing, to be safe
zero.pred2 <- nearZeroVar(padel.pp) # 25 predictors
zero.pred2.names <- colnames(padel.pp)[zero.pred2]
padel.pp <- padel.pp[ , -zero.pred2]

# Removing highly correlated predictors
too.high <- findCorrelation(cor(padel.pp), 0.95) # 0.95 mostly arbitrary
corr <- names(padel.pp)[too.high]
padel.pp <- padel.pp[ , -too.high]
padel.pp <- cbind(padel.split1, padel.pp)
colnames(padel.pp) <- str_replace(colnames(padel.pp), "-", ".")

# Creating External Validation Set ----------------------------------------

set.seed(4)
ext.val.ind <- sample(x = 1:nrow(padel.pp), 
                      size = round(0.15 * nrow(padel.pp)))
ext.val <- padel.pp[ext.val.ind, ]

# Data Organization and Saving --------------------------------------------

sprse.padel <- sparse.model.matrix(~., padel.pp)
mat.padel <- as.matrix(padel.pp)

saveRDS(padel.pp, "./padel.pp.RDS")
saveRDS(sprse.padel, "./sprse.padel.RDS")
saveRDS(mat.padel, "./mat.padel.RDS")

saveRDS(ext.val, "./external validation set.RDS")