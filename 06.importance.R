# Important note: Only workin with Beta-CD models (for now) due to accuracy
# of those models and time constraint

# Libraries and Packages --------------------------------------------------

library(caret)
library(Cubist)
library(data.table)
library(glmnet)
library(randomForest)
library(svm)
library(tidyverse)

# Loading Models ----------------------------------------------------------

df.raw <- readRDS("./padel.pp.new.RDS")
df <- df.raw %>% select(-guest,-host,-data.source) 
beta <- df %>% filter(beta > 0)

# setwd("~/SREP LAB/qsar")
cube <- readRDS("./models/cubist/cube.beta.RDS")
# glm <- readRDS("./models/glmnet/glm.all.sprse.RDS")
pls <- readRDS("./models/pls/pls.beta.RDS")
rf <- readRDS("./models/rforest/rf.beta.RDS")
svm <- readRDS("./models/svm/polysvm.beta.RDS")

# Functions ---------------------------------------------------------------

# Rescale from 0 to 1
range01 <- function(x) {(x-min(x))/(max(x)-min(x))}

# Rescale from any range
new.range <- function(x, newMax, newMin) { 
  (x - min(x))/(max(x)-min(x)) * (newMax - newMin) + newMin 
}

# VIP for PLS
# Function sourced from http://mevik.net/work/software/pls.html
VIP <- function(object) {
  if (object$method != "oscorespls")
    stop("Only implemented for orthogonal scores algorithm.  Refit with 'method = \"oscorespls\"'")
  if (nrow(object$Yloadings) > 1)
    stop("Only implemented for single-response models")
  
  SS <- c(object$Yloadings)^2 * colSums(object$scores^2)
  Wnorm2 <- colSums(object$loading.weights^2)
  SSW <- sweep(object$loading.weights^2, 2, SS / Wnorm2, "*")
  sqrt(nrow(SSW) * apply(SSW, 1, cumsum) / cumsum(SS))
}

sort.vip <- function(vip.df) {
  df <- vip.df[1, ] %>% gather()
  colnames(df)[2] <- row.names(vip.df)[1]
  for(r in 2:nrow(vip.df)) {
    long <- vip.df[r, ] %>% gather()
    colnames(long)[2] <- row.names(vip.df)[r]
    df <- inner_join(df, long, by = "key")
    }
  return(df)
}

# "Easy" Calculations -----------------------------------------------------

cube.imp <- varImp(cube) %>%
  mutate(desc = rownames(.)) %>%
  mutate(model = "cubist")
cube.imp <- cube.imp[order(-cube.imp$Overall), ]
rf.imp <- varImp(rf) %>%
  mutate(desc = rownames(.)) %>%
  mutate(model = "rforest")
rf.imp <- rf.imp[order(-rf.imp$Overall), ]

# Two-step Calculations ---------------------------------------------------

# Glmnet

# Needs to be trained with caret::train

glm <- train(beta[ , -1], beta[ , 1], 
             method = "glmnet", metric = "RMSE")
glm.imp <- varImp(glm)[[1]] %>% 
  mutate(desc = rownames(.)) %>%
  mutate(model = "glm")
glm.imp <- glm.imp[order(-glm.imp$Overall), ]

# SVM
svm <- train(beta[ , -1], beta[ , 1], 
             method = "svmPoly", metric = "RMSE")
svm.imp <- varImp(svm)[[1]] %>% 
  mutate(desc = rownames(.)) %>%
  mutate(model = "svm") 
svm.imp$Overall[is.na(svm.imp$Overall)] <- 0
svm.imp <- svm.imp[order(-svm.imp$Overall), ]

# PLS
# See note in "functions"
pls.vip <- VIP(pls) %>% data.frame() 
pls.imp <- sort.vip(pls.vip)
pls.imp$Overall <- rowMeans(pls.imp[ , -1])
pls.imp <- pls.imp[order(-pls.imp$Overall), ]

# Plot --------------------------------------------------------------------

imp1 <- row.names(cube.imp[1:50, ])
imp2 <- row.names(rf.imp[1:50, ])
imp3 <- row.names(glm.imp[1:50, ])
imp4 <- row.names(svm.imp[1:50, ])
imp5 <- row.names(pls.imp[1:50, ])

imps <- c(imp1, imp2, imp3, imp4, imp5) 
imps <- imps[duplicated(imps)]

cube.imp$Overall <- range01(cube.imp$Overall)
rf.imp$Overall <- range01(rf.imp$Overall)
glm.imp$Overall <- range01(glm.imp$Overall)
svm.imp$Overall <- range01(svm.imp$Overall)
pls.imp$Overall <- range01(pls.imp$Overall)

pls.temp <- pls.imp[ , -2:-12] %>%
  mutate(model = "pls") %>%
  dplyr::rename(., desc = key)
comp <- rbind(
  cube.imp, rf.imp, glm.imp, svm.imp, pls.temp
)

comp <- comp %>% filter(desc %in% imps) %>%
  dplyr::rename(., importance = Overall)

ggplot(comp, aes(x = model, y = desc, fill = importance)) + 
  geom_tile() + 
  scale_fill_gradient2(low = "palegreen1", high = "#AEA5DE", mid = "white", midpoint = 0.5) + 
  theme_bw() + 
  labs(title = "Descriptor Importance for Each Model on Beta-CD", x = "Model", 
       y = "Descriptor Variable", fill = "Importance")
