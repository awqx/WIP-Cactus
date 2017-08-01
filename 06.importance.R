# Libraries and Packages --------------------------------------------------

library(caret)
library(Cubist)
library(glmnet)
library(randomForest)
library(svm)
library(tidyverse)

# Loading Models ----------------------------------------------------------

df.raw <- readRDS("./padel.pp.new.RDS")
df <- df.raw %>% select(-guest,-host,-data.source) 

# setwd("~/SREP LAB/qsar")
cube <- readRDS("./models/cubist/cube.all.RDS")
# glm <- readRDS("./models/glmnet/glm.all.sprse.RDS")
rf <- readRDS("./models/rforest/rforest.all.RDS")
svm <- readRDS("./models/svm/polysvm.all.RDS")

# Functions ---------------------------------------------------------------

# Rescale from 0 to 1
range01 <- function(x) {(x-min(x))/(max(x)-min(x))}

# Rescale from any range
new.range <- function(x, newMax, newMin) { 
  (x - min(x))/(max(x)-min(x)) * (newMax - newMin) + newMin 
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

glm <- train(df[ , -1], df[ , 1], 
             method = "glmnet", metric = "RMSE")
glm.imp <- varImp(glm)[[1]] %>% 
  mutate(desc = rownames(.)) %>%
  mutate(model = "glm")
glm.imp <- glm.imp[order(-glm.imp$Overall), ]

# SVM
svm <- train(df[ , -1], df[ , 1], 
             method = "svmPoly", metric = "RMSE")
svm.imp <- varImp(svm)[[1]] %>% 
  mutate(desc = rownames(.)) %>%
  mutate(model = "svm") %>%
  mutate(Overall = range01(Overall))
svm.imp <- svm.imp[order(-svm.imp$Overall), ]

# Plot --------------------------------------------------------------------

imp1 <- row.names(cube.imp[1:50, ])
imp2 <- row.names(rf.imp[1:50, ])
imp3 <- row.names(glm.imp[1:50, ])
imp4 <- row.names(svm.imp[1:50, ])

imps <- c(imp1, imp2, imp3, imp4) 
imps <- imps[duplicated(imps)]

cube.imp$Overall <- range01(cube.imp$Overall)
rf.imp$Overall <- range01(rf.imp$Overall)
glm.imp$Overall <- range01(glm.imp$Overall)

comp <- rbind(
  cube.imp, rf.imp, glm.imp, svm.imp
)

comp <- comp %>% filter(desc %in% imps) %>%
  dplyr::rename(., importance = Overall)

ggplot(comp, aes(x = model, y = desc, fill = importance)) + 
  geom_tile() + 
  scale_fill_gradient2(low = "palegreen1", high = "#AEA5DE", mid = "white", midpoint = 0.5) + 
  theme_bw() + 
  labs(title = "Descriptor Importance for Each Model", x = "Model", 
       y = "Descriptor Variable", fill = "Importance")
