# Essentially just a cleaner version of 04.2
library(tidyverse)
library(stringr)
library(data.table)
library(ggplot2)
library(Matrix)
library(glmnet)
library(caret)
library(e1071)
# Change file path as needed
ri.padel <- read_csv("~/SREP LAB/Rekharsky and Inoue/Cactus/03-PaDEL-Descriptor and Rekharsky and Inoue.csv")

# Preprocessing -----------------------------------------------------------

# 1. Variance Cleaning
#     Checks for predictors that don't vary much; removes them
#     Note: rp stands for ri.padel and rpt stands for ri.padel.trans 
zero.pred <- nearZeroVar(ri.padel)
rp.no.zero <- ri.padel[, -zero.pred] # View(rp.no.zero)
rpt <- preProcess(rp.no.zero, na.remove = T, 
                  method = c("center", "scale"
                             # , "BoxCox" throws an error: Na/NaN/Inf in foreign function call
                  )) %>%
  predict(., ri.padel)
