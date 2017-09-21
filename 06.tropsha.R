# Models evaluated by standards from Tropsha 2010, DOI: 10.1002/minf201000061
# "Best Practices for QSAR Model Development, Validation, and Exploitation"


# Set working directory as needed
# setwd("~/SREP LAB/qsar")

# Libraries and Packages --------------------------------------------------

library(tidyverse)
library(caret)
library(stats)

# Functions ---------------------------------------------------------------

# Requirements: results is a datframe w/ calculated values labeled as pred 
# and observed values listed as obs
find.k <- function(results) {
  top <- sum(results$pred * results$obs)
  bottom <- sum(results$obs * results$obs)
  return (top/bottom)
}

find.q2.f1 <- function(results) {
  top <- (results$obs - results$pred)^2 %>% sum()
  bottom <- (results$obs - mean(results$pred))^2 %>% sum()
  return (1 - (top/bottom))
}

find.r2r20 <- function(results) {
  ybarr0 <- find.k(results) * results$obs
  top <- (results$obs - ybarr0)^2 %>% sum()
  bottom <- (results$obs - mean(results$obs))^2 %>% sum()
  r20 <- 1 - (top/bottom)
  r2 <- defaultSummary(results)[2]
  return((r2-r20)/r2)
}

# Importing Models' Results -----------------------------------------------

cube <- readRDS("./models/cubist/compiled.results.RDS")
svm  <- readRDS("./models/svm/polysvm.tst.results.RDS")

# Evaluating Standards ----------------------------------------------------

# 0.85 <= k <= 1.15
find.k(cube) # 0.985
find.k(svm) # 0.954

# Q2-F1 > 0.5
find.q2.f1(cube) # 0.659
find.q2.f1(svm) # 0.601

# (R2-R02)/R2 < 0.1
find.r2r20(cube) # -0.441
find.r2r20(svm) # -0.601