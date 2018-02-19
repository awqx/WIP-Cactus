# Libraries and Packages -------------------------------------------------

library(caret)
library(Cubist)
library(e1071)
library(glmnet)
library(kernlab)
library(Matrix)
library(pls)
library(stats)
library(tidyverse)

# Download Models ---------------------------------------------------------

# setwd("~/SREP LAB/qsar/models/")
cube.a <- readRDS("./models/cubist/cubist.alpha.RDS")
cube.b <- readRDS("./models/cubist/cubist.beta.RDS")
cube.c <- readRDS("./models/cubist/cubist.gamma.RDS")

glm.a <- readRDS("./models/glmnet/glm.alpha.RDS")
glm.b <- readRDS("./models/glmnet/glm.beta.RDS")
glm.c <- readRDS("./models/glmnet/glm.gamma.RDS")

rf.a <- readRDS("./models/rforest/rf.alpha.RDS")
rf.b <- readRDS("./models/rforest/rf.beta.RDS")
rf.c <- readRDS("./models/rforest/rf.gamma.RDS")

pls.a <- readRDS("./models/pls/pls.alpha.RDS")
pls.b <- readRDS("./models/pls/pls.beta.RDS")
pls.c <- readRDS("./models/pls/pls.gamma.RDS")

svm.a <- readRDS("./models/svm/polysvm.alpha.RDS")
svm.b <- readRDS("./models/svm/polysvm.beta.RDS")
svm.c <- readRDS("./models/svm/polysvm.gamma.RDS")

# Averaging Test Results --------------------------------------------------

cubist.tst <- readRDS("./models/cubist/cubist.results.RDS") %>% 
  select(pred, obs, cd.type) %>%
  rename(pred.cubist = pred)
glmnet.tst <- readRDS("./models/glmnet/glmnet.tst.results.RDS") %>% 
  select(pred, obs, cd.type) %>%
  rename(pred.glmnet = pred)
pls.tst <- readRDS("./models/pls/pls.results.RDS") %>% 
  select(pred, obs, cd.type) %>%
  rename(pred.pls = pred)
rforest.tst <- readRDS("./models/rforest/rf.results.RDS") %>% 
  select(pred, obs, cd.type) %>%
  rename(pred.rf = pred)
svm.tst <- readRDS("./models/svm/polysvm.tst.results.RDS") %>% 
  select(pred, obs, cd.type) %>%
  rename(pred.svm = pred)

# Currently, averaging the results has some errors with joining observed values
# and repeating some observations that shouldn't be duplicates

# tst.lst <- list(cubist.tst, glmnet.tst, pls.tst, rforest.tst, svm.tst)
# tst.results <- Reduce(full_join, tst.lst) %>% 
#   mutate(avg = rowMeans(subset(., select = c(pred.cubist, pred.glmnet, pred.pls, pred.rf, pred.svm))))
# avg.results <- tst.results %>% select(obs, cd.type, avg) %>%
#   rename(pred = avg) %>%
#   as.data.frame()
# defaultSummary(avg.results) # 0.744
# 
# ggplot(avg.results, aes(x = obs, y = pred, color = cd.type)) + 
#   geom_point()

# External Validation -----------------------------------------------------

ext.val <- readRDS("./external validation set new.RDS") %>%
  select(-guest:-host) %>%
  select(-data.source) 

ext.val.a <- ext.val %>% filter(alpha > 0) 
ext.val.b <- ext.val %>% filter(beta > 0)
ext.val.c <- ext.val %>% filter(gamma > 0) 

ev.a.glm <-  predict.glmnet(glm.a, as.matrix(ext.val.a[ , -1]), 
                        s = tail(glm.a$lambda, n = 1)) %>%
  cbind(as.matrix(ext.val.a[ , 1])) %>% data.frame() %>%
  dplyr::rename(pred.glm = X1, obs = V2) %>%
  mutate(model = "GLMNet")

ev.a.rf <-  predict(rf.a, ext.val.a[ , -1]) %>%
  cbind(ext.val.a[ , 1]) %>% data.frame() %>%
  dplyr::rename(., pred.rf = `.`, obs = V2) %>%
  mutate(model = "Random Forest")

ev.a <- cbind(ev.a.glm, ev.a.rf)

# Graphing ----------------------------------------------------------------

ggplot(ev.a.glm, aes(x = obs, y = pred.glm)) +
  geom_point() + 
  geom_abline(slope = 1, intercept = 0)
