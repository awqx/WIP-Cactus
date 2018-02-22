# Libraries and Packages -------------------------------------------------

library(caret)
library(Cubist)
library(data.table)
library(e1071)
library(glmnet)
library(kernlab)
library(Matrix)
library(pls)
library(randomForest)
library(stats)
library(tidyverse)

# Download Models ---------------------------------------------------------

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

# Averaging Results -------------------------------------------------------

# Training results
cubist.trn <- readRDS("./models/cubist/cubist.trn.results.RDS") %>% 
  select(pred, obs, cd.type) # %>%
  # rename(pred.cubist = pred)
glmnet.trn <- readRDS("./models/glmnet/glmnet.trn.results.RDS") %>% 
  select(pred, obs, cd.type) # %>%
  # rename(pred.glmnet = pred)
pls.trn <- readRDS("./models/pls/pls.trn.results.RDS") %>% 
  select(pred, obs, cd.type) # %>%
  # rename(pred.pls = pred)
rforest.trn <- readRDS("./models/rforest/rf.trn.results.RDS") %>% 
  select(pred, obs, cd.type) # %>%
  # rename(pred.rf = pred)
svm.trn <- readRDS("./models/svm/polysvm.trn.results.RDS") %>% 
  select(pred, obs, cd.type) # %>%
  # rename(pred.svm = pred)

trn.lst <- list(cubist.trn, glmnet.trn, pls.trn, 
                rforest.trn, svm.trn)
lapply(trn.lst, defaultSummary)

# Test results
cubist.tst <- readRDS("./models/cubist/cubist.results.RDS") %>% 
  select(pred, obs, cd.type) # %>%
  # rename(pred.cubist = pred)
glmnet.tst <- readRDS("./models/glmnet/glmnet.tst.results.RDS") %>% 
  select(pred, obs, cd.type) # %>%
  # rename(pred.glmnet = pred)
pls.tst <- readRDS("./models/pls/pls.results.RDS") %>% 
  select(pred, obs, cd.type) # %>%
  # rename(pred.pls = pred)
rforest.tst <- readRDS("./models/rforest/rf.results.RDS") %>% 
  select(pred, obs, cd.type) # %>%
  # rename(pred.rf = pred)
svm.tst <- readRDS("./models/svm/polysvm.tst.results.RDS") %>% 
  select(pred, obs, cd.type) # %>%
  # rename(pred.svm = pred)

tst.lst <- list(cubist.tst, glmnet.tst, pls.tst, 
                rforest.tst, svm.tst)
lapply(tst.lst, defaultSummary)

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

# ext.val <- readRDS("~/SREP LAB/qsar/external validation set new.RDS") %>%
#   select(-guest:-host) %>%
#   select(-data.source)

ext.val.a <- ext.val %>% filter(alpha > 0) 
ext.val.b <- ext.val %>% filter(beta > 0)
ext.val.c <- ext.val %>% filter(gamma > 0) 

# Alpha-CD

ctrl <- cubistControl(
  seed = 10, 
  sample = 75
)

ev.a.cube <- predict(cube.a, ext.val.a[ , -1]) %>%
  cbind(ext.val.a[ , 1]) %>%
  data.frame() %>%
  dplyr::rename(., pred = `.`, obs = V2) %>%
  mutate(model = "Cubist")

ev.a.glm <-  predict.glmnet(glm.a, as.matrix(ext.val.a[ , -1]), 
                        s = tail(glm.a$lambda, n = 1)) %>%
  cbind(as.matrix(ext.val.a[ , 1])) %>% data.frame() %>%
  dplyr::rename(pred = X1, obs = V2) %>%
  mutate(model = "GLMNet")

ev.a.pls <- predict(pls.a, ncomp = 4, newdata = ext.val.a[ , -1]) %>%
  cbind(ext.val.a[ , 1]) %>%
  data.frame() %>%
  dplyr::rename(., pred = `.`, obs = V2) %>%
  mutate(model = "PLS")

ev.a.rf <-predict(rf.a, ext.val.a[ , -1]) %>%
  cbind(ext.val.a[ , 1]) %>% data.frame() %>%
  dplyr::rename(., pred = `.`, obs = V2) %>%
  mutate(model = "Random Forest")

ev.a.svm <- predict(svm.a, ext.val.a[ , -1]) %>%
  cbind(ext.val.a[ , 1]) %>% data.frame() %>%
  dplyr::rename(., pred = `.`, obs = V2) %>%
  mutate(model = "SVM")

ev.a <- rbind(ev.a.cube, ev.a.glm, ev.a.pls, ev.a.rf, ev.a.svm)
ev.a.avg <- ev.a %>% data.table(key = "obs")
ev.a.avg <- ev.a.avg[ , list(pred = mean(pred)), by = "obs"]

defaultSummary(ev.a.avg) # 0.121
defaultSummary(ev.a.avg[ -1, ]) # 0.531

# Beta-CD

ev.b.cube <- predict(cube.b, ext.val.b[ , -1]) %>%
  cbind(ext.val.b[ , 1]) %>%
  data.frame() %>%
  dplyr::rename(., pred = `.`, obs = V2) %>%
  mutate(model = "Cubist")

ev.b.glm <-  predict.glmnet(glm.b, as.matrix(ext.val.b[ , -1]), 
                            s = tail(glm.b$lambda, n = 1)) %>%
  cbind(as.matrix(ext.val.b[ , 1])) %>% data.frame() %>%
  dplyr::rename(pred = X1, obs = V2) %>%
  mutate(model = "GLMNet")

ev.b.pls <- predict(pls.b, ncomp = 4, newdata = ext.val.b[ , -1]) %>%
  cbind(ext.val.b[ , 1]) %>%
  data.frame() %>%
  dplyr::rename(., pred = `.`, obs = V2) %>%
  mutate(model = "PLS")

ev.b.rf <-predict(rf.b, ext.val.b[ , -1]) %>%
  cbind(ext.val.b[ , 1]) %>% data.frame() %>%
  dplyr::rename(., pred = `.`, obs = V2) %>%
  mutate(model = "Random Forest")

ev.b.svm <- predict(svm.b, ext.val.b[ , -1]) %>%
  cbind(ext.val.b[ , 1]) %>% data.frame() %>%
  dplyr::rename(., pred = `.`, obs = V2) %>%
  mutate(model = "SVM")

ev.b <- rbind(ev.b.cube, ev.b.glm, ev.b.pls, ev.b.rf, ev.b.svm)
ev.b.avg <- ev.b %>% data.table(key = "obs")
ev.b.avg <- ev.b.avg[ , list(pred = mean(pred)), by = "obs"]

defaultSummary(ev.b.avg) # 0.803

# Gamma-CD
ev.c.cube <- predict(cube.c, ext.val.c[ , -1]) %>%
  cbind(ext.val.c[ , 1]) %>%
  data.frame() %>%
  dplyr::rename(., pred = `.`, obs = V2) %>%
  mutate(model = "Cubist")

ev.c.glm <-  predict.glmnet(glm.c, as.matrix(ext.val.c[ , -1]), 
                            s = tail(glm.c$lambda, n = 1)) %>%
  cbind(as.matrix(ext.val.c[ , 1])) %>% data.frame() %>%
  dplyr::rename(pred = X1, obs = V2) %>%
  mutate(model = "GLMNet")

ev.c.pls <- predict(pls.c, ncomp = 4, newdata = ext.val.c[ , -1]) %>%
  cbind(ext.val.c[ , 1]) %>%
  data.frame() %>%
  dplyr::rename(., pred = `.`, obs = V2) %>%
  mutate(model = "PLS")

ev.c.rf <-predict(rf.c, ext.val.c[ , -1]) %>%
  cbind(ext.val.c[ , 1]) %>% data.frame() %>%
  dplyr::rename(., pred = `.`, obs = V2) %>%
  mutate(model = "Random Forest")

ev.c.svm <- predict(svm.c, ext.val.c[ , -1]) %>%
  cbind(ext.val.c[ , 1]) %>% data.frame() %>%
  dplyr::rename(., pred = `.`, obs = V2) %>%
  mutate(model = "SVM")

ev.c <- rbind(ev.c.cube, ev.c.glm, ev.c.pls, ev.c.rf, ev.c.svm)
ev.c.avg <- ev.c %>% data.table(key = "obs")
ev.c.avg <- ev.c.avg[ , list(pred = mean(pred)), by = "obs"]

defaultSummary(ev.c.avg)

# Combining CD
temp.a <- ev.a.avg %>% mutate(cd.type = "alpha") 
temp.b <- ev.b.avg %>% mutate(cd.type = "beta")
temp.c <- ev.c.avg %>% mutate(cd.type = "gamma")

ev.abc <- rbind(temp.a, temp.b, temp.c)
defaultSummary(ev.abc) # 0.417
saveRDS(ev.abc, "./ext val results.RDS")

# A single outlier brings down the r-squared
# removing the first row of the alpha results...

# temp.a <- ev.a.avg[-1, ] %>% mutate(cd.type = "alpha") 
# 
# ev.abc <- rbind(temp.a, temp.b, temp.c)
# defaultSummary(ev.abc) # 0.723

# ...gets a much higher R-squared

# Graphing ----------------------------------------------------------------

ggplot(ev.a, aes(x = obs, y = pred, color = model)) +
  geom_point() + 
  geom_abline(slope = 1, intercept = 0) + 
  facet_grid(~model)
ggplot(ev.a.avg, aes(x = obs, y = pred)) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0)

ggplot(ev.b, aes(x = obs, y = pred, color = model)) +
  geom_point() + 
  geom_abline(slope = 1, intercept = 0) + 
  facet_grid(~model)
ggplot(ev.b.avg, aes(x = obs, y = pred)) +
  geom_point() + 
  geom_abline(slope = 1, intercept = 0) 

ggplot(ev.abc, aes(x = obs, y = pred, color = cd.type)) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0, color = "maroon") + 
  theme.2018 + 
  labs(x = "Experimental dG, kJ/mol", y = "Predicted dG, kJ/mol", 
       title = "QSPR ensemble results on external validation", 
       color = "CD Type") + 
  coord_fixed(xlim = c(-45, 5), ylim = c(-45, 5))
ggsave("./ext val graph.png", dpi = 600)  
