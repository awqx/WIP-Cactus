
# Download Models ---------------------------------------------------------

# setwd("~/SREP LAB/qsar/models/")
cube.alpha <- readRDS("./cubist/cube.alpha.RDS")
cube.beta <- readRDS("./cubist/cube.beta.RDS")
cube.gamma <- readRDS("./cubist/cube.gamma.RDS")

glm.alpha <- readRDS("./glmnet/glm.alpha.RDS")
glm.beta <- readRDS("./glmnet/glm.beta.RDS")
glm.gamma <- readRDS("./glmnet/glm.gamma.RDS")

rf.alpha <- readRDS("./rforest/rf.alpha.RDS")

ext.val <- readRDS("./external validation set new.RDS") %>%
  select(-guest:-host) %>%
  select(-data.source) 

ext.val.a <- ext.val %>% filter(alpha > 0) 
ext.val.b <- ext.val %>% filter(beta > 0)
ext.val.c <- ext.val %>% filter(gamma > 0) 

ev.a.glm <-  predict.glmnet(glm.alpha, as.matrix(ext.val.a[ , -1]), 
                        s = tail(glm.alpha$lambda, n = 1)) %>%
  cbind(as.matrix(ext.val.a[ , 1])) %>% data.frame() %>%
  dplyr::rename(pred = X1, obs.glm = V2)

ev.a.rf <-  predict(rf.alpha, ext.val.a[ , -1]) %>%
  cbind(ext.val.a[ , 1]) %>% data.frame() %>%
  dplyr::rename(., pred = `.`, obs.rf = V2)
