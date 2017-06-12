set.seed(999)
library(dplyr)
library(ggplot2)
library(Matrix)
library(glmnet)
library(caret)
library(e1071)

# Preprocessing -----------------------------------------------------------

#     1. Variance Cleaning
# Checks for predictors that don't vary much
# and then removes them.
# Note: rp stands for ri.padel and rpt stands for ri.padel.trans 
zero.pred <- nearZeroVar(ri.padel)
rp.no.zero <- ri.padel[, -zero.pred] # View(rp.no.zero)
rpt <- preProcess(rp.no.zero, na.remove = T, 
                  method = c("center", "scale"
                             # , "BoxCox" throws an error: Na/NaN/Inf in foreign function call
                  )) %>%
  predict(., ri.padel) 


# Beta Model (Lasso) ------------------------------------------------------

rpt.bp <-  rpt %>%
  filter(host == "Beta-CD") %>% 
  filter(!is.na(log.K)) %>%
  dplyr::select(., -X1:-Solvent.composition, -log.K.Uncertainty:-`bind.aff, kcal/mol`)

rpt.btherm <-  rpt %>%
  filter(host == "Beta-CD") %>%
  dplyr::select(., -X1:-log.K.Uncertainty, -`DelG.Uncertainty`:-`bind.aff, kcal/mol`)

rpt.b.sprse <- sparse.model.matrix(~., rpt.btherm)
set.seed(512)
# This sample is the only place where the seed is relevant
rpt.b.trn.ind <- sample(x = 1:nrow(rpt.b.sprse), size = round(0.75*nrow(rpt.b.sprse)))

rpt.b.trn <- rpt.b.sprse[rpt.b.trn.ind,  ]
rpt.b.tst <- rpt.b.sprse[-rpt.b.trn.ind, ]
# -1:-6, alternatively
rpt.b.glm.2 <- cv.glmnet(x = rpt.b.trn[,-1:-7], y = rpt.b.trn[, 6])
rpt.b.glm.3 <- cv.glmnet(x = rpt.b.trn[,-1], y = rpt.b.trn[, 2])

# Visualization of the model on test
predict(rpt.b.glm.2, rpt.b.tst[,-1:-7]) %>% 
  cbind(rpt.b.tst[,6]) %>% 
  data.frame() %>% 
  rename(predicted = X1, experimental = V2) %>%
  ggplot(., aes(x = experimental, y = predicted))+
  geom_point()+
  geom_abline(slope = 1, intercept = 0)
# Visualization of the model on training
predict(rpt.b.glm.2, rpt.b.trn[,-1:-7]) %>% 
  cbind(rpt.b.trn[,6]) %>% 
  data.frame() %>% 
  rename(predicted = X1, experimental = V2) %>%
  ggplot(., aes(x = experimental, y = predicted))+
  geom_point()+
  geom_abline(slope = 1, intercept = 0)

rpt.b.glm.3 <- cv.glmnet(x = rpt.b.trn[,-1:-2], y = rpt.b.trn[,2], nfolds = 68, grouped = F)
predict(rpt.b.glm.3, rpt.b.tst[,-1:-2]) %>% 
  cbind(rpt.b.tst[,2]) %>% 
  data.frame() %>% 
  rename(predicted = X1, experimental = V2) %>%
  ggplot(., aes(x = experimental, y = predicted))+
  geom_point()+
  geom_abline(slope = 1, intercept = 0)
predict(rpt.b.glm.3, rpt.b.trn[,-1:-2]) %>% 
  cbind(rpt.b.trn[,2]) %>% 
  data.frame() %>% 
  rename(predicted = X1, experimental = V2) %>%
  ggplot(., aes(x = experimental, y = predicted))+
  geom_point()+
  geom_abline(slope = 1, intercept = 0)

# Alpha Model (Lasso) -----------------------------------------------------
rpt.ap <-  rpt %>%
  filter(host == "Alpha-CD") %>% 
  filter(!is.na(log.K)) %>%
  dplyr::select(., -X1:-Solvent.composition, -log.K.Uncertainty:-`bind.aff, kcal/mol`)

rpt.atherm <-  rpt %>%
  filter(host == "Alpha-CD") %>%
  dplyr::select(X1:Solvent.composition, log.K.Uncertainty:`bind.aff, kcal/mol`)
rpt.a.sprse <- sparse.model.matrix(~., rpt.ap)

rpt.a.trn.ind <- sample(x = 1:nrow(rpt.a.sprse), size = round(0.75*nrow(rpt.a.sprse)))

rpt.a.trn <- rpt.a.sprse[rpt.a.trn.ind, ]
rpt.a.tst <- rpt.a.sprse[-rpt.a.trn.ind,]
# rpt.a.glm <- glmnet(x = rpt.a.trn[,-1:-6], y = rpt.a.trn[, 6])
rpt.a.glm.2 <- cv.glmnet(x = rpt.a.trn[,-1:-6], y = rpt.a.trn[, 6], nfolds = 48)

predict(rpt.glm.2, rpt.a.tst[,-1:-6]) %>% 
  cbind(rpt.a.tst[,6]) %>% 
  data.frame() %>% 
  rename(predicted = X1, experimental = V2) %>%
  ggplot(., aes(x = experimental, y = predicted))+
  geom_point()+
  geom_abline(slope = 1, intercept = 0)

# Gamma Model (Lasso) -----------------------------------------------------
rpt.cp <-  rpt %>%
  filter(host == "Gamma-CD") %>% 
  filter(!is.na(log.K)) %>%
  dplyr::select(., -X1:-Solvent.composition, -log.K.Uncertainty:-`bind.aff, kcal/mol`)

rpt.ctherm <-  rpt %>%
  filter(host == "Gamma-CD") %>%
  dplyr::select(X1:Solvent.composition, log.K.Uncertainty:`bind.aff, kcal/mol`)
rpt.c.sprse <- sparse.model.matrix(~., rpt.cp)

rpt.c.trn.ind <- sample(x = 1:nrow(rpt.c.sprse), size = round(0.75*nrow(rpt.c.sprse)))

rpt.c.trn <- rpt.c.sprse[rpt.c.trn.ind, ]
rpt.c.tst <- rpt.c.sprse[-rpt.c.trn.ind,]

rpt.c.glm <- cv.glmnet(x = rpt.c.sprse[,-1:-6], y = rpt.c.sprse[, 6])

predict(rpt.c.glm, rpt.c.tst[,-1:-6]) %>% 
  cbind(rpt.c.tst[,6]) %>% 
  data.frame() %>% 
  rename(predicted = X1, experimental = V2) %>%
  ggplot(., aes(x = experimental, y = predicted))+
  geom_point()+
  geom_abline(slope = 1, intercept = 0)


# All 3 CDs Model (Lasso) -------------------------------------------------

rpt.comp <- rpt %>% mutate(cd.type = ifelse(str_detect(host, "Alpha"), 1, ifelse(str_detect(host, "Beta"), 2, 3)))
rpt.p <-  rpt.comp %>%
  filter(!is.na(log.K)) %>%
  dplyr::select(., -X1:-Solvent.composition, -log.K.Uncertainty:-`bind.aff, kcal/mol`)

rpt.therm <-  rpt.comp %>%
  dplyr::select(., -X1:-log.K.Uncertainty, -`DelG.Uncertainty`:-`bind.aff, kcal/mol`)
rpt.sprse <- sparse.model.matrix(~., rpt.p)
rpt.sprse.dg <- sparse.model.matrix(~., rpt.therm)

rpt.trn.ind <- sample(x = 1:nrow(rpt.sprse), size = round(0.75*nrow(rpt.sprse)))

rpt.trn <- rpt.sprse[rpt.trn.ind, ]
rpt.tst <- rpt.sprse[-rpt.trn.ind,]
# For ALogP
rpt.glm <- cv.glmnet(x = rpt.trn[,-1:-6], y = rpt.trn[, 6])

    # Visualization on test set
predict(rpt.glm, rpt.tst[,-1:-6]) %>% 
  cbind(rpt.tst[,6], as.factor(rpt.tst[,2761])) %>% 
  data.frame() %>% 
  rename(predicted = X1, experimental = V2, cd.type = V3) %>%
  ggplot(., aes(x = experimental, y = predicted)) +
  geom_point(aes(color = factor(cd.type))) +
  geom_abline(slope = 1, intercept = 0) + 
  theme_bw() + 
  xlab("Experimental ALogP") +
  ylab("Predicted ALogP")

    # Visualization on training set
predict(rpt.glm, rpt.trn[,-1:-6]) %>% 
  cbind(rpt.trn[,6]) %>% 
  data.frame() %>% 
  rename(predicted = X1, experimental = V2) %>%
  ggplot(., aes(x = experimental, y = predicted))+
  geom_point()+
  geom_abline(slope = 1, intercept = 0)

# For Log.K
rpt.glm.lk <- cv.glmnet(x = rpt.trn[,-4], y = rpt.trn[, 4])

    # Visualization on test set
predict(rpt.glm.lk, rpt.tst[,-4]) %>% 
  cbind(rpt.tst[, 4], as.factor(rpt.tst[,2761])) %>% 
  data.frame() %>% 
  rename(predicted = X1, experimental = V2, cd.type = V3) %>%
  ggplot(., aes(x = experimental, y = predicted)) +
  geom_point(aes(color = factor(cd.type))) +
  geom_abline(slope = 1, intercept = 0) + 
  theme_bw() + 
  xlab("Experimental LogK") +
  ylab("Predicted LogK")

    # Visualization on training set
predict(rpt.glm.lk, rpt.trn[,-1:-6]) %>% 
  cbind(rpt.trn[, 4]) %>% 
  data.frame() %>% 
  rename(predicted = X1, experimental = V2) %>%
  ggplot(., aes(x = experimental, y = predicted))+
  geom_point()+
  geom_abline(slope = 1, intercept = 0)

# DelG
rpt.sprse.dg <- sparse.model.matrix(~., rpt.therm)
set.seed(222)
rpt.trn.ind.dg <- sample(x = 1:nrow(rpt.sprse.dg), size = round(0.75*nrow(rpt.sprse.dg)))

rpt.trn.dg <- rpt.sprse.dg[rpt.trn.ind, ]
rpt.tst.dg <- rpt.sprse.dg[-rpt.trn.ind,]
rpt.glm.dg <- cv.glmnet(x = rpt.trn.dg[,-2], y = rpt.trn.dg[, 2])
predict(rpt.glm.dg, rpt.trn.dg[,-2], s = "lambda.min") %>% 
  cbind(rpt.trn.dg[, 2]) %>% 
  data.frame() %>% 
  rename(predicted = X1, experimental = V2) %>%
  ggplot(., aes(x = experimental, y = predicted))+
  geom_point()+
  geom_abline(slope = 1, intercept = 0)
predict(rpt.glm.dg, rpt.tst.dg[,-2], s = "lambda.min") %>% 
  cbind(rpt.tst.dg[, 2]) %>% 
  data.frame() %>% 
  rename(predicted = X1, experimental = V2) %>%
  ggplot(., aes(x = experimental, y = predicted))+
  geom_point()+
  geom_abline(slope = 1, intercept = 0)

# Using "glmulti"
install.packages("glmulti")
library(glmulti)
View(rpt.therm)
    # In the interest of saving time, the df will be summarized to the 15
    # most important factors as found in random forest
set.seed(512)
rpt.trn.ind <- sample(x = 1:nrow(rpt.comp), size = round(0.75*nrow(rpt.comp)))

rpt.trn.df <- rpt.comp[rpt.trn.ind, ]
rpt.tst.df <- rpt.comp[-rpt.trn.ind,]

glm.x <- rpt.trn.df[ , c(row.names(rf.factors.df)[1:10], "DelG")]
glm.y <- rpt.trn.df[ , "DelG"]
dg.glmulti <- glmulti(y = "DelG", xr = colnames(glm.x), data = glm.x, method = "g", 
                      mutrate = 0.05, deltaB = 0.1, conseq = 2, confsetsize = 64)
summary.glmulti(dg.glmulti)$bestmodel

# Beta-CD (SVM) -----------------------------------------------------------

rpt.b.x <- rpt.b.trn[, -1:-2]
rpt.b.y <- rpt.b.trn[ , 2]
rpt.b.svm <- svm(rpt.b.x, rpt.b.y, kernel = "linear", cross = 68)

rpt.svm.b.trn <- predict(rpt.b.svm, rpt.b.x)%>% 
  cbind(rpt.b.y) %>% 
  data.frame()
colnames(rpt.svm.b.trn)[1] <- "predicted"
colnames(rpt.svm.b.trn)[2] <- "experimental"
# Visualization of training data set
rpt.svm.b.trn %>%
  ggplot(., aes(x = experimental, y = predicted))+
  geom_point()+
  geom_abline(slope = 1, intercept = 0)

rpt.svm.b.tst <- predict(rpt.b.svm, rpt.b.tst[,-1:-2]) %>% 
  cbind(rpt.b.tst[, 2]) %>% 
  data.frame()
colnames(rpt.svm.b.tst)[1] <- "predicted"
colnames(rpt.svm.b.tst)[2] <- "experimental"
# Visualization of test data set
rpt.svm.b.tst %>%
  ggplot(., aes(x = experimental, y = predicted))+
  geom_point()+
  geom_abline(slope = 1, intercept = 0)
# All 3 CDs (SVM) ---------------------------------------------------------
zero.pred <- nearZeroVar(rpt.comp)
rpt.comp.nz <- rpt.comp[, -zero.pred]
rpt.comp.t <- preProcess(rpt.comp.nz, na.remove = T, 
                         method = c("center", "scale"
                                    # , "BoxCox" throws an error: Na/NaN/Inf in foreign function call
                         )) %>%
  predict(., rpt.comp.nz)
rpt.therm <-  rpt.comp.t %>%
  dplyr::select(., -X1:-log.K.Uncertainty, -`DelG.Uncertainty`:-`bind.aff, kcal/mol`)
rpt.sprse.dg <- sparse.model.matrix(~., rpt.therm)

# DelG
rpt.trn.ind <- sample(x = 1:nrow(rpt.sprse.dg), size = round(0.8*nrow(rpt.sprse.dg)))
rpt.trn.dg <- rpt.sprse.dg[rpt.trn.ind, ]
rpt.tst.dg <- rpt.sprse.dg[-rpt.trn.ind,]
rpt.y <- rpt.trn.dg[,2]
rpt.x <- rpt.trn.dg[,-1:-2]
rpt.svm <- e1071::svm(x = rpt.x, y = rpt.y, kernel = "polynomial", cross=125)


rpt.svm.df.trn <- predict(rpt.svm, rpt.x)%>% 
  cbind(rpt.trn.dg[,2]) %>% 
  data.frame()
colnames(rpt.svm.df.trn)[1] <- "predicted"
colnames(rpt.svm.df.trn)[2] <- "experimental"
# Visualization of training data set
rpt.svm.df.trn %>%
  ggplot(., aes(x = experimental, y = predicted))+
  geom_point()+
  geom_abline(slope = 1, intercept = 0)

rpt.svm.df.tst <- predict(rpt.svm, rpt.tst.dg[,-1:-2]) %>% 
  cbind(rpt.tst.dg[,2]) %>% 
  data.frame()
colnames(rpt.svm.df.tst)[1] <- "predicted"
colnames(rpt.svm.df.tst)[2] <- "experimental"
# Visualization of test data set
rpt.svm.df.tst %>%
  ggplot(., aes(x = experimental, y = predicted))+
  geom_point()+
  geom_abline(slope = 1, intercept = 0)

 

# All 3 CDs (RF) ----------------------------------------------------------

install.packages("randomForest")
library(randomForest)
# Matrix for DelG
rpt.mat <- rpt.comp.t %>%
  dplyr::select(., -X1:-log.K.Uncertainty, -`DelG.Uncertainty`:-`bind.aff, kcal/mol`) %>%
  data.matrix()

set.seed(512)
# This sample is the only place where the seed is relevant
rpt.rf.trn.ind <- sample(x = 1:nrow(rpt.mat), size = round(0.75*nrow(rpt.mat)))

rpt.rf.trn <- rpt.mat[rpt.rf.trn.ind,  ]
rpt.rf.tst <- rpt.mat[-rpt.rf.trn.ind,  ] %>% as.data.frame()
rpt.rf.y <- rpt.rf.trn[,1] 
rpt.rf.x <- rpt.rf.trn[,-1] 
# rpt.rf <- randomForest(DelG ~., data = rpt.rf.df, ntree = 1000)
rpt.rf <- rfImpute(x = rpt.rf.x, y = rpt.rf.y, ntree = 100, na.action = na.omit)
rpt.rf.mod <- randomForest(y ~., data = rpt.rf, ntree = 100, na.action = na.omit)
rpt.rf.df <- predict(rpt.rf.mod, rpt.rf.tst[,-1]) %>% 
  cbind(rpt.rf.tst[,1]) %>% 
  data.frame() 
colnames(rpt.rf.df)[1] <- "predicted"
colnames(rpt.rf.df)[2] <- "experimental"
# Visualization of test data set
rpt.rf.df %>%
  ggplot(., aes(x = experimental, y = predicted))+
  geom_point()+
  geom_abline(slope = 1, intercept = 0)

rpt.rf.df <- predict(rpt.rf.mod, rpt.rf.trn[,-1]) %>% 
  cbind(rpt.rf.trn[,1]) %>% 
  data.frame() 
colnames(rpt.rf.df)[1] <- "predicted"
colnames(rpt.rf.df)[2] <- "experimental"
# Visualization of training data set
rpt.rf.df %>%
  ggplot(., aes(x = experimental, y = predicted))+
  geom_point()+
  geom_abline(slope = 1, intercept = 0)

# Analysis of important factors
rf.factors.df <- importance(rpt.rf.mod) %>% data.frame()
rf.factors.df <- rf.factors.df[order(-rf.factors.df$IncNodePurity), , drop = F] %>% data.frame()
rf.factors <- row.names(rf.factors.df)[1:10]

# colnames(rpt.therm) %in% rf.factors


