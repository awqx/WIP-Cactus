# Libraries and Packages --------------------------------------------------

library(caret)
library(e1071)
library(kernlab)
library(Matrix)
library(stats)
library(stringr)
library(tidyverse)

# Loading Data ------------------------------------------------------------

setwd("~/SREP LAB/qsar")
dir.create("./models")
dir.create("./models/svm")
df.raw <- readRDS("./padel.pp.new.RDS") 
df <- df.raw %>%
  select(., -guest:-host) %>%
  select(., -data.source)

sprse <- sparse.model.matrix(~., df)

set.seed(25)
trn.ind <- sample(x = 1:nrow(df), size = round(0.7 * nrow(df)))
trn <- df[trn.ind, ]
tst <- df[-trn.ind, ]

sprse <- sparse.model.matrix(~., df)
trn.ind <- sample(x = 1:nrow(sprse), size = round(0.7 * nrow(sprse)))

sprse.trn <- sprse[trn.ind, ]
sprse.trn.x <- sprse.trn[ , -1:-2]
sprse.trn.y <- sprse.trn[ , 2]

sprse.tst <- sprse[-trn.ind, ]
sprse.tst.x <- sprse.tst[ , -1:-2]
sprse.tst.y <- sprse.tst[ , 2]


# Polynomial Kernel -------------------------------------------------------

#     All Data ----------------------------------------------------

svm.all <- svm(x = sprse.trn.x, 
               y = sprse.trn.y, 
               coef0 = 1.5, 
               cost = 4096, 
               epsilon = 1, 
               kernel = "polynomial", 
               gamma = 0.03, 
               degree = 2)

svm.all.tst <- predict(svm.all, sprse.tst.x) %>%
  cbind(sprse.tst.y) %>%
  data.frame() %>%
  rename(., pred = `.`, obs = sprse.tst.y)
svm.all.trn <- predict(svm.all, sprse.trn.x) %>%
  cbind(sprse.trn.y) %>%
  data.frame() %>%
  rename(., pred = `.`, obs = sprse.trn.y)

defaultSummary(svm.all.tst)[2] # 0.529
defaultSummary(svm.all.trn)[2] # 0.969

#     Alpha -------------------------------------------------------

df.a <- filter(df, alpha > 0)
sprse.a <- sparse.model.matrix(~., df.a)

set.seed(24)
trn.ind <- sample(x = 1:nrow(sprse.a), 
                  size = round(0.7 * nrow(sprse.a)))
sprse.a.trn <- sprse.a[trn.ind, ]
sprse.a.tst <- sprse.a[-trn.ind, ]

svm.a <- svm(x = sprse.a.trn[ , -1:-2],
             y = sprse.a.trn[ , 2], 
             kernel = "polynomial", 
             degree = 2,
             cost = 4096, 
             gamma = 0.03, 
             epsilon = 1, 
             coef0 = 1.5)
svm.a.tst <- predict(svm.a, sprse.a.tst[ , -1:-2]) %>%
  cbind(sprse.a.tst[ , 2]) %>%
  data.frame() %>%
  rename(., pred = `.`, obs = V2)

defaultSummary(svm.a.tst)[2]

#     Beta --------------------------------------------------------
df.b <- filter(df, beta > 0)
sprse.b <- sparse.model.matrix(~., df.b)

set.seed(25)
trn.ind <- sample(x = 1:nrow(sprse.b), 
                  size = round(0.7 * nrow(sprse.b)))
sprse.b.trn <- sprse.b[trn.ind, ]
sprse.b.tst <- sprse.b[-trn.ind, ]

svm.b <- svm(x = sprse.b.trn[ , -1:-2],
             y = sprse.b.trn[ , 2], 
             kernel = "polynomial", 
             degree = 2,
             cost = 4096, 
             gamma = 0.03, 
             epsilon = 1, 
             coef0 = 1.5)
svm.b.tst <- predict(svm.b, sprse.b.tst[ , -1:-2]) %>%
  cbind(sprse.b.tst[ , 2]) %>%
  data.frame() %>%
  rename(., pred = `.`, obs = V2)

defaultSummary(svm.b.tst)[2]

#     Gamma -------------------------------------------------------

df.c <- filter(df, gamma > 0)
sprse.c <- sparse.model.matrix(~., df.c)

set.seed(12)
trn.ind <- sample(x = 1:nrow(sprse.c), 
                  size = round(0.7 * nrow(sprse.c)))
sprse.c.trn <- sprse.c[trn.ind, ]
sprse.c.tst <- sprse.c[-trn.ind, ]

svm.c <- svm(x = sprse.c.trn[ , -1:-2],
             y = sprse.c.trn[ , 2], 
             kernel = "polynomial", 
             degree = 2,
             cost = 4096, 
             gamma = 0.03, 
             epsilon = 1, 
             coef0 = 1.5)
svm.c.tst <- predict(svm.c, sprse.c.tst[ , -1:-2]) %>%
  cbind(sprse.c.tst[ , 2]) %>%
  data.frame() %>%
  rename(., pred = `.`, obs = V2)

defaultSummary(svm.c.tst)[2]

#     Compiled CDs --------------------------------------------------------

temp.a <- svm.a.tst %>% 
  mutate(cd.type = rep("alpha", length(svm.a.tst$pred)))
temp.b <- svm.b.tst %>% 
  mutate(cd.type = rep("beta", length(svm.b.tst$pred)))
temp.c <- svm.c.tst %>% 
  mutate(cd.type = rep("gamma", length(svm.c.tst$pred)))
svm.abc.tst <- rbind(temp.a, temp.b, temp.c, 
                     make.row.names = F) %>%
  mutate(residual = pred - obs)
defaultSummary(svm.abc.tst)

#####
# Radial Basis Kernel -----------------------------------------------------
#     All ------------------------------------------------------------
rbf.all <- svm(x = sprse.trn.x, 
               y = sprse.trn.y,
               cost = 2048, 
               epsilon = 0.5, 
               kernel = "radial", 
               gamma = 0.001)

rbf.all.tst <- predict(rbf.all, sprse.tst.x) %>%
  cbind(sprse.tst.y) %>%
  data.frame() %>%
  rename(., pred = `.`, obs = sprse.tst.y)
rbf.all.trn <- predict(rbf.all, sprse.trn.x) %>%
  cbind(sprse.trn.y) %>%
  data.frame() %>%
  rename(., pred = `.`, obs = sprse.trn.y)

defaultSummary(rbf.all.tst)[2]

#     Alpha ---------------------------------------------------------

rbf.a <- svm(x = sprse.a.trn[ , -1:-2],
             y = sprse.a.trn[ , 2], 
             kernel = "radial", 
             cost = 8, 
             gamma = 0.0005, 
             epsilon = 0.5)
rbf.a.tst <- predict(rbf.a, sprse.a.tst[ , -1:-2]) %>%
  cbind(sprse.a.tst[ , 2]) %>%
  data.frame() %>%
  rename(., pred = `.`, obs = V2)

defaultSummary(rbf.a.tst)[2]

#     Beta -----------------------------------------------------------

rbf.b <- svm(x = sprse.b.trn[ , -1:-2],
             y = sprse.b.trn[ , 2], 
             kernel = "radial", 
             cost = 8, 
             gamma = 0.0005, 
             epsilon = 0.5)
rbf.b.tst <- predict(rbf.b, sprse.b.tst[ , -1:-2]) %>%
  cbind(sprse.b.tst[ , 2]) %>%
  data.frame() %>%
  rename(., pred = `.`, obs = V2)

defaultSummary(rbf.b.tst)[2]

#     Gamma ----------------------------------------------------------

rbf.c <- svm(x = sprse.c.trn[ , -1:-2],
             y = sprse.c.trn[ , 2], 
             kernel = "radial", 
             cost = 8, 
             gamma = 0.0005, 
             epsilon = 0.5)
rbf.c.tst <- predict(rbf.c, sprse.c.tst[ , -1:-2]) %>%
  cbind(sprse.c.tst[ , 2]) %>%
  data.frame() %>%
  rename(., pred = `.`, obs = V2)

defaultSummary(rbf.c.tst)[2]

#     Compiled CDs --------------------------------------------------------
temp.a <- rbf.a.tst %>% 
  mutate(cd.type = rep("alpha", length(rbf.a.tst$pred)))
temp.b <- rbf.b.tst %>% 
  mutate(cd.type = rep("beta", length(rbf.b.tst$pred)))
temp.c <- rbf.c.tst %>% 
  mutate(cd.type = rep("gamma", length(rbf.c.tst$pred)))
rbf.abc.tst <- rbind(temp.a, temp.b, temp.c, 
                     make.row.names = F)
defaultSummary(rbf.abc.tst)

# Sigmoid Kernel -----------------------------------------------------
#     All ------------------------------------------------------------
sig.all <- svm(x = sprse.trn.x, 
               y = sprse.trn.y,
               cost = 1, 
               epsilon = 0.125, 
               kernel = "sigmoid", 
               gamma = 0.001)

sig.all.tst <- predict(sig.all, sprse.tst.x) %>%
  cbind(sprse.tst.y) %>%
  data.frame() %>%
  rename(., pred = `.`, obs = sprse.tst.y)
sig.all.trn <- predict(sig.all, sprse.trn.x) %>%
  cbind(sprse.trn.y) %>%
  data.frame() %>%
  rename(., pred = `.`, obs = sprse.trn.y)

defaultSummary(sig.all.tst)[2]
defaultSummary(sig.all.trn)
# Plots -------------------------------------------------------------------
# Notes: Outlier is 3-methylbenzoic acid
#     Polynomial ----------------------------------------------------------
#         Test Set --------------------------------------------------------
# SVM with all data points
ggplot(svm.all.tst, aes(x = obs, y = pred)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Experimental DelG, kJ/mol", y = "Predicted DelG, kJ/mol", 
       title = "Polynomial SVM - All Data Points") + 
  coord_fixed() + 
  theme_bw()
ggsave("./models/svm/2017-07-20 poly all data.png")

# Alpha-CD
ggplot(svm.a.tst, aes(x = obs, y = pred)) + 
  geom_point(fill = "salmon") + 
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Experimental DelG, kJ/mol", y = "Predicted DelG, kJ/mol", 
       title = "Polynomial SVM - Alpha CD") + 
  coord_fixed() + 
  theme_bw()
ggsave("./models/svm/2017-07-20 poly alpha cd.png")

# Beta-CD
ggplot(svm.b.tst, aes(x = obs, y = pred)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Experimental DelG, kJ/mol", y = "Predicted DelG, kJ/mol", 
       title = "Polynomial SVM - Beta CD") + 
  coord_fixed() + 
  theme_bw()
ggsave("./models/svm/2017-07-20 poly beta cd.png")

# Gamma-CD
ggplot(svm.c.tst, aes(x = obs, y = pred)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Experimental DelG, kJ/mol", y = "Predicted DelG, kJ/mol", 
       title = "Polynomial SVM - Gamma CD") + 
  coord_fixed() + 
  theme_bw()
ggsave("./models/svm/2017-07-20 poly gamma cd.png")

# All-CDs
ggplot(svm.abc.tst, aes(x = obs, y = pred, 
                        group = cd.type, 
                        color = cd.type)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Experimental DelG, kJ/mol", y = "Predicted DelG, kJ/mol", 
       title = "Polynomial SVM - Compiled CD Types", 
       color = "Cyclodextrin Type") + 
  coord_fixed() + 
  theme_bw()
ggsave("./models/svm/2017-07-20 poly compiled cds.png")

#    Residuals
ggplot(svm.abc.tst, aes(x = obs, y = residual, color = cd.type)) + 
  geom_point() + 
  geom_hline(yintercept = 0) + 
  labs(x = "Experimental DelG, kJ/mol", y = "Residual, kJ/mol", 
       title = "Polynomial SVM - Residuals of Compiled CD Types") + 
  coord_fixed() + 
  theme_bw()
ggsave("./models/svm/2017-07-20 poly compiled cds resid.png")

#         Training Set ----------------------------------------------------
# Very clearly overfitted
ggplot(svm.all.trn, aes(x = obs, y = pred)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Experimental DelG, kJ/mol", y = "Predicted DelG, kJ/mol", 
       title = "Polynomial SVM - All Data Points, Training") + 
  coord_fixed() + 
  theme_bw()
ggsave("./models/svm/2017-07-07 poly all data points trn.png")

#     RBF -----------------------------------------------------------------
ggplot(rbf.all.tst, aes(x = obs, y = pred)) + 
  geom_point(shape = 10, size = 2) + 
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Experimental DelG, kJ/mol", y = "Predicted DelG, kJ/mol", title = "Radial SVM - All Data Points") + 
  coord_fixed() + 
  theme_bw()
ggsave("./models/svm/2017-07-20 rbf all data points.png")

ggplot(rbf.a.tst, aes(x = obs, y = pred)) + 
  geom_point(shape = 10, size = 2) + 
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Experimental DelG, kJ/mol", y = "Predicted DelG, kJ/mol", 
       title = "Radial SVM - Alpha CD") + 
  coord_fixed() + 
  theme_bw()
ggsave("./models/svm/2017-07-20 rbf alpha cd.png")

ggplot(rbf.b.tst, aes(x = obs, y = pred)) + 
  geom_point(shape = 10, size = 2) + 
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Experimental DelG, kJ/mol", y = "Predicted DelG, kJ/mol", title = "Radial SVM - Beta CD") + 
  coord_fixed() + 
  theme_bw()
ggsave("./models/svm/2017-07-20 rbf beta cd.png")

ggplot(rbf.c.tst, aes(x = obs, y = pred)) + 
  geom_point(shape = 10, size = 2) + 
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Experimental DelG, kJ/mol", y = "Predicted DelG, kJ/mol", title = "Radial SVM - Gamma CD") + 
  coord_fixed() + 
  theme_bw()
ggsave("./models/svm/2017-07-20 rbf gamma cd.png")

ggplot(rbf.abc.tst, aes(x = obs, y = pred, color = cd.type)) + 
  geom_point(shape = 10, size = 2) + 
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Experimental DelG, kJ/mol", y = "Predicted DelG, kJ/mol", 
       title = "Radial SVM - Compiled CDs", 
       color = "Cyclodextrin Type") + 
  coord_fixed() + 
  theme_bw()
ggsave("./models/svm/2017-07-20 rbf compiled cd.png")
#####
# Model with Untransformed Data -------------------------------------------

df.raw <- readRDS("./padel.nopp.RDS") 
df <- df.raw %>%
  select(., -guest:-host) %>%
  select(., -data.source)
df <- lapply(df, as.numeric) %>% data.frame()
set.seed(25)
trn.ind <- sample(x = 1:nrow(df), size = round(0.7 * nrow(df)))
trn <- df[trn.ind, ]
tst <- df[-trn.ind, ]

sprse <- sparse.model.matrix(~., df)
trn.ind <- sample(x = 1:nrow(sprse), size = round(0.7 * nrow(sprse)))

sprse.trn <- sprse[trn.ind, ]
sprse.trn.x <- sprse.trn[ , -1:-2]
sprse.trn.y <- sprse.trn[ , 2]

sprse.tst <- sprse[-trn.ind, ]
sprse.tst.x <- sprse.tst[ , -1:-2]
sprse.tst.y <- sprse.tst[ , 2]

svm.all <- svm(x = trn[ , -1], 
               y = trn[ , 1], 
               coef0 = 1.5, 
               cost = 256, 
               epsilon = 1, 
               kernel = "polynomial", 
               gamma = 0.03, 
               degree = 2)

svm.all.tst <- predict(svm.all, tst[ , -1]) %>%
  cbind(tst[ , 1]) %>%
  data.frame() %>%
  rename(., pred = `.`, obs = V2)
svm.all.trn <- predict(svm.all, trn[ , -1]) %>%
  cbind(trn[ , 1]) %>%
  data.frame() %>%
  rename(., pred = `.`, obs = V2)

defaultSummary(svm.all.tst)[2] # 0.529
defaultSummary(svm.all.trn)[2] # 0.969
