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
df <- readRDS("./rpt.RDS") %>%
  select(-guest:-host) %>% select(-data.source)
set.seed(2)
trn.ind <- sample(x = 1:nrow(df), 
                  size = round(0.7 * nrow(df)))
df.trn <- df[trn.ind, ]
df.tst <- df[-trn.ind, ]

set.seed(1)
# sprse <- readRDS("./DelG.sparse.RDS")
sprse <- readRDS(".feature.select/sprse.ga.RDS")
trn.ind <- sample(x = 1:nrow(sprse), 
                  size = round(0.7 * nrow(sprse)))
sprse.trn <- sprse[trn.ind, ]
sprse.tst <- sprse[-trn.ind, ]

# Polynomial Kernel -------------------------------------------------------
#     All Data ----------------------------------------------------

svm.all <- svm(x = sprse.trn[ , -1:-2], 
               y = sprse.trn[ , 2], 
               coef0 = 2, 
               cost = 8, 
               epsilon = 0.1, 
               kernel = "polynomial", 
               gamma = 0.5, 
               degree = 2)

svm.all.tst <- predict(svm.all, sprse.tst[ , -1:-2]) %>%
  cbind(sprse.tst[ , 2]) %>%
  data.frame() %>%
  rename(., pred = `.`, obs = V2)
svm.all.trn <- predict(svm.all, sprse.trn[ , -1:-2]) %>%
  cbind(sprse.trn[ , 2]) %>%
  data.frame() %>%
  rename(., pred = `.`, obs = V2)

defaultSummary(svm.all.tst)[2] # 0.00971
defaultSummary(svm.all.trn)[2] # 0.985

#     Alpha -------------------------------------------------------
df.a <- filter(df, alpha > 0)
set.seed(1)

sprse.a <- sparse.model.matrix(~., df.a)
trn.ind <- sample(x = 1:nrow(sprse.a), 
                  size = round(0.7 * nrow(sprse.a)))
sprse.a.trn <- sprse.a[trn.ind, ]
sprse.a.tst <- sprse.a[-trn.ind, ]

svm.a <- svm(x = sprse.a.trn[ , -1:-2],
             y = sprse.a.trn[ , 2], 
             kernel = "polynomial", 
             degree = 2,
             cost = 1024, 
             gamma = 0.5, 
             epsilon = 0.1, 
             coef0 = 2)
svm.a.tst <- predict(svm.a, sprse.a.tst[ , -1:-2]) %>%
  cbind(sprse.a.tst[ , 2]) %>%
  data.frame() %>%
  rename(., pred = `.`, obs = V2)

defaultSummary(svm.a.tst)[2] # 0.651

#     Beta --------------------------------------------------------
df.b <- filter(df, beta > 0)
# 20 # 34 (outlier) # 54
set.seed(2)
trn.ind <- sample(x = 1:nrow(df.b), 
                  size = round(0.7 * nrow(df.b)))
df.b.trn <- df.b[trn.ind, ]
df.b.tst <- df.b[-trn.ind, ]

sprse.b <- sparse.model.matrix(~., df.b)
trn.ind <- sample(x = 1:nrow(sprse.b), 
                  size = round(0.7 * nrow(sprse.b)))
sprse.b.trn <- sprse.b[trn.ind, ]
sprse.b.tst <- sprse.b[-trn.ind, ]

svm.b <- svm(x = sprse.b.trn[ , -1:-2],
             y = sprse.b.trn[ , 2], 
             kernel = "polynomial", 
             degree = 2,
             cost = 1024, 
             gamma = 0.5, 
             epsilon = 0.1, 
             coef0 = 2)
svm.b.tst <- predict(svm.b, sprse.b.tst[ , -1:-2]) %>%
  cbind(sprse.b.tst[ , 2]) %>%
  data.frame() %>%
  rename(., pred = `.`, obs = V2)

defaultSummary(svm.b.tst)[2] # 0.617

#     Gamma -------------------------------------------------------

df.c <- filter(df, gamma > 0)
set.seed(3)
trn.ind <- sample(x = 1:nrow(df.c), size = round(0.7 * nrow(df.c)))
df.c.trn <- df.c[trn.ind, ]
df.c.tst <- df.c[-trn.ind, ]

sprse.c <- sparse.model.matrix(~., df.c)
trn.ind <- sample(x = 1:nrow(sprse.c), size = round(0.7 * nrow(sprse.c)))
sprse.c.trn <- sprse.c[trn.ind, ]
sprse.c.tst <- sprse.c[-trn.ind, ]

svm.c <- svm(x = sprse.c.trn[ , -1:-2],
             y = sprse.c.trn[ , 2], 
             kernel = "polynomial", 
             degree = 2,
             cost = 1024, 
             gamma = 0.5, 
             epsilon = 0.1, 
             coef0 = 2)
svm.c.tst <- predict(svm.c, sprse.c.tst[ , -1:-2]) %>%
  cbind(sprse.c.tst[ , 2]) %>%
  data.frame() %>%
  rename(., pred = `.`, obs = V2)

defaultSummary(svm.c.tst)[2] # 0.858


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
defaultSummary(svm.abc.tst) # 0.672

#####
# Radial Basis Kernel -----------------------------------------------------
#     All ------------------------------------------------------------
rbf.all <- svm(x = sprse.trn.x, 
               y = sprse.trn.y,
               cost = 8, 
               epsilon = 0.5, 
               kernel = "radial", 
               gamma = 0.0005)

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
                     make.row.names = F) %>%
  mutate(residual = pred - obs)
defaultSummary(rbf.abc.tst)

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
ggsave("./graphs/svm/2017-07-06 poly all data points.png")

# Alpha-CD
ggplot(svm.a.tst, aes(x = obs, y = pred)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Experimental DelG, kJ/mol", y = "Predicted DelG, kJ/mol", 
       title = "Polynomial SVM - Alpha CD") + 
  coord_fixed() + 
  theme_bw()
ggsave("./graphs/svm/2017-07-06 poly alpha cd.png")

# Beta-CD
ggplot(svm.b.tst, aes(x = obs, y = pred)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Experimental DelG, kJ/mol", y = "Predicted DelG, kJ/mol", 
       title = "Polynomial SVM - Beta CD") + 
  coord_fixed() + 
  theme_bw()
ggsave("./graphs/svm/2017-07-06 poly beta cd.png")

# Gamma-CD
ggplot(svm.c.tst, aes(x = obs, y = pred)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Experimental DelG, kJ/mol", y = "Predicted DelG, kJ/mol", 
       title = "Polynomial SVM - Gamma CD") + 
  coord_fixed() + 
  theme_bw()
ggsave("./graphs/svm/2017-07-06 poly gamma cd.png")

# All-CDs
ggplot(svm.abc.tst, aes(x = obs, y = pred, 
                        group = cd.type, 
                        color = cd.type)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Experimental DelG, kJ/mol", y = "Predicted DelG, kJ/mol", 
       title = "Polynomail SVM - Compiled CD Types", 
       color = "Cyclodextrin Type") + 
  coord_fixed() + 
  theme_bw()
ggsave("./graphs/svm/2017-07-11 poly compiled cds.png")

#    Residuals
df.guest <- readRDS("./DelG.df.RDS") %>%
  select(., guest, DelG) %>%
  rename(obs = DelG)
svm.abc.labels <- inner_join(df.guest, svm.abc.tst, by = "obs")
ggplot(svm.abc.tst, aes(x = obs, y = residual)) + 
  geom_point() + 
  geom_hline(yintercept = 0) + 
  labs(x = "Experimental DelG, kJ/mol", y = "Residual, kJ/mol", 
       title = "Polynomial SVM - Residuals of Compiled CD Types") + 
  coord_fixed() + 
  theme_bw()
ggsave("./graphs/svm/2017-07-11 poly compiled cds resid.png")
ggplot(svm.abc.labels, aes(x = obs, y = residual, label = guest)) + 
  geom_point() + 
  geom_text(check_overlap = T, hjust = 0, nudge_x = 0.5) + 
  geom_hline(yintercept = 0) + 
  labs(x = "Experimental DelG, kJ/mol", y = "Residual, kJ/mol", 
       title = "Polynomial SVM - Residuals of Compiled CD Types") + 
  coord_fixed() + 
  theme_bw()

#         Training Set ----------------------------------------------------
# Very clearly overfitted
ggplot(svm.all.trn, aes(x = obs, y = pred)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Experimental DelG, kJ/mol", y = "Predicted DelG, kJ/mol", 
       title = "Polynomial SVM - All Data Points, Training") + 
  coord_fixed() + 
  theme_bw()
ggsave("./graphs/svm/2017-07-07 poly all data points trn.png")

#     RBF -----------------------------------------------------------------
ggplot(rbf.all.tst, aes(x = obs, y = pred)) + 
  geom_point(shape = 10, size = 2) + 
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Experimental DelG, kJ/mol", y = "Predicted DelG, kJ/mol", title = "Radial SVM - All Data Points") + 
  coord_fixed() + 
  theme_bw()
ggsave("./graphs/svm/2017-07-11 rbf all data points.png")

ggplot(rbf.a.tst, aes(x = obs, y = pred)) + 
  geom_point(shape = 10, size = 2) + 
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Experimental DelG, kJ/mol", y = "Predicted DelG, kJ/mol", 
       title = "Radial SVM - Alpha CD") + 
  coord_fixed() + 
  theme_bw()
ggsave("./graphs/svm/2017-07-11 rbf alpha cd.png")

ggplot(rbf.b.tst, aes(x = obs, y = pred)) + 
  geom_point(shape = 10, size = 2) + 
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Experimental DelG, kJ/mol", y = "Predicted DelG, kJ/mol", title = "Radial SVM - Beta CD") + 
  coord_fixed() + 
  theme_bw()
ggsave("./graphs/svm/2017-07-11 rbf beta cd.png")

ggplot(rbf.c.tst, aes(x = obs, y = pred)) + 
  geom_point(shape = 10, size = 2) + 
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Experimental DelG, kJ/mol", y = "Predicted DelG, kJ/mol", title = "Radial SVM - Gamma CD") + 
  coord_fixed() + 
  theme_bw()
ggsave("./graphs/svm/2017-07-11 rbf gamma cd.png")

ggplot(rbf.abc.tst, aes(x = obs, y = pred, color = cd.type)) + 
  geom_point(shape = 10, size = 2) + 
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Experimental DelG, kJ/mol", y = "Predicted DelG, kJ/mol", 
       title = "Radial SVM - Compiled CDs", 
       color = "Cyclodextrin Type") + 
  coord_fixed() + 
  theme_bw()
ggsave("./graphs/svm/2017-07-11 rbf compiled cd.png")
defaultSummary(rbf.abc.tst)

ggplot(rbf.abc.tst, aes(x = obs, y = residual, color = cd.type)) + 
  geom_point(shape = 10, size = 2) + 
  geom_hline(yintercept = 0) + 
  coord_fixed() + 
  theme_bw()
ggsave("./graphs/svm/2017-07-11 rbf compiled cd resid.png")

# Suzuki ------------------------------------------------------------------

suz <- readRDS("./suzuki.clean.RDS") %>%
  select(-guest)
suz.sprse <- sparse.model.matrix(~., suz)
set.seed(1)
trn.ind <- sample(x = 1:nrow(suz.sprse), 
                  size = round(0.7 * nrow(suz.sprse)))
suz.sprse.trn <- suz.sprse[trn.ind, ]
suz.sprse.trn.x <- suz.sprse.trn[ , -1:-2]
suz.sprse.trn.y <- suz.sprse.trn[ , 2]
suz.sprse.tst <- suz.sprse[-trn.ind, ]
suz.sprse.tst.x <- suz.sprse.tst[ , -1:-2]
suz.sprse.tst.y <- suz.sprse.tst[ , 2]

suz.svm.all <- svm(x = suz.sprse.trn.x, 
               y = suz.sprse.trn.y, 
               coef0 = 2, 
               cost = 1024, 
               epsilon = 0.1, 
               kernel = "polynomial", 
               gamma = 0.5, 
               degree = 2)

suz.svm.all.tst <- predict(suz.svm.all, suz.sprse.tst.x) %>%
  cbind(suz.sprse.tst.y) %>%
  data.frame() %>%
  rename(., pred = `.`, obs = suz.sprse.tst.y)
suz.svm.all.trn <- predict(suz.svm.all, suz.sprse.trn.x) %>%
  cbind(suz.sprse.trn.y) %>%
  data.frame() %>%
  rename(., pred = `.`, obs = suz.sprse.trn.y)

defaultSummary(suz.svm.all.tst) # 2.84 # 0.680

ggplot(suz.svm.all.tst, aes(x = obs, y = pred)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Experimental DelG, kJ/mol", y = "Predicted DelG, kJ/mol", 
       title = "Polynomial suz.svm - All Data Points, Suzuki") + 
  coord_fixed() + 
  theme_bw()
