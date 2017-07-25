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
#####
# Model on Dataframe (not sparse) -----------------------------------------

trn.x <- trn[ , -1]
trn.y <- trn[ , 1]
tst.x <- tst[ , -1]
tst.y <- tst[ , 1]

svm.all <- svm(x = trn.x, 
               y = trn.y, 
               coef0 = 1.5, 
               cost = 4096, 
               epsilon = 1, 
               kernel = "polynomial", 
               gamma = 0.03, 
               degree = 2)

svm.all.tst <- predict(svm.all, tst.x) %>%
  cbind(tst.y) %>%
  data.frame() %>%
  rename(., pred = `.`, obs = tst.y)

svm.all.trn <- predict(svm.all, trn.x) %>%
  cbind(trn.y) %>%
  data.frame() %>%
  rename(., pred = `.`, obs = trn.y)

defaultSummary(svm.all.tst)
defaultSummary(svm.all.trn)

# Polynomial Kernel -------------------------------------------------------

#     Alpha-CD ------------------------------------------------------------

alpha <- df %>% filter(alpha > 0)
set.seed(10)
trn.ind <- sample(x = 1:nrow(alpha), size = round(0.7 * nrow(alpha)))
a.trn.x <- alpha[trn.ind, -1]
a.trn.y <- alpha[trn.ind, 1]

a.tst.x <- alpha[-trn.ind, -1]
a.tst.y <- alpha[-trn.ind, 1]

svm.alpha <- svm(x = a.trn.x, 
                 y = a.trn.y, 
                 coef0 = 2, 
                 cost = 2048, 
                 epsilon = 0.1, 
                 kernel = "polynomial", 
                 gamma = 0.03, 
                 degree = 2)

svm.a.tst <- predict(svm.alpha, a.tst.x) %>%
  cbind(a.tst.y) %>% data.frame() %>%
  rename(., pred = `.`, obs = a.tst.y)
svm.a.trn <- predict(svm.alpha, a.trn.x) %>%
  cbind(a.trn.y) %>% data.frame() %>%
  rename(., pred = `.`, obs = a.trn.y)
defaultSummary(svm.a.tst) # 0.607
defaultSummary(svm.a.trn) # 0.993

saveRDS(svm.alpha, "./models/svm/polysvm.alpha.RDS")

#     Beta-CD -------------------------------------------------------------

beta <- df %>% filter(beta > 0)
set.seed(10)
trn.ind <- sample(x = 1:nrow(beta), size = round(0.7 * nrow(beta)))
b.trn.x <- beta[trn.ind, -1]
b.trn.y <- beta[trn.ind, 1]

b.tst.x <- beta[-trn.ind, -1]
b.tst.y <- beta[-trn.ind, 1]

svm.beta <- svm(x = b.trn.x, 
                 y = b.trn.y, 
                 coef0 = 2, 
                 cost = 2048, 
                 epsilon = 0.1, 
                 kernel = "polynomial", 
                 gamma = 0.03, 
                 degree = 2)

svm.b.tst <- predict(svm.beta, b.tst.x) %>%
  cbind(b.tst.y) %>% data.frame() %>%
  rename(., pred = `.`, obs = b.tst.y)
svm.b.trn <- predict(svm.beta, b.trn.x) %>%
  cbind(b.trn.y) %>% data.frame() %>%
  rename(., pred = `.`, obs = b.trn.y)
defaultSummary(svm.b.tst) # 0.634
defaultSummary(svm.b.trn) # 0.998

saveRDS(svm.alpha, "./models/svm/polysvm.beta.RDS")

#     Gamma-CD ------------------------------------------------------------

gamma <- df %>% filter(gamma > 0)
set.seed(10)
trn.ind <- sample(x = 1:nrow(gamma), size = round(0.7 * nrow(gamma)))
c.trn.x <- gamma[trn.ind, -1]
c.trn.y <- gamma[trn.ind, 1]

c.tst.x <- gamma[-trn.ind, -1]
c.tst.y <- gamma[-trn.ind, 1]

svm.gamma <- svm(x = c.trn.x, 
                 y = c.trn.y, 
                 coef0 = 2, 
                 cost = 2048, 
                 epsilon = 0.1, 
                 kernel = "polynomial", 
                 gamma = 0.03, 
                 degree = 2)

svm.c.tst <- predict(svm.gamma, c.tst.x) %>%
  cbind(c.tst.y) %>% data.frame() %>%
  rename(., pred = `.`, obs = c.tst.y)
svm.c.trn <- predict(svm.gamma, c.trn.x) %>%
  cbind(c.trn.y) %>% data.frame() %>%
  rename(., pred = `.`, obs = c.trn.y)
defaultSummary(svm.c.tst) # 0.214
defaultSummary(svm.c.trn) # 0.99986

saveRDS(svm.gamma, "./models/svm/polysvm.gamma.RDS")

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

temp.a <- svm.a.trn %>% 
  mutate(cd.type = rep("alpha", length(svm.a.trn$pred)))
temp.b <- svm.b.trn %>% 
  mutate(cd.type = rep("beta", length(svm.b.trn$pred)))
temp.c <- svm.c.trn %>% 
  mutate(cd.type = rep("gamma", length(svm.c.trn$pred)))
svm.abc.trn <- rbind(temp.a, temp.b, temp.c, 
                     make.row.names = F) %>%
  mutate(residual = pred - obs)

defaultSummary(svm.abc.tst) # 0.617
defaultSummary(svm.abc.trn) # 0.997


# Graphs ------------------------------------------------------------------
dir.create("./graphs")
dir.create("./graphs/svm")

# ggplot_build(p)$data
ggplot(svm.a.tst, aes(x = pred, y = obs)) + 
  geom_point(color = "#F8766D") + 
  geom_abline(slope = 1, intercept = 0) + 
  theme_bw() + 
  coord_fixed() + 
  labs(title = "Polynomial SVM - Alpha-CD", 
       x = "Predicted DelG, kJ/mol", 
       y = "Experimental DelG, kJ/mol")
ggsave("./graphs/svm/2017-07-24 poly alpha on df.png")

ggplot(svm.b.tst, aes(x = pred, y = obs)) + 
  geom_point(color = "#00BA38") + 
  geom_abline(slope = 1, intercept = 0) + 
  theme_bw() + 
  coord_fixed() + 
  labs(title = "Polynomial SVM - Beta-CD", 
       x = "Predicted DelG, kJ/mol", 
       y = "Experimental DelG, kJ/mol")
ggsave("./graphs/svm/2017-07-24 poly beta on df.png")

ggplot(svm.c.tst, aes(x = pred, y = obs)) + 
  geom_point(color = "#619CFF") + 
  geom_abline(slope = 1, intercept = 0) + 
  theme_bw() + 
  coord_fixed() + 
  labs(title = "Polynomial SVM - Gamma-CD", 
       x = "Predicted DelG, kJ/mol", 
       y = "Experimental DelG, kJ/mol")
ggsave("./graphs/svm/2017-07-24 poly gamma on df.png")

ggplot(svm.abc.tst, aes(x = pred, y = obs, color = cd.type)) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0) + 
  theme_bw() + 
  coord_fixed() + 
  labs(title = "Polynomial SVM - Compiled CD Types", 
       x = "Predicted DelG, kJ/mol", 
       y = "Experimental DelG, kJ/mol", 
       color = "Cyclodextrin Type")
ggsave("./graphs/svm/2017-07-24 poly compiled on df.png")

ggplot(svm.all.tst, aes(x = pred, y = obs)) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0) + 
  theme_bw() + 
  coord_fixed() + 
  labs(title = "Polynomial SVM", 
       x = "Predicted DelG, kJ/mol", 
       y = "Experimental DelG, kJ/mol")
ggsave("./graphs/svm/2017-07-24 poly all on df.png")

# On training data
ggplot(svm.abc.trn, aes(x = pred, y = obs, color = cd.type)) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0) + 
  theme_bw() + 
  coord_fixed() + 
  labs(title = "Polynomial SVM - Training Data", 
       x = "Predicted DelG, kJ/mol", 
       y = "Experimental DelG, kJ/mol")
ggsave("./graphs/svm/2017-07-24 poly all on df trn.png")

#####

# External  Validation ----------------------------------------------------

#####

ext.val <- readRDS("./external validation set new.RDS") %>%
  select(-guest:-host) %>%
  select(-data.source) 

ext.val.a <- ext.val %>% filter(alpha > 0) 
ext.val.b <- ext.val %>% filter(beta > 0)
ext.val.c <- ext.val %>% filter(gamma > 0) 

ev.a <-  predict(svm.alpha, ext.val.a[ , -1]) %>%
  cbind(ext.val.a[ , 1]) %>% data.frame() %>%
  rename(., pred = `.`, obs = V2)

ev.b <-  predict(svm.beta, ext.val.b[ , -1]) %>%
  cbind(ext.val.b[ , 1]) %>% data.frame() %>%
  rename(., pred = `.`, obs = V2)

ev.c <-  predict(svm.gamma, ext.val.c[ , -1]) %>%
  cbind(ext.val.c[ , 1]) %>% data.frame() %>%
  rename(., pred = `.`, obs = V2)

temp.a <- ev.a %>% mutate(cd.type = "alpha")
temp.b <- ev.b %>% mutate(cd.type = "beta")
temp.c <- ev.c %>% mutate(cd.type = "gamma")
# temp.a <- temp.a[-2, ]
ev.abc <- rbind(temp.a, temp.b, temp.c) %>%
  mutate(resid = pred - obs)

defaultSummary(ev.abc) # 0.358 normally, 0.595 without outlier
ggplot(ev.abc, aes(x = pred, y = obs, color = cd.type)) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0) + 
  coord_fixed() + 
  theme_bw() + 
  labs(title = "Polynomial SVM - External Validation", 
       x = "Predicted DelG, kJ/mol", y = "Experimental DelG, kJ/mol", 
       color = "Cyclodextrin Type")
ggsave("./graphs/svm/2017-07-25 polysvm extval.png")
