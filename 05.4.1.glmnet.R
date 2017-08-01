# Libraries ---------------------------------------------------------------

library(caret)
library(glmnet)
library(tidyverse)

# Loading Data ------------------------------------------------------------

# setwd("~/SREP LAB/qsar")
dir.create("./graphs/glmnet")

df.raw <- readRDS("./padel.pp.new.RDS")
df <- df.raw %>% select(-guest:-host) %>%
  select(-data.source)
# mat <- sparse.model.matrix(~., df)
mat <- as.matrix(df)

set.seed(25)
trn.ind <- sample(x = 1:nrow(mat), 
                  size = round(0.7 * nrow(mat)))
trn.x <- mat[trn.ind, -1]
trn.y <- mat[trn.ind, 1]
tst.x <- mat[-trn.ind, -1]
tst.y <- mat[-trn.ind, 1]

# All Predictors ----------------------------------------------------------

glm.all <- glmnet(x = trn.x, y = trn.y, 
                  dfmax = 32, 
                  alpha = 1, 
                  family = "mgaussian")
glm.df <- predict.glmnet(glm.all, tst.x, 
                         s = tail(glm.all$lambda, n = 1)) %>%
  cbind(tst.y) %>% data.frame() %>%
  dplyr::rename(pred = X1, obs = tst.y) %>%
  mutate(resid = pred - obs)

defaultSummary(glm.df) # 0.505


# Separating Cyclodextrins ------------------------------------------------

#     Alpha-CD ------------------------------------------------------------

alpha <- df %>% filter(alpha > 0) %>%
  as.matrix()
set.seed(25)
trn.ind <- sample(x = 1:nrow(alpha), 
                  size = round(0.7 * nrow(alpha)))
a.trn.x <- alpha[trn.ind, -1]
a.trn.y <- alpha[trn.ind, 1]
a.tst.x <- alpha[-trn.ind, -1]
a.tst.y <- alpha[-trn.ind, 1]

glm.alpha <- glmnet(x = a.trn.x, y = a.trn.y, 
                    dfmax = 32, 
                    alpha = 1, 
                    family = "mgaussian")
glm.alpha.df <- predict.glmnet(glm.alpha, a.tst.x, 
                         s = tail(glm.alpha$lambda, n = 1)) %>%
  cbind(a.tst.y) %>% data.frame() %>%
  dplyr::rename(., pred = X1, obs = a.tst.y)  %>%
  mutate(resid = pred - obs)

defaultSummary(glm.alpha.df) # 0.473

#     Beta-CD -------------------------------------------------------------

beta <- df %>% filter(beta > 0) %>%
  as.matrix()
set.seed(25)
trn.ind <- sample(x = 1:nrow(beta), 
                  size = round(0.7 * nrow(beta)))
b.trn.x <- beta[trn.ind, -1]
b.trn.y <- beta[trn.ind, 1]
b.tst.x <- beta[-trn.ind, -1]
b.tst.y <- beta[-trn.ind, 1]

glm.beta <- glmnet(x = b.trn.x, y = b.trn.y, 
                    dfmax = 32, 
                    alpha = 1, 
                    family = "mgaussian")
glm.beta.df <- predict.glmnet(glm.beta, b.tst.x, 
                               s = tail(glm.beta$lambda, n = 1)) %>%
  cbind(b.tst.y) %>% data.frame() %>%
  dplyr::rename(pred = X1, obs = b.tst.y) %>%
  mutate(resid = pred - obs)

defaultSummary(glm.beta.df) # 0.718

#     Gamma-CD ------------------------------------------------------------

gamma <- df %>% filter(gamma > 0) %>%
  as.matrix()
set.seed(25)
trn.ind <- sample(x = 1:nrow(gamma), 
                  size = round(0.7 * nrow(gamma)))
c.trn.x <- gamma[trn.ind, -1]
c.trn.y <- gamma[trn.ind, 1]
c.tst.x <- gamma[-trn.ind, -1]
c.tst.y <- gamma[-trn.ind, 1]

glm.gamma <- glmnet(x = c.trn.x, y = c.trn.y, 
                   dfmax = 32, 
                   alpha = 1, 
                   family = "mgaussian")
glm.gamma.df <- predict.glmnet(glm.gamma, c.tst.x, 
                              s = tail(glm.gamma$lambda, n = 1)) %>%
  cbind(c.tst.y) %>% data.frame() %>%
  dplyr::rename(pred = X1, obs = c.tst.y) %>%
  mutate(resid = pred - obs)

defaultSummary(glm.gamma.df) # 0.345

#     Compiled CDs --------------------------------------------------------

temp.a <- glm.alpha.df %>% 
  mutate(cd.type = rep("alpha", length(glm.alpha.df$pred)))
temp.b <- glm.beta.df %>% 
  mutate(cd.type = rep("beta", length(glm.beta.df$pred)))
temp.c <- glm.gamma.df %>% 
  mutate(cd.type = rep("gamma", length(glm.gamma.df$pred)))
glm.abc.tst <- rbind(temp.a, temp.b, temp.c, 
                     make.row.names = F) %>%
  mutate(residual = pred - obs)
defaultSummary(glm.abc.tst)

# Graphs ------------------------------------------------------------------

# all datapoints
ggplot(glm.df, aes(x = obs, y = pred)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) + 
  theme_bw() + 
  coord_fixed(xlim = c(-45, 5), ylim = c(-45, 5)) + 
  labs(x = "Observed DelG, kJ/mol", y = "Predicted DelG, kJ/mol", 
       title = "GLMnet on All Datapoints")
ggsave("./graphs/glmnet/2017-07-28 glmnet all datapoints.png")

ggplot(glm.alpha.df, aes(x = obs, y = pred)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) + 
  coord_fixed(xlim = c(-45, 5), ylim = c(-45, 5)) + 
  theme_bw() + 
  labs(x = "Observed DelG, kJ/mol", y = "Predicted DelG, kJ/mol", 
       title = "GLMnet on Alpha-CD")
ggsave("./graphs/glmnet/2017-07-28 glmnet alpha-cd.png")

ggplot(glm.beta.df, aes(x = obs, y = pred)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) + 
  coord_fixed(xlim = c(-45, 5), ylim = c(-45, 5)) + 
  theme_bw() + 
  labs(x = "Observed DelG, kJ/mol", y = "Predicted DelG, kJ/mol", 
       title = "GLMnet on Beta-CD")
ggsave("./graphs/glmnet/2017-07-28 glmnet beta-cd.png")

ggplot(glm.gamma.df, aes(x = obs, y = pred)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) + 
  # coord_fixed(xlim = c(-45, 5), ylim = c(-45, 5)) + 
  theme_bw() + 
  labs(x = "Observed DelG, kJ/mol", y = "Predicted DelG, kJ/mol", 
       title = "GLMnet on Gamma-CD")
ggsave("./graphs/glmnet/2017-07-28 glmnet gamma-cd.png")

ggplot(glm.abc.tst, aes(x = obs, y = pred, color = cd.type)) + 
  geom_point() + 
  coord_fixed(xlim = c(-45, 5), ylim = c(-45, 5)) + 
  geom_abline(intercept = 0, slope = 1) + 
  theme_bw() + 
  labs(x = "Observed DelG, kJ/mol", y = "Predicted DelG, kJ/mol", 
       title = "GLMnet on Separate CDs", color = "Cyclodextrin")
ggsave("./graphs/glmnet/2017-07-28 glmnet compiled.png")

ggplot(glm.abc.tst, aes(x = obs, y = residual, color = cd.type)) + 
  geom_point() + 
  coord_fixed(xlim = c(-45, 5), ylim = c(-45, 5)) + 
  geom_hline(yintercept = 0) + 
  theme_bw() + 
  labs(y = "Residuals, kJ/mol", x = "Expeirmental DelG, kJ/mol", 
       title = "GLMnet Residuals on Separate CDs", 
       color = "Cyclodextrin")
ggsave("./graphs/glmnet/2017-07-28 glmnet compiled resid.png")

#####

# Modeling on Regular Matrix  ---------------------------------------------

#####

df.raw <- readRDS("./padel.pp.new.RDS")
df <- df.raw %>% select(-guest:-host) %>%
  select(-data.source)
mat <- as.matrix(df)

set.seed(144)
trn.ind <- sample(x = 1:nrow(mat), 
                  size = round(0.7 * nrow(mat)))
trn.x <- mat[trn.ind, -1]
trn.y <- mat[trn.ind, 1]
tst.x <- mat[-trn.ind, -1]
tst.y <- mat[-trn.ind, 1]

dir.create("./models/glmnet")

# All Predictors ----------------------------------------------------------

glm.all <- glmnet(x = trn.x, y = trn.y, 
                  dfmax = 32, 
                  alpha = 1, 
                  family = "mgaussian")
glm.df <- predict.glmnet(glm.all, tst.x, 
                         s = tail(glm.all$lambda, n = 1)) %>%
  cbind(tst.y) %>% data.frame() %>%
  dplyr::rename(pred = X1, obs = tst.y)
ggplot(glm.df, aes(x = obs, y = pred)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) + 
  coord_fixed(xlim = c(-45, 5), ylim = c(-45, 5)) + 
  theme_bw()
defaultSummary(glm.df) # 0.635

# Separating Cyclodextrins ------------------------------------------------

#     Alpha-CD ------------------------------------------------------------

alpha <- df %>% filter(alpha > 0) %>%
  as.matrix()
set.seed(35)
trn.ind <- sample(x = 1:nrow(alpha), 
                  size = round(0.7 * nrow(alpha)))
a.trn.x <- alpha[trn.ind, -1]
a.trn.y <- alpha[trn.ind, 1]
a.tst.x <- alpha[-trn.ind, -1]
a.tst.y <- alpha[-trn.ind, 1]

glm.alpha <- glmnet(x = a.trn.x, y = a.trn.y, 
                    dfmax = 32, 
                    alpha = 1, 
                    family = "mgaussian")
glm.alpha.df <- predict.glmnet(glm.alpha, a.tst.x, 
                               s = tail(glm.alpha$lambda, n = 1)) %>%
  cbind(a.tst.y) %>% data.frame() %>%
  dplyr::rename(pred = X1, obs = a.tst.y)
glm.alpha.trn <- predict.glmnet(glm.alpha, a.trn.x, 
                               s = tail(glm.alpha$lambda, n = 1)) %>%
  cbind(a.trn.y) %>% data.frame() %>%
  dplyr::rename(pred = X1, obs = a.trn.y)

defaultSummary(glm.alpha.df) # 0.629
saveRDS(glm.alpha, "./models/glmnet/glm.alpha.RDS")

#     Beta-CD -------------------------------------------------------------

beta <- df %>% filter(beta > 0) %>%
  as.matrix()
set.seed(25)
trn.ind <- sample(x = 1:nrow(beta), 
                  size = round(0.7 * nrow(beta)))
b.trn.x <- beta[trn.ind, -1]
b.trn.y <- beta[trn.ind, 1]
b.tst.x <- beta[-trn.ind, -1]
b.tst.y <- beta[-trn.ind, 1]

glm.beta <- glmnet(x = b.trn.x, y = b.trn.y, 
                   dfmax = 32, 
                   alpha = 1, 
                   family = "mgaussian")
glm.beta.df <- predict.glmnet(glm.beta, b.tst.x, 
                              s = tail(glm.beta$lambda, n = 1)) %>%
  cbind(b.tst.y) %>% data.frame() %>%
  dplyr::rename(pred = X1, obs = b.tst.y)
glm.beta.trn <- predict.glmnet(glm.beta, b.trn.x, 
                              s = tail(glm.beta$lambda, n = 1)) %>%
  cbind(b.trn.y) %>% data.frame() %>%
  dplyr::rename(pred = X1, obs = b.trn.y)

defaultSummary(glm.beta.df) # 0.718
saveRDS(glm.beta, "./models/glmnet/glm.beta.RDS")

#     Gamma-CD ------------------------------------------------------------

gamma <- df %>% filter(gamma > 0) %>%
  as.matrix()
set.seed(12)
trn.ind <- sample(x = 1:nrow(gamma), 
                  size = round(0.7 * nrow(gamma)))
c.trn.x <- gamma[trn.ind, -1]
c.trn.y <- gamma[trn.ind, 1]
c.tst.x <- gamma[-trn.ind, -1]
c.tst.y <- gamma[-trn.ind, 1]

glm.gamma <- glmnet(x = c.trn.x, y = c.trn.y, 
                    dfmax = 32, 
                    alpha = 1, 
                    family = "mgaussian")
glm.gamma.df <- predict.glmnet(glm.gamma, c.tst.x, 
                               s = tail(glm.gamma$lambda, n = 1)) %>%
  cbind(c.tst.y) %>% data.frame() %>%
  dplyr::rename(pred = X1, obs = c.tst.y)
glm.gamma.trn <- predict.glmnet(glm.gamma, c.trn.x, 
                               s = tail(glm.gamma$lambda, n = 1)) %>%
  cbind(c.trn.y) %>% data.frame() %>%
  dplyr::rename(pred = X1, obs = c.trn.y)

defaultSummary(glm.gamma.df) # 0.304
saveRDS(glm.gamma, "./models/glmnet/glm.gamma.RDS")

#     Compiled CDs --------------------------------------------------------

temp.a <- glm.alpha.df %>% 
  mutate(cd.type = rep("alpha", length(glm.alpha.df$pred)))
temp.b <- glm.beta.df %>% 
  mutate(cd.type = rep("beta", length(glm.beta.df$pred)))
temp.c <- glm.gamma.df %>% 
  mutate(cd.type = rep("gamma", length(glm.gamma.df$pred)))
glm.abc.tst <- rbind(temp.a, temp.b, temp.c, 
                     make.row.names = F) %>%
  mutate(residual = pred - obs)
defaultSummary(glm.abc.tst)

temp.a <- glm.alpha.trn %>% 
  mutate(cd.type = rep("alpha", length(glm.alpha.trn$pred)))
temp.b <- glm.beta.trn %>% 
  mutate(cd.type = rep("beta", length(glm.beta.trn$pred)))
temp.c <- glm.gamma.trn %>% 
  mutate(cd.type = rep("gamma", length(glm.gamma.trn$pred)))
glm.abc.trn <- rbind(temp.a, temp.b, temp.c, 
                     make.row.names = F) %>%
  mutate(residual = pred - obs)
defaultSummary(glm.abc.trn)

# Graphs ------------------------------------------------------------------

# all datapoints
ggplot(glm.df, aes(x = obs, y = pred)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) + 
  theme_bw() + 
  coord_fixed(xlim = c(-45, 5), ylim = c(-45, 5)) + 
  labs(x = "Observed DelG, kJ/mol", y = "Predicted DelG, kJ/mol", 
       title = "GLMnet on All Datapoints")
ggsave("./graphs/glmnet/2017-07-28 glmnet all datapoints df.png")

ggplot(glm.alpha.df, aes(x = obs, y = pred)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) + 
  coord_fixed(xlim = c(-45, 5), ylim = c(-45, 5)) + 
  theme_bw() + 
  labs(x = "Observed DelG, kJ/mol", y = "Predicted DelG, kJ/mol", 
       title = "GLMnet on Alpha-CD")
ggsave("./graphs/glmnet/2017-07-28 glmnet alpha-cd df.png")

ggplot(glm.beta.df, aes(x = obs, y = pred)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) + 
  coord_fixed(xlim = c(-45, 5), ylim = c(-45, 5)) + 
  theme_bw() + 
  labs(x = "Observed DelG, kJ/mol", y = "Predicted DelG, kJ/mol", 
       title = "GLMnet on Beta-CD")
ggsave("./graphs/glmnet/2017-07-28 glmnet beta-cd df.png")

ggplot(glm.gamma.df, aes(x = obs, y = pred)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) + 
  # coord_fixed(xlim = c(-45, 5), ylim = c(-45, 5)) + 
  theme_bw() + 
  labs(x = "Observed DelG, kJ/mol", y = "Predicted DelG, kJ/mol", 
       title = "GLMnet on Gamma-CD")
ggsave("./graphs/glmnet/2017-07-28 glmnet gamma-cd df.png")

ggplot(glm.abc.tst, aes(x = obs, y = pred, color = cd.type)) + 
  geom_point() + 
  coord_fixed(xlim = c(-45, 5), ylim = c(-45, 5)) + 
  geom_abline(intercept = 0, slope = 1) + 
  theme_bw() + 
  labs(x = "Observed DelG, kJ/mol", y = "Predicted DelG, kJ/mol", 
       title = "GLMnet on Separate CDs", color = "Cyclodextrin")
ggsave("./graphs/glmnet/2017-07-28 glmnet compiled df.png")

ggplot(glm.abc.tst, aes(x = obs, y = residual, color = cd.type)) + 
  geom_point() + 
  geom_hline(yintercept = 0) + 
  theme_bw() + 
  labs(y = "Residuals, kJ/mol", x = "Observed DelG, kJ/mol", 
       title = "GLMnet Residuals on Separate CDs", 
       color = "Cyclodextrin")
ggsave("./graphs/glmnet/2017-07-28 glmnet compiled resid df.png")

ggplot(glm.abc.trn, aes(x = obs, y = pred, color = cd.type)) + 
  geom_point() + 
  coord_fixed(xlim = c(-45, 5), ylim = c(-45, 5)) + 
  geom_abline(intercept = 0, slope = 1) + 
  theme_bw() + 
  labs(x = "Observed DelG, kJ/mol", y = "Predicted DelG, kJ/mol", 
       title = "GLMnet on Separate CDs, Training", color = "Cyclodextrin")
ggsave("./graphs/glmnet/2017-07-28 glmnet compiled trn.png")


# External Validation -----------------------------------------------------

ext.val <- readRDS("./external validation set new.RDS") %>%
  select(-guest:-host) %>%
  select(-data.source) 

ext.val.a <- ext.val %>% filter(alpha > 0) %>% as.matrix()
ext.val.b <- ext.val %>% filter(beta > 0) %>% as.matrix()
ext.val.c <- ext.val %>% filter(gamma > 0) 

ev.a <-  predict.glmnet(glm.alpha, ext.val.a[ , -1], 
                        s = tail(glm.alpha$lambda, n = 1)) %>%
  cbind(ext.val.a[ , 1]) %>% data.frame() %>%
  dplyr::rename(pred = X1, obs = V2)

ev.b <-  predict.glmnet(glm.beta, ext.val.b[ , -1], 
                        s = tail(glm.beta$lambda, n = 1)) %>%
  cbind(ext.val.b[ , 1]) %>% data.frame() %>%
  dplyr::rename(pred = X1, obs = V2)

ext.val.temp <- ext.val.c
ext.val.temp[ , 1] <- NULL
ext.val.temp <- as.matrix(ext.val.temp)

ev.c <-  predict.glmnet(glm.gamma, ext.val.temp, 
                        s = tail(glm.gamma$lambda, n = 1)) %>%
  cbind(ext.val.c[ , 1]) %>% data.frame() %>%
  dplyr::rename(pred = X1, obs = V2)

temp.a <- ev.a %>% mutate(cd.type = "alpha")
temp.b <- ev.b %>% mutate(cd.type = "beta")
temp.c <- ev.c %>% mutate(cd.type = "gamma")
# temp.a <- temp.a[-2, ]
ev.abc <- rbind(temp.a, temp.b, temp.c)

defaultSummary(ev.abc) # 0.429 normally, 0.649 without outlier
ggplot(ev.abc, aes(x = obs, y = pred, color = cd.type)) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0) + 
  coord_fixed(xlim = c(-45, 5), ylim = c(-45, 5)) + 
  theme_bw() + 
  labs(x = "Observed DelG, kJ/mol", y = "Predicted DelG, kJ/mol", 
       title = "GLMnet External Validation", color = "Cyclodextrin")
ggsave("./graphs/glmnet/2017-07-28 glm extval.png")
