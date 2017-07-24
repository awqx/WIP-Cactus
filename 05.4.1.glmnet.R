# Libraries ---------------------------------------------------------------

library(caret)
library(ggplot2)
library(glmnet)

# Loading Data ------------------------------------------------------------

setwd("~/SREP LAB/qsar")
dir.create("./graphs/glmnet")

df.raw <- readRDS("./padel.pp.new.RDS")
df <- df.raw %>% select(-guest:-host) %>%
  select(-data.source)
mat <- sparse.model.matrix(~., df)

set.seed(25)
trn.ind <- sample(x = 1:nrow(mat), 
                  size = round(0.7 * nrow(mat)))
trn.x <- mat[trn.ind, -1:-2]
trn.y <- mat[trn.ind, 2]
tst.x <- mat[-trn.ind, -1:-2]
tst.y <- mat[-trn.ind, 2]

# All Predictors ----------------------------------------------------------

glm.all <- glmnet(x = trn.x, y = trn.y, 
                  dfmax = 32, 
                  alpha = 1, 
                  family = "mgaussian")
glm.df <- predict.glmnet(glm.all, tst.x, 
                         s = tail(glm.all$lambda, n = 1)) %>%
  cbind(tst.y) %>% data.frame() %>%
  rename(pred = X1, obs = tst.y)
ggplot(glm.df, aes(x = pred, y = obs)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) + 
  coord_fixed() + 
  theme_bw()
defaultSummary(glm.df)


# Separating Cyclodextrins ------------------------------------------------

#     Alpha-CD ------------------------------------------------------------

alpha <- df %>% filter(alpha > 0) %>%
  sparse.model.matrix(~., .)
set.seed(25)
trn.ind <- sample(x = 1:nrow(alpha), 
                  size = round(0.7 * nrow(alpha)))
a.trn.x <- alpha[trn.ind, -1:-2]
a.trn.y <- alpha[trn.ind, 2]
a.tst.x <- alpha[-trn.ind, -1:-2]
a.tst.y <- alpha[-trn.ind, 2]

glm.alpha <- glmnet(x = a.trn.x, y = a.trn.y, 
                    dfmax = 32, 
                    alpha = 1, 
                    family = "mgaussian")
glm.alpha.df <- predict.glmnet(glm.alpha, a.tst.x, 
                         s = tail(glm.alpha$lambda, n = 1)) %>%
  cbind(a.tst.y) %>% data.frame() %>%
  rename(pred = X1, obs = a.tst.y)

defaultSummary(glm.alpha.df)

#     Beta-CD -------------------------------------------------------------

beta <- df %>% filter(beta > 0) %>%
  sparse.model.matrix(~., .)
set.seed(25)
trn.ind <- sample(x = 1:nrow(beta), 
                  size = round(0.7 * nrow(beta)))
b.trn.x <- beta[trn.ind, -1:-2]
b.trn.y <- beta[trn.ind, 2]
b.tst.x <- beta[-trn.ind, -1:-2]
b.tst.y <- beta[-trn.ind, 2]

glm.beta <- glmnet(x = b.trn.x, y = b.trn.y, 
                    dfmax = 32, 
                    alpha = 1, 
                    family = "mgaussian")
glm.beta.df <- predict.glmnet(glm.beta, b.tst.x, 
                               s = tail(glm.beta$lambda, n = 1)) %>%
  cbind(b.tst.y) %>% data.frame() %>%
  rename(pred = X1, obs = b.tst.y)

defaultSummary(glm.beta.df)

#     Gamma-CD ------------------------------------------------------------

gamma <- df %>% filter(gamma > 0) %>%
  sparse.model.matrix(~., .)
set.seed(12)
trn.ind <- sample(x = 1:nrow(gamma), 
                  size = round(0.7 * nrow(gamma)))
c.trn.x <- gamma[trn.ind, -1:-2]
c.trn.y <- gamma[trn.ind, 2]
c.tst.x <- gamma[-trn.ind, -1:-2]
c.tst.y <- gamma[-trn.ind, 2]

glm.gamma <- glmnet(x = c.trn.x, y = c.trn.y, 
                   dfmax = 32, 
                   alpha = 1, 
                   family = "mgaussian")
glm.gamma.df <- predict.glmnet(glm.gamma, c.tst.x, 
                              s = tail(glm.gamma$lambda, n = 1)) %>%
  cbind(c.tst.y) %>% data.frame() %>%
  rename(pred = X1, obs = c.tst.y)

defaultSummary(glm.gamma.df)

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
ggplot(glm.df, aes(x = pred, y = obs)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) + 
  theme_bw() + 
  coord_fixed() + 
  labs(x = "Predicted DelG, kJ/mol", y = "Expeirmental DelG, kJ/mol", 
       title = "GLMnet on All Datapoints")
ggsave("./graphs/glmnet/2017-07-21 glmnet all datapoints.png")

ggplot(glm.alpha.df, aes(x = pred, y = obs)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) + 
  coord_fixed() + 
  theme_bw() + 
  labs(x = "Predicted DelG, kJ/mol", y = "Expeirmental DelG, kJ/mol", 
       title = "GLMnet on Alpha-CD")
ggsave("./graphs/glmnet/2017-07-21 glmnet alpha-cd.png")

ggplot(glm.beta.df, aes(x = pred, y = obs)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) + 
  coord_fixed() + 
  theme_bw() + 
  labs(x = "Predicted DelG, kJ/mol", y = "Expeirmental DelG, kJ/mol", 
       title = "GLMnet on Beta-CD")
ggsave("./graphs/glmnet/2017-07-21 glmnet beta-cd.png")

ggplot(glm.gamma.df, aes(x = pred, y = obs)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) + 
  # coord_fixed() + 
  theme_bw() + 
  labs(x = "Predicted DelG, kJ/mol", y = "Expeirmental DelG, kJ/mol", 
       title = "GLMnet on Gamma-CD")
ggsave("./graphs/glmnet/2017-07-21 glmnet gamma-cd.png")

ggplot(glm.abc.tst, aes(x = pred, y = obs, color = cd.type)) + 
  geom_point() + 
  coord_fixed() + 
  geom_abline(intercept = 0, slope = 1) + 
  theme_bw() + 
  labs(x = "Predicted DelG, kJ/mol", y = "Expeirmental DelG, kJ/mol", 
       title = "GLMnet on Separate CDs", color = "Cyclodextrin")
ggsave("./graphs/glmnet/2017-07-21 glmnet compiled.png")

ggplot(glm.abc.tst, aes(x = obs, y = residual, color = cd.type)) + 
  geom_point() + 
  coord_fixed() + 
  geom_hline(yintercept = 0) + 
  theme_bw() + 
  labs(y = "Residuals, kJ/mol", x = "Expeirmental DelG, kJ/mol", 
       title = "GLMnet Residuals on Separate CDs", 
       color = "Cyclodextrin")
ggsave("./graphs/glmnet/2017-07-21 glmnet compiled resid.png")

#####

# Modeling on Regular Matrix  ---------------------------------------------

#####

df.raw <- readRDS("./padel.pp.new.RDS")
df <- df.raw %>% select(-guest:-host) %>%
  select(-data.source)
mat <- as.matrix(df)

set.seed(25)
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
  rename(pred = X1, obs = tst.y)
ggplot(glm.df, aes(x = pred, y = obs)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) + 
  coord_fixed() + 
  theme_bw()
defaultSummary(glm.df)

# Separating Cyclodextrins ------------------------------------------------

#     Alpha-CD ------------------------------------------------------------

alpha <- df %>% filter(alpha > 0) %>%
  sparse.model.matrix(~., .)
set.seed(25)
trn.ind <- sample(x = 1:nrow(alpha), 
                  size = round(0.7 * nrow(alpha)))
a.trn.x <- alpha[trn.ind, -1:-2]
a.trn.y <- alpha[trn.ind, 2]
a.tst.x <- alpha[-trn.ind, -1:-2]
a.tst.y <- alpha[-trn.ind, 2]

glm.alpha <- glmnet(x = a.trn.x, y = a.trn.y, 
                    dfmax = 32, 
                    alpha = 1, 
                    family = "mgaussian")
glm.alpha.df <- predict.glmnet(glm.alpha, a.tst.x, 
                               s = tail(glm.alpha$lambda, n = 1)) %>%
  cbind(a.tst.y) %>% data.frame() %>%
  rename(pred = X1, obs = a.tst.y)

defaultSummary(glm.alpha.df)
saveRDS(glm.alpha, "./models/glmnet/glm alpha.RDS")

#     Beta-CD -------------------------------------------------------------

beta <- df %>% filter(beta > 0) %>%
  sparse.model.matrix(~., .)
set.seed(25)
trn.ind <- sample(x = 1:nrow(beta), 
                  size = round(0.7 * nrow(beta)))
b.trn.x <- beta[trn.ind, -1:-2]
b.trn.y <- beta[trn.ind, 2]
b.tst.x <- beta[-trn.ind, -1:-2]
b.tst.y <- beta[-trn.ind, 2]

glm.beta <- glmnet(x = b.trn.x, y = b.trn.y, 
                   dfmax = 32, 
                   alpha = 1, 
                   family = "mgaussian")
glm.beta.df <- predict.glmnet(glm.beta, b.tst.x, 
                              s = tail(glm.beta$lambda, n = 1)) %>%
  cbind(b.tst.y) %>% data.frame() %>%
  rename(pred = X1, obs = b.tst.y)

defaultSummary(glm.beta.df)
saveRDS(glm.beta, "./models/glmnet/glm beta.RDS")

#     Gamma-CD ------------------------------------------------------------

gamma <- df %>% filter(gamma > 0) %>%
  sparse.model.matrix(~., .)
set.seed(12)
trn.ind <- sample(x = 1:nrow(gamma), 
                  size = round(0.7 * nrow(gamma)))
c.trn.x <- gamma[trn.ind, -1:-2]
c.trn.y <- gamma[trn.ind, 2]
c.tst.x <- gamma[-trn.ind, -1:-2]
c.tst.y <- gamma[-trn.ind, 2]

glm.gamma <- glmnet(x = c.trn.x, y = c.trn.y, 
                    dfmax = 32, 
                    alpha = 1, 
                    family = "mgaussian")
glm.gamma.df <- predict.glmnet(glm.gamma, c.tst.x, 
                               s = tail(glm.gamma$lambda, n = 1)) %>%
  cbind(c.tst.y) %>% data.frame() %>%
  rename(pred = X1, obs = c.tst.y)

defaultSummary(glm.gamma.df)
saveRDS(glm.gamma, "./models/glmnet/glm gamma.RDS")

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
ggplot(glm.df, aes(x = pred, y = obs)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) + 
  theme_bw() + 
  coord_fixed() + 
  labs(x = "Predicted DelG, kJ/mol", y = "Expeirmental DelG, kJ/mol", 
       title = "GLMnet on All Datapoints")
ggsave("./graphs/glmnet/2017-07-21 glmnet all datapoints df.png")

ggplot(glm.alpha.df, aes(x = pred, y = obs)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) + 
  coord_fixed() + 
  theme_bw() + 
  labs(x = "Predicted DelG, kJ/mol", y = "Expeirmental DelG, kJ/mol", 
       title = "GLMnet on Alpha-CD")
ggsave("./graphs/glmnet/2017-07-21 glmnet alpha-cd df.png")

ggplot(glm.beta.df, aes(x = pred, y = obs)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) + 
  coord_fixed() + 
  theme_bw() + 
  labs(x = "Predicted DelG, kJ/mol", y = "Expeirmental DelG, kJ/mol", 
       title = "GLMnet on Beta-CD")
ggsave("./graphs/glmnet/2017-07-21 glmnet beta-cd df.png")

ggplot(glm.gamma.df, aes(x = pred, y = obs)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) + 
  # coord_fixed() + 
  theme_bw() + 
  labs(x = "Predicted DelG, kJ/mol", y = "Expeirmental DelG, kJ/mol", 
       title = "GLMnet on Gamma-CD")
ggsave("./graphs/glmnet/2017-07-21 glmnet gamma-cd df.png")

ggplot(glm.abc.tst, aes(x = pred, y = obs, color = cd.type)) + 
  geom_point() + 
  coord_fixed() + 
  geom_abline(intercept = 0, slope = 1) + 
  theme_bw() + 
  labs(x = "Predicted DelG, kJ/mol", y = "Expeirmental DelG, kJ/mol", 
       title = "GLMnet on Separate CDs", color = "Cyclodextrin")
ggsave("./graphs/glmnet/2017-07-21 glmnet compiled df.png")

ggplot(glm.abc.tst, aes(x = obs, y = residual, color = cd.type)) + 
  geom_point() + 
  coord_fixed() + 
  geom_hline(yintercept = 0) + 
  theme_bw() + 
  labs(y = "Residuals, kJ/mol", x = "Expeirmental DelG, kJ/mol", 
       title = "GLMnet Residuals on Separate CDs", 
       color = "Cyclodextrin")
ggsave("./graphs/glmnet/2017-07-21 glmnet compiled resid df.png")
