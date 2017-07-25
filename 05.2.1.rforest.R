# Libraries and Packages --------------------------------------------------

library(caret)
library(randomForest)
library(tidyverse)

# Data Organization -------------------------------------------------------

# setwd("~/SREP LAB/qsar")
df <- readRDS("./padel.pp.new.RDS") %>%
  select(-guest:-host) %>%
  select(-data.source)
mat <- df %>% as.matrix()

set.seed(2)
trn.ind <- sample(x = 1:nrow(mat), size = round(0.7 * nrow(mat)))
trn <- mat[trn.ind, ]
trn.x <- mat[trn.ind, -1]
trn.y <- mat[trn.ind, 1]
tst.x <- mat[-trn.ind, -1]
tst.y <- mat[-trn.ind, 1]

# Random Forest Model -----------------------------------------------------

rf <- randomForest(x = trn.x, y = trn.y, 
                   ntree = 500, na.action = na.omit, 
                   mtry = 700, nodesize = 1, 
                   importance = T)
rf.tst.df <- predict(rf, tst.x) %>%
  cbind(tst.y) %>% 
  data.frame() %>%
  rename(., pred = `.`, obs = tst.y) %>%
  mutate(tst.resid = obs - pred)

# Training Data
rf.trn.df <- predict(rf, trn.x) %>%
  cbind(trn.y) %>% 
  data.frame() %>%
  rename(., pred = `.`, obs = trn.y) %>%
  mutate(trn.resid = obs - pred)

defaultSummary(rf.tst.df) # 0.680
defaultSummary(rf.trn.df) # 0.950

#     Alpha-CD ------------------------------------------------------------

alpha <- df %>% filter(alpha > 0) %>% as.matrix()

set.seed(2)
trn.ind <- sample(x = 1:nrow(alpha), size = round(0.7 * nrow(alpha)))
a.trn.x <- alpha[trn.ind, -1]
a.trn.y <- alpha[trn.ind, 1]
a.tst.x <- alpha[-trn.ind, -1]
a.tst.y <- alpha[-trn.ind, 1]

rf.alpha <- randomForest(x = a.trn.x, y = a.trn.y,
  ntree = 500, mtry = 700,
  na.action = na.omit,
  nodesize = 1,  importance = T
)

rf.tst.a <- predict(rf.alpha, a.tst.x) %>%
  cbind(a.tst.y) %>% 
  data.frame() %>%
  rename(., pred = `.`, obs = a.tst.y) %>%
  mutate(tst.resid = obs - pred)

# Training Data
rf.trn.a <- predict(rf.alpha, a.trn.x) %>%
  cbind(a.trn.y) %>% 
  data.frame() %>%
  rename(., pred = `.`, obs = a.trn.y) %>%
  mutate(trn.resid = obs - pred)

defaultSummary(rf.tst.a) # 0.541
defaultSummary(rf.trn.a) # 0.957

#     Beta-CD -------------------------------------------------------------

beta <- df %>% filter(beta > 0) %>% as.matrix()

set.seed(2)
trn.ind <- sample(x = 1:nrow(beta), size = round(0.7 * nrow(beta)))
b.trn.x <- beta[trn.ind, -1]
b.trn.y <- beta[trn.ind, 1]
b.tst.x <- beta[-trn.ind, -1]
b.tst.y <- beta[-trn.ind, 1]

rf.beta <- randomForest(x = b.trn.x, y = b.trn.y,
                         ntree = 500, mtry = 700,
                         na.bction = na.omit,
                         nodesize = 1,  importance = T
)

rf.tst.b <- predict(rf.beta, b.tst.x) %>%
  cbind(b.tst.y) %>% 
  data.frame() %>%
  rename(., pred = `.`, obs = b.tst.y) %>%
  mutate(tst.resid = obs - pred)

# Training Data
rf.trn.b <- predict(rf.beta, b.trn.x) %>%
  cbind(b.trn.y) %>% 
  data.frame() %>%
  rename(., pred = `.`, obs = b.trn.y) %>%
  mutate(trn.resid = obs - pred)

defaultSummary(rf.tst.b) # 0.750
defaultSummary(rf.trn.b) # 0.959

#     Gamma-CD ------------------------------------------------------------

gamma <- df %>% filter(gamma > 0) %>% as.matrix()

set.seed(2)
trn.ind <- sample(x = 1:nrow(gamma), size = round(0.7 * nrow(gamma)))
c.trn.x <- gamma[trn.ind, -1]
c.trn.y <- gamma[trn.ind, 1]
c.tst.x <- gamma[-trn.ind, -1]
c.tst.y <- gamma[-trn.ind, 1]

rf.gamma <- randomForest(x = c.trn.x, y = c.trn.y,
                        ntree = 500, mtry = 700,
                        na.cction = na.omit,
                        nodesize = 1,  importance = T
)

rf.tst.c <- predict(rf.gamma, c.tst.x) %>%
  cbind(c.tst.y) %>% 
  data.frame() %>%
  rename(., pred = `.`, obs = c.tst.y) %>%
  mutate(tst.resid = obs - pred)

# Training Data
rf.trn.c <- predict(rf.gamma, c.trn.x) %>%
  cbind(c.trn.y) %>% 
  data.frame() %>%
  rename(., pred = `.`, obs = c.trn.y) %>%
  mutate(trn.resid = obs - pred)

defaultSummary(rf.tst.c) # 0.634
defaultSummary(rf.trn.c) # 0.980

#     Compiled CDs --------------------------------------------------------

temp.a <- rf.tst.a %>% 
  mutate(cd.type = "alpha")
temp.b <- rf.tst.b %>% 
  mutate(cd.type = "beta")
temp.c <- rf.tst.c %>% 
  mutate(cd.type = "gamma")
rf.abc.tst <- rbind(temp.a, temp.b, temp.c, 
                     make.row.names = F) 

temp.a <- rf.trn.a %>% 
  mutate(cd.type = "alpha")
temp.b <- rf.trn.b %>% 
  mutate(cd.type = "beta")
temp.c <- rf.trn.c %>% 
  mutate(cd.type = "gamma")
rf.abc.trn <- rbind(temp.a, temp.b, temp.c, 
                    make.row.names = F)

defaultSummary(rf.abc.tst) # 0.700
defaultSummary(rf.abc.trn) # 0.960

# Graphs ------------------------------------------------------------------

dir.create("./graphs")
dir.create("./graphs/rforest")

#     Test Data -----------------------------------------------------------

# All Datapoints
rf.tst.df %>% ggplot(., aes(x = obs, y = pred)) +
  geom_point() + 
  geom_abline(slope = 1, intercept = 0) + 
  coord_fixed() + 
  theme_bw() + 
  labs(x = "Experimental DelG Observations, kJ/mol", 
       y = "Predicted DelG, kJ/mol", 
       title = "Random Forest on Testing Data")
ggsave("./graphs/rforest/2017-07-24 rf test.png")

ggplot(rf.tst.a, aes(x = obs, y = pred)) +
  geom_point() + 
  geom_abline(slope = 1, intercept = 0) + 
  coord_fixed() + 
  theme_bw() + 
  labs(x = "Experimental DelG Observations, kJ/mol", 
       y = "Predicted DelG, kJ/mol", 
       title = "Random Forest - Alpha CD")
ggsave("./graphs/rforest/2017-07-24 rf alphacd tst.png")

ggplot(rf.tst.b, aes(x = obs, y = pred)) +
  geom_point() + 
  geom_abline(slope = 1, intercept = 0) + 
  coord_fixed() + 
  theme_bw() + 
  labs(x = "Experimental DelG Observations, kJ/mol", 
       y = "Predicted DelG, kJ/mol", 
       title = "Random Forest - Beta CD")
ggsave("./graphs/rforest/2017-07-24 rf betacd tst.png")

ggplot(rf.tst.c, aes(x = obs, y = pred)) +
  geom_point() + 
  geom_abline(slope = 1, intercept = 0) + 
  coord_fixed() + 
  theme_bw() + 
  labs(x = "Experimental DelG Observations, kJ/mol", 
       y = "Predicted DelG, kJ/mol", 
       title = "Random Forest - Gamma CD")
ggsave("./graphs/rforest/2017-07-24 rf gammacd tst.png")

ggplot(rf.abc.tst, aes(x = obs, y = pred, color = cd.type)) +
  geom_point() + 
  geom_abline(slope = 1, intercept = 0) + 
  coord_fixed() + 
  theme_bw() + 
  labs(x = "Experimental DelG Observations, kJ/mol", 
       y = "Predicted DelG, kJ/mol", 
       title = "Random Forest - Compiled CD")
ggsave("./graphs/rforest/2017-07-24 rf compiled cd tst.png")


#     Training Data -------------------------------------------------------

rf.trn.df %>% ggplot(., aes(x = obs, y = pred)) +
  geom_point() + 
  geom_abline(slope = 1, intercept = 0) + 
  coord_fixed() + 
  theme_bw() + 
  labs(x = "Experimental DelG Observations, kJ/mol", 
       y = "Predicted DelG, kJ/mol", 
       title = "Random Forest on Training Data")
ggsave("./graphs/rforest/2017-07-24 rf trn.png")

ggplot(rf.abc.trn, aes(x = obs, y = pred, color = cd.type)) +
  geom_point() + 
  geom_abline(slope = 1, intercept = 0) + 
  coord_fixed() + 
  theme_bw() + 
  labs(x = "Experimental DelG Observations, kJ/mol", 
       y = "Predicted DelG, kJ/mol", 
       title = "Random Forest - Compiled CD, Training Data")
ggsave("./graphs/rforest/2017-07-24 rf compiled cd trn.png")

#     Residuals -----------------------------------------------------------

ggplot(rf.abc.tst, aes(x = obs, y = tst.resid, color = cd.type)) + 
  geom_point() + 
  geom_hline(yintercept = 0) + 
  coord_fixed() + 
  theme_bw() + 
  labs(x = "Experimental DelG Observations, kJ/mol", 
       y = "Residuals, kJ/mol", 
       title = "Random Forest Residuals on Testing Data")
ggsave("./graphs/rforest/2017-07-24 rf tst resid.png")

rf.abc.trn %>% ggplot(., aes(x = obs, y = trn.resid, color = cd.type)) + 
  geom_point() + 
  geom_hline(yintercept = 0) + 
  coord_fixed() + 
  theme_bw() + 
  labs(x = "Experimental DelG Observations, kJ/mol", 
       y = "Residuals, kJ/mol", 
       title = "Random Forest Residuals on Training Data")
ggsave("./graphs/rforest/2017-07-13 rforest trn resid.png")
