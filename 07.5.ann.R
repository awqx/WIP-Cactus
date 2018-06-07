# Libraries and Packages --------------------------------------------------

library(caret)
library(neuralnet)
library(tidyverse)

# Data --------------------------------------------------------------------

# setwd("~/SREP LAB/qsar")
df.raw <- readRDS("./padel.pp.new.RDS")
df <- df.raw %>% select(-guest:-host) %>%
  select(-data.source) %>%
  as.matrix()

set.seed(12)
trn.ind <- sample(x = 1:nrow(df), 
                  size = round(0.7 * nrow(df)))
trn <- df[trn.ind, ]
tst <- df[-trn.ind, ]

formula <- as.formula(paste("DelG", paste(colnames(df)[2:712], collapse = " + "), sep = " ~ "))

# Neural Network ----------------------------------------------------------

#     All CD --------------------------------------------------------------

nn.all <- neuralnet(formula, data = trn, hidden = 7)
nn.all.tst <- neuralnet::compute(nn.all, tst[ , -1]) 
nn.all.tst <- nn.all.tst$net.result %>%
  cbind(tst[ , 1]) %>%
  data.frame() %>%
  rename(pred = X1, obs = X2) 
defaultSummary(nn.all.tst) # 0.491

#     Alpha-CD ------------------------------------------------------------

alpha <- df.raw %>% select(-guest:-host) %>%
  select(-data.source) %>%
  filter(alpha > 0) %>%
  as.matrix()

set.seed(12)
trn.ind <- sample(x = 1:nrow(alpha), 
                  size = round(0.7 * nrow(alpha)))
tst.alpha <- alpha[trn.ind, ]
trn.alpha <- alpha[-trn.ind, ]

nn.alpha <- neuralnet(formula, data = trn.alpha, hidden = 5)
nn.alpha.tst <- neuralnet::compute(nn.alpha, tst.alpha[ , -1]) 
nn.alpha.tst <- nn.alpha.tst$net.result %>%
  cbind(tst.alpha[ , 1]) %>%
  data.frame() %>%
  rename(pred = X1, obs = X2) 
nn.alpha.trn <- neuralnet::compute(nn.alpha, trn.alpha[ , -1]) 
nn.alpha.trn <- nn.alpha.trn$net.result %>%
  cbind(trn.alpha[ , 1]) %>%
  data.frame() %>%
  rename(pred = X1, obs = X2) 

defaultSummary(nn.alpha.tst)
defaultSummary(nn.alpha.trn)

#     Beta-CD -------------------------------------------------------------

beta <- df.raw %>% select(-guest:-host) %>%
  select(-data.source) %>%
  filter(beta > 0) %>%
  as.matrix()

set.seed(12)
trn.ind <- sample(x = 1:nrow(beta), 
                  size = round(0.7 * nrow(beta)))
tst.beta <- beta[trn.ind, ]
trn.beta <- beta[-trn.ind, ]

nn.beta <- neuralnet(formula, data = trn.beta, hidden = 5)
nn.beta.tst <- neuralnet::compute(nn.beta, tst.beta[ , -1]) 
nn.beta.tst <- nn.beta.tst$net.result %>%
  cbind(tst.beta[ , 1]) %>%
  data.frame() %>%
  rename(pred = X1, obs = X2) 
nn.beta.trn <- neuralnet::compute(nn.beta, trn.beta[ , -1]) 
nn.beta.trn <- nn.beta.trn$net.result %>%
  cbind(trn.beta[ , 1]) %>%
  data.frame() %>%
  rename(pred = X1, obs = X2) 

defaultSummary(nn.beta.tst)
defaultSummary(nn.beta.trn)

#     Gamma-CD ------------------------------------------------------------

gamma <- df.raw %>% select(-guest:-host) %>%
  select(-data.source) %>%
  filter(gamma > 0) %>%
  as.matrix()

set.seed(12)
trn.ind <- sample(x = 1:nrow(gamma), 
                  size = round(0.7 * nrow(gamma)))
tst.gamma <- gamma[trn.ind, ]
trn.gamma <- gamma[-trn.ind, ]

nn.gamma <- neuralnet(formula, data = trn.gamma, hidden = 7)
nn.gamma.tst <- neuralnet::compute(nn.gamma, tst.gamma[ , -1]) 
nn.gamma.tst <- nn.gamma.tst$net.result %>%
  cbind(tst.gamma[ , 1]) %>%
  data.frame() %>%
  rename(pred = X1, obs = X2) 
nn.gamma.trn <- neuralnet::compute(nn.gamma, trn.gamma[ , -1]) 
nn.gamma.trn <- nn.gamma.trn$net.result %>%
  cbind(trn.gamma[ , 1]) %>%
  data.frame() %>%
  rename(pred = X1, obs = X2) 

defaultSummary(nn.gamma.tst)
defaultSummary(nn.gamma.trn)

#     Compile -------------------------------------------------------------

nn.abc.tst <- rbind(mutate(nn.alpha.tst, cd.type = "alpha"), 
                    mutate(nn.beta.tst, cd.type = "beta"), 
                    mutate(nn.gamma.tst, cd.type = "gamma"))
nn.abc.trn <- rbind(mutate(nn.alpha.trn, cd.type = "alpha"), 
                    mutate(nn.beta.trn, cd.type = "beta"), 
                    mutate(nn.gamma.trn, cd.type = "gamma"))

defaultSummary(nn.abc.tst) # 0.360
defaultSummary(nn.abc.trn) # 0.640

nn.ev.a <- compute

# External Validation -----------------------------------------------------

ext.val <- readRDS("./external validation set new.RDS") %>%
  select(-guest:-host) %>%
  select(-data.source) 

ext.val.a <- ext.val %>% filter(alpha > 0) %>% as.matrix()
ext.val.b <- ext.val %>% filter(beta > 0) %>% as.matrix()
ext.val.c <- ext.val %>% filter(gamma > 0) %>% as.matrix()

nn.ev.a <- neuralnet::compute(nn.alpha, ext.val.a[ , -1]) 
nn.ev.a <- nn.ev.a$net.result %>%
  cbind(ext.val.a[ , 1]) %>%
  data.frame() %>%
  rename(pred = X1, obs = X2) 

nn.ev.b <- neuralnet::compute(nn.beta, ext.val.b[ , -1]) 
nn.ev.b <- nn.ev.b$net.result %>%
  cbind(ext.val.b[ , 1]) %>%
  data.frame() %>%
  rename(pred = X1, obs = X2) 

temp <- rbind(ext.val.b, ext.val.c)
nn.ev.c <- neuralnet::compute(nn.gamma, temp[ , -1]) 
nn.ev.c <- nn.ev.c$net.result %>%
  cbind(temp[ , 1]) %>%
  data.frame() %>%
  rename(pred = X1, obs = X2) 
nn.ev.c <- nn.ev.c[58, ]

nn.abc.ev <- rbind(mutate(nn.ev.a, cd.type = "alpha"), 
                   mutate(nn.ev.b, cd.type = "beta"), 
                   mutate(nn.ev.c, cd.type = "gamma"))
defaultSummary(nn.abc.ev)

# Graphs ------------------------------------------------------------------

dir.create("./graphs/ann")

ggplot(nn.abc.tst, aes(x = obs, y = pred, color = cd.type)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) + 
  coord_fixed(xlim = c(-45, 5), ylim = c(-45, 5)) + 
  theme_bw() + 
  labs(x = "Observed Delg, kJ/mol", y = "Predicted DelG, kJ/mol", 
       title = "Artificial Neural Network", 
       color = "Cyclodextrin")
ggsave("./graphs/ann/ann tst.png")

ggplot(nn.abc.trn, aes(x = obs, y = pred, color = cd.type)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) + 
  coord_fixed(xlim = c(-45, 5), ylim = c(-45, 5)) + 
  theme_bw() + 
  labs(x = "Observed Delg, kJ/mol", y = "Predicted DelG, kJ/mol", 
       title = "Artificial Neural Network, Training", 
       color = "Cyclodextrin")
ggsave("./graphs/ann/ann trn.png")

ggplot(nn.abc.ev, aes(x = obs, y = pred, color = cd.type)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) + 
  coord_fixed(xlim = c(-45, 5), ylim = c(-45, 5)) + 
  theme_bw() + 
  labs(x = "Observed Delg, kJ/mol", y = "Predicted DelG, kJ/mol", 
       title = "Artificial Neural Network, External Validation", 
       color = "Cyclodextrin")
ggsave("./graphs/ann/ann extval.png")
