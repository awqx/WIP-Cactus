# Libraries and Packages --------------------------------------------------

library(caret)
library(randomForest)
library(tidyverse)
# Data Organization -------------------------------------------------------

# setwd("~/SREP LAB/qsar")
rpt <- readRDS("./rpt.RDS")
mat.dg <- rpt %>% 
  dplyr::select(., -X1:-log.K.Uncertainty,
                -DelG.Uncertainty:-`bind.aff, kcal/mol`)

set.seed(2)
trn.ind <- sample(x = 1:nrow(mat.dg), size = round(0.8 * nrow(mat.dg)))
trn <- mat.dg[trn.ind, ]
tst <- mat.dg[-trn.ind, ]

# Random Forest Model -----------------------------------------------------

rf.x <- trn[ , -1]
rf.y <- trn[ , 1]
rf <- randomForest(DelG ~., data = trn, 
                   ntree = 500, na.action = na.omit, 
                   mtry = 700, nodesize = 1, 
                   importance = T)
rf.tst.df <- predict(rf, tst[ , -1]) %>%
  cbind(tst[ , 1]) %>% 
  data.frame()
colnames(rf.tst.df)[1] <- "pred"
colnames(rf.tst.df)[2] <- "obs"

rf.tst.df <- rf.tst.df %>%
  mutate(tst.resid = obs - pred)

defaultSummary(rf.tst.df)

# Training Data
rf.trn.df <- predict(rf, trn[ , -1]) %>%
  cbind(trn[ , 1]) %>% 
  data.frame()
colnames(rf.trn.df)[1] <- "pred"
colnames(rf.trn.df)[2] <- "obs"

rf.trn.df <- rf.trn.df %>%
  mutate(trn.resid = obs - pred)

defaultSummary(rf.trn.df)

# Graphs ------------------------------------------------------------------

dir.create("./graphs")
dir.create("./graphs/rforest")
# Test data
rf.tst.df %>% ggplot(., aes(x = obs, y = pred)) +
  geom_point() + 
  geom_abline(slope = 1, intercept = 0) + 
  coord_fixed() + 
  theme_bw() + 
  labs(x = "Experimental DelG Observations, kJ/mol", 
       y = "Predicted DelG, kJ/mol", 
       title = "Random Forest on Testing Data")
ggsave("./graphs/rforest/2017-07-13 rforest test.png")

# Training data
rf.trn.df %>% ggplot(., aes(x = obs, y = pred)) +
  geom_point() + 
  geom_abline(slope = 1, intercept = 0) + 
  coord_fixed() + 
  theme_bw() + 
  labs(x = "Experimental DelG Observations, kJ/mol", 
       y = "Predicted DelG, kJ/mol", 
       title = "Random Forest on Training Data")
ggsave("./graphs/rforest/2017-07-13 rforest trn.png")

# Residuals
ggplot(rf.tst.df, aes(x = obs, y = tst.resid)) + 
  geom_point() + 
  geom_hline(yintercept = 0) + 
  coord_fixed() + 
  theme_bw() + 
  labs(x = "Experimental DelG Observations, kJ/mol", 
       y = "Residuals, kJ/mol", 
       title = "Random Forest Residuals on Testing Data")
ggsave("./graphs/rforest/2017-07-13 rforest tst resid.png")

rf.trn.df %>% ggplot(., aes(x = obs, y = trn.resid)) + 
  geom_point() + 
  geom_hline(yintercept = 0) + 
  coord_fixed() + 
  theme_bw() + 
  labs(x = "Experimental DelG Observations, kJ/mol", 
       y = "Residuals, kJ/mol", 
       title = "Random Forest Residuals on Training Data")
ggsave("./graphs/rforest/2017-07-13 rforest trn resid.png")
