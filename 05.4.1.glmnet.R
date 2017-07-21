# Libraries ---------------------------------------------------------------

library(caret)
library(ggplot2)
library(glmnet)

# Loading Data ------------------------------------------------------------

setwd("~/SREP LAB/qsar")

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


#     Beta-CD -------------------------------------------------------------


#     Gamma-CD ------------------------------------------------------------


#     Compiled CDs --------------------------------------------------------


# Graphs ------------------------------------------------------------------

ggplot(glm.df, aes(x = pred, y = obs)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) + 
  theme_bw()
