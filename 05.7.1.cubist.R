# Libraries and Packages --------------------------------------------------

library(caret)
library(Cubist)
library(tidyverse)

# Loading Data ------------------------------------------------------------

setwd("~/SREP LAB/qsar")
dir.create("./models/cubist")

df.raw <- readRDS("./padel.pp.new.RDS")
df <- df.raw %>% select(-guest:-host) %>%
  select(-data.source)
mat <- as.matrix(df)

set.seed(12)
trn.ind <- sample(x = 1:nrow(mat), 
                  size = round(0.7 * nrow(mat)))
trn.x <- mat[trn.ind, -1]
trn.y <- mat[trn.ind, 1]
tst.x <- mat[-trn.ind, -1]
tst.y <- mat[-trn.ind, 1]


# All Predictors ----------------------------------------------------------
ctrl <- cubistControl(
  seed = 12, 
  sample = 90
)
cube <- cubist(trn.x, trn.y, 
               control = ctrl, 
               committees = 90)
all.tst <- predict(cube, tst.x) %>%
  cbind(tst.y) %>%
  data.frame() %>%
  rename(., pred = `.`, obs = tst.y)

defaultSummary(all.tst) # 0.531

# Separate Cyclodextrin ---------------------------------------------------

# Alpha-CD ----------------------------------------------------------------

alpha <- df %>% filter(alpha > 0) %>% as.matrix()
set.seed(12)
trn.ind <- sample(x = 1:nrow(alpha), 
                  size = round(0.7 * nrow(alpha)))
a.trn.x <- alpha[trn.ind, -1]
a.trn.y <- alpha[trn.ind, 1]
a.tst.x <- alpha[-trn.ind, -1]
a.tst.y <- alpha[-trn.ind, 1]

cube.alpha <- cubist(a.trn.x, a.trn.y, 
               control = ctrl, 
               committees = 90)
alpha.tst <- predict(cube.alpha, a.tst.x) %>%
  cbind(a.tst.y) %>%
  data.frame() %>%
  rename(., pred = `.`, obs = a.tst.y)

defaultSummary(alpha.tst) # 0.683

# Beta-CD ----------------------------------------------------------------

beta <- df %>% filter(beta > 0) %>% as.matrix()
set.seed(13)
trn.ind <- sample(x = 1:nrow(beta), 
                  size = round(0.7 * nrow(beta)))
b.trn.x <- beta[trn.ind, -1]
b.trn.y <- beta[trn.ind, 1]
b.tst.x <- beta[-trn.ind, -1]
b.tst.y <- beta[-trn.ind, 1]

cube.beta <- cubist(b.trn.x, b.trn.y, 
                     control = ctrl, 
                     committees = 90)
beta.tst <- predict(cube.beta, b.tst.x) %>%
  cbind(b.tst.y) %>%
  data.frame() %>%
  rename(., pred = `.`, obs = b.tst.y)

defaultSummary(beta.tst) # 0.799

# Gamma-CD ----------------------------------------------------------------

gamma <- df %>% filter(gamma > 0) %>% as.matrix()
set.seed(13)
trn.ind <- sample(x = 1:nrow(gamma), 
                  size = round(0.7 * nrow(gamma)))
c.trn.x <- gamma[trn.ind, -1]
c.trn.y <- gamma[trn.ind, 1]
c.tst.x <- gamma[-trn.ind, -1]
c.tst.y <- gamma[-trn.ind, 1]

cube.gamma <- cubist(c.trn.x, c.trn.y, 
                    control = ctrl, 
                    committees = 90)
gamma.tst <- predict(cube.gamma, c.tst.x) %>%
  cbind(c.tst.y) %>%
  data.frame() %>%
  rename(., pred = `.`, obs = c.tst.y)

defaultSummary(gamma.tst) # 0.799


#     Compiled CDs --------------------------------------------------------

temp.a <- alpha.tst %>% 
  mutate(cd.type = "alpha")
temp.b <- beta.tst %>% 
  mutate(cd.type = "beta")
temp.c <- gamma.tst %>% 
  mutate(cd.type = "gamma")
cube.abc.tst <- rbind(temp.a, temp.b, temp.c, 
                     make.row.names = F) %>%
  mutate(residual = pred - obs)
defaultSummary(cube.abc.tst) # 0.738

# Graphs ------------------------------------------------------------------

ggplot(all.tst, aes(x = pred, y = obs)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) + 
  theme_bw() + 
  coord_fixed() + 
  labs(x = "Predicted DelG, kJ/mol", y = "Experimental DelG, kJ/mol", 
       title = "Cubist - All Data Points") 
ggsave("./models/cubist/2017-07-21 cubist all data.png")

ggplot(alpha.tst, aes(x = pred, y = obs)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) + 
  theme_bw() + 
  coord_fixed() + 
  labs(x = "Predicted DelG, kJ/mol", y = "Experimental DelG, kJ/mol", 
       title = "Cubist - Alpha-CD") 
ggsave("./models/cubist/2017-07-21 cubist alpha-cd.png")

ggplot(beta.tst, aes(x = pred, y = obs)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) + 
  theme_bw() + 
  coord_fixed() + 
  labs(x = "Predicted DelG, kJ/mol", y = "Experimental DelG, kJ/mol", 
       title = "Cubist - Beta-CD") 
ggsave("./models/cubist/2017-07-21 cubist beta-cd.png")

ggplot(gamma.tst, aes(x = pred, y = obs)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) + 
  theme_bw() + 
  coord_fixed() + 
  labs(x = "Predicted DelG, kJ/mol", y = "Experimental DelG, kJ/mol", 
       title = "Cubist - Gamma-CD") 
ggsave("./models/cubist/2017-07-21 cubist gamma-cd.png")

ggplot(cube.abc.tst, aes(x = pred, y = obs, color = cd.type)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) + 
  theme_bw() + 
  coord_fixed() + 
  labs(x = "Predicted DelG, kJ/mol", y = "Experimental DelG, kJ/mol", 
       title = "Cubist - Compiled CDs", color = "Cyclodextrin") 
ggsave("./models/cubist/2017-07-21 cubist gamma-cd.png")
