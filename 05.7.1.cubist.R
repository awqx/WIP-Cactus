# Libraries and Packages --------------------------------------------------

library(caret)
library(Cubist)
library(tidyverse)

# Loading Data ------------------------------------------------------------

# setwd("~/SREP LAB/qsar")
df.raw <- readRDS("./padel.pp.new.RDS")
df <- df.raw %>% select(-guest:-host) %>%
  select(-data.source)
mat <- as.matrix(df)

set.seed(10)
trn.ind <- sample(x = 1:nrow(mat), 
                  size = round(0.7 * nrow(mat)))
trn.x <- mat[trn.ind, -1]
trn.y <- mat[trn.ind, 1]
tst.x <- mat[-trn.ind, -1]
tst.y <- mat[-trn.ind, 1]


# All Predictors ----------------------------------------------------------
ctrl <- cubistControl(
  seed = 10, 
  sample = 75
)

# Took around 40 seconds to create model
cube <- cubist(trn.x, trn.y, 
               control = ctrl, 
               committees = 75)

all.tst <- predict(cube, tst.x) %>%
  cbind(tst.y) %>%
  data.frame() %>%
  dplyr::rename(., pred = `.`, obs = tst.y)

defaultSummary(all.tst) # 0.530
saveRDS(cube, "./models/cubist/cube.all.RDS")

# Separate Cyclodextrin ---------------------------------------------------

#     Alpha-CD ------------------------------------------------------------

alpha <- df %>% filter(alpha > 0) %>% as.matrix()
set.seed(10)
trn.ind <- sample(x = 1:nrow(alpha), 
                  size = round(0.7 * nrow(alpha)))
a.trn.x <- alpha[trn.ind, -1]
a.trn.y <- alpha[trn.ind, 1]
a.tst.x <- alpha[-trn.ind, -1]
a.tst.y <- alpha[-trn.ind, 1]

cube.a <- cubist(a.trn.x, a.trn.y, 
               control = ctrl, 
               committees = 90)
cube.a.tst <- predict(cube.a, a.tst.x) %>%
  cbind(a.tst.y) %>%
  data.frame() %>%
  dplyr::rename(., pred = `.`, obs = a.tst.y)
cube.a.trn <- predict(cube.a, a.trn.x) %>%
  cbind(a.trn.y) %>%
  data.frame() %>%
  dplyr::rename(., pred = `.`, obs = a.trn.y)

defaultSummary(cube.a.tst) # 0.601
defaultSummary(cube.a.trn) # 0.997

#     Beta-CD ------------------------------------------------------------

beta <- df %>% filter(beta > 0) %>% as.matrix()
set.seed(10)
trn.ind <- sample(x = 1:nrow(beta), 
                  size = round(0.7 * nrow(beta)))
b.trn.x <- beta[trn.ind, -1]
b.trn.y <- beta[trn.ind, 1]
b.tst.x <- beta[-trn.ind, -1]
b.tst.y <- beta[-trn.ind, 1]

cube.b <- cubist(b.trn.x, b.trn.y, 
                     control = ctrl, 
                     committees = 90)
cube.b.tst <- predict(cube.b, b.tst.x) %>%
  cbind(b.tst.y) %>%
  data.frame() %>%
  dplyr::rename(., pred = `.`, obs = b.tst.y)
cube.b.trn <- predict(cube.b, b.trn.x) %>%
  cbind(b.trn.y) %>%
  data.frame() %>%
  dplyr::rename(., pred = `.`, obs = b.trn.y)

defaultSummary(cube.b.tst) # 0.701
defaultSummary(cube.b.trn) # 0.999

# Gamma-CD ----------------------------------------------------------------

gamma <- df %>% filter(gamma > 0) %>% as.matrix()
set.seed(10)
trn.ind <- sample(x = 1:nrow(gamma), 
                  size = round(0.7 * nrow(gamma)))
c.trn.x <- gamma[trn.ind, -1]
c.trn.y <- gamma[trn.ind, 1]
c.tst.x <- gamma[-trn.ind, -1]
c.tst.y <- gamma[-trn.ind, 1]

cube.c <- cubist(c.trn.x, c.trn.y, 
                    control = ctrl, 
                    committees = 70)
cube.c.tst <- predict(cube.c, c.tst.x) %>%
  cbind(c.tst.y) %>%
  data.frame() %>%
  dplyr::rename(., pred = `.`, obs = c.tst.y)
cube.c.trn <- predict(cube.c, c.trn.x) %>%
  cbind(c.trn.y) %>%
  data.frame() %>%
  dplyr::rename(., pred = `.`, obs = c.trn.y)

defaultSummary(cube.c.tst) # 0.332
defaultSummary(cube.c.trn) # 1.000

#     Compiled CDs --------------------------------------------------------

temp.a <- cube.a.tst %>% 
  mutate(cd.type = "alpha")
temp.b <- cube.b.tst %>% 
  mutate(cd.type = "beta")
temp.c <- cube.c.tst %>% 
  mutate(cd.type = "gamma")
cube.abc.tst <- rbind(temp.a, temp.b, temp.c, 
                     make.row.names = F) %>%
  mutate(residual = pred - obs)
defaultSummary(cube.abc.tst) # 0.622

temp.a <- cube.a.trn %>% 
  mutate(cd.type = "alpha")
temp.b <- cube.b.trn %>% 
  mutate(cd.type = "beta")
temp.c <- cube.c.trn %>% 
  mutate(cd.type = "gamma")
cube.abc.trn <- rbind(temp.a, temp.b, temp.c, 
                      make.row.names = F) %>%
  mutate(residual = pred - obs)
defaultSummary(cube.abc.trn) # 0.998

# Saving Models -----------------------------------------------------------

dir.create("./models/cubist")
saveRDS(cube.a, "./models/cubist/cubist.alpha.RDS")
saveRDS(cube.b, "./models/cubist/cubist.beta.RDS")
saveRDS(cube.c, "./models/cubist/cubist.gamma.RDS")

saveRDS(cube.abc.tst, "./models/cubist/cubist.results.RDS")
saveRDS(cube.abc.trn, "./models/cubist/cubist.trn.results.RDS")

# Graphs ------------------------------------------------------------------

dir.create("./graphs/cubist")

ggplot(all.tst, aes(x = obs, y = pred)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) + 
  theme_bw() + 
  coord_fixed(xlim = c(-45, 5), ylim = c(-45, 5)) + 
  labs(y = "Predicted DelG, kJ/mol", x = "Observed DelG, kJ/mol", 
       title = "Cubist - All Data Points") 
ggsave("./graphs/cubist/2017-07-28 cubist all data.png")

ggplot(alpha.tst, aes(y = pred, x = obs)) + 
  geom_point(color = "#F8766D") + 
  geom_abline(intercept = 0, slope = 1) + 
  theme_bw() + 
  coord_fixed(xlim = c(-45, 5), ylim = c(-45, 5)) + 
  labs(y = "Predicted DelG, kJ/mol", x = "Observed DelG, kJ/mol", 
       title = "Cubist - Alpha-CD") 
ggsave("./graphs/cubist/2017-07-28 cubist alpha-cd.png")

ggplot(beta.tst, aes(y = pred, x = obs)) + 
  geom_point(color = "#00BA38") + 
  geom_abline(intercept = 0, slope = 1) + 
  theme_bw() + 
  coord_fixed(xlim = c(-45, 5), ylim = c(-45, 5)) + 
  labs(y = "Predicted DelG, kJ/mol", x = "Observed DelG, kJ/mol", 
       title = "Cubist - Beta-CD") 
ggsave("./graphs/cubist/2017-07-28 cubist beta-cd.png")

ggplot(gamma.tst, aes(y = pred, x = obs)) + 
  geom_point(color = "#619CFF" ) + 
  geom_abline(intercept = 0, slope = 1) + 
  theme_bw() + 
  coord_fixed(xlim = c(-45, 5), ylim = c(-45, 5)) + 
  labs(y = "Predicted DelG, kJ/mol", x = "Observed DelG, kJ/mol", 
       title = "Cubist - Gamma-CD") 
ggsave("./graphs/cubist/2017-07-28 cubist gamma-cd.png")

ggplot(cube.abc.tst, aes(y = pred, x = obs, color = cd.type)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) + 
  theme_bw() + 
  coord_fixed(xlim = c(-45, 5), ylim = c(-45, 5)) + 
  labs(y = "Predicted DelG, kJ/mol", x = "Observed DelG, kJ/mol", 
       title = "Cubist - Compiled CDs", color = "Cyclodextrin") 
ggsave("./graphs/cubist/2017-07-28 cubist compiled cd.png")

ggplot(cube.abc.tst, aes(x = obs, y = residual, color = cd.type)) + 
  geom_point() + 
  geom_hline(yintercept = 0) + 
  theme_bw() + 
  coord_fixed() + 
  labs(y = "Residuals, kJ/mol", x = "Observed DelG, kJ/mol", 
       title = "Cubist - Compiled CDs Residuals", color = "Cyclodextrin")
ggsave("./graphs/cubist/2017-07-28 cubist compiled cd resid.png")

ggplot(cube.abc.trn, aes(y = pred, x = obs, color = cd.type)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) + 
  theme_bw() + 
  coord_fixed(xlim = c(-45, 5), ylim = c(-45, 5)) + 
  labs(y = "Predicted DelG, kJ/mol", x = "Observed DelG, kJ/mol", 
       title = "Cubist - Compiled CDs, Training Data", color = "Cyclodextrin") 
ggsave("./graphs/cubist/2017-07-28 cubist compiled trn cd.png")

#####

# External  Validation ----------------------------------------------------

#####

ext.val <- readRDS("./external validation set new.RDS") %>%
  select(-guest:-host) %>%
  select(-data.source) 

ext.val.a <- ext.val %>% filter(alpha > 0) 
ext.val.b <- ext.val %>% filter(beta > 0)
ext.val.c <- ext.val %>% filter(gamma > 0) 

ev.a <-  predict(cube.alpha, ext.val.a[ , -1]) %>%
  cbind(ext.val.a[ , 1]) %>% data.frame() %>%
  dplyr::rename(., pred = `.`, obs = V2)

ev.b <-  predict(cube.beta, ext.val.b[ , -1]) %>%
  cbind(ext.val.b[ , 1]) %>% data.frame() %>%
  dplyr::rename(., pred = `.`, obs = V2)

ev.c <-  predict(cube.gamma, ext.val.c[ , -1]) %>%
  cbind(ext.val.c[ , 1]) %>% data.frame() %>%
  dplyr::rename(., pred = `.`, obs = V2)

temp.a <- ev.a %>% mutate(cd.type = "alpha")
temp.b <- ev.b %>% mutate(cd.type = "beta")
temp.c <- ev.c %>% mutate(cd.type = "gamma")
# temp.a <- temp.a[-2, ]
ev.abc <- rbind(temp.a, temp.b, temp.c) %>%
  mutate(resid = pred - obs)
# Prearing copy for compilation
ev.cube <- ev.abc %>% select(-resid) %>%
  mutate(Model = "Cubist")

defaultSummary(ev.abc) # 0.384 normally, 0.658 without outlier
ggplot(ev.abc, aes(y = pred, x = obs, color = cd.type)) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0) + 
  coord_fixed(xlim = c(-45, 5), ylim = c(-45, 5)) + 
  theme_bw() + 
  labs(title = "Cubist - External Validation", 
       x = "Predicted DelG, kJ/mol", y = "Observed DelG, kJ/mol", 
       color = "Cyclodextrin")
ggsave("./graphs/cubist/2017-07-28 cubist extval.png")

