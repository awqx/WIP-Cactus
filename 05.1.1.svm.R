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

set.seed(25)
trn.ind <- sample(x = 1:nrow(df), size = round(0.7 * nrow(df)))
trn <- df[trn.ind, ]
tst <- df[-trn.ind, ]
trn.x <- trn[ , -1]
trn.y <- trn[ , 1]
tst.x <- tst[ , -1]
tst.y <- tst[ , 1]

# Polynomial Kernel -------------------------------------------------------
#     All Data ------------------------------------------------------------

svm.all <- svm(x = trn.x, 
               y = trn.y, 
               coef0 = 1.5, 
               cost = 4096, 
               epsilon = 0.1, 
               kernel = "polynomial", 
               gamma = 0.03, 
               degree = 2)

svm.all.tst <- predict(svm.all, tst.x) %>%
  cbind(tst.y) %>%
  data.frame() %>%
  rename(., pred = `.`, obs = tst.y) %>%
  mutate(resid = obs - pred)

svm.all.trn <- predict(svm.all, trn.x) %>%
  cbind(trn.y) %>%
  data.frame() %>%
  rename(., pred = `.`, obs = trn.y) %>%
  mutate(resid = obs - pred)

defaultSummary(svm.all.tst) # 0.470
defaultSummary(svm.all.trn) # 0.990

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
  rename(., pred = `.`, obs = a.tst.y) %>%
  mutate(resid = obs - pred)

svm.a.trn <- predict(svm.alpha, a.trn.x) %>%
  cbind(a.trn.y) %>% data.frame() %>%
  rename(., pred = `.`, obs = a.trn.y)%>%
  mutate(resid = obs - pred)

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
  rename(., pred = `.`, obs = b.tst.y) %>%
  mutate(resid = obs - pred)
svm.b.trn <- predict(svm.beta, b.trn.x) %>%
  cbind(b.trn.y) %>% data.frame() %>%
  rename(., pred = `.`, obs = b.trn.y) %>%
  mutate(resid = obs - pred)
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
  rename(., pred = `.`, obs = c.tst.y) %>%
  mutate(resid = obs - pred)

svm.c.trn <- predict(svm.gamma, c.trn.x) %>%
  cbind(c.trn.y) %>% data.frame() %>%
  rename(., pred = `.`, obs = c.trn.y)%>%
  mutate(resid = obs - pred)

defaultSummary(svm.c.tst) # 0.214
defaultSummary(svm.c.trn) # 0.99986

saveRDS(svm.gamma, "./models/svm/polysvm.gamma.RDS")

#     Compiled CDs --------------------------------------------------------

temp.a <- svm.a.tst %>% mutate(cd.type = "alpha")
temp.b <- svm.b.tst %>% mutate(cd.type = "beta")
temp.c <- svm.c.tst %>% mutate(cd.type = "gamma")
poly.abc.tst <- rbind(temp.a, temp.b, temp.c, make.row.names = F)

temp.a <- svm.a.trn %>% mutate(cd.type = "alpha")
temp.b <- svm.b.trn %>% mutate(cd.type = "beta")
temp.c <- svm.c.trn %>% mutate(cd.type = "gamma")
poly.abc.trn <- rbind(temp.a, temp.b, temp.c, make.row.names = F) 

defaultSummary(svm.abc.tst) # 0.614
defaultSummary(svm.abc.trn) # 0.997

saveRDS(poly.abc.tst, "./models/svm/polysvm.tst.results.RDS")
saveRDS(poly.abc.trn, "./models/svm/polysvm.trn.results.RDS")

#####
# Radial Basis Kernel -----------------------------------------------------
#     All ------------------------------------------------------------

rbf.all <- svm(x = trn.x, 
               y = trn.y,
               cost = 4, 
               epsilon = 0.5, 
               kernel = "radial", 
               gamma = 0.001)

rbf.all.tst <- predict(rbf.all, tst.x) %>%
  cbind(tst.y) %>%
  data.frame() %>%
  rename(., pred = `.`, obs = tst.y) %>%
  mutate(resid = obs - pred)

rbf.all.trn <- predict(rbf.all, trn.x) %>%
  cbind(trn.y) %>%
  data.frame() %>%
  rename(., pred = `.`, obs = trn.y) %>%
  mutate(resid = obs - pred)

defaultSummary(rbf.all.tst) # 0.567
defaultSummary(rbf.all.trn) # 0.791

#     Alpha ---------------------------------------------------------

rbf.alpha <- svm(x = a.trn.x,
             y = a.trn.y, 
             kernel = "radial", 
             cost = 4, 
             gamma = 0.001, 
             epsilon = 0.5)
rbf.a.tst <- predict(rbf.alpha, a.tst.x) %>%
  cbind(a.tst.y) %>%
  data.frame() %>%
  rename(., pred = `.`, obs = a.tst.y) %>%
  mutate(resid = obs - pred)

rbf.a.trn <- predict(rbf.alpha, a.trn.x) %>%
  cbind(a.trn.y) %>%
  data.frame() %>%
  rename(., pred = `.`, obs = a.trn.y) %>%
  mutate(resid = obs - pred)

defaultSummary(rbf.a.tst) # 0.564
defaultSummary(rbf.a.trn) # 0.783

#     Beta -----------------------------------------------------------

rbf.beta <- svm(x = b.trn.x,
                 y = b.trn.y, 
                 kernel = "radial", 
                 cost = 4, 
                 gamma = 0.001, 
                 epsilon = 0.5)

rbf.b.tst <- predict(rbf.beta, b.tst.x) %>%
  cbind(b.tst.y) %>%
  data.frame() %>%
  rename(., pred = `.`, obs = b.tst.y) %>%
  mutate(resid = obs - pred)

rbf.b.trn <- predict(rbf.beta, b.trn.x) %>%
  cbind(b.trn.y) %>%
  data.frame() %>%
  rename(., pred = `.`, obs = b.trn.y) %>%
  mutate(resid = obs - pred)

defaultSummary(rbf.b.tst) # 0.735
defaultSummary(rbf.b.trn) # 0.7996
#     Gamma ----------------------------------------------------------

rbf.gamma <- svm(x = c.trn.x,
                y = c.trn.y, 
                kernel = "radial", 
                cost = 4, 
                gamma = 0.001, 
                epsilon = 0.5)

rbf.c.tst <- predict(rbf.gamma, c.tst.x) %>%
  cbind(c.tst.y) %>%
  data.frame() %>%
  rename(., pred = `.`, obs = c.tst.y) %>%
  mutate(resid = obs - pred)

rbf.c.trn <- predict(rbf.gamma, c.trn.x) %>%
  cbind(c.trn.y) %>%
  data.frame() %>%
  rename(., pred = `.`, obs = c.trn.y) %>%
  mutate(resid = obs - pred)

defaultSummary(rbf.c.tst) # 0.279
defaultSummary(rbf.c.trn) # 0.905

#     Compiled CDs --------------------------------------------------------

temp.a <- rbf.a.tst %>% mutate(cd.type = "alpha")
temp.b <- rbf.b.tst %>% mutate(cd.type = "beta")
temp.c <- rbf.c.tst %>% mutate(cd.type = "gamma")
rbf.abc.tst <- rbind(temp.a, temp.b, temp.c, make.row.names = F)

temp.a <- rbf.a.trn %>% mutate(cd.type = "alpha")
temp.b <- rbf.b.trn %>% mutate(cd.type = "beta")
temp.c <- rbf.c.trn %>% mutate(cd.type = "gamma")
rbf.abc.trn <- rbind(temp.a, temp.b, temp.c, make.row.names = F)

defaultSummary(rbf.abc.tst) # 0.625
defaultSummary(rbf.abc.trn) # 0.806

saveRDS(rbf.alpha, "./models/svm/rbf.alpha.RDS")
saveRDS(rbf.beta, "./models/svm/rbf.beta.RDS")
saveRDS(rbf.gamma, "./models/svm/rbf.gamma.RDS")

# Sigmoid Kernel -----------------------------------------------------

#     All ------------------------------------------------------------

sig.all <- svm(x = trn.x, 
               y = trn.y,
               cost = 1, 
               epsilon = 0.125, 
               kernel = "sigmoid", 
               gamma = 0.001)

sig.all.tst <- predict(sig.all, tst.x) %>%
  cbind(tst.y) %>%
  data.frame() %>%
  rename(., pred = `.`, obs = tst.y) %>%
  mutate(resid = obs - pred)

sig.all.trn <- predict(sig.all, trn.x) %>%
  cbind(trn.y) %>%
  data.frame() %>%
  rename(., pred = `.`, obs = trn.y) %>%
  mutate(resid = obs - pred)

defaultSummary(sig.all.tst) # 0.470
defaultSummary(sig.all.trn) # 0.474

#     Alpha ---------------------------------------------------------

sig.alpha <- svm(x = a.trn.x,
                 y = a.trn.y, 
                 kernel = "sigmoid", 
                 cost = 1, 
                 gamma = 0.001, 
                 epsilon = 0.125)
sig.a.tst <- predict(sig.alpha, a.tst.x) %>%
  cbind(a.tst.y) %>%
  data.frame() %>%
  rename(., pred = `.`, obs = a.tst.y) %>%
  mutate(resid = obs - pred)

sig.a.trn <- predict(sig.alpha, a.trn.x) %>%
  cbind(a.trn.y) %>%
  data.frame() %>%
  rename(., pred = `.`, obs = a.trn.y) %>%
  mutate(resid = obs - pred)

defaultSummary(sig.a.tst) # 0.420
defaultSummary(sig.a.trn) # 0.496

#     Beta -----------------------------------------------------------

sig.beta <- svm(x = b.trn.x,
                y = b.trn.y, 
                kernel = "sigmoid", 
                cost = 1, 
                gamma = 0.001, 
                epsilon = 0.125)

sig.b.tst <- predict(sig.beta, b.tst.x) %>%
  cbind(b.tst.y) %>%
  data.frame() %>%
  rename(., pred = `.`, obs = b.tst.y) %>%
  mutate(resid = obs - pred)

sig.b.trn <- predict(sig.beta, b.trn.x) %>%
  cbind(b.trn.y) %>%
  data.frame() %>%
  rename(., pred = `.`, obs = b.trn.y) %>%
  mutate(resid = obs - pred)

defaultSummary(sig.b.tst) # 0.6997
defaultSummary(sig.b.trn) # 0.635

#     Gamma ----------------------------------------------------------

sig.gamma <- svm(x = c.trn.x,
                 y = c.trn.y, 
                 kernel = "sigmoid", 
                 cost = 1, 
                 gamma = 0.001, 
                 epsilon = 0.125)

sig.c.tst <- predict(sig.gamma, c.tst.x) %>%
  cbind(c.tst.y) %>%
  data.frame() %>%
  rename(., pred = `.`, obs = c.tst.y) %>%
  mutate(resid = obs - pred)

sig.c.trn <- predict(sig.gamma, c.trn.x) %>%
  cbind(c.trn.y) %>%
  data.frame() %>%
  rename(., pred = `.`, obs = c.trn.y) %>%
  mutate(resid = obs - pred)

defaultSummary(sig.c.tst) # 0.222
defaultSummary(sig.c.trn) # 0.737

#     Compiled CDs --------------------------------------------------------

temp.a <- sig.a.tst %>% mutate(cd.type = "alpha")
temp.b <- sig.b.tst %>% mutate(cd.type = "beta")
temp.c <- sig.c.tst %>% mutate(cd.type = "gamma")
sig.abc.tst <- rbind(temp.a, temp.b, temp.c, make.row.names = F)

temp.a <- sig.a.trn %>% mutate(cd.type = "alpha")
temp.b <- sig.b.trn %>% mutate(cd.type = "beta")
temp.c <- sig.c.trn %>% mutate(cd.type = "gamma")
sig.abc.trn <- rbind(temp.a, temp.b, temp.c, make.row.names = F)

defaultSummary(sig.abc.tst) # 0.494
defaultSummary(sig.abc.trn) # 0.575

# Not saving the model because it's the worst of the SVMs
# saveRDS(sig.alpha, "./models/svm/sig.alpha.RDS")
# saveRDS(sig.beta, "./models/svm/sig.beta.RDS")
# saveRDS(sig.gamma, "./models/svm/sig.gamma.RDS")

# Plots -------------------------------------------------------------------

# Notes: Outlier is 3-methylbenzoic acid
# dir.create("./graphs)
dir.create("./graphs/svm")

#     Polynomial ----------------------------------------------------------
#         Test Set --------------------------------------------------------

ggplot(svm.all.tst, aes(x = obs, y = pred)) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0) + 
  theme_bw() + 
  coord_fixed(xlim = c(-45, 5), ylim = c(-45, 5)) + 
  labs(title = "Polynomial SVM", 
       x = "Observed DelG, kJ/mol", 
       y = "Predicted DelG, kJ/mol")
ggsave("./graphs/svm/2017-07-28 polysvm all.png")

# ggplot_build(p)$data
ggplot(svm.a.tst, aes(x = obs, y = pred)) + 
  geom_point(color = "#F8766D") + 
  geom_abline(slope = 1, intercept = 0) + 
  theme_bw() + 
  coord_fixed(xlim = c(-45, 5), ylim = c(-45, 5)) + 
  labs(title = "Polynomial SVM - Alpha-CD", 
       x = "Observed DelG, kJ/mol", 
       y = "Predicted DelG, kJ/mol")
ggsave("./graphs/svm/2017-07-28 polysvm alpha.png")

ggplot(svm.b.tst, aes(x = obs, y = pred)) + 
  geom_point(color = "#00BA38") + 
  geom_abline(slope = 1, intercept = 0) + 
  theme_bw() + 
  coord_fixed(xlim = c(-45, 5), ylim = c(-45, 5)) + 
  labs(title = "Polynomial SVM - Beta-CD", 
       x = "Observed DelG, kJ/mol", 
       y = "Predicted DelG, kJ/mol")
ggsave("./graphs/svm/2017-07-28 polysvm beta.png")

ggplot(svm.c.tst, aes(x = obs, y = pred)) + 
  geom_point(color = "#619CFF") + 
  geom_abline(slope = 1, intercept = 0) + 
  theme_bw() + 
  coord_fixed(xlim = c(-45, 5), ylim = c(-45, 5)) + 
  labs(title = "Polynomial SVM - Gamma-CD", 
       x = "Observed DelG, kJ/mol", 
       y = "Predicted DelG, kJ/mol")
ggsave("./graphs/svm/2017-07-28 polysvm gamma.png")

ggplot(poly.abc.tst, aes(x = obs, y = pred, color = cd.type)) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0) + 
  theme_bw() + 
  coord_fixed(xlim = c(-45, 5), ylim = c(-45, 5)) + 
  labs(title = "Polynomial SVM - Compiled CD Types", 
       x = "Observed DelG, kJ/mol", 
       y = "Predicted DelG, kJ/mol", 
       color = "Cyclodextrin")
ggsave("./graphs/svm/2017-07-28 polysvm compiled.png")

#         Training Set ----------------------------------------------------

ggplot(svm.all.trn, aes(x = obs, y = pred)) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0) + 
  theme_bw() + 
  coord_fixed(xlim = c(-45, 5), ylim = c(-45, 5)) + 
  labs(title = "Polynomial SVM - Training", 
       x = "Observed DelG, kJ/mol", 
       y = "Predicted DelG, kJ/mol")
ggsave("./graphs/svm/2017-07-27 polysvm all trn.png")

# ggplot_build(p)$data
ggplot(svm.a.trn, aes(x = obs, y = pred)) + 
  geom_point(color = "#F8766D") + 
  geom_abline(slope = 1, intercept = 0) + 
  theme_bw() + 
  coord_fixed(xlim = c(-45, 5), ylim = c(-45, 5)) + 
  labs(title = "Polynomial SVM - Alpha-CD, Training", 
       x = "Observed DelG, kJ/mol", 
       y = "Predicted DelG, kJ/mol")
ggsave("./graphs/svm/2017-07-27 polysvm alpha trn.png")

ggplot(svm.b.trn, aes(x = obs, y = pred)) + 
  geom_point(color = "#00BA38") + 
  geom_abline(slope = 1, intercept = 0) + 
  theme_bw() + 
  coord_fixed(xlim = c(-45, 5), ylim = c(-45, 5)) + 
  labs(title = "Polynomial SVM - Beta-CD, Training", 
       x = "Observed DelG, kJ/mol", 
       y = "Predicted DelG, kJ/mol")
ggsave("./graphs/svm/2017-07-27 polysvm beta trn.png")

ggplot(svm.c.trn, aes(x = obs, y = pred)) + 
  geom_point(color = "#619CFF") + 
  geom_abline(slope = 1, intercept = 0) + 
  theme_bw() + 
  coord_fixed(xlim = c(-45, 5), ylim = c(-45, 5)) + 
  labs(title = "Polynomial SVM - Gamma-CD, Training", 
       x = "Observed DelG, kJ/mol", 
       y = "Predicted DelG, kJ/mol")
ggsave("./graphs/svm/2017-07-27 polysvm gamma trn.png")

ggplot(poly.abc.trn, aes(x = obs, y = pred, color = cd.type)) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0) + 
  theme_bw() + 
  coord_fixed(xlim = c(-45, 5), ylim = c(-45, 5)) + 
  labs(title = "Polynomial SVM - Compiled CD Types, Training", 
       x = "Observed DelG, kJ/mol", 
       y = "Predicted DelG, kJ/mol", 
       color = "Cyclodextrin")
ggsave("./graphs/svm/2017-07-28 polysvm compiled trn.png")

#     RBF -----------------------------------------------------------------
#         Test Set --------------------------------------------------------

ggplot(rbf.all.tst, aes(x = obs, y = pred)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Observed DelG, kJ/mol", y = "Predicted DelG, kJ/mol", 
       title = "Radial SVM - All Data Points") + 
  coord_fixed(xlim = c(-45, 5), ylim = c(-45, 5)) + 
  theme_bw()
ggsave("./graphs/svm/2017-07-27 rbf all data points.png")

ggplot(rbf.a.tst, aes(x = obs, y = pred)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Observed DelG, kJ/mol", y = "Predicted DelG, kJ/mol", 
       title = "Radial SVM - Alpha CD") + 
  coord_fixed(xlim = c(-45, 5), ylim = c(-45, 5)) + 
  theme_bw()
ggsave("./graphs/svm/2017-07-27 rbf alpha cd.png")

ggplot(rbf.b.tst, aes(x = obs, y = pred)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Observed DelG, kJ/mol", y = "Predicted DelG, kJ/mol", 
       title = "Radial SVM - Beta CD") + 
  coord_fixed(xlim = c(-45, 5), ylim = c(-45, 5)) + 
  theme_bw()
ggsave("./graphs/svm/2017-07-27 rbf beta cd.png")

ggplot(rbf.c.tst, aes(x = obs, y = pred)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Observed DelG, kJ/mol", y = "Predicted DelG, kJ/mol", 
       title = "Radial SVM - Gamma CD") + 
  coord_fixed(xlim = c(-45, 5), ylim = c(-45, 5)) + 
  theme_bw()
ggsave("./graphs/svm/2017-07-27 rbf gamma cd.png")

ggplot(rbf.abc.tst, aes(x = obs, y = pred, color = cd.type)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Observed DelG, kJ/mol", y = "Predicted DelG, kJ/mol", 
       title = "Radial SVM - Compiled CDs", 
       color = "Cyclodextrin") + 
  coord_fixed(xlim = c(-45, 5), ylim = c(-45, 5)) + 
  theme_bw()
ggsave("./graphs/svm/2017-07-27 rbf compiled cd.png")

#         Training Set ----------------------------------------------------

ggplot(rbf.all.trn, aes(x = obs, y = pred)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Observed DelG, kJ/mol", y = "Predicted DelG, kJ/mol", 
       title = "Radial SVM - All Data, Training") + 
  coord_fixed(xlim = c(-45, 5), ylim = c(-45, 5)) + 
  theme_bw()
ggsave("./graphs/svm/2017-07-27 rbf all data points trn.png")

ggplot(rbf.a.trn, aes(x = obs, y = pred)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Observed DelG, kJ/mol", y = "Predicted DelG, kJ/mol", 
       title = "Radial SVM - Alpha CD, Training") + 
  coord_fixed(xlim = c(-45, 5), ylim = c(-45, 5)) + 
  theme_bw()
ggsave("./graphs/svm/2017-07-27 rbf alpha cdtrn.png")

ggplot(rbf.b.trn, aes(x = obs, y = pred)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Observed DelG, kJ/mol", y = "Predicted DelG, kJ/mol", 
       title = "Radial SVM - Beta CD, Training") + 
  coord_fixed(xlim = c(-45, 5), ylim = c(-45, 5)) + 
  theme_bw()
ggsave("./graphs/svm/2017-07-27 rbf beta cd trn.png")

ggplot(rbf.c.trn, aes(x = obs, y = pred)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Observed DelG, kJ/mol", y = "Predicted DelG, kJ/mol", 
       title = "Radial SVM - Gamma CD, Training") + 
  coord_fixed(xlim = c(-45, 5), ylim = c(-45, 5)) + 
  theme_bw()
ggsave("./graphs/svm/2017-07-27 rbf gamma cd trn.png")

ggplot(rbf.abc.trn, aes(x = obs, y = pred, color = cd.type)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Observed DelG, kJ/mol", y = "Predicted DelG, kJ/mol", 
       title = "Radial SVM - Compiled CDs, Training", 
       color = "Cyclodextrin") + 
  coord_fixed(xlim = c(-45, 5), ylim = c(-45, 5)) + 
  theme_bw()
ggsave("./graphs/svm/2017-07-27 rbf compiled cd trn.png")

#     Sigmoid -------------------------------------------------------------
#         Test Set --------------------------------------------------------

ggplot(sig.all.tst, aes(x = obs, y = pred)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Observed DelG, kJ/mol", y = "Predicted DelG, kJ/mol", 
       title = "Sigmoid SVM - All Data Points") + 
  coord_fixed(xlim = c(-45, 5), ylim = c(-45, 5)) + 
  theme_bw()
ggsave("./graphs/svm/2017-07-27 sig all data points.png")

ggplot(sig.a.tst, aes(x = obs, y = pred)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Observed DelG, kJ/mol", y = "Predicted DelG, kJ/mol", 
       title = "Sigmoid SVM - Alpha CD") + 
  coord_fixed(xlim = c(-45, 5), ylim = c(-45, 5)) + 
  theme_bw()
ggsave("./graphs/svm/2017-07-27 sig alpha cd.png")

ggplot(sig.b.tst, aes(x = obs, y = pred)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Observed DelG, kJ/mol", y = "Predicted DelG, kJ/mol", 
       title = "Sigmoid SVM - Beta CD") + 
  coord_fixed(xlim = c(-45, 5), ylim = c(-45, 5)) + 
  theme_bw()
ggsave("./graphs/svm/2017-07-27 sig beta cd.png")

ggplot(sig.c.tst, aes(x = obs, y = pred)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Observed DelG, kJ/mol", y = "Predicted DelG, kJ/mol", 
       title = "Sigmoid SVM - Gamma CD") + 
  coord_fixed(xlim = c(-45, 5), ylim = c(-45, 5)) + 
  theme_bw()
ggsave("./graphs/svm/2017-07-27 sig gamma cd.png")

ggplot(sig.abc.tst, aes(x = obs, y = pred, color = cd.type)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Observed DelG, kJ/mol", y = "Predicted DelG, kJ/mol", 
       title = "Sigmoid SVM - Compiled CDs", 
       color = "Cyclodextrin") + 
  coord_fixed(xlim = c(-45, 5), ylim = c(-45, 5)) + 
  theme_bw()
ggsave("./graphs/svm/2017-07-27 sig compiled cd.png")

#         Training Set ----------------------------------------------------

ggplot(sig.all.trn, aes(x = obs, y = pred)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Observed DelG, kJ/mol", y = "Predicted DelG, kJ/mol", 
       title = "Sigmoid SVM - All Data, Training") + 
  coord_fixed(xlim = c(-45, 5), ylim = c(-45, 5)) + 
  theme_bw()
ggsave("./graphs/svm/2017-07-27 sig all data points trn.png")

ggplot(sig.a.trn, aes(x = obs, y = pred)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Observed DelG, kJ/mol", y = "Predicted DelG, kJ/mol", 
       title = "Sigmoid SVM - Alpha CD, Training") + 
  coord_fixed(xlim = c(-45, 5), ylim = c(-45, 5)) + 
  theme_bw()
ggsave("./graphs/svm/2017-07-27 sig alpha cdtrn.png")

ggplot(sig.b.trn, aes(x = obs, y = pred)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Observed DelG, kJ/mol", y = "Predicted DelG, kJ/mol", 
       title = "Sigmoid SVM - Beta CD, Training") + 
  coord_fixed(xlim = c(-45, 5), ylim = c(-45, 5)) + 
  theme_bw()
ggsave("./graphs/svm/2017-07-27 sig beta cd trn.png")

ggplot(sig.c.trn, aes(x = obs, y = pred)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Observed DelG, kJ/mol", y = "Predicted DelG, kJ/mol", 
       title = "Sigmoid SVM - Gamma CD, Training") + 
  coord_fixed(xlim = c(-45, 5), ylim = c(-45, 5)) + 
  theme_bw()
ggsave("./graphs/svm/2017-07-27 sig gamma cd trn.png")

ggplot(sig.abc.trn, aes(x = obs, y = pred, color = cd.type)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Observed DelG, kJ/mol", y = "Predicted DelG, kJ/mol", 
       title = "Sigmoid SVM - Compiled CDs, Training", 
       color = "Cyclodextrin") + 
  coord_fixed(xlim = c(-45, 5), ylim = c(-45, 5)) + 
  theme_bw()
ggsave("./graphs/svm/2017-07-27 sig compiled cd trn.png")

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
ev.abc.poly <- rbind(temp.a, temp.b, temp.c) %>%
  mutate(resid = pred - obs)

ev.a <-  predict(rbf.alpha, ext.val.a[ , -1]) %>%
  cbind(ext.val.a[ , 1]) %>% data.frame() %>%
  rename(., pred = `.`, obs = V2)

ev.b <-  predict(rbf.beta, ext.val.b[ , -1]) %>%
  cbind(ext.val.b[ , 1]) %>% data.frame() %>%
  rename(., pred = `.`, obs = V2)

ev.c <-  predict(rbf.gamma, ext.val.c[ , -1]) %>%
  cbind(ext.val.c[ , 1]) %>% data.frame() %>%
  rename(., pred = `.`, obs = V2)

temp.a <- ev.a %>% mutate(cd.type = "alpha")
temp.b <- ev.b %>% mutate(cd.type = "beta")
temp.c <- ev.c %>% mutate(cd.type = "gamma")
# temp.a <- temp.a[-2, ]
ev.abc.rbf <- rbind(temp.a, temp.b, temp.c) %>%
  mutate(resid = pred - obs)

ev.a <-  predict(sig.alpha, ext.val.a[ , -1]) %>%
  cbind(ext.val.a[ , 1]) %>% data.frame() %>%
  rename(., pred = `.`, obs = V2)

ev.b <-  predict(sig.beta, ext.val.b[ , -1]) %>%
  cbind(ext.val.b[ , 1]) %>% data.frame() %>%
  rename(., pred = `.`, obs = V2)

ev.c <-  predict(sig.gamma, ext.val.c[ , -1]) %>%
  cbind(ext.val.c[ , 1]) %>% data.frame() %>%
  rename(., pred = `.`, obs = V2)

temp.a <- ev.a %>% mutate(cd.type = "alpha")
temp.b <- ev.b %>% mutate(cd.type = "beta")
temp.c <- ev.c %>% mutate(cd.type = "gamma")
# temp.a <- temp.a[-2, ]
ev.abc.sig <- rbind(temp.a, temp.b, temp.c) %>%
  mutate(resid = pred - obs)

#     Plots ---------------------------------------------------------------

defaultSummary(ev.abc.poly) # 0.358 normally, 0.595 without outlier
ggplot(ev.abc.poly, aes(x = obs, y = pred, color = cd.type)) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0) + 
  coord_fixed(xlim = c(-45, 5), ylim = c(-45, 5)) + 
  theme_bw() + 
  labs(title = "Polynomial SVM - External Validation", 
       x = "Observed DelG, kJ/mol", y = "Predicted DelG, kJ/mol", 
       color = "Cyclodextrin")
ggsave("./graphs/svm/2017-07-28 polysvm extval.png")

defaultSummary(ev.abc.rbf)
ggplot(ev.abc.rbf, aes(x = pred, y = obs, color = cd.type)) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0) + 
  coord_fixed(xlim = c(-45, 5), ylim = c(-45, 5)) + 
  theme_bw() + 
  labs(title = "Radial SVM - External Validation", 
       x = "Observed DelG, kJ/mol", y = "Predicted DelG, kJ/mol", 
       color = "Cyclodextrin")
ggsave("./graphs/svm/2017-07-28 rbfsvm extval.png")

defaultSummary(ev.abc.sig)
ggplot(ev.abc.sig, aes(x = pred, y = obs, color = cd.type)) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0) + 
  coord_fixed(xlim = c(-45, 5), ylim = c(-45, 5)) + 
  theme_bw() + 
  labs(title = "Radial SVM - External Validation", 
       x = "Observed DelG, kJ/mol", y = "Predicted DelG, kJ/mol", 
       color = "Cyclodextrin")
ggsave("./graphs/svm/2017-07-28 sigsvm extval.png")