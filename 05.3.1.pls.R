# Libraries and Packages --------------------------------------------------

library(pls)
library(tidyverse)

# Loading Data ------------------------------------------------------------

# setwd("~/SREP LAB/qsar")
# dir.create("./models/pls")
df.raw <- readRDS("./padel.pp.new.RDS")
df <- df.raw %>% select(-guest, -host, -data.source)

set.seed(10)
trn.ind <- sample(x = 1:nrow(df), 
                  size = round(0.8 * nrow(df)))
trn <- df[trn.ind, ]
tst <- df[-trn.ind, ]

# All Predictors ----------------------------------------------------------

# Creating a basic PLS model
pls <- plsr(DelG ~ ., 
                ncomp = 11, 
                data = trn, 
                validation = "LOO", 
                method = "oscorespls")
summary(pls)
# Plotting number of components vs. RMSEP
plot(RMSEP(pls), legendpos = "topright") # 3 comps best
# Plotting one of the validation sets
# plot(pls, ncomp = 3, asp = 1, line = T)
# Test data prediction plot
# predict(pls.mod, ncomp = 3, newdata = pls.tst) %>%
#   cbind(pls.tst[ , 1]) %>%
#   data.frame() %>%
#   ggplot(., aes(x = ., y = V2)) + 
#   geom_point() + 
#   geom_abline(slope = 1, intercept = 0) + 
#   coord_fixed()
# Testing R-squared
pls.tst <- predict(pls, ncomp = 11, newdata = tst) %>%
  cbind(tst[ , 1]) %>%
  data.frame() %>%
  dplyr::rename(., pred = `.`, obs = V2)
defaultSummary(pls.tst) # R-squared = 0.654

pls.trn <- predict(pls, ncomp = 11, newdata = trn) %>%
  cbind(trn[ , 1]) %>%
  data.frame() %>%
  dplyr::rename(., pred = `.`, obs = V2)
defaultSummary(pls.trn) # R-squared = 0.758

saveRDS(pls, "./models/pls/pls.all.RDS")

# Separate CD -------------------------------------------------------------

#     Alpha-CD ------------------------------------------------------------

alpha <- df %>% filter(alpha > 0)
set.seed(10)
trn.ind <- sample(x = 1:nrow(alpha), size = round(0.7 * nrow(alpha)))
a.trn <- alpha[trn.ind, ]
a.tst <- alpha[-trn.ind, ]

pls.alpha <- plsr(DelG ~ ., 
                 ncomp = 4, 
                 data = a.trn, 
                 validation = "LOO", 
                 method = "oscorespls", 
                 seed = 10)

pls.tst.a <- predict(pls.alpha, ncomp = 4, newdata = a.tst) %>%
  cbind(a.tst[ , 1]) %>%
  data.frame() %>%
  dplyr::rename(., pred = `.`, obs = V2)
defaultSummary(pls.tst.a) # R-squared = 0.566

#     Beta-CD -------------------------------------------------------------

beta <- df %>% filter(beta > 0)
set.seed(10)
trn.ind <- sample(x = 1:nrow(beta), size = round(0.7 * nrow(beta)))
b.trn <- beta[trn.ind, ]
b.tst <- beta[-trn.ind, ]

pls.beta <- plsr(DelG ~ ., 
            ncomp = 4, 
            data = b.trn, 
            validation = "LOO", 
            method = "oscorespls", 
            seed = 10)

pls.tst.b <- predict(pls.beta, ncomp = 4, newdata = b.tst) %>%
  cbind(b.tst[ , 1]) %>%
  data.frame() %>%
  dplyr::rename(., pred = `.`, obs = V2)
defaultSummary(pls.tst.b) # R-squared = 0.755

saveRDS(pls.beta, "./models/pls/pls.beta.RDS")

# Gamma-CD ----------------------------------------------------------------

gamma <- df %>% filter(gamma > 0)
set.seed(10)
trn.ind <- sample(x = 1:nrow(gamma), size = round(0.7 * nrow(gamma)))
c.trn <- gamma[trn.ind, ]
c.tst <- gamma[-trn.ind, ]

pls.gamma <- plsr(DelG ~ ., 
                 ncomp = 5, 
                 data = c.trn, 
                 validation = "LOO", 
                 method = "oscorespls", 
                 seed = 10)

pls.tst.c <- predict(pls.gamma, ncomp = 5, newdata = c.tst) %>%
  cbind(c.tst[ , 1]) %>%
  data.frame() %>%
  dplyr::rename(., pred = `.`, obs = V2)
defaultSummary(pls.tst.c) # R-squared = 0.319

#     Compiled-CD ---------------------------------------------------------

temp.a <- pls.tst.a %>% mutate(cd.type = "alpha")
temp.b <- pls.tst.b %>% mutate(cd.type = "beta")
temp.c <- pls.tst.c %>% mutate(cd.type = "gamma")
pls.abc.tst <- rbind(temp.a, temp.b, temp.c, 
                     make.row.names = F) %>%
  mutate(residual = pred - obs)
defaultSummary(glm.abc.tst) # 0.708

ggplot(pls.abc.tst, aes(x = obs, y = pred, color = cd.type)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1)

# Saving Models -----------------------------------------------------------

saveRDS(pls.alpha, "./models/pls/pls.alpha.RDS")
saveRDS(pls.beta, "./models/pls/pls.beta2.RDS")
saveRDS(pls.gamma, "./models/pls/pls.gamma.RDS")

saveRDS(pls.abc.tst, "./models/pls/pls.results.RDS")

# External Validation -----------------------------------------------------

ext.val <- readRDS("./external validation set new.RDS") %>%
  select(-guest:-host) %>%
  select(-data.source) 

ext.val.a <- ext.val %>% filter(alpha > 0) 
ext.val.b <- ext.val %>% filter(beta > 0)
ext.val.c <- ext.val %>% filter(gamma > 0) 

ev.a <-  predict(pls.alpha, ncomp = 5, newdata = ext.val.a) %>%
  cbind(ext.val.a[ , 1]) %>% data.frame() %>%
  dplyr::rename(., pred = `.`, obs = V2)

ev.b <-  predict(pls.beta, ncomp = 5, newdata = ext.val.b) %>%
  cbind(ext.val.b[ , 1]) %>% data.frame() %>%
  dplyr::rename(., pred = `.`, obs = V2)

ev.c <-  predict(pls.gamma, ncomp = 5, newdata = ext.val.c) %>%
  cbind(ext.val.c[ , 1]) %>% data.frame() %>%
  dplyr::rename(., pred = `.`, obs = V2)

temp.a <- ev.a %>% mutate(cd.type = "alpha")
temp.b <- ev.b %>% mutate(cd.type = "beta")
temp.c <- ev.c %>% mutate(cd.type = "gamma")
# temp.a <- temp.a[-2, ]
ev.abc <- rbind(temp.a, temp.b, temp.c) %>%
  mutate(resid = pred - obs)
# Prearing copy for compilation
ev.pls <- ev.abc %>% select(-resid) %>%
  mutate(Model = "PLS")
