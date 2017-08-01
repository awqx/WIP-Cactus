# Libraries and Packages --------------------------------------------------

library(pls)
library(tidyverse)

# Loading Data ------------------------------------------------------------

# setwd("~/SREP LAB/qsar")
dir.create("./models/pls")
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


