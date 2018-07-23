# MARS - Multivariate Area Regression 
dir.create("./models/mars")

# Libraries ---------------------------------------------------------------

if(!require("pacman"))
  install.packages("pacman")
library(pacman)
p_load(caret, tidyverse, earth)
source("./07.model.functions.R")

# Functions ---------------------------------------------------------------


# Test --------------------------------------------------------------------

trn.alpha <- readRDS("./pre-process/alpha/1/pp.RDS") %>%
  select(., -guest)
features <- readRDS("./feature.selection/alpha.vars.RDS")
colnames(trn.alpha) <- str_replace(colnames(trn.alpha), "-", ".")
trn.alpha <- trn.alpha %>% select(., DelG, features) 
trn.alpha.x <- select(trn.alpha, -DelG) 
trn.alpha.y <- trn.alpha$DelG

# Tuning parameters
# degree: 1-10, int
# penalty: 0+ (preferably < 10)
# thresh: 0-0.1
# minspan: 0-number of cases, int
# fast.k: 0, 5, 10, 20         
mars <- earth(x = trn.alpha.x, 
              y = trn.alpha.y, 
              degree = 2, 
              penalty = 2, 
              thresh = 0.01, 
              fast.k = 5)
temp <- predict(mars, trn.alpha.x) %>% 
  data.frame(obs = trn.alpha.y, .) %>%
  rename(pred = trn.alpha.y)

ggplot(temp, aes(x = obs, y = pred)) + 
  geom_point() + 
  theme_light() + 
  coord_fixed() + 
  geom_abline(slope = 1, intercept = 0)
     
  