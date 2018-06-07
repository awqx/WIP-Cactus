# Genetic algorithm feature selection from caret (Kuhn) used 
# Less descriptors potentially leads to less noise, which means less
# overfitting

# Libraries and packages --------------------------------------------------
install.packages("devtools")
devtools::install_github('topepo/caret/pkg/caret')
library(caret)
library(doParallel)
library(tidyverse)

# Data --------------------------------------------------------------------

set.seed(1001)
data <- readRDS("./data/padel.pp.RDS") %>% select(-data.source)
reduced.desc <- colnames(data) %>% str_detect(., "ATS")
data <- data[ , !reduced.desc]

alpha <- data %>% filter(host == "alpha")
trn <- sample(x = 1:nrow(alpha), size = round(0.8 * nrow(alpha)))
alpha.x <- alpha[trn , -1:-3]
alpha.y <- alpha[trn , 3] 
alpha.x.tst <- alpha[-trn, -1:-3]
alpha.y.tst <- alpha[-trn, 3]
beta <- data %>% filter(host == "beta")
trn <- sample(x = 1:nrow(beta), size = round(0.8 * nrow(beta)))
beta.x <- beta[trn , -1:-3]
beta.y <- beta[trn , 3]
beta.x.tst <- beta[-trn, -1:-3]
beta.y.tst <- beta[-trn, 3]
gamma <- data %>% filter(host == "gamma")
trn <- sample(x = 1:nrow(gamma), size = round(0.8 * nrow(gamma)))
gamma.x <- gamma[trn , -1:-3]
gamma.y <- gamma[trn , 3]
gamma.x.tst <- gamma[-trn, -1:-3]
gamma.y.tst <- gamma[-trn, 3]

dir.create("./feature-selection")
saveRDS(alpha.x, "./feature-selection/alpha.x.RDS")
saveRDS(alpha.y, "./feature-selection/alpha.y.RDS")
saveRDS(beta.x, "./feature-selection/beta.x.RDS")
saveRDS(beta.y, "./feature-selection/beta.y.RDS")
saveRDS(gamma.x, "./feature-selection/gamma.x.RDS")
saveRDS(gamma.y, "./feature-selection/gamma.y.RDS")
saveRDS(alpha.x.tst, "./feature-selection/alpha.x.tst.RDS")
saveRDS(alpha.y.tst, "./feature-selection/alpha.y.tst.RDS")
saveRDS(beta.x.tst, "./feature-selection/beta.x.tst.RDS")
saveRDS(beta.y.tst, "./feature-selection/beta.y.tst.RDS")
saveRDS(gamma.x.tst, "./feature-selection/gamma.x.tst.RDS")
saveRDS(gamma.y.tst, "./feature-selection/gamma.y.tst.RDS")


# GAFS --------------------------------------------------------------------

registerDoParallel(4)
getDoParWorkers()
# From previous investigation, GLMNet works best w/ lambda = 0.362 (alpha), 
# 0.320 (beta), 0.121 (gamma), and alpha = 1 (all)
ga.ctrl <- gafsControl(functions = caretGA, # Assess fitness with caret::train
                       method = "cv",       # 10 fold cross validation
                       genParallel=TRUE,    # Use parallel programming
                       allowParallel = TRUE)
set.seed(1001)
lev <- c("PS", "WS")
system.time(
gafs.a <- gafs(method = "glmnet",
              x = alpha.x, y = as.numeric(alpha.y), 
              lambda = 0.362,
              iters = 100, 
              popSize = 20, 
              gafsControl = ga.ctrl) 
)
train(method = "glmnet", lambda = 0.362,
      x = alpha.x, y = as.numeric(alpha.y))
summary(gafs.a)

system.time(
  gafs.b <- gafs(method = "glmnet",
                 x = beta.x, y = as.numeric(beta.y), 
                 lambda = 0.32,
                 iters = 100, 
                 popSize = 20, 
                 gafsControl = ga.ctrl)
)
