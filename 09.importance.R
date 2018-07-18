# Important note: Only workin with Beta-CD models (for now) due to accuracy
# of those models and time constraint

# Libraries and Packages --------------------------------------------------

library(caret)
library(Cubist)
library(data.table)
library(e1071)
library(glmnet)
library(randomForest)
library(tidyverse)

# Functions ---------------------------------------------------------------

# Rescale from 0 to 1
range01 <- function(x) {(x-min(x))/(max(x)-min(x))}

# Rescale from any range
new.range <- function(x, newMax, newMin) { 
  (x - min(x))/(max(x)-min(x)) * (newMax - newMin) + newMin 
}

# VIP for PLS
# Function sourced from http://mevik.net/work/software/pls.html
VIP <- function(object) {
  if (object$method != "oscorespls")
    stop("Only implemented for orthogonal scores algorithm.  Refit with 'method = \"oscorespls\"'")
  if (nrow(object$Yloadings) > 1)
    stop("Only implemented for single-response models")
  
  SS <- c(object$Yloadings)^2 * colSums(object$scores^2)
  Wnorm2 <- colSums(object$loading.weights^2)
  SSW <- sweep(object$loading.weights^2, 2, SS / Wnorm2, "*")
  sqrt(nrow(SSW) * apply(SSW, 1, cumsum) / cumsum(SS))
}

sort.vip <- function(vip.df) {
  df <- vip.df[1, ] %>% gather()
  colnames(df)[2] <- row.names(vip.df)[1]
  for(r in 2:nrow(vip.df)) {
    long <- vip.df[r, ] %>% gather()
    colnames(long)[2] <- row.names(vip.df)[r]
    df <- inner_join(df, long, by = "key")
  }
  return(df)
}
# Alpha-CD -----------------------------------------------------------------

#     Loading Models ------------------------------------------------------

# setwd("~/SREP LAB/qsar")
cube <- readRDS("./models/alpha/cube.RDS")[[2]]
glm  <- readRDS("./models/alpha/glmnet.RDS")[[2]]
pls  <- readRDS("./models/alpha/pls.RDS")[[2]]
rf   <- readRDS("./models/alpha/rf.RDS")[[2]]
polysvm <- readRDS("./models/alpha/polysvm.RDS")[[2]]
rbfsvm  <- readRDS("./models/alpha/rbfsvm.RDS")[[2]]

#     Variable importance -----------------------------------------------------

# Cubist
cube.imp <- varImp(cube) %>%
  mutate(desc = rownames(.)) %>%
  mutate(model = "cubist")
cube.imp <- cube.imp[order(-cube.imp$Overall), ]

# Random Forest
rf.imp <- importance(rf) %>%
  as.data.frame() %>%
  mutate(desc = rownames(.)) %>%
  mutate(model = "rforest") %>%
  rename(Overall = IncNodePurity)
rf.imp <- rf.imp[order(-rf.imp$Overall), ]

# PLS
# See note in "functions"
pls.vip <- VIP(pls) %>% data.frame() 
pls.imp <- sort.vip(pls.vip)
pls.imp <- pls.imp %>% select(key, `Comp 8`) %>%
  .[order(-pls.imp$`Comp 8`), ] %>%
  mutate(model = "pls") %>%
  rename(Overall = `Comp 8`, desc = key)

# Attempting to get some sort of result by training in caret
# SVM (Polynomial)
trn.alpha <- readRDS("./pre-process/alpha/3/pp.RDS") %>%
  select(., -guest)
features <- readRDS("./feature.selection/alpha.vars.RDS")
colnames(trn.alpha) <- str_replace(colnames(trn.alpha), "-", ".")
trn.alpha <- trn.alpha %>% select(., DelG, features) 
trn.alpha.x <- select(trn.alpha, -DelG) 
trn.alpha.y <- trn.alpha$DelG

polysvm.caret <- train(trn.alpha.x, trn.alpha.y,
                       degree = 3, 
                       method = "svmPoly", 
                       metric = "RMSE")
polysvm.imp <- varImp(polysvm.caret)[[1]] %>% 
  mutate(desc = rownames(.)) %>%
  mutate(model = "polysvm") 
polysvm.imp <- polysvm.imp[order(-polysvm.imp$Overall), ]

# SVM (Radial)
trn.alpha <- readRDS("./pre-process/alpha/3/pp.RDS") %>%
  select(., -guest)
features <- readRDS("./feature.selection/alpha.vars.RDS")
colnames(trn.alpha) <- str_replace(colnames(trn.alpha), "-", ".")
trn.alpha <- trn.alpha %>% select(., DelG, features) 
trn.alpha.x <- select(trn.alpha, -DelG) 
trn.alpha.y <- trn.alpha$DelG

rbfsvm.caret <- train(trn.alpha.x, trn.alpha.y,
                      method = "svmRadial", 
                      metric = "RMSE")
rbfsvm.imp <- varImp(rbfsvm.caret)[[1]] %>% 
  mutate(desc = rownames(.)) %>%
  mutate(model = "rbfsvm") 
rbfsvm.imp <- rbfsvm.imp[order(-rbfsvm.imp$Overall), ]


# Glmnet
trn.alpha <- readRDS("./pre-process/alpha/3/pp.RDS") %>%
  select(., -guest)
features <- readRDS("./feature.selection/alpha.vars.RDS")
colnames(trn.alpha) <- str_replace(colnames(trn.alpha), "-", ".")
trn.alpha <- trn.alpha %>% select(., DelG, features) 
trn.alpha.x <- select(trn.alpha, -DelG) 
trn.alpha.y <- trn.alpha$DelG

glm.caret <- train(trn.alpha.x, trn.alpha.y,
                   lambda = tail(glm$lambda, n = 1),
                   method = "glmnet", metric = "RMSE")
glm.imp <- varImp(glm.caret, lambda = tail(glm$lambda, n = 1))[[1]] %>% 
  mutate(desc = rownames(.)) %>%
  mutate(model = "glmnet")
glm.imp <- glm.imp[order(-glm.imp$Overall), ]

alpha.imp <- rbind(
  cube.imp %>% mutate(Overall = range01(Overall)), 
  glm.imp  %>% mutate(Overall = range01(Overall)), 
  pls.imp  %>% mutate(Overall = range01(Overall)), 
  rf.imp   %>% mutate(Overall = range01(Overall)),
  polysvm.imp  %>% mutate(Overall = range01(Overall)), 
  rbfsvm.imp   %>% mutate(Overall = range01(Overall))
) %>% rename(importance = Overall) %>%
  mutate(host = "alpha")

# Beta-CD -----------------------------------------------------------------

#     Loading Models ------------------------------------------------------

# setwd("~/SREP LAB/qsar")
cube <- readRDS("./models/beta/cube.RDS")[[2]]
glm  <- readRDS("./models/beta/glmnet.RDS")[[2]]
pls  <- readRDS("./models/beta/pls.RDS")[[2]]
rf   <- readRDS("./models/beta/rf.RDS")[[2]]
polysvm <- readRDS("./models/beta/polysvm.RDS")[[2]]
rbfsvm  <- readRDS("./models/beta/rbfsvm.RDS")[[2]]

#     Variable importance -----------------------------------------------------

# Cubist
cube.imp <- varImp(cube) %>%
  mutate(desc = rownames(.)) %>%
  mutate(model = "cubist")
cube.imp <- cube.imp[order(-cube.imp$Overall), ]

# Random Forest
rf.imp <- importance(rf) %>%
  as.data.frame() %>%
  mutate(desc = rownames(.)) %>%
  mutate(model = "rforest") %>%
  rename(Overall = IncNodePurity)
rf.imp <- rf.imp[order(-rf.imp$Overall), ]

# PLS
# See note in "functions"
pls.vip <- VIP(pls) %>% data.frame() 
pls.imp <- sort.vip(pls.vip)
pls.imp <- pls.imp %>% select(key, `Comp 25`) %>%
  .[order(-pls.imp$`Comp 25`), ] %>%
  mutate(model = "pls") %>%
  rename(Overall = `Comp 25`, desc = key)

# Attempting to get some sort of result by training in caret
# SVM (Polynomial)
trn.beta <- readRDS("./pre-process/beta/2/pp.RDS") %>%
  select(., -guest)
features <- readRDS("./feature.selection/beta.vars.RDS")
colnames(trn.beta) <- str_replace(colnames(trn.beta), "-", ".")
trn.beta <- trn.beta %>% select(., DelG, features) 
trn.beta.x <- select(trn.beta, -DelG) 
trn.beta.y <- trn.beta$DelG

polysvm.caret <- train(trn.beta.x, trn.beta.y,
                       degree = 3, 
                       method = "svmPoly", 
                       metric = "RMSE")
polysvm.imp <- varImp(polysvm.caret)[[1]] %>% 
  mutate(desc = rownames(.)) %>%
  mutate(model = "polysvm") 
polysvm.imp <- polysvm.imp[order(-polysvm.imp$Overall), ]

# SVM (Radial)
trn.beta <- readRDS("./pre-process/beta/3/pp.RDS") %>%
  select(., -guest)
features <- readRDS("./feature.selection/beta.vars.RDS")
colnames(trn.beta) <- str_replace(colnames(trn.beta), "-", ".")
trn.beta <- trn.beta %>% select(., DelG, features) 
trn.beta.x <- select(trn.beta, -DelG) 
trn.beta.y <- trn.beta$DelG

rbfsvm.caret <- train(trn.beta.x, trn.beta.y,
                      method = "svmRadial", 
                      metric = "RMSE")
rbfsvm.imp <- varImp(rbfsvm.caret)[[1]] %>% 
  mutate(desc = rownames(.)) %>%
  mutate(model = "rbfsvm") 
rbfsvm.imp <- rbfsvm.imp[order(-rbfsvm.imp$Overall), ]


# Glmnet
trn.beta <- readRDS("./pre-process/beta/5/pp.RDS") %>%
  select(., -guest)
features <- readRDS("./feature.selection/beta.vars.RDS")
colnames(trn.beta) <- str_replace(colnames(trn.beta), "-", ".")
trn.beta <- trn.beta %>% select(., DelG, features) 
trn.beta.x <- select(trn.beta, -DelG) 
trn.beta.y <- trn.beta$DelG

glm.caret <- train(trn.beta.x, trn.beta.y,
                   lambda = tail(glm$lambda, n = 1),
                   method = "glmnet", metric = "RMSE")
glm.imp <- varImp(glm.caret, lambda = tail(glm$lambda, n = 1))[[1]] %>% 
  mutate(desc = rownames(.)) %>%
  mutate(model = "glmnet")
glm.imp <- glm.imp[order(-glm.imp$Overall), ]

# Compiling everything
beta.imp <- rbind(
  cube.imp %>% mutate(Overall = range01(Overall)), 
  glm.imp  %>% mutate(Overall = range01(Overall)), 
  pls.imp  %>% mutate(Overall = range01(Overall)), 
  rf.imp   %>% mutate(Overall = range01(Overall)),
  polysvm.imp  %>% mutate(Overall = range01(Overall)), 
  rbfsvm.imp   %>% mutate(Overall = range01(Overall))
) %>% rename(importance = Overall) %>% 
  mutate(host = "beta")

# Data compiling ----------------------------------------------------------

comb.desc <- inner_join(alpha.imp, beta.imp, by = 'desc') %>%
  .$desc %>% unique()
comb.imp <- rbind(alpha.imp %>% filter(desc %in% comb.desc), 
                  beta.imp %>% filter(desc %in% comb.desc))

# Creating wide tables to visualize averages
alpha.dt <- alpha.imp %>% select(-host) %>% 
  mutate(desc = as.factor(desc)) %>%
  mutate(model = as.factor(model)) %>%
  data.table()
alpha.dt <- dcast.data.table(alpha.dt, desc ~ model, 
                             value.var = "importance")
alpha.dt$cubist <-  replace(alpha.dt$cubist, is.na(alpha.dt$cubist), 0)
alpha.dt <- alpha.dt %>% mutate(overall = 
                                  rowMeans(select(., cubist:rforest))) %>%
  data.frame()

beta.dt <- beta.imp %>% select(-host) %>% 
  mutate(desc = as.factor(desc)) %>%
  mutate(model = as.factor(model)) %>%
  data.table()
beta.dt <- dcast.data.table(beta.dt, desc ~ model, 
                            value.var = "importance")
beta.dt$cubist <-  replace(beta.dt$cubist, is.na(beta.dt$cubist), 0)
beta.dt$rforest <-  replace(beta.dt$rforest, is.na(beta.dt$rforest), 0)
beta.dt <- beta.dt %>% mutate(overall = 
                                rowMeans(select(., cubist:rforest))) %>%
  data.frame()

# Adding a column to the original importance table of 'ensemble'
# (basically the average importance)
alpha.ens <- alpha.dt %>% select(desc, overall) %>%
  mutate(model = 'ensemble', host = 'alpha') %>%
  rename(importance = overall) 
alpha.imp <- rbind(alpha.imp, alpha.ens) %>%
  mutate(desc = as.factor(desc), model = as.factor(model))

beta.ens <- beta.dt %>% select(desc, overall) %>%
  mutate(model = 'ensemble', host = 'beta') %>%
  rename(importance = overall) 
beta.imp <- rbind(beta.imp, beta.ens) %>%
  mutate(desc = as.factor(desc), model = as.factor(model))

dir.create('var.imp')
saveRDS(alpha.imp, 'var.imp/alpha.varimp.RDS')
saveRDS(beta.imp, 'var.imp/beta.varimp.RDS')
saveRDS(comb.imp, 'var.imp/comb.varimp.RDS')
saveRDS(alpha.dt, 'var.imp/alpha.wide.RDS')
saveRDS(beta.dt, 'var.imp/beta.wide.RDS')

# Graphs ------------------------------------------------------------------

alpha.order <- alpha.imp %>% filter(model == 'ensemble') %>%
  .[order(.$importance), ] %>% .$desc
alpha.imp$desc <- factor(alpha.imp$desc, levels = alpha.order)
alpha.imp$model <- factor(alpha.imp$model, 
                          levels = c('cubist', 'glmnet', 'pls', 
                                     'polysvm', 'rbfsvm', 'rforest', 
                                     'ensemble'))
alpha.temp <- alpha.imp[!is.na(alpha.imp$desc), ]
ggplot(alpha.temp, aes(x = model, y = desc, fill = importance)) + 
  theme.plos +
  theme(text = element_text(size=11), 
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10)) +
  geom_tile() + 
  scale_fill_gradient2(high = "#f8766d", low = "#619cff",
                       mid = "white", midpoint = 0.5,
                       name = "Importance") +
  # scale_fill_gradientn(colors = terrain.colors(8), 
  #                      name = "Importance")
  scale_x_discrete(labels = c("Cubist", "GLMNet", "PLS", "Poly-SVM",
                              "RBF-SVM", "RF", "Ensemble")) +
  labs(x = "Model", y = "Descriptor Variable", 
       title = "Alpha-CD") +
  coord_fixed(ratio = 0.32)
ggsave("./graphs/alpha.varimp.png", scale = 1, dpi = 600)


beta.order <- beta.imp %>% filter(model == 'ensemble') %>%
  .[order(.$importance), ] %>% .$desc
beta.imp$desc <- factor(beta.imp$desc, levels = beta.order)
beta.imp$model <- factor(beta.imp$model, 
                          levels = c('cubist', 'glmnet', 'pls', 
                                     'polysvm', 'rbfsvm', 'rforest', 
                                     'ensemble'))
beta.temp <- beta.imp[!is.na(beta.imp$desc), ]
ggplot(beta.temp, aes(x = model, y = desc, fill = importance)) + 
  theme.plos +
  theme(text = element_text(size=11), 
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10)) +
  geom_tile() + 
  scale_fill_gradient2(high = "#f8766d", low = "#619cff",
                       mid = "white", midpoint = 0.5) + 
  scale_x_discrete(labels = c("Cubist", "GLMNet", "PLS", "Poly-SVM",
                              "RBF-SVM", "RF", "Ensemble")) +
  labs(x = "Model", y = "Descriptor Variable", fill = "Importance", 
       title = "Beta-CD") +
  coord_fixed(ratio = 0.36)
ggsave("./graphs/beta.varimp.png", scale = 1, dpi = 600)

# Variables shared between alpha- and beta-CD
# Not super helpful, because they share only 7 variables
ggplot(comb.imp, aes(x = model, y = desc, fill = importance)) + 
  theme.plos +
  theme(text=element_text(size=12)) +
  geom_tile() + 
  scale_fill_gradient2(high = "#f8766d", low = "#619cff",
                       mid = "white", midpoint = 0.5) + 
  scale_x_discrete(labels = c("Cubist", "GLMNet", "PLS", "Poly-SVM", 
                              "RBF-SVM", "Random forest")) + 
  labs(x = "Model", y = "Descriptor Variable", fill = "Importance") +
  coord_fixed(ratio = 0.7) + 
  facet_grid(.~host)
ggsave("./graphs/comb.varimp.png", scale = 1, dpi = 600)
