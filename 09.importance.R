# Important note: Only workin with Beta-CD models (for now) due to accuracy
# of those models and time constraint
dir.create('var.imp')
# Libraries and Packages --------------------------------------------------

if(!require('pacman')) { 
  install.packages('pacman')
  library(pacman)
} else
  library(pacman)

# pacman::p_load
p_load(caret, Cubist, e1071, earth, kernlab, randomForest, tidyverse)

# Functions ---------------------------------------------------------------

# Rescale from 0 to 1
range01 <- function(x) {(x-min(x))/(max(x)-min(x))}

# Rescale from any range
new.range <- function(x, newMax, newMin) { 
  (x - min(x))/(max(x)-min(x)) * (newMax - newMin) + newMin 
}

# VIP for PLS
# Function sourced from http://mevik.net/work/software/pls.html
# VIP <- function(object) {
#   if (object$method != 'oscorespls')
#     stop('Only implemented for orthogonal scores algorithm.  Refit with 'method = \'oscorespls\''')
#   if (nrow(object$Yloadings) > 1)
#     stop('Only implemented for single-response models')
#   
#   SS <- c(object$Yloadings)^2 * colSums(object$scores^2)
#   Wnorm2 <- colSums(object$loading.weights^2)
#   SSW <- sweep(object$loading.weights^2, 2, SS / Wnorm2, '*')
#   sqrt(nrow(SSW) * apply(SSW, 1, cumsum) / cumsum(SS))
# }

# sort.vip <- function(vip.df) {
#   df <- vip.df[1, ] %>% gather()
#   colnames(df)[2] <- row.names(vip.df)[1]
#   for(r in 2:nrow(vip.df)) {
#     long <- vip.df[r, ] %>% gather()
#     colnames(long)[2] <- row.names(vip.df)[r]
#     df <- inner_join(df, long, by = 'key')
#   }
#   return(df)
# }

read.models <- function(read.dir) {
  name.models <- list.files('models/alpha') %>% 
    str_remove('.RDS$')
  path.models <- paste0(read.dir, list.files(read.dir))
  
  for(n in 1:length(name.models))
    assign( # Creates variable names from characters
      name.models[n], 
      readRDS(path.models[n])[[2]], 
      envir = .GlobalEnv
      )
}
# Alpha-CD -----------------------------------------------------------------

read.models('models/alpha/')
# name.models <- list.files('models/alpha') %>% str_remove('.RDS$')
# path.models <- paste0('models/alpha/', list.files('models/alpha'))
# for(n in 1:length(name.models))
#   assign(name.models[n], 
#          readRDS(path.models[n])[[2]])

#     Variable importance -----------------------------------------------------

# Cubist
cube.imp <- varImp(cube) %>%
  mutate(desc = rownames(.)) %>%
  rename(importance = Overall) %>%
  mutate(model = 'Cubist')
cube.imp <- cube.imp[order(-cube.imp$importance), ]

# GBM (gradient boosted)
gbm.imp <- summary.gbm(gbm) %>% data.frame() %>%
  rename(desc = var, importance = rel.inf) %>%
  mutate(model = 'GBM') %>%
  remove_rownames() # A cosmetic change
# MARS (built with earth)
mars.evimp <- evimp(mars, trim = F)
mars.imp <- data.frame(
  desc = row.names(mars.evimp) %>% 
    str_remove('-unused'), # A quirk of the evimp function
  # importance can be evaluated as either GCV ([ , 4]) or RSS
  importance = mars.evimp[ , 6]) %>% # This grabs the RSS value 
  mutate(model = 'MARS') %>%
  remove_rownames() 

# Random Forest
rf.imp <- importance(rf) %>%
  as.data.frame() %>%
  mutate(desc = rownames(.)) %>%
  mutate(model = 'RF') %>%
  rename(importance = IncNodePurity)
rf.imp <- rf.imp[order(-rf.imp$importance), ]

# Attempting to get some sort of result by training in caret
# SVM (Polynomial)
trn.alpha <- readRDS('pre-process/alpha/6/pp.RDS') %>%
  select(., -guest)
features <- readRDS('feature.selection/alpha.vars.RDS')
colnames(trn.alpha) <- str_replace(colnames(trn.alpha), '-', '.')
trn.alpha <- trn.alpha %>% select(., DelG, features) 
trn.alpha.x <- select(trn.alpha, -DelG) 
trn.alpha.y <- trn.alpha$DelG

polysvm.caret <- train(trn.alpha.x, trn.alpha.y,
                       degree = 3, 
                       method = 'svmPoly', 
                       metric = 'RMSE')
polysvm.imp <- varImp(polysvm.caret)[[1]] %>% 
  mutate(desc = rownames(.), model = 'SVM-Poly') %>%
  rename(importance = Overall) 
polysvm.imp <- polysvm.imp[order(-polysvm.imp$importance), ]

# SVM (Radial)
trn.alpha <- readRDS('pre-process/alpha/10/pp.RDS') %>%
  select(., -guest)
features <- readRDS('feature.selection/alpha.vars.RDS')
colnames(trn.alpha) <- str_replace(colnames(trn.alpha), '-', '.')
trn.alpha <- trn.alpha %>% select(., DelG, features) 
trn.alpha.x <- select(trn.alpha, -DelG) 
trn.alpha.y <- trn.alpha$DelG

rbfsvm.caret <- train(trn.alpha.x, trn.alpha.y,
                      method = 'svmRadial', 
                      metric = 'RMSE')
rbfsvm.imp <- varImp(rbfsvm.caret)[[1]] %>% 
  mutate(desc = rownames(.), model = 'SVM-RBF') %>%
  rename(importance = Overall)
rbfsvm.imp <- rbfsvm.imp[order(-rbfsvm.imp$importance), ]

#     Compiling/handling data ---------------------------------------------

# Simple rbind of the importance data frames
alpha.imp <- 
  do.call(rbind,
          lapply(list(cube.imp, gbm.imp, mars.imp,
                      rf.imp, polysvm.imp, rbfsvm.imp),
                 function(x)
                   mutate(x, importance = range01(importance))
                 )) %>%
  mutate(host = 'Alpha-CD')

# Creating a wide table to visualize averages
alpha.dt <- alpha.imp %>% select(-host) %>% 
  mutate(desc = as.factor(desc), model = as.factor(model)) %>%
  data.table()
alpha.dt <- dcast.data.table(alpha.dt, desc ~ model, 
                             value.var = 'importance')
alpha.dt <- alpha.dt %>%
  mutate(overall =
           rowMeans(select(., `Cubist`:`SVM-RBF`))) %>%
  data.frame()

# Adding a column to the original importance table of 'ensemble'
# (basically the average importance)
alpha.ens <- alpha.dt %>% select(desc, overall) %>%
  mutate(model = 'Ensemble', host = 'alpha') %>%
  rename(importance = overall) 
alpha.imp <- rbind(alpha.imp, alpha.ens) %>%
  mutate(desc = as.factor(desc), model = as.factor(model))

# Creating a vector of the most importance variables
alpha.order <- alpha.imp %>% 
  filter(model == 'Ensemble') %>%
  # descending = T orders the variables correctly in the dataframe, 
  # but reverses the order when plotted
  .[order(.$importance), ] %>% 
  .$desc
# Coercing the factor levels 
alpha.imp$desc <- factor(alpha.imp$desc, levels = alpha.order)
# Re-ordering the factor levels of 'model' so ensemble is last
print(levels(alpha.imp$model)) # check the names
alpha.imp$model <- factor(alpha.imp$model, 
                          levels = c('Cubist', 'GBM', 'MARS', 'RF',
                                     'SVM-Poly', 'SVM-RBF', 'Ensemble'))
alpha.temp <- alpha.imp[!is.na(alpha.imp$desc), ] # Spot check
# Plotting a tile graph
graph.alpha <- 
  ggplot(alpha.temp, aes(x = model, y = desc, fill = importance)) + 
    # theme.plos +
    theme_bw() +
    theme(text = element_text(size=11), 
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 10)) +
    geom_tile() + 
    scale_fill_gradient2(high = '#f8766d', # red
                         low = '#619cff', # blue
                         mid = 'white', 
                         midpoint = 0.5, 
                         name = 'Importance') +
    # scale_fill_gradientn(colors = terrain.colors(8), 
    #                      name = 'Importance')
    scale_x_discrete(labels = c('Cubist', 'GBM', 'MARS', 'RF',
                                'SVM-Poly', 'SVM-RBF', 'Ensemble')) +
    labs(x = 'Model', y = 'Predictor', 
         title = 'Variable importance for alpha-CD models') 
print(graph.alpha)
ggsave('graphs/2018-08 alpha.varimp.png', scale = 1, dpi = 600)

saveRDS(alpha.imp, 'var.imp/alpha.varimp.RDS')
saveRDS(alpha.dt, 'var.imp/alpha.wide.RDS')

# Beta-CD -----------------------------------------------------------------

# Explantory comments can be found in the section 'Alpha-CD'
read.models('models/beta/')

#     Variable importance -----------------------------------------------------

# Cubist
cube.imp <- varImp(cube) %>%
  mutate(desc = rownames(.)) %>%
  rename(importance = Overall) %>%
  mutate(model = 'Cubist')
cube.imp <- cube.imp[order(-cube.imp$importance), ]

# GBM (gradient boosted)
gbm.imp <- summary.gbm(gbm) %>% data.frame() %>%
  rename(desc = var, importance = rel.inf) %>%
  mutate(model = 'GBM') %>% 
  remove_rownames() 

# MARS (built with earth)
mars.evimp <- evimp(mars, trim = F)
mars.imp <- data.frame(
  desc = row.names(mars.evimp) %>% 
    str_remove('-unused'), 
  importance = mars.evimp[ , 6]) %>% 
  mutate(model = 'MARS') %>%
  remove_rownames() 

# Random Forest
rf.imp <- importance(rf) %>%
  as.data.frame() %>%
  mutate(desc = rownames(.)) %>%
  mutate(model = 'RF') %>%
  rename(importance = IncNodePurity)
rf.imp <- rf.imp[order(-rf.imp$importance), ]

# SVM (Polynomial)
trn.beta <- readRDS('pre-process/beta/2/pp.RDS') %>%
  select(., -guest)
features <- readRDS('feature.selection/beta.vars.RDS')
colnames(trn.beta) <- str_replace(colnames(trn.beta), '-', '.')
trn.beta <- trn.beta %>% select(., DelG, features) 
trn.beta.x <- select(trn.beta, -DelG) 
trn.beta.y <- trn.beta$DelG

polysvm.caret <- train(trn.beta.x, trn.beta.y,
                       degree = 3, 
                       method = 'svmPoly', 
                       metric = 'RMSE')
polysvm.imp <- varImp(polysvm.caret)[[1]] %>% 
  mutate(desc = rownames(.), model = 'SVM-Poly') %>%
  rename(importance = Overall) 
polysvm.imp <- polysvm.imp[order(-polysvm.imp$importance), ]

# SVM (Radial)
trn.beta <- readRDS('pre-process/beta/7/pp.RDS') %>%
  select(., -guest)
features <- readRDS('feature.selection/beta.vars.RDS')
colnames(trn.beta) <- str_replace(colnames(trn.beta), '-', '.')
trn.beta <- trn.beta %>% select(., DelG, features) 
trn.beta.x <- select(trn.beta, -DelG) 
trn.beta.y <- trn.beta$DelG

rbfsvm.caret <- train(trn.beta.x, trn.beta.y,
                      method = 'svmRadial', 
                      metric = 'RMSE')
rbfsvm.imp <- varImp(rbfsvm.caret)[[1]] %>% 
  mutate(desc = rownames(.), model = 'SVM-RBF') %>%
  rename(importance = Overall)
rbfsvm.imp <- rbfsvm.imp[order(-rbfsvm.imp$importance), ]

#     Compiling/handling data ---------------------------------------------

# compiling models
beta.imp <- 
  do.call(rbind,
          lapply(list(cube.imp, gbm.imp, mars.imp,
                      rf.imp, polysvm.imp, rbfsvm.imp),
                 function(x)
                   mutate(x, importance = range01(importance))
          )) %>%
  mutate(host = 'Beta-CD')

# Creating a wide data.table
beta.dt <- beta.imp %>% select(-host) %>% 
  mutate(desc = as.factor(desc), model = as.factor(model)) %>%
  data.table()
beta.dt <- dcast.data.table(beta.dt, desc ~ model, 
                             value.var = 'importance')
beta.dt <- beta.dt %>%
  mutate(overall =
           rowMeans(select(., `Cubist`:`SVM-RBF`))) %>%
  data.frame()

# Creating an ensemble column
beta.ens <- beta.dt %>% select(desc, overall) %>%
  mutate(model = 'Ensemble', host = 'Beta-CD') %>%
  rename(importance = overall) 
beta.imp <- rbind(beta.imp, beta.ens) %>%
  mutate(desc = as.factor(desc), model = as.factor(model))

# Ordering by importance overall
beta.order <- beta.imp %>% 
  filter(model == 'Ensemble') %>%
  .[order(.$importance), ] %>% 
  .$desc
beta.imp$desc <- factor(beta.imp$desc, levels = beta.order)
beta.imp$model <- factor(beta.imp$model, 
                          levels = c('Cubist', 'GBM', 'MARS', 'RF',
                                     'SVM-Poly', 'SVM-RBF', 'Ensemble'))
beta.temp <- beta.imp[!is.na(beta.imp$desc), ] 

# Tile graph
graph.beta <- 
  ggplot(beta.temp, aes(x = model, y = desc, fill = importance)) + 
  # theme.plos +
    theme_bw() +
    theme(text = element_text(size=11), 
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 10)) +
    geom_tile() + 
    scale_fill_gradient2(high = '#f8766d', # red
                         low = '#619cff', # blue
                         mid = 'white', 
                         midpoint = 0.5, 
                         name = 'Importance') +
    # scale_fill_gradientn(colors = terrain.colors(8), 
    #                      name = 'Importance')
    scale_x_discrete(labels = c('Cubist', 'GBM', 'MARS', 'RF',
                                'SVM-Poly', 'SVM-RBF', 'Ensemble')) +
    labs(x = 'Model', y = 'Predictor', 
         title = 'Variable importance for beta-CD models') 
print(graph.beta)
ggsave('graphs/2018-08 beta.varimp.png', scale = 1, dpi = 600)

saveRDS(beta.imp, 'var.imp/beta.varimp.RDS')
saveRDS(beta.dt, 'var.imp/beta.wide.RDS')

# Other graphs ------------------------------------------------------------

# Saving the graphs as 'A' and 'B' to join them together later
graph.a <- 
  graph.alpha + 
    labs(title = 'A') + 
    scale_fill_gradient2(
      high = '#f8766d', low = '#619cff', 
      mid = 'white', midpoint = 0.5,
      guide = F # removing the legend to make a cleaner join
      )
ggsave('graphs/a.varimp.png', graph.a, 
       height = 6.5, width = 5.5, dpi = 600)
graph.b <- 
  graph.beta + 
    labs(title = 'B', y = NULL, x = NULL)
ggsave('graphs/b.varimp.png', graph.b, 
       height = 6.5, width = 6.5, dpi = 600)

# Obsolete QSARs ----------------------------------------------------------

# # PLS
# # See note in 'functions'
# pls.vip <- VIP(pls) %>% data.frame() 
# pls.imp <- sort.vip(pls.vip)
# pls.imp <- pls.imp %>% select(key, `Comp 8`) %>%
#   .[order(-pls.imp$`Comp 8`), ] %>%
#   mutate(model = 'pls') %>%
#   rename(Overall = `Comp 8`, desc = key)

# Glmnet
# trn.alpha <- readRDS('pre-process/alpha/3/pp.RDS') %>%
#   select(., -guest)
# features <- readRDS('feature.selection/alpha.vars.RDS')
# colnames(trn.alpha) <- str_replace(colnames(trn.alpha), '-', '.')
# trn.alpha <- trn.alpha %>% select(., DelG, features) 
# trn.alpha.x <- select(trn.alpha, -DelG) 
# trn.alpha.y <- trn.alpha$DelG
# 
# glm.caret <- train(trn.alpha.x, trn.alpha.y,
#                    lambda = tail(glm$lambda, n = 1),
#                    method = 'glmnet', metric = 'RMSE')
# glm.imp <- varImp(glm.caret, lambda = tail(glm$lambda, n = 1))[[1]] %>% 
#   mutate(desc = rownames(.)) %>%
#   mutate(model = 'glmnet')
# glm.imp <- glm.imp[order(-glm.imp$Overall), ]
