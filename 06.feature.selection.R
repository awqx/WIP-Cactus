source('data.handling.R')
dir.create('feature.selection')
dir.create('feature.selection/alpha')
dir.create('feature.selection/beta')
dir.create('feature.selection/gamma')

# Recursive Feature Elimination (RFE) -------------------------------------

# registerDoParallel(4)
# getDoParWorkers() # Confirming 4 cores

rfe.combos <- expand.grid(c("alpha", "beta", "gamma"), c(1:10)) %>%
  rename(cd.type = Var1, num = Var2)
cd.combos <- rfe.combos$cd.type
num.combos <- rfe.combos$num

# If Error: worker initialization failed occurs, re-run mapply
mapply(FUN = save.rfe, 
       cd.type = cd.combos, num = num.combos, 
       pp.dir = "./pre-process", write.dir = "./feature.selection")

# Creating a vector of variable names
rfe.alpha <- paste0("rfe", c(1:10), ".alpha")
rfe.beta <- paste0("rfe", c(1:10), ".beta")
rfe.gamma <- paste0("rfe", c(1:10), ".gamma")

# Vector of locations of all the files
rfe.alpha.files <- paste0("./feature.selection/alpha/", list.files("./feature.selection/alpha"))
read.rfe(rfe.alpha, rfe.alpha.files)

rfe.beta.files <- paste0("./feature.selection/beta/", list.files("./feature.selection/beta"))
read.rfe(rfe.beta, rfe.beta.files)

rfe.gamma.files <- paste0("./feature.selection/gamma/", list.files("./feature.selection/gamma"))
read.rfe(rfe.gamma, rfe.gamma.files)

#     Analyzing patterns --------------------------------------------------

# Alpha-CD
rfe.alpha.vars <- list(
  rfe1.alpha,
  rfe2.alpha,
  rfe3.alpha,
  rfe4.alpha,
  rfe5.alpha, 
  rfe6.alpha, 
  rfe7.alpha, 
  rfe8.alpha, 
  rfe9.alpha, 
  rfe10.alpha
)

pred.alpha <- unlist(lapply(FUN = predictors, rfe.alpha.vars))
pred.alpha.uq <- unique(pred.alpha)
pred.alpha.pattern <- paste0("^", pred.alpha.uq) 
pred.alpha.pattern <- paste0(pred.alpha.pattern, "$")
count.alpha <- lapply(FUN = str_count, X = pred.alpha.pattern, string = pred.alpha) %>%
  lapply(FUN = sum, X = .) %>% unlist()
varimp.alpha <- data.frame(pred.alpha.uq, count.alpha) %>%
  rename(predictor = pred.alpha.uq, frequency = count.alpha) %>%
  mutate(predictor = as.character(predictor)) %>%
  .[order(.$frequency, decreasing = T), ]

# Limiting to the variables that appeared in 100% models
alpha.vars <- varimp.alpha %>% filter(frequency == 10 ) %>% .$predictor

saveRDS(varimp.alpha, "./feature.selection/varimp.alpha.RDS")
saveRDS(alpha.vars, "./feature.selection/alpha.vars.RDS")

# Beta-CD
rfe.beta.vars <- list(
  rfe1.beta,
  rfe2.beta,
  rfe3.beta,
  rfe4.beta,
  rfe5.beta, 
  rfe6.beta, 
  rfe7.beta, 
  rfe8.beta, 
  rfe9.beta, 
  rfe10.beta
)

pred.beta <- unlist(lapply(FUN = predictors, rfe.beta.vars))
pred.beta.uq <- unique(pred.beta)
pred.beta.pattern <- paste0("^", pred.beta.uq) 
pred.beta.pattern <- paste0(pred.beta.pattern, "$")
count.beta <- lapply(FUN = str_count, X = pred.beta.pattern, string = pred.beta) %>%
  lapply(FUN = sum, X = .) %>% unlist()
varimp.beta <- data.frame(pred.beta.uq, count.beta) %>%
  rename(predictor = pred.beta.uq, frequency = count.beta) %>%
  mutate(predictor = as.character(predictor)) %>%
  .[order(.$frequency, decreasing = T), ]

beta.vars <- varimp.beta %>% filter(frequency == 10) %>% .$predictor

saveRDS(varimp.beta, "./feature.selection/varimp.beta.RDS")
saveRDS(beta.vars, "./feature.selection/beta.vars.RDS")

# Gamma-CD
rfe.gamma.vars <- list(
  rfe1.gamma,
  rfe2.gamma,
  rfe3.gamma,
  rfe4.gamma,
  rfe5.gamma, 
  rfe6.gamma, 
  rfe7.gamma, 
  rfe8.gamma, 
  rfe9.gamma, 
  rfe10.gamma
)

pred.gamma <- unlist(lapply(FUN = predictors, rfe.gamma.vars))
pred.gamma.uq <- unique(pred.gamma)
pred.gamma.pattern <- paste0("^", pred.gamma.uq) 
pred.gamma.pattern <- paste0(pred.gamma.pattern, "$")
count.gamma <- lapply(FUN = str_count, X = pred.gamma.pattern, string = pred.gamma) %>%
  lapply(FUN = sum, X = .) %>% unlist()
varimp.gamma <- data.frame(pred.gamma.uq, count.gamma) %>%
  rename(predictor = pred.gamma.uq, frequency = count.gamma) %>%
  mutate(predictor = as.character(predictor)) %>%
  .[order(.$frequency, decreasing = T), ]

gamma.vars <- varimp.gamma %>% filter(frequency > 9) %>% .$predictor

# Large number of predictors indicates over-fittingm unfortunately
saveRDS(varimp.gamma, "./feature.selection/varimp.gamma.RDS")
saveRDS(gamma.vars, "./feature.selection/gamma.vars.RDS")
