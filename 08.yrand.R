# Models evaluated by standards from Tropsha 2010, DOI: 10.1002/minf201000061
# "Best Practices for QSAR Model Development, Validation, and Exploitation"
# Standards include y-randomization. 

# According Rucker, Rucker, and Meringer in 'y-Randomization and its Variants in
# QSAR/QSOR' (DOI: 10.1021/ci700157b), there are 5 modes of exploring the
# viability of the model using different forms of permutation and randomization.
# For the sake of time, I'll only be using Mode 1: scrambling (permuting) the
# observed DelG (y) values while keeping descriptors the same

dir.create('yrand')
source('data.handling.R')
source('modeling.R')

# Functions ---------------------------------------------------------------

# df: data.frame  contains the observed DelG (in col1) and descriptors (num)
# n : int         number of repetitions

# Returns a data.frame with permuted y (first column) values
permy <- function(df) {
  new.rows <- sample(nrow(df), replace = F)
  y.orig <- df[ , 1]
  y.new <- c()
  for(i in 1:nrow(df))
    y.new[i] <- y.orig[new.rows[i]]
  df.new <- data.frame(y.new, df[ , -1])
  colnames(df.new) <- colnames(df)
  return(df.new)
}

# Returns a list of permy data.frames
permy.ntimes <- function(df, n) {
  results <- list()
  for(i in 1:n) {
    results[[i]] <- permy(df)
  }
  return(results)
}

# name: char    the name of the dir to save the files
save.perms <- function(permlist, dir) {
  # Fixing the directory names
  if(!str_detect(dir, '/S'))
    dir <- paste0(dir, '/')
  for (i in 1:length(permlist)) {
    dir.create(paste0(dir, i)) # Create folders
    saveRDS(permlist[[i]], paste0(dir, i, '/data.RDS'))
  }
}

tempdf <- data.frame(y = 1:20, x = 1:20)

# Alpha-CD ----------------------------------------------------------------

dir.create('yrand/alpha')

#     Organizing data -----------------------------------------------------

# 1. Initial data ====
  # Merging the external validation data and the modeling data because this
  # is the easiest way to get the data without outliers
alpha1 <- readRDS('ext.validation/alpha.RDS') %>%
  select(., -guest, -host, -data.source)
alpha2 <- readRDS('model.data/alpha.md.RDS') %>%
  select(., -guest, -host, -data.source)
alpha.info <- rbind(alpha1, alpha2)

  # 25 is the number of iterations used in the Rucker et al. paper
set.seed(101)
alpha.perm <- permy.ntimes(alpha.info, 25)
save.perms(alpha.perm, 'yrand/alpha')

# 2. Creating external validation sets ====
alpha.path <- paste0('yrand/alpha/', 1:25)
  # Sorting out the modeling data
alpha.md <- mapply(ev.split,
                   alpha.perm, alpha.path, 
                   SIMPLIFY = F)

# 3. Creating train-test splits  ====
    # Only 5 each for time's sake
alpha.tt.path <- paste0(alpha.path, '/tt/')
lapply(alpha.tt.path, dir.create)
mapply(tt.split, 
       alpha.md, alpha.tt.path, 
       MoreArgs = list(times = 5))

# 4. Pre-processing ====
set.seed(101)
alpha.pp.path <- paste0(alpha.path, '/pp/')
lapply(alpha.pp.path, dir.create)
mapply(pp.split, 
       alpha.tt.path, alpha.pp.path)

# 5. Feature selection ====
    # Using recursive feature elimination
set.seed(101)
rfe.yrand('yrand/alpha/', 25, 5)
    # Though I've since removed the function's ability to print the list of
    # variables to the console, I can attest that the number of variables
    # selected as important varies greatly among the different trials, 
    # ranging from 2 to 200+. 

# 6. Model building ====
# Create directories to store the results and models
lapply(c(paste0(alpha.path, '/models'), 
         paste0(alpha.path, '/results')),
       dir.create)

# Using lapply for the build.[model] functions to avoid an addition for loop
# X in the function indicates the ntrial
  # Cubist
lapply(X = c(1:25), 
       FUN = build.cube, 
       host = 'alpha', 
       nsplit = 5, 
       seed = 101)
  # GBM
lapply(X = c(1:25), 
       FUN = build.gbm, 
       host = 'alpha', 
       nsplit = 5, 
       seed = 101)

  # GLMNet
lapply(X = c(1:25), 
       FUN = build.glm, 
       host = 'alpha', 
       nsplit = 5, 
       seed = 101)

  # MARS
lapply(X = c(1:25), 
       FUN = build.mars, 
       host = 'alpha', 
       nsplit = 5, 
       seed = 101)

  # Random forest
lapply(X = c(1:25), 
       FUN = build.rf, 
       host = 'alpha', 
       nsplit = 5, 
       seed = 101)
  # SVM
    # Polynomial
lapply(X = c(1:25), 
       FUN = build.polysvm, 
       host = 'alpha', 
       nsplit = 5, 
       seed = 101)
    # Radial 
lapply(X = c(1:25), 
       FUN = build.rbfsvm, 
       host = 'alpha', 
       nsplit = 5, 
       seed = 101)
    # Sigmoid
lapply(X = c(1:25), 
       FUN = build.sigsvm, 
       host = 'alpha', 
       nsplit = 5, 
       seed = 101)

  # PLS 
lapply(X = c(1:25), 
       FUN = build.pls, 
       host = 'alpha', 
       nsplit = 5, 
       seed = 101)



# 7. Statistics -----------------------------------------------------------

# Get q2 mean and standard deviation
models <- c('cubist', 'gbm', 'glm', 'mars', 'pls', 'rf', 
            'svm-poly', 'svm-rbf', 'svm-sig')
q2.alpha <- do.call(
  rbind, 
  lapply(
    models, 
    get.q2.sd, 
    host = 'alpha', 
    ntrial = 25
  )
)
saveRDS(q2.alpha, 'yrand/alpha/q2.results.RDS')

# Beta-CD -----------------------------------------------------------------

dir.create('yrand/beta')

#     Organizing data -----------------------------------------------------

# 1. Initial data ====
# Merging the external validation data and the modeling data because this
# is the easiest way to get the data without outliers
beta1 <- readRDS('ext.validation/beta.RDS') %>%
  select(., -guest, -host, -data.source)
beta2 <- readRDS('model.data/beta.md.RDS') %>%
  select(., -guest, -host, -data.source)
beta.info <- rbind(beta1, beta2)

# 25 is the number of iterations used in the Rucker et al. paper
set.seed(101)
beta.perm <- permy.ntimes(beta.info, 25)
save.perms(beta.perm, 'yrand/beta')

# 2. Creating external validation sets ====
beta.path <- paste0('yrand/beta/', 1:25)
# Sorting out the modeling data
beta.md <- mapply(ev.split,
                   beta.perm, beta.path, 
                   SIMPLIFY = F)

# 3. Creating train-test splits  ====
# Only 5 each for time's sake
beta.tt.path <- paste0(beta.path, '/tt/')
lapply(beta.tt.path, dir.create)
mapply(tt.split, 
       beta.md, beta.tt.path, 
       MoreArgs = list(times = 5))

# 4. Pre-processing ====
set.seed(101)
beta.pp.path <- paste0(beta.path, '/pp/')
lapply(beta.pp.path, dir.create)
mapply(pp.split, 
       beta.tt.path, beta.pp.path)

# 5. Feature selection ====
# Using recursive feature elimination
set.seed(101)
rfe.yrand('yrand/beta/', 25, 5)
# Though I've since removed the function's ability to print the list of
# variables to the console, I can attest that the number of variables
# selected as important varies greatly among the different trials, 
# ranging from 2 to 200+. 

# 6. Model building ====
# Create directories to store the results and models
lapply(c(paste0(beta.path, '/models'), 
         paste0(beta.path, '/results')),
       dir.create)

# Using lapply for the build.[model] functions to avoid an addition for loop
# X in the function indicates the ntrial
# Cubist
lapply(X = c(1:12), 
       FUN = build.cube, 
       host = 'beta', 
       nsplit = 5, 
       seed = 101)
# Had to re-do 13 b/c code needed to be fixed
# Would re-run everything, but unfortunately cubist takes very long
build.cube('beta', 13, 5, 101)
lapply(X = c(14:25), 
       FUN = build.cube, 
       host = 'beta', 
       nsplit = 5, 
       seed = 101)

  # GBM
lapply(X = c(1:12), 
       FUN = build.gbm, 
       host = 'beta', 
       nsplit = 5, 
       seed = 101)
# Doesn't appear GBM is compatible w/ single variable
lapply(X = c(14:25), 
       FUN = build.gbm, 
       host = 'beta', 
       nsplit = 5, 
       seed = 101)

  # GLMNet
lapply(X = c(1:12), # not trial 13 compatible
       FUN = build.glm, 
       host = 'beta', 
       nsplit = 5, 
       seed = 101)
lapply(X = c(14:25), 
       FUN = build.glm, 
       host = 'beta', 
       nsplit = 5, 
       seed = 101)

  # MARS
lapply(X = c(1:25), 
       FUN = build.mars, 
       host = 'beta', 
       nsplit = 5, 
       seed = 101)

  # Random forest
lapply(X = c(1:25), 
       FUN = build.rf, 
       host = 'beta', 
       nsplit = 5, 
       seed = 101)
# SVM
# Polynomial
lapply(X = c(1:25), 
       FUN = build.polysvm, 
       host = 'beta', 
       nsplit = 5, 
       seed = 101)
# Radial 
lapply(X = c(1:25), 
       FUN = build.rbfsvm, 
       host = 'beta', 
       nsplit = 5, 
       seed = 101)
# Sigmoid
lapply(X = c(1:25), 
       FUN = build.sigsvm, 
       host = 'beta', 
       nsplit = 5, 
       seed = 101)

# PLS 
lapply(X = c(1:25), 
       FUN = build.pls, 
       host = 'beta', 
       nsplit = 5, 
       seed = 101)

# 7. Statistics -----------------------------------------------------------

skip13.models <- c('gbm', 'glm')
models <- models[!models %in% skip13.models]
# Q2
q2.beta1 <- do.call(
  rbind, 
  lapply(
    skip13.models, 
    get.q2.sd, 
    host = 'beta', 
    ntrial = 25, 
    skip13 = T
  )
)
q2.beta2 <- do.call(
  rbind, 
  lapply(
    models, 
    get.q2.sd, 
    host = 'beta', 
    ntrial = 25
  )
)

q2.beta <- rbind(q2.beta1, q2.beta2)
saveRDS(q2.beta, 'yrand/beta/q2.results.RDS')
