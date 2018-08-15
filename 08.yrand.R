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

# 6. Tuning ====
# Create directories to store the results and models
lapply(c(paste0(alpha.path, '/models'), 
         paste0(alpha.path, '/results')),
       dir.create)

# Using lapply for the build.[model] functions to avoid an addition for loop
# X in the function indicates the ntrial
  # GLMNet

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


