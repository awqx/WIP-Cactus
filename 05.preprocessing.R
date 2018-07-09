dir.create("./pre-process")

# Libraries and Packages --------------------------------------------------

library(caret)
library(data.table)
library(Matrix)
library(stringr)
library(tidyverse)

# Functions ---------------------------------------------------------------

split.train.test <- function(times, data, info, path) {
  for(n in 1:times) {
    # Initialize a sample for reference for maxDissim
    data.init <- sample_n(data, 20)
    
    # Generate representative training set
    trn.ind <- maxDissim(a = data.init, b = data, 
                         n = round(nrow(data)*0.75), 
                         na.rm = T)
    tst.data <- data[-trn.ind, ]
    trn.data <- data[trn.ind, ]
    
    # Adding the guest molecule name and correct DelG
    tst <- cbind(info[-trn.ind, ], tst.data)
    trn <- cbind(info[trn.ind, ], trn.data)
    
    # Create logical file path name
    tst.path <- paste0(path, "tst", n, ".RDS")
    trn.path <- paste0(path, "trn", n, ".RDS")
    saveRDS(tst, tst.path)
    saveRDS(trn, trn.path)
    
    message("Split ", n, " completed")
  }
}

# Splitting data by folds instead of maxDissim

# data: data.frame with first column being guest, second column being DelG or
# the response 
# path: string of the target directory, ending with a backslash
split.train.test <- function(times, data, path) {
  folds <- createFolds(data[ , 2], k = times)
  for(i in 1:times) {
    tst.ind <- folds[[i]]
    tst <- data[tst.ind, ]
    trn <- data[-tst.ind, ]
    
    tst.path <- paste0(path, "tst", i, ".RDS")
    trn.path <- paste0(path, "trn", i, ".RDS")
    saveRDS(tst, tst.path)
    saveRDS(trn, trn.path)
  }
}

# Performing pre-processing on the files created by split.train.test
preprocess.splits <- function(filepath, writepath) {
  files <- list.files(path = filepath) %>% 
    str_extract(., "trn[[:print:]]+") %>% 
    .[!is.na(.)]
  for(n in 1:length(files)) {
    dirpath <- paste0(writepath, n)
    dir.create(dirpath)
    data <- readRDS(paste0(filepath, files[n]))
    info <- data %>% dplyr::select(guest, DelG)
    desc <- data %>% dplyr::select(-guest, -DelG)
    
    # Replacing inconvenient values
    desc <- do.call(data.frame, lapply(desc, 
                                       function(x)
                                         replace(x, is.infinite(x), NA)))
    
    # Initial pre-processing
    pp.settings <- preProcess(desc, na.remove = T, 
                              method = c("knnImpute", "center", "scale"), 
                              verbose = F)
    desc.pp <- predict(pp.settings, desc)
  
    # Removing predictors with near-zero variance
    zero.pred <- nearZeroVar(desc.pp)
    desc.pp <- desc.pp[ , -zero.pred]
    
    # Removing highly correlated predictors
    print(colnames(desc.pp))
    too.high <- findCorrelation(cor(desc.pp, use = "pairwise.complete.obs"), 0.95) # 0.95 mostly arbitrary
    desc.pp <- desc.pp[ , -too.high]
    
    desc.pp <- cbind(info, desc.pp)
    saveRDS(pp.settings, paste0(dirpath, "/pp.settings.RDS"))
    saveRDS(desc.pp, paste0(dirpath, "/pp.RDS"))
    saveRDS(zero.pred, paste0(dirpath, "/zero.pred.RDS"))
    saveRDS(too.high, paste0(dirpath, "/high.cor.RDS"))
    
    message("Pre-processing of split ", n, " completed.")
  }
}

# Cleaning for Outliers ---------------------------------------------------

dir.create("./pre-process/outliers")
alpha.info <- readRDS("./descriptors/alpha.padel.RDS")
beta.info <- readRDS("./descriptors/beta.padel.RDS")
gamma.info <- readRDS("./descriptors/gamma.padel.RDS")

# Using the method described by Roy 2015: Determining Applicability Domain of
# QSAR Models
source("./10.0.ad.functions.R")
alpha <- alpha.info %>% select(., -host, -DelG, -data.source)
# Technically not necessary for the statistical analysis, but it makes
# things easier
pp.settings <- preProcess(alpha, na.remove = T, 
                          method = c("center", "scale"), 
                          verbose = F)
alpha.scaled <- predict(pp.settings, alpha)
# Removing these to save time for the function, especially since they 
# definitely won't be used for the actual model building
zero.pred <- nearZeroVar(alpha.scaled)
alpha.scaled <- alpha.scaled[ , -zero.pred]
alpha.ad <- domain.num(alpha.scaled)
alpha.outliers <- alpha.ad %>% filter(domain == "outside") %>% .$guest

beta <- beta.info %>% select(., -host, -DelG, -data.source)
pp.settings <- preProcess(beta, na.remove = T, 
                          method = c("center", "scale"), 
                          verbose = F)
beta.scaled <- predict(pp.settings, beta)
zero.pred <- nearZeroVar(beta.scaled)
beta.scaled <- beta.scaled[ , -zero.pred]
beta.ad <- domain.num(beta.scaled)
beta.outliers <- beta.ad %>% filter(domain == "outside") %>% .$guest

gamma <- gamma.info %>% select(., -host, -DelG, -data.source)
pp.settings <- preProcess(gamma, na.remove = T, 
                          method = c("center", "scale"), 
                          verbose = F)
gamma.scaled <- predict(pp.settings, gamma)
zero.pred <- nearZeroVar(gamma.scaled)
gamma.scaled <- gamma.scaled[ , -zero.pred]
gamma.ad <- domain.num(gamma.scaled)
gamma.outliers <- gamma.ad %>% filter(domain == "outside") %>% .$guest

saveRDS(alpha.outliers, "./pre-process/outliers/alpha.RDS")
saveRDS(beta.outliers, "./pre-process/outliers/beta.RDS")
saveRDS(gamma.outliers, "./pre-process/outliers/gamma.RDS")

# Splitting Data ----------------------------------------------------------

#     External validation vs. modeling ------------------------------------

# Citing Tropsha, Best proactices for QSAR (DOI:10.1002/minf.201000061) an
# external validation (constituting 10-20% of the original data) should be set
# aside completely. The paper cites that this particular set should be randomly
# selected... probably indicating a simple random sample

# Of course, separate external validation sets should be created for each 
# cyclodextrin type (as they construct completely different models)
alpha <- readRDS("./descriptors/alpha.padel.RDS")
alpha.outliers <- readRDS('./pre-process/outliers/alpha.RDS')
alpha <- alpha %>% filter(!guest %in% alpha.outliers)
beta <- readRDS("./descriptors/beta.padel.RDS")
beta.outliers <- readRDS('./pre-process/outliers/beta.RDS')
beta <- beta %>% filter(!guest %in% beta.outliers)
gamma <- readRDS("./descriptors/gamma.padel.RDS")
gamma.outliers <- readRDS('./pre-process/outliers/gamma.RDS')
gamma <- gamma %>% filter(!guest %in% gamma.outliers)

dir.create("./ext.validation")

set.seed(101) # for reproducibility
alpha.ev <- sample_frac(alpha, size = 0.15)
beta.ev <- sample_frac(beta, size = 0.15)
gamma.ev <- sample_frac(gamma, size = 0.15)

saveRDS(alpha.ev, "./ext.validation/alpha.RDS")
saveRDS(beta.ev, "./ext.validation/beta.RDS")
saveRDS(gamma.ev, "./ext.validation/gamma.RDS")

# The remaining dataset (modeling data, according to Tropsha) should be kept
# separate. The relevant suffix is .md (modeling data)

dir.create("./model.data")
alpha.md <- alpha[!row.names(alpha) %in% row.names(alpha.ev), ]
beta.md <- beta[!row.names(beta) %in% row.names(beta.ev), ]
gamma.md <- gamma[!row.names(gamma) %in% row.names(gamma.ev), ]

saveRDS(alpha.md, "./model.data/alpha.md.RDS")
saveRDS(beta.md, "./model.data/beta.md.RDS")
saveRDS(gamma.md, "./model.data/gamma.md.RDS")

#     Test vs. train ------------------------------------------------------

# Multiple combinations of test and train sets should be created in order
# to fully validate the models. No specification was made in the paper
# as to how many different splits should be created, exactly, so I
# decided (arbitrarily) that 10 sets would be created.

# Though Sphere Exclusion modeling may be preferred, there is no R package
# that handles that algorith, so caret::createFolds was used here.

# Folds are created on Gibbs free energy, not the structures. This is 
# a weakness of the process that should be rectified later. It may be
# possible to create a Sphere Exclusion model using the instructions provided
# by Tropsha, but it would probably not be practical given the large
# number of descriptors. 

# Loading data
alpha <- readRDS("./model.data/alpha.md.RDS") %>% 
  select(-host, -data.source)
# %>% 
#   dplyr::select(., -guest:-data.source)
# alpha.info <- readRDS("./model.data/alpha.md.RDS") %>%
#   dplyr::select(guest, DelG)
beta <- readRDS("./model.data/beta.md.RDS") %>% 
  select(-host, -data.source)
gamma <- readRDS("./model.data/gamma.md.RDS") %>% 
  select(-host, -data.source, -guest.charge) %>%
  as.data.frame()

dir.create("./model.data/alpha")
dir.create("./model.data/beta")
dir.create("./model.data/gamma")
set.seed(101)
split.train.test(10, alpha, "./model.data/alpha/")
split.train.test(10, beta, "./model.data/beta/")
split.train.test(10, gamma, "./model.data/gamma/")

# Pre-processing and cleaning ---------------------------------------------

dir.create("./pre-process")
dir.create("./pre-process/alpha")
dir.create("./pre-process/beta")
dir.create("./pre-process/gamma")
preprocess.splits <- function(filepath, writepath) {
  files <- list.files(path = filepath) %>% 
    str_extract(., "trn[[:print:]]+") %>% 
    .[!is.na(.)]
  for(n in 1:length(files)) {
    dirpath <- paste0(writepath, n)
    dir.create(dirpath)
    data <- readRDS(paste0(filepath, files[n]))
    info <- data %>% dplyr::select(guest, DelG)
    desc <- data %>% dplyr::select(-guest, -DelG)
    
    # Replacing inconvenient values
    desc <- do.call(data.frame, lapply(desc, 
                                       function(x)
                                         replace(x, is.infinite(x), NA)))
    
    # Initial pre-processing
    pp.settings <- preProcess(desc, na.remove = T, 
                              method = c("knnImpute", "center", "scale"), 
                              verbose = F)
    desc.pp <- predict(pp.settings, desc)
    
    # Removing predictors with near-zero variance
    zero.pred <- nearZeroVar(desc.pp)
    desc.pp <- desc.pp[ , -zero.pred]
    
    # Removing highly correlated predictors
    desc.pp <- desc.pp[ , sapply(desc.pp, is.numeric)]
    too.high <- findCorrelation(cor(desc.pp, use = "pairwise.complete.obs"), 0.95) # 0.95 mostly arbitrary
    desc.pp <- desc.pp[ , -too.high]
    
    desc.pp <- cbind(info, desc.pp)
    saveRDS(pp.settings, paste0(dirpath, "/pp.settings.RDS"))
    saveRDS(desc.pp, paste0(dirpath, "/pp.RDS"))
    saveRDS(zero.pred, paste0(dirpath, "/zero.pred.RDS"))
    saveRDS(too.high, paste0(dirpath, "/high.cor.RDS"))
    
    message("Pre-processing of split ", n, " completed.")
  }
}
preprocess.splits(filepath = "./model.data/alpha/", 
                  writepath = "./pre-process/alpha/")
preprocess.splits(filepath = "./model.data/beta/", 
                  writepath = "./pre-process/beta/")
preprocess.splits(filepath = "./model.data/gamma/", 
                  writepath = "./pre-process/gamma/")

# # Suzuki Only -------------------------------------------------------------
# 
# suz <- readRDS("./descriptors/suz.all.padel.RDS")
# suz.nz <- suz[ , -zero.pred]
# suz.cd <- suz.nz  %>% 
#   mutate(alpha = ifelse(str_detect(host, "alpha"), 1, 0)) %>%
#   mutate(beta = ifelse(str_detect(host, "beta"), 1, 0)) %>%
#   mutate(gamma = ifelse(str_detect(host, "gamma"), 1, 0))
# 
# suz.split <- suz.cd %>% dplyr::select(., -guest:-data.source)
# suz.temp <- do.call(data.frame,
#                     lapply(suz.split, function(x) replace(x, is.infinite(x),NA)))
# suz.pp <- predict(pp.settings, suz.temp)
# suz.pp <- suz.pp[ , -zero.pred2]
# suz.pp <- suz.pp[ , -too.high]
# suz.pp <- cbind(suz.nz[ , 1:4], suz.pp)
# colnames(suz.pp) <- str_replace(colnames(suz.pp), "-", ".")
# 
# #     External Validation -------------------------------------------------
# 
# set.seed(4)
# ext.val.ind <- sample(x = 1:nrow(suz.pp), 
#                       size = round(0.15 * nrow(suz.pp)))
# ext.val <- suz.pp[ext.val.ind, ]
# suz.pp.all <- suz.pp
# suz.pp <- suz.pp[-ext.val.ind, ]
# 
# #     Saving --------------------------------------------------------------
# 
# saveRDS(suz.pp, "./data/suz.pp.RDS")
# saveRDS(ext.val, "./data/suz.extval.RDS")
