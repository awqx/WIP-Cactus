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
    data.init <- sample_n(data, 4)
    
    # Generate representative test set
    tst.ind <- maxDissim(a = data.init, b = data, 
                         n = round(nrow(data)*0.2), 
                         na.rm = T)
    tst.data <- data[tst.ind, ]
    trn.data <- data[-tst.ind, ]
    
    # Adding the guest molecule name and correct DelG
    tst <- cbind(info[tst.ind, ], tst.data)
    trn <- cbind(info[-tst.ind, ], trn.data)
    
    # Create logical file path name
    tst.path <- paste0(path, "tst", n, ".RDS")
    trn.path <- paste0(path, "trn", n, ".RDS")
    saveRDS(tst, tst.path)
    saveRDS(trn, trn.path)
    
    message("Split ", n, " completed")
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
    
    # Removing predictors with near-zero variance
    zero.pred <- nearZeroVar(desc)
    desc <- desc[ , -zero.pred]
    
    # Replacing inconvenient values
    desc <- do.call(data.frame, lapply(desc, 
                                       function(x)
                                         replace(x, is.infinite(x), NA)))
    
    # Initial pre-processing
    pp.settings <- preProcess(desc, na.remove = T, 
                              method = c("knnImpute", "center", "scale"))
    desc.pp <- predict(pp.settings, desc)
    
    # More nearZeroVar analysis
    zero.pred2 <- nearZeroVar(desc.pp)
    if (length(zero.pred2) >= 1) {
      desc.pp <- desc.pp[ , -zero.pred2]
      saveRDS(zero.pred2, paste0(dirpath, "/zero.pred2.RDS"))
    }
    
    # Removing highly correlated predictors
    too.high <- findCorrelation(cor(desc.pp), 0.95) # 0.95 mostly arbitrary
    desc.pp <- desc.pp[ , -too.high]
    
    desc.pp <- cbind(info, desc.pp)
    saveRDS(pp.settings, paste0(dirpath, "/pp.settings.RDS"))
    saveRDS(desc.pp, paste0(dirpath, "/pp.RDS"))
    saveRDS(zero.pred, paste0(dirpath, "/zero.pred.RDS"))
    saveRDS(too.high, paste0(dirpath, "/high.cor.RDS"))
    
    message("Pre-processing of split ", n, " completed.")
  }
}

# Splitting Data ----------------------------------------------------------

#     External validation vs. modeling ------------------------------------

# Citing Tropsha, Best proactices for QSAR (DOI:10.1002/minf.201000061) an
# external validation (constituting 10-20% of the original data) should be set
# aside completely. The paper cites that this particular set should be randomly
# selected... probably indicating a simple random sample

# Of course, separate external validation sets should be created for each 
# cyclodextrin type (as they construct completely different models)
 
dir.create("./ext.validation")

set.seed(101) # for reproducibility
alpha <- readRDS("./descriptors/alpha.padel.RDS")
alpha.ev <- sample_frac(alpha, size = 0.15)
beta <- readRDS("./descriptors/beta.padel.RDS")
beta.ev <- sample_frac(beta, size = 0.15)

saveRDS(alpha.ev, "./ext.validation/alpha.RDS")
saveRDS(beta.ev, "./ext.validation/beta.RDS")

# The remaining dataset (modeling data, according to Tropsha) should be kept
# separate. The relevant suffix is .md (modeling data)

dir.create("./model.data")
alpha.md <- alpha[!row.names(alpha) %in% row.names(alpha.ev), ]
beta.md <- beta[!row.names(beta) %in% row.names(beta.ev), ]

saveRDS(alpha.md, "./model.data/alpha.md.RDS")
saveRDS(beta.md, "./model.data/beta.md.RDS")

#     Test vs. train ------------------------------------------------------

# Multiple combinations of test and train sets should be created in order
# to fully validate the models. No specification was made in the paper
# as to how many different splits should be created, exactly, so I
# decided (arbitrarily) that 10 sets would be created

# Though Sphere Exclusion modeling may be preferred, there is no R package
# that handles that algorith, so caret::maxDissim was used here

# Loading data
alpha <- readRDS("./model.data/alpha.md.RDS") %>% 
  dplyr::select(., -guest:-data.source)
alpha.info <- readRDS("./model.data/alpha.md.RDS") %>%
  dplyr::select(guest, DelG)
beta <- readRDS("./model.data/beta.md.RDS") %>%
  dplyr::select(-guest:-data.source)
beta.info <- readRDS("./model.data/beta.md.RDS") %>%
  dplyr::select(guest, DelG)

# Splitting data
set.seed(101)
dir.create("./model.data/alpha")
dir.create("./model.data/beta")
split.train.test(10, alpha, alpha.info, "./model.data/alpha/")
split.train.test(10, beta, beta.info, "./model.data/beta/")

# Pre-processing and cleaning ---------------------------------------------

dir.create("./pre-process")
dir.create("./pre-process/alpha")
dir.create("./pre-process/beta")

preprocess.splits(filepath = "./model.data/alpha/", 
                  writepath = "./pre-process/alpha/")
preprocess.splits(filepath = "./model.data/beta/", 
                  writepath = "./pre-process/beta/")

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
