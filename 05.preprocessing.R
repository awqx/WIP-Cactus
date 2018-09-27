dir.create("./pre-process")
source('data.handling.R')

# Cleaning for Outliers ---------------------------------------------------

dir.create("./pre-process/outliers")
alpha.info <- readRDS("./descriptors/alpha.padel.RDS")
beta.info <- readRDS("./descriptors/beta.padel.RDS")
gamma.info <- readRDS("./descriptors/gamma.padel.RDS")

# Structural outliers
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

# Activity outliers
alpha <- readRDS("./descriptors/alpha.padel.RDS") %>%
  select(-host)
alpha.act.out <- find.activity.outlier(alpha, 2.5) 
alpha.act.out <- alpha[alpha.act.out, "guest"]
alpha.outliers <- c(alpha.outliers, alpha.act.out)

beta <- readRDS("./descriptors/beta.padel.RDS") %>%
  select(-host)
beta.act.out <- find.activity.outlier(beta, 2.5) 
beta.act.out <- beta[beta.act.out, "guest"]
beta.outliers <- c(beta.outliers, beta.act.out)

gamma <- readRDS("./descriptors/gamma.padel.RDS") %>%
  select(-host) %>% data.frame()
gamma.act.out <- find.activity.outlier(gamma, 2.5) 
gamma.act.out <- gamma[gamma.act.out, "guest"]
gamma.outliers <- c(gamma.outliers, gamma.act.out)

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
readRDS('results/alpha/cube.RDS')[-6, ] %>% defaultSummary()
readRDS('results/alpha/glmnet.RDS') %>% defaultSummary()
readRDS('results/alpha/pls.RDS') %>% defaultSummary()
readRDS('results/alpha/polysvm.RDS') %>% defaultSummary()
readRDS('results/alpha/rbfsvm.RDS') %>% defaultSummary()
readRDS('results/alpha/rf.RDS') %>% defaultSummary()

readRDS('results/beta/cube.RDS')[-6, ] %>% defaultSummary()
readRDS('results/beta/glmnet.RDS') %>% defaultSummary()
readRDS('results/beta/pls.RDS') %>% defaultSummary()
readRDS('results/beta/polysvm.RDS') %>% defaultSummary()
readRDS('results/beta/rbfsvm.RDS') %>% defaultSummary()
readRDS('results/beta/rf.RDS') %>% defaultSummary()