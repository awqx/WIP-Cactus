setwd("~/SREP LAB/qsar")
source("./model.code/tuning.functions.R")

# Data Organization -------------------------------------------------------

rpt <- readRDS("./rpt.RDS")
mat.dg <- rpt %>% 
  dplyr::select(., -X1:-log.K.Uncertainty,
                -DelG.Uncertainty:-`bind.aff, kcal/mol`)

set.seed(25)
trn.ind <- sample(x = 1:nrow(mat.dg), size = round(0.7 * nrow(mat.dg)))
trn <- mat.dg[trn.ind, ]
tst <- mat.dg[-trn.ind, ]

dir.create("./tuning/cubist")

cube.tst <- cubist(x = trn[ , -1], y = trn[ , 1])
predict(cube.tst, tst[ , -1]) %>%
  cbind(tst[ , 1])

# Tuning ------------------------------------------------------------------


#     Committees ----------------------------------------------------------

cube.cmte <- do.call(
  rbind,
  lapply(
    2 ^ (2:9),
    FUN = tune.cubist.cmte,
    data = trn,
    nfolds = 10, 
    seed = 1
  )
)
