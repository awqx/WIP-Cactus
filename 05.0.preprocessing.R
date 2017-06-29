library(caret)
library(Matrix)
library(stringr)
library(tidyverse)

# Pre-processing and cleaning ---------------------------------------------

ri.padel <- read_csv("~/SREP LAB/Rekharsky and Inoue/Cactus/03-PaDEL-Descriptor and Rekharsky and Inoue.csv")
# Removing predictors with near zero variance
zero.pred <- nearZeroVar(ri.padel)
rp.no.zero <- ri.padel[ , -zero.pred]

# Binning alpha, beta, and gamma
rp.cd <- rp.no.zero %>% 
  mutate(alpha = ifelse(str_detect(host, "Alpha"), 1, 0)) %>%
  mutate(beta = ifelse(str_detect(host, "Beta"), 1, 0)) %>%
  mutate(gamma = ifelse(str_detect(host, "Gamma"), 1, 0))

# Separating responses and characters from the predictors
rp.split1 <- rp.cd %>% dplyr::select(., X1:`bind.aff, kcal/mol`)
rp.split2 <- rp.cd %>% dplyr::select(., -X1:-`bind.aff, kcal/mol`)

rpt <- preProcess(rp.split2, na.remove = T, 
                  method = c("knnImpute", "center", "scale")) %>%
  predict(., rp.split2)
rpt <- preProcess(rpt, method = "BoxCox") %>%
  predict(., rpt)
zero.pred2 <- nearZeroVar(rpt)
rpt <- rpt[ , -zero.pred2]
rpt <- cbind(rp.split1, rpt)


rpt.pred <- rpt[ , -1:-17]
too.high <- findCorrelation(cor(rpt.pred), 0.95)
corr.pred <- names(rpt.pred)[too.high]
rpt.pred <- rpt.pred[ , -too.high]
rpt <- cbind(rpt[ , 1:17], rpt.pred)
saveRDS(rpt, "./rpt.RDS")

# Data Organization -------------------------------------------------------
# Delta G
df.dg <- dplyr::select(rpt, -X1:-log.K.Uncertainty, -DelG.Uncertainty:-`bind.aff, kcal/mol`)
sparse.dg <- sparse.model.matrix(~., df.dg)
mat.dg <- data.matrix(df.dg)
saveRDS(df.dg, "./DelG.df.RDS")
saveRDS(sparse.dg, "./DelG.sparse.RDS")
saveRDS(mat.dg, "./DelG.matrix.RDS")

# Log K
df.logk <- dplyr::select(rpt, -X1:-Solvent.composition, -log.K.Uncertainty:-`bind.aff, kcal/mol`)
sparse.logk <- sparse.model.matrix(~., df.logk)
mat.logk <- data.matrix(df.logk)
saveRDS(df.logk, "./LogK.df.RDS")
saveRDS(sparse.logk, "./LogK.sparse.RDS")
saveRDS(mat.logk, "./LogK.matrix.RDS")

# set.seed(512)
# trn.ind <- sample(x = 1:nrow(sparse.dg), size = round(0.8 * nrow(sparse.dg)))
# sparse.dg.trn <- sparse.dg[trn.ind, ]
# sparse.dg.tst <- sparse.dg[-trn.ind, ]