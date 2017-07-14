library(caret)
library(data.table)
library(Matrix)
library(RCurl)
library(stringr)
library(tidyverse)
library(XML)

# Pre-processing and cleaning ---------------------------------------------

ri.padel <- read_csv("~/SREP LAB/Rekharsky and Inoue/Cactus/03-PaDEL-Descriptor and Rekharsky and Inoue.csv")
# Removing predictors with near zero variance
zero.pred <- nearZeroVar(ri.padel)
zero.pred.names <- colnames(ri.padel)[zero.pred]
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
zero.pred2.names <- colnames(rpt)[zero.pred2]
rpt <- rpt[ , -zero.pred2]
rpt <- cbind(rp.split1, rpt)

# Extracting the predictors, removing highly correlated ones
rpt.pred <- rpt[ , -1:-17]
too.high <- findCorrelation(cor(rpt.pred), 0.95) # 0.95 mostly arbitrary
corr.pred <- names(rpt.pred)[too.high]
rpt.pred <- rpt.pred[ , -too.high]
rpt <- cbind(rpt[ , 1:17], rpt.pred)
saveRDS(rpt, "./rpt.RDS")

# Suzuki Dataset ----------------------------------------------------------
# Basically repeat same process as RI dataset
url   <- "http://pubs.acs.org/doi/full/10.1021/ci010295f"
file  <-
  getURL(url, .opts = curlOptions(followlocation = TRUE, cookiefile = "nosuchfile"))
suzuki.list <-
  readHTMLTable(
    file,
    header = T,
    as.data.frame = T,
    stringAsFactors = F
  )
index11 <- suzuki.list %>% 
  lapply(names) %>% 
  lapply(length) == 11 %>% as.vector()
suzuki.raw <- rbindlist(suzuki.list[index11], fill = T) %>% 
  lapply(as.character) %>% 
  as.data.frame(stringsAsFactors = F)

# Column names got shifted to row one, select relevant columns 
suzuki <- suzuki.raw[-1, c(2, 4, 8)]
# only relevant columns are guest, and DelG obsd
colnames(suzuki) <- c("guest", "DelG.a", "DelG.b")
suzuki <- suzuki %>% # Replacing weird characters
  mutate(DelG.a = str_replace(DelG.a, "Ã¢\u0088\u0092", "-")) %>%
  mutate(DelG.a = str_replace(DelG.a, "(Ã¢|Ã)", "")) %>%
  mutate(DelG.b = str_replace(DelG.b, "Ã¢\u0088\u0092", "-")) %>%
  mutate(DelG.b = str_replace(DelG.b, "(Ã¢|Ã)", "")) %>%
  mutate(guest = str_replace(guest, "Ã¢\u0080\u0089", "-")) 
# Cleaning out rows that are just header names
h.ind <- suzuki$guest %>% str_detect("guest")
suzuki <- suzuki[!h.ind, ] 
suzuki[ , 2:3] <- sapply(suzuki[, 2:3], as.numeric) 
# Separating the two CD types 
suzuki.dg.alpha <- suzuki %>% 
  filter(!is.na(DelG.a)) %>% 
  select(guest, DelG.a) %>% 
  rename(DelG = DelG.a)
suzuki.dg.beta <- suzuki %>%  
  select(guest, DelG.b) %>%
  rename(DelG = DelG.b)
# Binding together the results from PaDEL
# I lost the code that I used to download this; will find later
suzuki.beta <- read_csv("./2017 06 16 Suzuki beta.csv") %>% 
  mutate(alpha = c(rep(0, 185))) %>%
  mutate(beta = c(rep(1, 185))) %>%
  mutate(gamma = c(rep(0, 185))) %>%
  rename(guest = Name) %>%
  inner_join(suzuki.dg.beta, by = "guest", .)
suzuki.alpha <- read_csv("./2017 06 16 Suzuki alpha comp.csv") %>%
  mutate(alpha = c(rep(1, 55))) %>%
  mutate(beta = c(rep(0, 55))) %>%
  mutate(gamma = c(rep(0, 55))) %>%
  rename(guest = Name) %>%
  inner_join(suzuki.dg.alpha, by = "guest", .)
suzuki.padel <- rbind(suzuki.alpha, suzuki.beta)
colnames(suzuki.padel) <- colnames(suzuki.padel) %>% str_replace("-", "\\.")
saveRDS(suzuki.padel, "./suzuki.padel.RDS")

sp.no.zero <- suzuki.padel[ , !colnames(suzuki.padel) %in% zero.pred.names]
st <- preProcess(sp.no.zero[ , -1:-2], na.remove = T, 
                 method = c("knnImpute", "center", "scale")) %>%
  predict(., sp.no.zero[ , -1:-2])

st <- st[ , !colnames(st) %in% zero.pred2.names]
st.nocor <- st[ , !colnames(st) %in% corr.pred]
suzuki.transformed <- cbind(sp.no.zero[ , 1:2], st.nocor)
saveRDS(suzuki.transformed, "./suzuki.clean.RDS")

# Attaching Data ----------------------------------------------------------

# Sorting the columns that are in both data sets
rpt.shared <- rpt[ , colnames(rpt) %in% colnames(suzuki.transformed)]
df <- rbind(rpt.shared, suzuki.transformed)

# Creating External Validation Set ----------------------------------------
set.seed(4)
ext.val.ind <- sample(x = 1:nrow(df), 
                      size = round(0.1 * nrow(df)))
ext.val <- df[ext.val.ind, ]
saveRDS(ext.val, "./external validation set.RDS")

# Data Organization -------------------------------------------------------
# Delta G
# df.dg <- dplyr::select(rpt, -X1:-log.K.Uncertainty, -DelG.Uncertainty:-`bind.aff, kcal/mol`)
# sparse.dg <- sparse.model.matrix(~., df.dg)
# mat.dg <- data.matrix(df.dg)
# saveRDS(df.dg, "./DelG.df.RDS")
# saveRDS(sparse.dg, "./DelG.sparse.RDS")
# saveRDS(mat.dg, "./DelG.matrix.RDS")
sparse.df <- sparse.model.matrix(~., df)
mat.df <- as.matrix(df)
saveRDS(df, "./DelG.df.RDS")
saveRDS(sparse.df, "./DelG.sprse.RDS")
saveRDS(mat.df, "./DelG.mat.RDS")


# Log K
# df.logk <- dplyr::select(rpt, -X1:-Solvent.composition, -log.K.Uncertainty:-`bind.aff, kcal/mol`)
# sparse.logk <- sparse.model.matrix(~., df.logk)
# mat.logk <- data.matrix(df.logk)
# saveRDS(df.logk, "./LogK.df.RDS")
# saveRDS(sparse.logk, "./LogK.sparse.RDS")
# saveRDS(mat.logk, "./LogK.matrix.RDS")

# set.seed(512)
# trn.ind <- sample(x = 1:nrow(sparse.dg), size = round(0.8 * nrow(sparse.dg)))
# sparse.dg.trn <- sparse.dg[trn.ind, ]
# sparse.dg.tst <- sparse.dg[-trn.ind, ]