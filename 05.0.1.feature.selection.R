library(caret)
library(devtools)
library(doParallel)
library(pls)
library(pso)
install_github("jakobbossek/acotsp")
library(acotsp)

setwd("~/SREP LAB/qsar")
dir.create("./feature.select")

# Data Organization -------------------------------------------------------

set.seed(48)
df <- readRDS("./padel.pp.RDS")
mat <- readRDS("./mat.padel.RDS")
trn.ind <- sample(x = 1:nrow(mat), size = round(0.8 * nrow(mat)))
mat.trn <- mat[trn.ind, ]
mat.tst <- mat[-trn.ind, ]
mat.x <- mat.trn[ , -1]
mat.y <- mat.trn[ , 1]
df.trn <- df[trn.ind, ]
df.tst <- df[-trn.ind, ]

# Genetic Algorithm -------------------------------------------------------
# rf.ga3 <- readRDS("./genetic alg.RDS")
set.seed(10)
registerDoParallel(4)
# getDoParWorkers()
ga.ctrl <- gafsControl(functions = rfGA, # Assess fitness with RF
                       method = "cv",    # 10 fold cross validation
                       genParallel=TRUE, # Use parallel programming
                       allowParallel = TRUE)
set.seed(10)
lev <- c("PS", "WS")
rf.ga3 <- gafs(x = mat.x, y = mat.y, 
               iters = 100, 
               popSize = 20, 
               levels = lev, 
               gafsControl = ga.ctrl) # Started at 9:47:40 AM. Ended at 23:00:00 PM
summary(rf.ga3)
ga.final <- rf.ga3$ga$final

df.ga <- df.dg[colnames(df.dg) %in% ga.final]
df.ga <- cbind(df.dg$DelG, df.ga)
colnames(df.ga)[1] <- "DelG"
saveRDS(rf.ga3, "./genetic alg.RDS")
saveRDS(df.ga, "./genetic alg results.RDS") 
saveRDS(df.ga, "./GAFS.RDS")
# mat.ga <- data.matrix(df.ga)
# sparse.ga <- sparse.model.matrix(~., df.ga)
# sparse.data <- sparse.ga
# 
# set.seed(256)
# trn.ind <- sample(x = 1:nrow(sparse.data), size = round(0.8 * nrow(sparse.data))) 
# sparse.trn <- sparse.data[trn.ind, ]
# sparse.tst <- sparse.data[-trn.ind, ]
# trn.ind <- sample(x = 1:nrow(df.ga), size = round(0.8 * nrow(df.ga))) 
# df.trn <- df.ga[trn.ind, ]
# df.tst <- df.ga[-trn.ind, ]
# df.x <- df.ga[ , -1]
# df.y <- df.ga[ , 1]
# df.trn.x <- df.trn[ , -1]
# df.trn.y <- df.trn[ , 1]
# df.tst.x <- df.tst[ , -1]
# df.tst.y <- df.tst[ , 1]
# sprse.x <- sparse.data[ , -1:-2]
# sprse.y <- sparse.data[ , 2]
# sprse.trn.x <- sparse.trn[ , -1:-2]
# sprse.trn.y <- sparse.trn[ , 2]
# sprse.tst.x <- sparse.tst[ , -1:-2]
# sprse.tst.y <- sparse.tst[ , 2]

# Ant Colony --------------------------------------------------------------
# net.dg <- makeNetwork(mat.dg)
# ant.cont <- makeACOTSPControl(n.ants = 100L, n.elite = 50L, prp.prob = 0.5)
# runACOTSP(net.dg, ant.cont)


# Particle Swarm Optimization ---------------------------------------------


# PLS ---------------------------------------------------------------------

pls.mod <- plsr(DelG ~ ., 
                ncomp = 10, 
                data = df.trn, 
                validation = "LOO", 
                method = "oscorespls")
summary(pls.mod)
plot(RMSEP(pls.mod), legendpos = "topright") 

#     VIP Analysis --------------------------------------------------------
# Function sourced from http://mevik.net/work/software/pls.html
VIP <- function(object) {
  if (object$method != "oscorespls")
    stop("Only implemented for orthogonal scores algorithm.  Refit with 'method = \"oscorespls\"'")
  if (nrow(object$Yloadings) > 1)
    stop("Only implemented for single-response models")
  
  SS <- c(object$Yloadings)^2 * colSums(object$scores^2)
  Wnorm2 <- colSums(object$loading.weights^2)
  SSW <- sweep(object$loading.weights^2, 2, SS / Wnorm2, "*")
  sqrt(nrow(SSW) * apply(SSW, 1, cumsum) / cumsum(SS))
}

sort.vip <- function(vip.df, row.num) {
  long <- vip.df[row.num, ] %>% gather()
  return(long[long[ , 2] > 1.2, ])
}
pls.vip <- VIP(pls.mod) %>% data.frame()
pls.vip1 <- sort.vip(pls.vip, 1)
pls.vip2 <- sort.vip(pls.vip, 2)
pls.vip3 <- sort.vip(pls.vip, 3)
pls.vip4 <- sort.vip(pls.vip, 4)
pls.vip5 <- sort.vip(pls.vip, 5)
pls.vip6 <- sort.vip(pls.vip, 6)
pls.vip7 <- sort.vip(pls.vip, 7)
pls.vip8 <- sort.vip(pls.vip, 8)
pls.vip9 <- sort.vip(pls.vip, 9)
pls.vip10 <- sort.vip(pls.vip, 10)

highest.vip <- inner_join(pls.vip1, pls.vip2, by = "key") %>%
  inner_join(., pls.vip3, by = "key") %>%
  inner_join(., pls.vip4, by = "key") %>%
  inner_join(., pls.vip5, by = "key") %>%
  inner_join(., pls.vip6, by = "key") %>%
  inner_join(., pls.vip7, by = "key") %>%
  inner_join(., pls.vip8, by = "key") %>%
  inner_join(., pls.vip9, by = "key") %>%
  inner_join(., pls.vip10, by = "key") 
highest.vip.desc <- highest.vip$key