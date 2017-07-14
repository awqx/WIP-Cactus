library(caret)
library(devtools)
library(doParallel)
library(pls)
library(pso)
install_github("jakobbossek/acotsp")
library(acotsp)


set.seed(48)
df.dg <- readRDS("./DelG.df.RDS")
mat.dg <- readRDS("./DelG.matrix.RDS")
mat.trn.ind <- sample(x = 1:nrow(mat.dg), size = round(0.8 * nrow(mat.dg)))
mat.trn <- mat.dg[mat.trn.ind, ]
mat.tst <- mat.dg[-mat.trn.ind, ]
mat.x <- mat.trn[ , -1]
mat.y <- mat.trn[ , 1]

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


