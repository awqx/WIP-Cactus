install.packages("caret")
install.packages("AppliedPredictiveModeling")
install.packages("lars")
install.packages("pls")
install.packages("elasticnet")
install.packages("broom")
library(MASS)
library(caret)
library(AppliedPredictiveModeling)
library(lars)
library(pls)
library(elasticnet)
library(dplyr)
library(stringr)
library(Amelia)
library(broom)

library(irlba)

data(solubility)
ls(pattern = "^solT")
training <- solTrainXtrans # box-cox transformed
training$solubility <- solTrainY
sample(names(training), 10)
lm.allpred <- lm(solubility ~ ., data = training)
summary(lm.allpred)
# redict test set
lm.testpred <- predict(lm.allpred, solTestXtrans)
# evaluation of model
lm.eval <- data.frame(obs = solTestY, pred = lm.testpred)
defaultSummary(lm.eval)
# plot predicted and observed values
plot(lm.eval, main = "Solubility Values", ylab = "Predicted Values", 
     xlab = "True Observed Valuese")
lines(c(-10, 1), c(-10, 1), col = "red")

# ===========================================================

gamma.chem.mm2 <- read.csv(paste0(chem.path, "/GammaID.mm2.csv"))
names(gamma.ds)[names(gamma.ds) == "formula"] <- "Name"
gamma.chem.mm2 <- full_join(gamma.chem.mm2, gamma.ds)
gamma.chem.mm2 <- gamma.chem.mm2[-c(59,60), ]
lm.gammapred <- lm(binding.affinity ~ ., data = gamma.chem.mm2)
summary(lm.gammapred)


#gamma.mmff94 <- read.csv(paste0(chem.path, "/GammaID.mmff942.csv"), header = T, stringsAsFactors = F)
#names(gamma.ds)[names(gamma.ds) == "formula"] <- "Name"
gamma.mmff94 <- full_join(gamma.mmff94, gamma.ds)
gamma.mmff94 <- gamma.mmff94[-c(59,60), ]
#lm.gammapred1 <- lm(binding.affinity ~ ., data = gamma.mmff94)
#summary(lm.gammapred)

gamma.mmff94 <- read.csv(paste0(chem.path, "/GammaID.mmff942.csv"), header = T, stringsAsFactors = F)
names(gamma.ds)[names(gamma.ds) == "formula"] <- "Name"
gamma.mmff94 <- full_join(gamma.mmff94, gamma.ds)

benzene.row <- gamma.mmff94[22, ]
bad.col <- which(is.na(benzene.row))
bad.col <- colnames(benzene.row[ , bad.col])
bad.col.v <- as.vector(bad.col)
gamma.m94.clean <- gamma.mmff94[ , -bad.col]
missmap(gamma.m94.clean)
not.missing <- gamma.mmff94[!is.na(gamma.mmff94$TDB1i), ]
not.missing.chem <- not.missing$Name
not.missing.values <- gamma.mmff94[which(!is.na(gamma.mmff94$TDB1i)), bad.col]
gamma94 <- gamma.m94.clean
gamma94$Name <- NULL
gamma94Train <- gamma94[1:40, ]
gamma94Test <- gamma94[41:58, ]
lm.gammapred <- lm(binding.affinity ~ ., data = gamma94Train)
summary(lm.gammapred)
tidy.gamma <- tidy(lm.gammapred)
gamma.pred <- predict(lm.gammapred, gamma94Test, level = 0.95)
#gamma.mmff94 <- gamma.mmff94[-c(59,60), ]
#gamma.mmff94$Name <- NULL
#lm.gammapred1 <- lm(binding.affinity ~ ., data = gamma.mmff94)
#summary(lm.gammapred)
#amelia.mmff94 <- gamma.mmff94
#missmap(amelia.mmff94)

#gamma.clean <- min.aff.gamma$formula %>%
#str_replace(pattern = "^[[:alnum:]]+_", "") %>%
#  str_replace(pattern = "_ghemical_[[:graph:]]+", "")

# prcomp_irlba(alpha.chem.noname.matrix.small, 5)

# rstep.alpha <- stepwise(alpha.y~., data = alpha.chem.noname) 


# ===========
# IRLBA Attempt
library(irlba)
alpha.chem.noname.matrix <- matrix(as.numeric(unlist(alpha.chem.noname)),nrow=nrow(alpha.chem.noname))
alpha.chem.noname.matrix.small <- alpha.chem.noname.matrix[ , c(1:400)]
alpha.irlba <- irlba(alpha.chem.noname.matrix.small, 8, verbose = TRUE)
alpha.p1 <- prcomp_irlba(alpha.chem.noname.matrix.small, n = 5)
alpha.eigen <- partial_eigen(alpha.chem.noname.matrix.small, n = 5)

alpha.small.noaff <- alpha.chem.noname.matrix.small[ , -1]
alpha.irlba.noaff <- irlba(alpha.small.noaff, 8, verbose = TRUE)
alpha.p1.noaff <- prcomp_irlba(alpha.small.noaff, n = 5)

alpha.log <- log(alpha.small.noaff[ , 1:200])

alpha.matrix.250 <- alpha.chem.noname.matrix[ , 1:250]
alpha.irlba.250 <- irlba(alpha.matrix.250, 8, verbose = T)
alpha.pc.250 <- prcomp_irlba(alpha.matrix.250, n = 5)
summary(alpha.pc.250)

alpha.matrix.300 <- alpha.chem.noname.matrix[ , 1:300]
alpha.irlba.300 <- irlba(alpha.matrix.300, 8, verbose = T)
alpha.pc.300 <- prcomp_irlba(alpha.matrix.300, n = 5)
summary(alpha.pc.300)

alpha.matrix.400 <- alpha.chem.noname.matrix[ , 1:400]
alpha.irlba.400 <- irlba(alpha.matrix.400, 8, verbose = T)
alpha.pc.400 <- prcomp_irlba(alpha.matrix.400, n = 5)
summary(alpha.pc.400)

alpha.matrix.500 <- alpha.chem.noname.matrix[ , 1:500]
alpha.irlba.500 <- irlba(alpha.matrix.500, 8, verbose = T)
alpha.pc.500 <- prcomp_irlba(alpha.matrix.500, n = 5)
summary(alpha.pc.500)

alpha.matrix.450 <- alpha.chem.noname.matrix[ , 1:450]
alpha.irlba.450 <- irlba(alpha.matrix.450, 8, verbose = T)
alpha.pc.450 <- prcomp_irlba(alpha.matrix.450, n = 5)
summary(alpha.pc.450)

alpha.matrix.475 <- alpha.chem.noname.matrix[ , 1:475]
alpha.irlba.475 <- irlba(alpha.matrix.475, 8, verbose = T)
alpha.pc.475 <- prcomp_irlba(alpha.matrix.475, n = 5)
summary(alpha.pc.475)

alpha.stats.pc <- prcomp(alpha.chem.noname.matrix.small[ ,-c(12)], center = T, scale. = T)
autoplot(alpha.stats.pc, alpha.chem.noname.matrix.small)

