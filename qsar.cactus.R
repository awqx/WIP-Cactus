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
