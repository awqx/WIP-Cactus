# I'm not even sure you need all these packages but can't hurt, right?
install.packages("stringr")
install.packages("dplyr")
install.packages("Amelia")
install.packages("caret")
install.packages("AppliedPredictiveModeling")
install.packages("lars")
install.packages("pls")
install.packages("elasticnet")
install.packages("broom")
install.packages("irlba")
library(stringr)
library(dplyr)
library(ggplot2)
library(Amelia)
library(MASS)
library(caret)
library(AppliedPredictiveModeling)
library(lars)
library(pls)
library(elasticnet)
library(broom)
library(stats)
library(irlba)

# ========================================================================
#                  Rekharsky and Inoue Dataset
#                     Initial Downloading
# ========================================================================
# download the cactus.functions.R file from GitHub, put it into this folder
dir.create("~/Cactus")
source(file = "~/Cactus/cactus.functions.R")
folder  <- "~/Cactus/"
dataset <-
  readRDS(file = "~/Cactus/RI.rds") 
# Edited note: apparently I can't find the URL to download the dataset successfully, 
# so "dataset" has to be read from a pre-prepared rds file

# Fixing chemical names that affect URL at download-cactus stage
# Typos
dataset[48, "guest"] <- "biebrich scarlet"
dataset[1081, "guest"] <- "4-[(4-hydroxyphenyl)azo]benzoate"
dataset[23, "guest"] <- "3-(aminomethyl)-proxyl"
dataset[82, "guest"] <- "3-carbamoyl-proxyl"
# Cactus uses alternate identifying name for the following:
dataset[42, "guest"] <- "1-o-benzyl-rac-glycerol"
dataset[126, "guest"] <-
  "cis-diammine(1,1-cyclobutanedicarboxylato)platinum(II)"
dataset[c(142, 143, 727, 728), "guest"] <-
  "[(1R,2S)-1-hydroxy-1-phenylpropan-2-yl]-methylazanium"
dataset[c(144, 145, 729, 730), "guest"] <-
  "[(1S,2R)-1-hydroxy-1-phenylpropan-2-yl]-methylazanium"

data.clean <- dataset[dataset$T.K == 298|dataset$T.K == 303,]
data.clean <- data.clean[!str_detect(data.clean$ref, "b|c|g|i|j|m"),]
# Unicode detects the characters alpha and beta, which won't show otherwise
# This varies from computer to computer
# Note to self -- create error handling for this
# Alpha == \u03b1, beta == \u03b2
alpha.guest.clean <- unique(data.clean$guest[data.clean$host == "1\u03b1"])
beta.guest.clean  <- unique(data.clean$guest[data.clean$host == "1\u03b2"])
gamma.guest.clean <- unique(data.clean$guest[data.clean$host == "1γ"]) # Gamma works w/o unicode

# ========================================================================
#                          Cactus Downloading
# ========================================================================

#--------------------------------------
#          Download Gamma
#--------------------------------------

gamma.dest <- create.host.dir(folder, "GammaCD")

gamma.results <-
  do.call(
    rbind,
    lapply(
      gamma.guest,
      download.cactus.results,
      path = gamma.dest,
      chemical.format = "SDF"
    )
  )

all.results <- create.host.dir(folder, "Download Results")

saveRDS(gamma.results, file = paste0(all.results, "/gamma.results.RDS"))

success.gamma <- gamma.results[gamma.results$downloaded == "yes", ]

gamma.fail <- gamma.results[gamma.results$downloaded == "no", ]

gamma.empty <- filter.filesize(
  path = gamma.dest,
  pattern = "SDF",
  size = 100,
  logical = "<"
)

saveRDS(gamma.empty, file = paste0(all.results, "/empty_gamma_sdfiles.RDS"))

file.remove(gamma.empty$filepath)

# This step adjusts names of molecules that failed to download
# Into a form more likely to be used by Cactus
gamma.empty$molecule <- str_replace(string      = gamma.empty$molecule,
                                    pattern     = "\\α",
                                    replacement = "alpha")

# Still don't know what to do with these
gamma.empty.results <-
  do.call(
    rbind,
    lapply(
      gamma.empty$molecule,
      download.cactus.results,
      path = gamma.dest,
      chemical.format = "SDF"
    )
  )

gamma.comb <- rbind(gamma.empty.results, gamma.results)

#--------------------------------------
#         Downloading Beta
#--------------------------------------
# Expect this step to take a while, depending on speed of internet and computer
beta.dest <- create.host.dir(folder, "BetaCD")

beta.results <-
  do.call(
    rbind,
    lapply(
      beta.guest,
      download.cactus.results,
      path = beta.dest,
      chemical.format = "SDF"
    )
  )

beta.empty <-
  filter.filesize(
    path = beta.dest,
    pattern = "SDF",
    size = 100,
    logical = "<"
  )
success.beta <- beta.results[beta.results$downloaded == "yes", ]
file.remove(beta.empty$filepath)

saveRDS(beta.results, file = paste0(beta.dest, "/beta.results.RDS"))

saveRDS(beta.empty, file = paste0(beta.dest, "/empty_beta_sdfiles.RDS"))

# Trying to modify failed downloads to work properly
beta.empty$molecule <-
  str_replace(string = beta.empty$molecule,
              pattern = "\\β",
              replacement = "beta")

beta.empty$molecule <-
  str_replace(string = beta.empty$molecule,
              pattern = "\\α",
              replacement = "alpha")

beta.empty$molecule <-
  str_replace(
    string = beta.empty$molecule,
    pattern =  "\\(\\-*±*[0-9A-Z]*,*
    [0-9A-Z]*\\)\\-\\(*\\+*\\-*\\)*\\-*|nor(?!t)|
    \\([a-z]*\\,*\\s*[a-z]*(I[0-9]\\-)*\\)"
    ,
    replacement = ""
    )

beta.dest.2 <- create.host.dir(paste0(beta.dest), "/2nd.Try")

beta.results.2 <- do.call(
  rbind,
  lapply(
    beta.empty$molecule,
    download.cactus.results,
    path = beta.dest.2,
    chemical.format = "SDF"
  )
)
# Nothing changed. Maybe remove this step?
#--------------------------------------
#          Download Alpha
#--------------------------------------
alpha.dest <- create.host.dir(folder, "AlphaCD")

results.alpha <-
  do.call(
    rbind,
    lapply(
      alpha.guest,
      download.cactus.results,
      path = alpha.dest,
      chemical.format = "SDF"
    )
  )

success.alpha <- results.alpha[results.alpha$downloaded == "yes", ]

fail.alpha <- results.alpha[results.alpha$downloaded == "no", ]

empty.alpha <- filter.filesize(
  path = alpha.dest,
  pattern = "SDF",
  size = 100,
  logical = "<"
)

file.remove(empty.alpha$filepath)

saveRDS(results.alpha, file = paste0(alpha.dest, "/results.alpha.RDS"))

saveRDS(empty.alpha, file = paste0(alpha.dest, "/empty_alpha_sdfiles.RDS"))
# Not sure if the following code does anything, so I'm commenting it out for now
# empty.alpha$molecule <-
# str_replace(string = empty.alpha$molecule,
#              pattern =  "\\(\\-*\\±*[0-9A-Z]*\\,*[0-9A-Z]*\\)\\-\\(*\\+*\\-*\\)*\\-*|nor(?!t)|\\([a-z]*\\,*\\s*[a-z]*(I[0-9]\\-)*\\)",
#              replacement = "")

# empty.alpha$molecule <-
#  str_replace(string = empty.alpha$molecule,
#              pattern = "\\β",
#              replacement = "beta")

# empty.alpha$molecule <-
#  str_replace(string = empty.alpha$molecule,
#              pattern = "\\α",
#              replacement = "alpha")

# download_cactus(
#  guest = empty.alpha$molecule,
#  host = "AlphaCD",
#  path = folder,
#  chemical.format = "SDF"
#)

# alpha.dest.2 <- create.host.dir(paste0(alpha.dest), "/2nd.Try")

# results.alpha.2 <- do.call(
#  rbind,
#  lapply(
#    empty.alpha$molecule,
#    download.cactus.results,
#    path = alpha.dest.2,
#    chemical.format = "SDF"
#  )
# )
#=========================================================================
# Combining all the SDF files into one file
# For Gamma
all.sdf <- create.host.dir(folder, "All.SDF")

gamma.success <- as.character(success.gamma$guest)
all.gamma <- do.call(bind_rows, lapply(gamma.success, fix.sdf, path = gamma.dest))
write.table(
  all.gamma,
  file = paste0(all.sdf, "/allgamma.SDF"),
  append = TRUE,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

# For Beta
beta.success <- as.character(success.beta$guest)
e.dansyl <- fix.sdf("e-dansyl-l-lysine", beta.dest)
d.valero <- fix.sdf("d-valerolactam", beta.dest)
all.beta <- do.call(bind_rows, lapply(beta.success, fix.sdf, path = beta.dest))
all.beta <- rbind(e.dansyl, d.valero, all.beta)
write.table(
  all.beta,
  file = paste0(all.sdf, "/allbeta.SDF"),
  append = TRUE,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)
# For Alpha
alpha.success <- as.character(success.alpha$guest)
a.methyl1 <- fix.sdf("(+)-a-methylbenzylamine", alpha.dest)
a.methyl2 <- fix.sdf("(-)-a-methylbenzylamine", alpha.dest)
d.valero.a <- fix.sdf("d-valerolactam", alpha.dest)
all.alpha <- do.call(bind_rows, lapply(alpha.success, fix.sdf, path = alpha.dest))
all.alpha <- rbind(a.methyl1, a.methyl2, d.valero.a, all.alpha)
write.table(
  all.alpha,
  file = paste0(all.sdf, "/allalpha.SDF"),
  append = TRUE,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

# Creating a small sample
set.seed(4)
sample.aff.dir <- create.host.dir(folder, "Affinity/Sample")
common.chem.names <- Reduce(intersect, list(alpha.names, beta.names, gamma.names))
sample.chem.list <- sample(common.chem.names, 3)

# Sample should include cinnarizine, (E)-stilbene, and 4-nitrophenol
sample.files <- c(paste0(gamma.dest, "/", sample.chem.list, ".SDF"))
sample.chem <- lapply(sample.files, read.table, header = FALSE, sep = "\t")
sample.chem <- bind_rows(sample.chem)
write.table(
  sample.chem,
  file = paste0(sample.aff.dir, "/sample.SDF"),
  append = TRUE,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)






#==============================================================================
#                              Missing Data
#==============================================================================
# The alpha, beta, and gamma guests from Rekharsky and Inoue will serve as a 
# "full" data set
all.guests <- unique(c(alpha.guest, beta.guest, gamma.guest))
blank.space <- rep(NA, 555)
thing <- data.frame(all.guests, blank.space)
thing2 <- data.frame(all.guests.clean.na, blank.space)
thing2$all.guests <- all.guests.clean.na
thing2 <- thing2[ , c(3, 2, 1)]
thing3 <- full_join(thing, thing2, by = "all.guests")
thing3 <- thing3[-c(556:574), -c(2,3)]
guest.success <- unique(c(alpha.success, beta.success, gamma.success))
guest.success.df <- data.frame(guest.success)
guest.success.df$all.guests <- guest.success
thing3 <- full_join(thing3, guest.success.df)
missmap(thing3, rank.order = F, col =  c("seagreen1", "slateblue"), main = "Data Loss per Step")


#===============================================================================
#                       Predictive Modeling -- Linear Model
#===============================================================================
# Results from PaDEL Chem Descriptor Lab
# For Alpha
# Organizing a data table on experimental solubility
binding.aff.clean <- unlist(Map(convert.kj.kcal, data.clean$DelG))
data.clean$binding.affinity <- binding.aff.clean
alpha.exp.aff <- data.clean[data.clean$host == "1\u03b1",c("guest", "binding.affinity")]
names(alpha.exp.aff)[names(alpha.exp.aff) == "guest"] <- "Name"
# Importing the results
alpha.results.folder <- paste0(folder, "PaDEL Results/")
alpha.chem <- read.csv(paste0(alpha.results.folder, "alpha.mmff94.csv"), header = T, stringsAsFactors = F)
alpha.chem <- full_join(alpha.exp.aff, alpha.chem)
alpha.chem.clean <- alpha.chem[!is.na(alpha.chem$Name),]
alpha.chem.clean <- alpha.chem.clean[!is.na(alpha.chem.clean$binding.affinity), ]
alpha.chem.clean <- alpha.chem.clean[!is.na(alpha.chem.clean$nAcid), ]
alpha.chem.noname <- alpha.chem.clean
alpha.chem.noname$Name <- NULL 
# random sample of clean data for training
set.seed(1)
alpha.sample.size <- nrow(alpha.chem.clean) * 3 / 4
alpha.train.names <- sample(alpha.chem.clean$Name, alpha.sample.size, replace = F)
alpha.train <- alpha.chem.clean[alpha.chem.clean$Name %in% alpha.train.names, ]
alpha.train$Name <- NULL
alpha.test <- alpha.chem.clean[!alpha.chem.clean$Name %in% alpha.train.names, ]
alpha.test$Name <- NULL
lm.alpha.pred <- lm(binding.affinity ~ ., data = alpha.train)
summary(lm.alpha.pred) 
tidy(lm.alpha.pred) 
alpha.pred <- predict.lm(lm.alpha.pred, alpha.test, level = 0.95)
alpha.results <- alpha.chem.clean[!alpha.chem.clean$Name %in% alpha.train.names, ]
alpha.results <- alpha.results[ , colnames(alpha.results) == "Name" | colnames(alpha.results) == "binding.affinity"]
alpha.eval <- data.frame(alpha.results, pred = alpSha.pred)
alpha.eval <- alpha.eval[!is.na(alpha.eval$pred), ]
alpha.eval.noperc <- alpha.eval
alpha.eval.noperc$percent <- NULL
defaultSummary(alpha.eval)
plot(alpha.eval.noperc, main = "Solubility Values", ylab = "Predicted Values", xlab = "True Observed Values")
lines(c(-10, 1), c(-10, 1), col = "red")
summary(alpha.eval)
boxcox(lm., family="yjPower", plotit = T)

alpha.int <- sapply(alpha.train, is.integer)
alpha.train.int <- alpha.train[ , alpha.int]
alpha.dup.col <- sapply(alpha.train.int, n_distinct) == 1
alpha.dup <- alpha.train.int[ , alpha.dup.col]
alpha.dup.table <- data.frame(sapply(alpha.train.int, n_distinct))
alpha.dup.table <- summarise_each(alpha.train.int, n_distinct)
ggplot(alpha.dup.table, aes())
gather(alpha.dup.table) %>% 
  arrange(., value) %>% 
  ggplot(., aes(reorder(key, value), y = value)) + 
    geom_bar(stat="identity")

# Percent Yields
alpha.eval$percent <- alpha.eval$pred / alpha.eval$alpha.results
ggplot(alpha.eval, aes(x = row.names(alpha.eval), y = alpha.eval$percent)) +
  geom_point() +
  scale_y_log10()

ggplot(dataset, aes(x = dataset$DelG)) +
  geom_histogram() +
  xlab("DelG")

# irlba
alpha.chem.noname.matrix <- data.matrix(alpha.chem.noname)
irlba(alpha.chem.noname.matrix, nv = 5, maxit = 1000, work = nv + 7, reorth = TRUE,
      tol = 1e-05, v = NULL, right_only = FALSE, verbose = FALSE, scale)
