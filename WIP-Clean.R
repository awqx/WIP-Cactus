install.packages("stringr")
install.packages("dplyr")
install.packages("Amelia")
# install.packages("caret")
# install.packages("AppliedPredictiveModeling")
# install.packages("lars")
# install.packages("pls")
# install.packages("elasticnet")
install.packages("broom")
install.packages("irlba")
library(stringr)
library(dplyr)
library(ggplot2)
library(Amelia)
# library(MASS)
# library(caret)
# library(AppliedPredictiveModeling)
# library(lars)
# library(pls)
# library(elasticnet)
library(broom)
library(stats)
library(irlba)

# I'm not even sure you need all these packages but can't hurt, right?

# Note to self: add code that creates directories

# ========================================================================
#                  Rekharsky and Inoue Dataset
#                     Initial Downloading
# ========================================================================
source(file = "~/SREP LAB/Rekharsky and Inoue/Cactus/read_cactus.R")
folder  <- "~/SREP LAB/Rekharsky and Inoue/Cactus"
dataset <-
  readRDS(file = "~/SREP LAB/Rekharsky and Inoue/Cactus/RI.rds") 
# Note to self: find original WIP and URL where table can be downloaded

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
# Cactus is a chemical identifier resolver and can be used
# to find and download SD files used by PyRx
# But there's sometimes hiccups in downloading, so I used a tryCatch function
# Solution to function inspired by hfty's answer to
# http://stackoverflow.com/questions/33423725/
# using-trycatch-to-populate-a-data-frame-inside-a-loop-nicely

download.cactus.results <- function(guest, path, chemical.format) {
  report <- tryCatch({
    destfile       <- paste0(path, "/", guest, ".SDF")
    # Chemical format must be parsed to match all the outputs from NCI cactus
    Guest.URL      <- unlist(lapply(guest, URLencode, reserved = T))
    URL            <- paste0(
      "https://cactus.nci.nih.gov/chemical/structure/",
      Guest.URL,
      "/",
      chemical.format
    )
    Map(download.file, url = URL, destfile = destfile)
    data.frame(
      guest = guest,
      downloaded = "yes",
      warning = "no",
      error = "no"
    )
  },
  warning = function(warn) {
    message("Warning: either URL error or already existing directory.")
    destfile       <- paste0(path, "/", guest, ".SDF")
    Guest.URL      <- unlist(lapply(guest, URLencode, reserved = T))
    URL            <- paste0(
      "https://cactus.nci.nih.gov/chemical/structure/",
      Guest.URL, "/", chemical.format
    )
    Map(download.file, url = URL, destfile = destfile)
    data.frame(
      guest = guest,
      downloaded = "yes",
      warning = "yes",
      error = "no"
    )
  },
  error = function(err) {
    message("An error occurred")
    data.frame(
      guest = guest,
      downloaded = "no",
      warning = "yes",
      error = "yes"
    )
    
  },
  finally = {
    message("Chemical name processed")
  })
  return(report)
}

# Find files that are too small to be SDF
filter.filesize <- function(path, pattern, size, logical) {
  filenames <- list.files(path = path, pattern = pattern)
  location  <- dir(path = path,
                   full.names = T,
                   recursive = T)
  data      <- file.info(paste0(path, "/", filenames))[1]
  index     <- eval(call(logical, data, size))
  return(data.frame(
    filepath = location[index],
    molecule = filenames[index],
    stringsAsFactors = F
  ))
}
# Create a directory to store files

create.host.dir <- function(path, host) {
  host.directory <- paste0(path, host)
  dir.create(path = host.directory)
  return(host.directory)
}
# Read SDF files and assign a chemical name -- this makes things easier later
fix.sdf <- function(guest, path){
  report <- tryCatch({
    file <- paste0(path, "/", guest, ".SDF")
    raw.sdf <- read.table(file = file, header = F, sep = "\t", stringsAsFactors = F)
    raw.sdf[1, ] <- paste0(raw.sdf[1,], " ", guest)
    return(raw.sdf)
  }, 
  warning = function(warn){
    message("A warning occurred:")
    message(last.warning)
    file <- paste0(path, "/", guest, ".SDF")
    raw.sdf <- read.table(file = file, header = F, sep = "\t", stringsAsFactors = F)
    raw.sdf[1, ] <- paste0(raw.sdf[1,], " ", guest)
    return(raw.sdf)
  }, 
  error = function(err){
    message("An error occurred")
  })
  return(report)
}

#--------------------------------------
#          Download Gamma
#--------------------------------------

gamma.dest <- create.host.dir(folder, "GammaCD")

results.gamma <-
  do.call(
    rbind,
    lapply(
      gamma.guest,
      download.cactus.results,
      path = gamma.dest,
      chemical.format = "SDF"
    )
  )


results.all <- create.host.dir(folder, "Download Results")

saveRDS(results.gamma, file = paste0(results.all, "/results.gamma.RDS"))

success.gamma <- results.gamma[results.gamma$downloaded == "yes", ]

fail.gamma <- results.gamma[results.gamma$downloaded == "no", ]

empty.gamma <- filter_filesize(
  path = gamma.dest,
  pattern = "SDF",
  size = 100,
  logical = "<"
)

saveRDS(empty.gamma, file = paste0(results.all, "/empty_gamma_sdfiles.RDS"))

file.remove(empty.gamma$filepath)

# This step adjusts names of molecules that failed to download
# Into a form more likely to be used by Cactus
empty.gamma$molecule <- str_replace(string      = empty.gamma$molecule,
                                    pattern     = "\\α",
                                    replacement = "alpha")

# Still don't know what to do with these
results.empty.gamma <-
  do.call(
    rbind,
    lapply(
      empty.gamma$molecule,
      download.cactus.results,
      path = gamma.dest,
      chemical.format = "SDF"
    )
  )

gamma.comb <- rbind(results.empty.gamma, results.gamma)

#--------------------------------------
#         Downloading Beta
#--------------------------------------
# Expect this step to take a while, depending on speed of internet and computer
beta.dest <- create.host.dir(folder, "BetaCD")

results.beta <-
  do.call(
    rbind,
    lapply(
      beta.guest,
      download.cactus.results,
      path = beta.dest,
      chemical.format = "SDF"
    )
  )

empty.beta <-
  filter_filesize(
    path = beta.dest,
    pattern = "SDF",
    size = 100,
    logical = "<"
  )
success.beta <- results.beta[results.beta$downloaded == "yes", ]
file.remove(empty.beta$filepath)

saveRDS(results.beta, file = paste0(beta.dest, "/results.beta.RDS"))

saveRDS(empty.beta, file = paste0(beta.dest, "/empty_beta_sdfiles.RDS"))

# Trying to modify failed downloads to work properly
empty.beta$molecule <-
  str_replace(string = empty.beta$molecule,
              pattern = "\\β",
              replacement = "beta")

empty.beta$molecule <-
  str_replace(string = empty.beta$molecule,
              pattern = "\\α",
              replacement = "alpha")

empty.beta$molecule <-
  str_replace(
    string = empty.beta$molecule,
    pattern =  "\\(\\-*±*[0-9A-Z]*,*
    [0-9A-Z]*\\)\\-\\(*\\+*\\-*\\)*\\-*|nor(?!t)|
    \\([a-z]*\\,*\\s*[a-z]*(I[0-9]\\-)*\\)"
    ,
    replacement = ""
    )

beta.dest.2 <- create.host.dir(paste0(beta.dest), "/2nd.Try")

results.beta.2 <- do.call(
  rbind,
  lapply(
    empty.beta$molecule,
    download.cactus.results,
    path = beta.dest.2,
    chemical.format = "SDF"
  )
)
# Nothing changed. Maybe remove this step?
#--------------------------------------
#          Download Alpha
#--------------------------------------
# Sometimes the code generates a table of just errors, even
# if the moleules download correctly
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

empty.alpha <- filter_filesize(
  path = alpha.dest,
  pattern = "SDF",
  size = 100,
  logical = "<"
)

file.remove(empty.alpha$filepath)

saveRDS(results.alpha, file = paste0(alpha.dest, "/results.alpha.RDS"))

saveRDS(empty.alpha, file = paste0(alpha.dest, "/empty_alpha_sdfiles.RDS"))

empty.alpha$molecule <-
  str_replace(string = empty.alpha$molecule,
              pattern =  "\\(\\-*\\±*[0-9A-Z]*\\,*[0-9A-Z]*\\)\\-\\(*\\+*\\-*\\)*\\-*|nor(?!t)|\\([a-z]*\\,*\\s*[a-z]*(I[0-9]\\-)*\\)",
              replacement = "")

empty.alpha$molecule <-
  str_replace(string = empty.alpha$molecule,
              pattern = "\\β",
              replacement = "beta")

empty.alpha$molecule <-
  str_replace(string = empty.alpha$molecule,
              pattern = "\\α",
              replacement = "alpha")

download_cactus(
  guest = empty.alpha$molecule,
  host = "AlphaCD",
  path = folder,
  chemical.format = "SDF"
)

alpha.dest.2 <- create.host.dir(paste0(alpha.dest), "/2nd.Try")

results.alpha.2 <- do.call(
  rbind,
  lapply(
    empty.alpha$molecule,
    download.cactus.results,
    path = alpha.dest.2,
    chemical.format = "SDF"
  )
)

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
sample.aff.dir <- create.host.dir(folder, "Affinity/Sample")

common.chem.names <- Reduce(intersect, list(alpha.names, beta.names, gamma.names))
set.seed(4)
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
plot(alpha.eval.noperc, main = "Prediction via Linear Model", ylab = "Predicted Values", xlab = "True Observed Values")
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


#==========
# Sort
data.with.aff <- dataset
binding.aff.all <- unlist(Map(convert.kj.kcal, dataset$DelG))
data.with.aff$binding.affinity <- binding.aff.all
View(data.with.aff)
alpha.greatest.aff <- data.with.aff[data.with.aff$host == "1\u03b1", ] %>% 
  arrange (binding.affinity) %>%
  head(10)
beta.greatest.aff <- data.with.aff[data.with.aff$host == "1\u03b2", ] %>% 
  arrange (binding.affinity) %>%
  head(10)
gamma.greatest.aff <- data.with.aff[data.with.aff$host == "1γ", ] %>% 
  arrange (binding.affinity) %>%
  head(10)


#
# find molecules with most solvents
guest.solvents <- data.with.aff %>% 
  group_by(host, guest) %>%
  mutate(number.repeated = n_distinct(log.K)) %>%
  filter(row_number() > 1)
  # filter(row_number() > 1)

#==================================================================
#                    Combined Table
#===============================================================
# Creates a table with all the info from RI, cactus, and PaDel

# View(data.clean)
comb.table <- data.clean[str_detect(data.clean$host, "^1\u03b1|^1\u03b2|^1γ"),]
comb.table <- comb.table[comb.table$pH <= 8 & comb.table$pH >= 5, ]
comb.table <- comb.table[str_detect(comb.table$solvent, "^H2O$|^D2O$"),]
comb.table <- comb.table[is.na(comb.table$solvent.specs), ]
names(comb.table)[names(comb.table) == "binding.affinity"] <- "bind.aff, kcal/mol"

# -------- Testing the duplicates

# Splitting the table into A, B, C
comb.table.gamma <-
  comb.table[str_detect(comb.table$host, "1\u03b1"), ]
comb.table.beta  <-
  comb.table[str_detect(comb.table$host, "1\u03b2"), ]
comb.table.gamma <- 
  comb.table[str_detect(comb.table$host, "1γ"), ]

# Finding duplicates
dup.alpha <-
  comb.table.alpha[duplicated(comb.table.alpha$guest) |
                     duplicated(comb.table.alpha$guest, fromLast = T), ]
dup.beta  <-
  comb.table.beta[duplicated(comb.table.beta$guest) |
                    duplicated(comb.table.beta$guest, fromLast = T), ]
dup.gamma <-
  comb.table.gamma[duplicated(comb.table.gamma$guest) |
                     duplicated(comb.table.gamma$guest, fromLast = T), ]

# Graphing to see trends better
ggplot(dup.alpha, aes(x = dup.alpha$guest, y = dup.alpha$DelG)) +
  geom_point(stat = "identity",
             alpha = 0.4,
             color = "red") +
  scale_y_continuous(limits = c(-30, 0)) +
  geom_boxplot()

ggplot(dup.beta, aes(x = dup.beta$guest, y = dup.beta$DelG)) +
  geom_point(stat = "identity",
             alpha = 0.4,
             color = "blue") +
  scale_y_continuous(limits = c(-30, 0)) +
  geom_boxplot(alpha = 0.2)
ggplot(dup.gamma, aes(x = dup.gamma$guest, y = dup.gamma$DelG)) +
  geom_point(stat = "identity",
             alpha = 0.4,
             color = "green") +
  scale_y_continuous(limits = c(-30, 0)) +
  geom_boxplot()


# dply --> duplicate --> group_by(guest) --> apply max, allpy min --> take difference --> calc measure of centrality and dispersion
# build a function, which one of the group _by is > greater than the average
# boxplot on rage vcxbgf

# Finding Deviation
alpha.summary <- dup.alpha %>%
  group_by(guest) %>% summarize(
    minDelG = min(DelG),
    maxDelG = max(DelG),
    range = maxDelG - minDelG,
    meanDelG = mean(DelG)
  )

beta.summary <- dup.beta %>%
  group_by(guest) %>% summarize(
    minDelG = min(DelG),
    maxDelG = max(DelG),
    range = maxDelG - minDelG, 
    meanDelG = mean(DelG)
  )

gamma.summary <- dup.gamma %>%
  group_by(guest) %>% summarize(
    minDelG = min(DelG),
    maxDelG = max(DelG),
    range = maxDelG - minDelG,
    meanDelG = mean(DelG)
  )

# analyze.duplicates <- function(duplicate, group, col) {
#   summary <- duplicate %>%
#     group_by_(group) %>% summarize_(
#       min = min(col), 
#       max = max(col), 
#       mean  = mean(col),
#       median = median(col)
#     )
#   return (summary)
# }

# alpha.summary <- analyze.duplicates(dup.alpha, "guest", "DelG")
# alpha.summary

alpha.summary.condensed <- alpha.summary %>% summarize(
  min.range = min(range), 
  max.range = max(range), 
  mean.range = mean(range), 
  stand.dev = sd(range), 
  cutoff = mean.range + 3 * stand.dev
)

beta.summary.condensed <- beta.summary %>% summarize(
  min.range = min(range), 
  max.range = max(range), 
  mean.range = mean(range), 
  stand.dev = sd(range), 
  cutoff = mean.range + 3 * stand.dev
)

gamma.summary.condensed <- gamma.summary %>% summarize(
  min.range = min(range), 
  max.range = max(range), 
  mean.range = mean(range), 
  stand.dev = sd(range), 
  cutoff = mean.range + 3 * stand.dev
)

ggplot(alpha.summary, aes(x = guest, y = range)) + 
  geom_point(stat = "identity",
             alpha = 0.4,
             color = "red")

ggplot(beta.summary, aes(x = guest, y = range)) + 
  geom_point(stat = "identity",
             alpha = 0.4,
             color = "blue")

ggplot(gamma.summary, aes(x = guest, y = range)) + 
  geom_point(stat = "identity",
             alpha = 0.4,
             color = "green")


alpha.beta.summary <- rbind(alpha.summary, beta.summary)
alpha.beta.summary$host <- c(rep("alpha", 25), rep("beta", 13))
cutoff <- c(alpha.summary.condensed$cutoff, beta.summary.condensed$cutoff)
range.cutoff <- data.frame(c("alpha", "beta"), cutoff)

ggplot(alpha.beta.summary, aes(x = guest, y = range)) + 
  geom_point(stat = "identity") +
  facet_grid(.~host) + 
  theme(axis.text.x=element_text(angle=90, hjust=1)) + 
  geom_hline(data = range.cutoff, aes(yintercept = cutoff))



# deciding to cut off the outlier
# create a redundant table...jsut in case
dup.alpha.pre <- dup.alpha
dup.alpha <- dup.alpha[-65, ]

# Reducing duplicates to averages
# everything numeric becomes a mean
# characters are collapsed with commas as the division
# NaNs occur when two NAs are averaged
dup.test <- dup.alpha
dup.alpha.clean <- dup.alpha %>%
  group_by(guest) %>%
  summarize(
    solvent = paste0(solvent, collapse = ", "), 
    solvent.specs = paste0(solvent.specs, collapse = ", "),
    Solvent.composition = paste0(Solvent.composition, collapse = ", "), 
    pH = mean(pH), 
    T.K = mean(T.K), 
    log.K = mean(log.K, na.rm = T), 
    log.K.Uncertainty = mean(log.K.Uncertainty, na.rm = T), 
    DelG = mean(DelG), 
    DelG.Uncertainty = mean(DelG.Uncertainty, na.rm = T), 
    DelH = mean(DelH, na.rm = T), 
    DelH.Uncertainty = mean(DelH.Uncertainty, na.rm = T), 
    TDelS = mean(TDelS), 
    TDelS.Uncertainty = mean(TDelS.Uncertainty, na.rm = T), 
    methoda = paste0(methoda, collapse = ", "),
    ref = paste0(ref, collapse = ", "),
    `bind.aff, kcal/mol` = mean(`bind.aff, kcal/mol`)
  )

dup.beta.clean <- dup.beta %>%
  group_by(guest) %>%
  summarize(
    solvent = paste0(solvent, collapse = ", "), 
    solvent.specs = paste0(solvent.specs, collapse = ", "),
    Solvent.composition = paste0(Solvent.composition, collapse = ", "), 
    pH = mean(pH), 
    T.K = mean(T.K), 
    log.K = mean(log.K, na.rm = T), 
    log.K.Uncertainty = mean(log.K.Uncertainty, na.rm = T), 
    DelG = mean(DelG), 
    DelG.Uncertainty = mean(DelG.Uncertainty, na.rm = T), 
    DelH = mean(DelH, na.rm = T), 
    DelH.Uncertainty = mean(DelH.Uncertainty, na.rm = T), 
    TDelS = mean(TDelS), 
    TDelS.Uncertainty = mean(TDelS.Uncertainty, na.rm = T), 
    methoda = paste0(methoda, collapse = ", "),
    ref = paste0(ref, collapse = ", "),
    `bind.aff, kcal/mol` = mean(`bind.aff, kcal/mol`)
  )

dup.gamma.clean <- dup.gamma %>%
  group_by(guest) %>%
  summarize(
    solvent = paste0(solvent, collapse = ", "), 
    solvent.specs = paste0(solvent.specs, collapse = ", "),
    Solvent.composition = paste0(Solvent.composition, collapse = ", "), 
    pH = mean(pH), 
    T.K = mean(T.K), 
    log.K = mean(log.K, na.rm = T), 
    log.K.Uncertainty = mean(log.K.Uncertainty, na.rm = T), 
    DelG = mean(DelG), 
    DelG.Uncertainty = mean(DelG.Uncertainty, na.rm = T), 
    DelH = mean(DelH, na.rm = T), 
    DelH.Uncertainty = mean(DelH.Uncertainty, na.rm = T), 
    TDelS = mean(TDelS), 
    TDelS.Uncertainty = mean(TDelS.Uncertainty, na.rm = T), 
    methoda = paste0(methoda, collapse = ", "),
    ref = paste0(ref, collapse = ", "),
    `bind.aff, kcal/mol` = mean(`bind.aff, kcal/mol`)
  )

# remove duplicates, add cleaned duplicates
cta.test <- comb.table.alpha[-c(which (comb.table.alpha$guest %in% dup.alpha$guest)), ]
dup.alpha.clean <- cbind(host = "1\u03b1", dup.alpha.clean)
comb.table.alpha <- rbind(cta.test, dup.alpha.clean)

ctb.temp <- comb.table.beta[-c(which (comb.table.beta$guest %in% dup.beta$guest)), ]
dup.beta.clean <- cbind(host = "1\u03b2", dup.beta.clean)
comb.table.beta <- rbind(ctb.temp, dup.beta.clean)

ctg.temp <- comb.table.gamma[-c(which (comb.table.gamma$guest %in% dup.gamma$guest)), ]
dup.gamma.clean <- cbind(host = "1\u03b2", dup.gamma.clean)
comb.table.gamma <- rbind(ctg.temp, dup.gamma.clean)

# re-downloading cactus files for convenience
comb.alpha.guest <- comb.table.alpha$guest
comb.beta.guest <- comb.table.beta$guest
comb.gamma.guest <- comb.table.gamma$guest

alpha.comb.dest <- create.host.dir(folder, "AlphaClean")
beta.comb.dest <- create.host.dir(folder, "BetaClean")
gamma.comb.dest <- create.host.dir(folder, "GammaClean")

# results: 2 failures: cis-1,2,3,4-tetraphenylcyclobutane and trans-1,2,3,4-tetraphenylcyclobutane
results.gamma.comb <-
  do.call(
    rbind,
    lapply(
      comb.gamma.guest,
      download.cactus.results,
      path = gamma.comb.dest,
      chemical.format = "SDF"
    )
  )
# l-trp
results.alpha.comb <-
  do.call(
    rbind,
    lapply(
      comb.alpha.guest,
      download.cactus.results,
      path = alpha.comb.dest,
      chemical.format = "SDF"
    )
  )

results.beta.comb <-
  do.call(
    rbind,
    lapply(
      comb.beta.guest,
      download.cactus.results,
      path = beta.comb.dest,
      chemical.format = "SDF"
    )
  )
# Which files don't process through PaDEL
# proadifen
# resorcinol
# adiphenine
# cis-diammine(1,1-cyclobutanedicarboxylato)platinum(II)
# perchlorate
# iodide
# perchloric acid



# Removing molecules that didn't download from the tables

# beta.cactus.fail <- results.beta.comb[results.beta.comb$downloaded == "no", results.beta.comb$guest] 


padel.filepath <- "~/Cactus/PaDEL Results/"
alpha.file <- paste0(padel.filepath, "alpha.csv")
beta.file <- paste0(padel.filepath, "beta.csv")
gamma.file <- paste0(padel.filepath, "gamma.csv")

alpha.padel <- read.csv(alpha.file)
beta.padel <- read.csv(beta.file)
gamma.padel <- read.csv(gamma.file)

names(alpha.padel)[names(alpha.padel) == "Name"] <- "guest"
names(beta.padel)[names(beta.padel) == "Name"] <- "guest"
names(gamma.padel)[names(gamma.padel) == "Name"] <- "guest"

alpha.all <- inner_join(comb.table.alpha, alpha.padel)
beta.all <- inner_join(comb.table.beta, beta.padel)
gamma.all <- inner_join(comb.table.gamma, gamma.padel)

all.clean <- rbind(alpha.all, beta.all, gamma.all)
all.clean <- str_replace_all(all.clean$host, "1\u03b1", "alpha CD")
all.clean <- str_replace(all.clean$host, "1γ", "gamma CD")
write.csv(all.clean, file = "cyclodextrin.csv")

alpha.guest.clean <- unique(data.clean$guest[data.clean$host == "1\u03b1"])
beta.guest.clean  <- unique(data.clean$guest[data.clean$host == "1\u03b2"])

