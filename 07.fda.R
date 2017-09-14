# Libraries and Packages --------------------------------------------------

library(stringr)
library(tidyverse)

# Functions ---------------------------------------------------------------

# Finds standard deviation for a single descriptor
# Requires a vector or single column; returns num
find.sd.desc <- function(data) {
  sd <- (data - mean(data)) ^ 2 %>% sum()
  sd <- sqrt(sd / (length(data) - 1))
  return(sd)
}

# Determines whether a chemical is within the applicability domain
# Requires a df or matrix and an index
# Can be used in do.call, rbind, lapply sequence
determine.domain <- function(index, data){
  results <- c(rep(0.0, length(data)  - 1))
  for (i in 2:length(data)) {
    sd <- find.sd.desc(data[ , i])
    results[i - 1] <- (data[index, i] - mean(data[ , i])) / sd
  }
}

# Precondition: first column of data is guests, rest is descriptors
# initial.standardize works on the data the model was trained on
initial.standardize <- function(data) {
  df <- data[ , -1]
  guest <- data[ , 1]
  for (c in 1:length(df)) {
    sd <- find.sd.desc(df[ , c])
    for (r in 1:nrow(df)) {
      ski <- abs(df[r, c] - mean(df[ , c])) / sd 
      df[r, c] <- ski
    }
    # message(paste("Column ", c, " completed."))
  }
  return(cbind(guest, df))
}

# Standardize works on new data
# sd.list should be retrieved from the data the model was trained on
standardize <- function(data, sd.list) {
  df <- data[ , -1]
  guest <- data[ , 1]
  for (c in 1:length(df)) {
    sd <- sd.list[c]
    for (r in 1:nrow(df)) {
      ski <- abs(df[r, c] - mean(df[ , c])) / sd 
      df[r, c] <- ski
    }
    # message(paste("Column ", c, " completed."))
  }
  return(cbind(guest, df))
}

# Precondition: data is the result of standardize
# guest is in first col, rest are numeric
domain <- function(data) {
  guest <- data[ , 1]
  results <- c(rep(NA, nrow(data)))
  df <- data[ , -1]
  for(r in 1:nrow(df)) {
    if (max(df[r, ] <= 3)) {
      results[r] <- "inside"
    } else if (min(df[r, ] > 3)){
      results[r] <- "outside"
    } else {
      newSk <- mean(as.numeric(df[r, ])) + 1.28 * find.sd.desc(as.numeric(df[r, ]))
      if (newSk <= 3)
        results[r] <- "inside"
      else
        results[r] <- "outside"
    }
  }
  return(data.frame(guest, results))
}

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

# Data --------------------------------------------------------------------

# Downloaded from http://zinc15.docking.org/substances/subsets/fda/, 
# analyzed via PaDEL-Descriptor
dir.create("./molecules/fdaSDFs")

fda <- read.csv("./molecules/fda/fda.csv")
fda.guest <- fda[ , 2]

guest <- fda.guest %>% as.vector()
guest1 <- fda.guest[1:100] %>% as.vector()

cactus1 <- 
  do.call(
    rbind,
    lapply(
      guest1,
      download.cactus.results,
      path = "./molecules/fdaSDFs",
      chemical.format = "SDF"
    )
  )

fda.results <- 
  do.call(
    rbind,
    lapply(
      guest,
      download.cactus.results,
      path = "./molecules/fdaSDFv2",
      chemical.format = "SDF"
    )
  )


# PaDEL-Descriptor --------------------------------------------------------

# Loading raw data
fda.padel1 <- read.csv("./molecules/fda/fda desc 0-3 kb.csv")
fda.padel2 <- read.csv("./molecules/fda/fda desc 4-5 kb.csv")
# fda.padel3 <- read.csv("./molecules/fda/fda desc 6-8 kb.csv")
# fda.padel4 <- read.csv("./molecules/fda/fda desc 9+ kb.csv")
fda.padel <- rbind(fda.padel1, fda.padel2)
fda.padel <- rbind(fda.padel1, fda.padel2, fda.padel3, fda.padel3)
fda.padel <- rbind(fda.padel, fda.padel, fda.padel)

# fda.alpha <- fda.padel %>% mutate(alpha = 1) %>%
#   mutate(beta = 0) %>% mutate(gamma = 0)
# fda.beta <- fda.padel %>% mutate(alpha = 0) %>%
#   mutate(beta = 1) %>% mutate(gamma = 0)
# fda.gamma <- fda.padel %>% mutate(alpha = 0) %>%
#   mutate(beta = 0) %>% mutate(gamma = 1)

pp.settings <- readRDS("./preprocess.settings.RDS")
zero.pred <- readRDS("./zero.pred.RDS") %>% str_replace(., "-", ".")
zero.pred2 <- readRDS("./zero.pred2.RDS")[1]
too.high <- readRDS("./high.cor.RDS") %>% str_replace(., "-", ".")

fda.desc <- fda.padel[ , !colnames(fda.padel) %in% zero.pred]

fda.desc <- fda.desc[ , -1] 
fda.desc <- fda.desc %>% mutate(alpha = c(rep(1, nrow(fda.desc)/3), rep(0, nrow(fda.desc)/3 * 2))) %>%
  mutate(beta = c(rep(0, nrow(fda.desc)/3), rep(1, nrow(fda.desc)/3), rep(0,nrow(fda.desc)/3))) %>%
  mutate(gamma = c(rep(0,nrow(fda.desc)/3 * 2), rep(1, nrow(fda.desc)/3)))

fda.pp <- predict(pp.settings, fda.desc)

fda.pp <- fda.pp[ , !colnames(fda.pp) %in% too.high]
fda.pp <- fda.pp[ , !colnames(fda.pp) %in% zero.pred2]

fda.data <- cbind(fda.padel[ , 1], fda.pp)
colnames(fda.data)[1] <- "guest"

#     Applicability Domain Analysis ---------------------------------------

sdevs <- readRDS("./domain/trn.sdevs.RDS")
fda.stand <- standardize(fda.data, sdevs)
fda.dmn <- domain(fda.stand)

# Models ------------------------------------------------------------------

#     SVM -----------------------------------------------------------------

svm.alpha <- readRDS("./models/svm/polysvm.alpha.RDS")
svm.beta <- readRDS("./models/svm/polysvm.beta.RDS")
svm.gamma <- readRDS("./models/svm/polysvm.gamma.RDS")

fda.alpha <- predict(svm.alpha, fda.pp)[1:830]
fda.beta <- predict(svm.beta, fda.pp)[831:1660]
fda.gamma <- predict(svm.gamma, fda.pp)[1661:2490]

fda.svm <- c(fda.alpha, fda.beta, fda.gamma)
guest <- fda.padel$Name %>% as.vector()

fda.svm <- data.frame(guest, fda.svm) %>%
  dplyr::rename(., pred = fda.svm) %>%
  mutate(cd.type = c(rep("alpha", 830), rep("beta", 830), rep("gamma", 830)))

ggplot(fda.svm, aes(x = guest, y = pred, color = cd.type, shape = cd.type)) +
  geom_point(alpha = 0.6) + theme_bw() + 
  scale_x_discrete(breaks = NULL) + 
  theme(
    axis.ticks.x=element_blank(),
    axis.line = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
  ) + 
  geom_hline(yintercept = -24, linetype = "dotted", color = "red", size = 1) +
  geom_hline(yintercept = -20, linetype = "dotted", color = "orange", size = 1) +
  labs(title = "Polynomial SVM - FDA Approved Drugs", 
       y = "Predicted DelG, kJ/mol", 
       color = "Cyclodextrin", shape = "Cyclodextrin")


#     Cubist --------------------------------------------------------------

cube.alpha <- readRDS("./models/cubist/cube.alpha.RDS")
cube.beta <- readRDS("./models/cubist/cube.beta.RDS")
cube.gamma <- readRDS("./models/cubist/cube.gamma.RDS")

fda.alpha <- predict(cube.alpha, fda.pp)[1:830]
fda.beta <- predict(cube.beta, fda.pp)[831:1660]
fda.gamma <- predict(cube.gamma, fda.pp)[1661:2490]

fda.cube <- c(fda.alpha, fda.beta, fda.gamma)
fda.cube <- data.frame(guest, fda.cube) %>%
  dplyr::rename(., pred = fda.cube) %>%
  mutate(cd.type = c(rep("alpha", 830), rep("beta", 830), rep("gamma", 830)))

ggplot(fda.cube, aes(x = guest, y = pred, color = cd.type, shape = cd.type)) +
  geom_point(alpha = 0.6) + theme_bw() + 
  scale_x_discrete(breaks = NULL) + 
  theme(
    axis.ticks.x=element_blank(),
    axis.line = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
  ) + 
  geom_hline(yintercept = -24, linetype = "dotted", color = "red", size = 1) +
  geom_hline(yintercept = -20, linetype = "dotted", color = "orange", size = 1) +
  labs(title = "Cubist - FDA Approved Drugs", 
       y = "Predicted DelG, kJ/mol", 
       color = "Cyclodextrin", shape = "Cyclodextrin")

#     GLMnet --------------------------------------------------------------

glm.alpha <- readRDS("./models/glmnet/glm.alpha.RDS")
glm.beta <- readRDS("./models/glmnet/glm.beta.RDS")
glm.gamma <- readRDS("./models/glmnet/glm.gamma.RDS")

fda.mat <- as.matrix(fda.pp)
fda.alpha <- predict.glmnet(glm.alpha, fda.mat, 
                            s = tail(glm.alpha$lambda, n = 1))[1:830]
fda.beta <- predict.glmnet(glm.beta, fda.mat, 
                           s = tail(glm.beta$lambda, n = 1))[831:1660]
fda.gamma <- predict.glmnet(glm.gamma, fda.mat, 
                            s = tail(glm.gamma$lambda, n = 1))[1661:2490]

fda.glm <- c(fda.alpha, fda.beta, fda.gamma)
fda.glm <- data.frame(guest, fda.glm) %>%
  dplyr::rename(., pred = fda.glm) %>%
  mutate(cd.type = c(rep("alpha", 830), rep("beta", 830), rep("gamma", 830)))
ggplot(fda.glm, aes(x = guest, y = pred, color = cd.type, shape = cd.type)) +
  geom_point(alpha = 0.6) + theme_bw() + 
  scale_x_discrete(breaks = NULL) + 
  theme(
    axis.ticks.x=element_blank(),
    axis.line = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
  ) + 
  geom_hline(yintercept = -24, linetype = "dotted", color = "red", size = 1) +
  geom_hline(yintercept = -20, linetype = "dotted", color = "orange", size = 1) +
  labs(title = "GLMnet - FDA Approved Drugs", 
       y = "Predicted DelG, kJ/mol", 
       color = "Cyclodextrin", shape = "Cyclodextrin")

#     Glmnet and SVM ------------------------------------------------------

glm.svm.fda <- rbind(fda.svm %>% mutate(Model = "SVM"), 
                     fda.glm %>% mutate(Model = "GLMnet"))
glm.svm.nogamma <- glm.svm.fda[glm.svm.fda$cd.type != "gamma", ]

ggplot(glm.svm.fda, aes(x = guest, y = pred, color = cd.type, shape = cd.type)) +
  geom_point(alpha = 0.6) + theme_bw() + 
  scale_x_discrete(breaks = NULL) + 
  theme(
    axis.ticks.x=element_blank(),
    axis.line = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
  ) + 
  facet_wrap(~Model) + 
  geom_hline(yintercept = -24, linetype = "dotted", color = "red", size = 1) +
  geom_hline(yintercept = -20, linetype = "dotted", color = "orange", size = 1) +
  labs(title = "GLMnet and SVM - FDA Approved Drugs", 
       y = "Predicted DelG, kJ/mol", 
       color = "Cyclodextrin", shape = "Cyclodextrin")

ggplot(glm.svm.nogamma, aes(x = guest, y = pred, color = cd.type, shape = cd.type)) +
  geom_point(alpha = 0.6) + theme_bw() + 
  scale_x_discrete(breaks = NULL) + 
  theme(
    axis.ticks.x=element_blank(),
    axis.line = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
  ) + 
  facet_wrap(~Model) + 
  geom_hline(yintercept = -24, linetype = "dotted", color = "red", size = 1) +
  geom_hline(yintercept = -20, linetype = "dotted", color = "orange", size = 1) +
  labs(title = "GLMnet and SVM - FDA Approved Drugs, No Gamma-CD", 
       y = "Predicted DelG, kJ/mol", 
       color = "Cyclodextrin", shape = "Cyclodextrin")
