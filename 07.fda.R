# Libraries and Packages --------------------------------------------------

library(caret)
library(data.table)
library(e1071)
library(Matrix)
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
fda.padel1 <- read.csv("./fda/fda desc 0-3 kb.csv")
fda.padel2 <- read.csv("./fda/fda desc 4-5 kb.csv")
fda.padel3 <- read.csv("./fda/fda desc 6-8 kb.csv")
fda.padel4 <- read.csv("./fda/fda desc 9+ kb.csv")
fda.padel <- rbind(fda.padel1, fda.padel2, fda.padel3, fda.padel3)
# fda.padel <- rbind(fda.padel, fda.padel, fda.padel)

# fda.alpha <- fda.padel %>% mutate(alpha = 1) %>%
#   mutate(beta = 0) %>% mutate(gamma = 0)
# fda.beta <- fda.padel %>% mutate(alpha = 0) %>%
#   mutate(beta = 1) %>% mutate(gamma = 0)
# fda.gamma <- fda.padel %>% mutate(alpha = 0) %>%
#   mutate(beta = 0) %>% mutate(gamma = 1)

pp.settings <- readRDS("./preprocess/preprocess.settings.RDS")
zero.pred <- readRDS("./preprocess/zero.pred.RDS") %>% str_replace(., "-", ".")
zero.pred2 <- readRDS("./preprocess/zero.pred2.RDS")[1]
too.high <- readRDS("./preprocess/high.cor.RDS") %>% str_replace(., "-", ".")

fda.desc <- fda.padel[ , !colnames(fda.padel) %in% zero.pred]

fda.desc <- fda.desc[ , -1] 
fda.desc <- rbind(fda.desc %>% mutate(alpha = 1, beta = 0, gamma = 0), 
                  fda.desc %>% mutate(alpha = 0, beta = 1, gamma = 0), 
                  fda.desc %>% mutate(alpha = 0, beta = 0, gamma = 1))
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

svm.a <- readRDS("./models/svm/polysvm.alpha.RDS")
svm.b <- readRDS("./models/svm/polysvm.beta.RDS")
svm.c <- readRDS("./models/svm/polysvm.gamma.RDS")

fda.svm.a <- predict(svm.a, fda.pp %>% filter(alpha > 0)) %>%
  as.data.frame() %>% mutate(cd.type = "alpha")
fda.svm.b <- predict(svm.b, fda.pp %>% filter(beta > 0)) %>%
  as.data.frame() %>% mutate(cd.type = "beta")
fda.svm.c <- predict(svm.c, fda.pp %>% filter(gamma > 0)) %>%
  as.data.frame() %>% mutate(cd.type = "gamma")

guest <- fda.padel$Name %>% as.vector() %>% rep(., 3)
fda.svm <- data.frame(guest, rbind(fda.svm.a, fda.svm.b, fda.svm.c)) %>%
  rename(., pred = `.`)

# Removing outliers
fda.svm <- fda.svm[-which.max(fda.svm$pred), ]
fda.svm <- fda.svm[-which.max(fda.svm$pred), ]
fda.svm <- fda.svm[-which.min(fda.svm$pred), ]
fda.svm <- fda.svm[-which.min(fda.svm$pred), ]

ggplot(fda.svm, aes(x = guest, y = pred, color = cd.type)) +
  geom_point(alpha = 0.6) + theme_bw() + 
  scale_x_discrete(breaks = NULL) + 
  theme.2018 +
  # geom_hline(yintercept = -24, linetype = "dotted", color = "red", size = 1) +
  # geom_hline(yintercept = -20, linetype = "dotted", color = "orange", size = 1) +
  labs(title = "Polynomial SVM - FDA Approved Drugs", 
       y = "Predicted DelG, kJ/mol", 
       color = "CD Type", shape = "Cyclodextrin", 
       x = NULL)


#     Cubist --------------------------------------------------------------

cube.a <- readRDS("./models/cubist/cubist.alpha.RDS")
cube.b <- readRDS("./models/cubist/cubist.beta.RDS")
cube.c <- readRDS("./models/cubist/cubist.gamma.RDS")

fda.cube.a <- predict(cube.a, fda.pp %>% filter(alpha > 0)) %>%
  as.data.frame() %>% mutate(cd.type = "alpha")
fda.cube.b <- predict(cube.b, fda.pp %>% filter(beta > 0)) %>%
  as.data.frame() %>% mutate(cd.type = "beta")
fda.cube.c <- predict(cube.c, fda.pp %>% filter(gamma > 0)) %>%
  as.data.frame() %>% mutate(cd.type = "gamma")

fda.cube <- data.frame(guest, rbind(fda.cube.a, fda.cube.b, fda.cube.c)) %>%
  rename(., pred = `.`)

ggplot(fda.cube, aes(x = guest, y = pred, color = cd.type)) +
  geom_point(alpha = 0.6) + theme_bw() + 
  scale_x_discrete(breaks = NULL) + 
  theme.2018 +
  # geom_hline(yintercept = -24, linetype = "dotted", color = "red", size = 1) +
  # geom_hline(yintercept = -20, linetype = "dotted", color = "orange", size = 1) +
  labs(title = "Cubist - FDA Approved Drugs", 
       y = "Predicted DelG, kJ/mol", 
       color = "CD Type", shape = "Cyclodextrin", 
       x = NULL)


#     GLMnet --------------------------------------------------------------

glm.a <- readRDS("./models/glmnet/glm.alpha.RDS")
glm.b <- readRDS("./models/glmnet/glm.beta.RDS")
glm.c <- readRDS("./models/glmnet/glm.gamma.RDS")

fda.a <- fda.pp %>% filter(alpha > 0) %>% as.matrix()
fda.b <- fda.pp %>% filter(beta > 0) %>% as.matrix()
fda.c <- fda.pp %>% filter(gamma > 0) %>% as.matrix()
fda.glm.a <- predict.glmnet(glm.a, fda.a, s = tail(glm.a$lambda, n = 1)) %>%
  as.data.frame() %>% mutate(cd.type = "alpha")
fda.glm.b <- predict.glmnet(glm.b, fda.b, s = tail(glm.b$lambda, n = 1)) %>%
  as.data.frame() %>% mutate(cd.type = "beta")
fda.glm.c <- predict.glmnet(glm.c, fda.c, s = tail(glm.c$lambda, n = 1)) %>%
  as.data.frame() %>% mutate(cd.type = "gamma")

fda.glm <- data.frame(guest, rbind(fda.glm.a, fda.glm.b, fda.glm.c)) %>%
  rename(., pred = X1)

ggplot(fda.glm, aes(x = guest, y = pred, color = cd.type)) +
  geom_point(alpha = 0.6) + theme_bw() + 
  scale_x_discrete(breaks = NULL) + 
  theme.2018 +
  # geom_hline(yintercept = -24, linetype = "dotted", color = "red", size = 1) +
  # geom_hline(yintercept = -20, linetype = "dotted", color = "orange", size = 1) +
  labs(title = "GLMNet FDA Approved Drugs", 
       y = "Predicted DelG, kJ/mol", 
       color = "CD Type", shape = "Cyclodextrin", 
       x = NULL)

#     PLS -----------------------------------------------------------------

pls.a <- readRDS("./models/pls/pls.alpha.RDS")
pls.b <- readRDS("./models/pls/pls.beta.RDS")
pls.c <- readRDS("./models/pls/pls.gamma.RDS")

fda.pls.a <- predict(pls.a, fda.pp %>% filter(alpha > 0), ncomp = 4) %>%
  as.data.frame() %>% mutate(cd.type = "alpha") %>%
  rename(pred = `DelG.4 comps`)
fda.pls.b <- predict(pls.b, fda.pp %>% filter(beta > 0), ncomp = 4) %>%
  as.data.frame() %>% mutate(cd.type = "beta") %>%
  rename(pred = `DelG.4 comps`)
fda.pls.c <- predict(pls.c, fda.pp %>% filter(gamma > 0), ncomp = 5) %>%
  as.data.frame() %>% mutate(cd.type = "gamma") %>%
  rename(pred = `DelG.5 comps`)

fda.pls <- data.frame(guest, rbind(fda.pls.a, fda.pls.b, fda.pls.c))

ggplot(fda.pls, aes(x = guest, y = pred, color = cd.type)) +
  geom_point(alpha = 0.6) + theme_bw() + 
  scale_x_discrete(breaks = NULL) + 
  theme.2018 +
  # geom_hline(yintercept = -24, linetype = "dotted", color = "red", size = 1) +
  # geom_hline(yintercept = -20, linetype = "dotted", color = "orange", size = 1) +
  labs(title = "Cubist - FDA Approved Drugs", 
       y = "Predicted DelG, kJ/mol", 
       color = "CD Type", shape = "Cyclodextrin", 
       x = NULL)

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

# Compiling ---------------------------------------------------------------

fda.comb <- rbind(fda.cube %>% mutate(model = "Cubist"), 
                  fda.glm %>% mutate(model = "GLMNet"), 
                  fda.svm %>% mutate(model = "SVM"))
ggplot(fda.comb, aes(x = guest, y = pred, color = cd.type)) +
  geom_point() + 
  theme.2018 + 
  facet_wrap(~model) + 
  labs(title = "QSPR predictions on FDA-approved drugs", 
       x = NULL, 
       y = "Predicted dG, kJ/mol", 
       color = "CD Type") + 
  ylim(-60, 25)  +
  scale_x_discrete(breaks = NULL) + 
  theme(
    axis.ticks.x=element_blank(),
    axis.line = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
  ) 
ggsave("./graphs/fda.png", dpi = 600)
