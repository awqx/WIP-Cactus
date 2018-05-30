# This script should be sourced for everything after retrieving PaDEL-Descriptor results.
# It includes Cactus, pre-processing and models for Cubist, GLMNet, PLS, random forest, and SVM

library(caret)
library(e1071)
library(glmnet)
library(RCurl)
library(stringr)
library(tidyverse)
library(XML)

# Cactus ------------------------------------------------------------------

make.regex <- function(string) {
  new.string <- str_replace(string, pattern = "\\(", 
                            replacement = "\\\\(")
  new.string <- str_replace(new.string, pattern = "\\)", 
                            replacement = "\\\\)")
  # new.string <- str_replace(new.string, pattern = "\\Î²|\\B\\s+H\\+|\u03b2", # alternatives: \\Î²|\\B\\s+H\\+
  #                           replacement = "\\(Î²|\u03b2\\)")
  return(new.string)
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

# Pre-processing ----------------------------------------------------------

pp.settings <- readRDS("./pre-process/pp.settings.RDS")
zero.pred <- readRDS("./pre-process/zero.pred.RDS") %>% str_replace(., "-", ".")
zero.pred2 <- readRDS("./pre-process/zero.pred2.RDS")[1]
too.high <- readRDS("./pre-process/high.cor.RDS") %>% str_replace(., "-", ".")

# Applicability domain ----------------------------------------------------

# Finds standard deviation for a single descriptor
# Requires a vector or single column; returns num
find.sd.desc <- function(data) {
  sd <- (data - mean(data, na.rm = T)) ^ 2 %>% sum()
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
standardize.withSDs <- function(data, sd.list) {
  df <- data[ , -1]
  guest <- data[ , 1]
  for (c in 1:length(df)) {
    sd <- sd.list[c]
    for (r in 1:nrow(df)) {
      ski <- abs(df[r, c] - mean(df[ , c])) / sd 
      df[r, c] <- ski
    }
    if(c %% 50 == 0) message(paste("Column ", c, " completed."))
  }
  return(cbind(guest, df))
}

# Precondition: data is the result of standardize
# standard deviation has been centered to be 1

domain.num <- function(data) {
  newSk <- c(rep(0, nrow(data)))
  if (class(data[, 1]) != "numeric") {
    guest <- data[, 1]
    data <- data[, -1]
    result <-
      apply(data, 1, function(x)
        mean(as.numeric(x), na.rm = T) + 1.28 * find.sd.desc(as.numeric(x))) %>%
      as.data.frame()
    result <- cbind(guest, result)
  } else
    result <- apply(data, 1, function(x)
      mean(as.numeric(x), na.rm = T) + 1.28 * find.sd.desc(as.numeric(x))) %>%
      as.data.frame()
  colnames(result)[1] <- "newSk"
  return(result %>% 
           mutate(domain = ifelse(result$newSk > 3, "outside", "inside")))
}

# Removes descriptors w/ very little variation (<= 2 unique values) in an 
# attempt to make applicability domain a little more useful
remove.binary <- function(data) {
  binary <- sapply(data, unique) %>% sapply(., length) %>% data.frame()
  binary.pred <- which(binary < 3)
  return(data[ , -binary.pred])
}

# SVM ---------------------------------------------------------------------

svm.a <- readRDS("./models/svm/polysvm.alpha.RDS")
svm.b <- readRDS("./models/svm/polysvm.beta.RDS")
svm.c <- readRDS("./models/svm/polysvm.gamma.RDS")

# GLMNet ------------------------------------------------------------------

glm.a <- readRDS("./models/glmnet/glm.alpha.RDS")
glm.b <- readRDS("./models/glmnet/glm.beta.RDS")
glm.c <- readRDS("./models/glmnet/glm.gamma.RDS")

