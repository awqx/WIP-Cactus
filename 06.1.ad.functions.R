# Function based off of Roy 2015: Determining Applicability Domain of QSAR Models
# Creating a separate file that allows for convenient sourcing

# Libraries and Packages --------------------------------------------------

library(caret)
library(stringr)
library(tidyverse)

# Functions ---------------------------------------------------------------

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