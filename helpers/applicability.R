# Function based off of Roy 2015: Determining Applicability Domain of QSAR Models
# Creating a separate file that allows for convenient sourcing

# Libraries and Packages --------------------------------------------------

if (!require(pacman)) install.packages("pacman")
require(pacman)
p_load(caret, stringr, tidyverse)

# Functions ---------------------------------------------------------------

# Finds standard deviation for a single descriptor
# Requires a vector or single column; returns num
find_sd <- function(desc) {
  desc <- desc[!is.na(desc)]
  sd <- sum((desc - mean(desc, na.rm = T)) ^ 2)
  sqrt(sd / length(desc)) 
}

# Determines whether a chemical is within the applicability domain
# Requires a df or matrix and an index
# Can be used in do.call, rbind, lapply sequence
determine_domain <- function(index, data){
  results <- c(rep(0.0, length(data)  - 1))
  for (i in 2:length(data)) {
    sd <- find_sd(data[ , i])
    results[i - 1] <- (data[index, i] - mean(data[ , i])) / sd
  }
}

# Precondition: first column of data is guests, rest is descriptors
# initial.standardize works on the data the model was trained on
initial_standardize <- function(data) {
  df <- data
  # check if first column is guests
  if (class(data[ , 1]) == "character") {
    guest <- df[ , 1]
    df <- df[ , -1]
  }
  for (c in 1:length(df)) {
    sd <- find_sd(df[ , c])
    for (r in 1:nrow(df)) {
     ski <- abs(df[r, c] - mean(df[ , c])) / sd 
     df[r, c] <- ski
    }
    # message(paste("Column ", c, " completed."))
  }
  if (exists("guest")) cbind(guest, df)
  df
}

# Standardize works on new data with a list of standard deviations
# sd.list should be retrieved from the data the model was trained on
use_sds <- function(data, sd_list) {
  df <- data
  # check if first column is guests
  if (class(data[ , 1]) == "character") {
    guest <- df[ , 1]
    df <- df[ , -1]
  }
  for (c in 1:length(df)) {
    sd <- sd_list[c]
    for (r in 1:nrow(df)) {
      ski <- abs(df[r, c] - mean(df[ , c])) / sd 
      df[r, c] <- ski
    }
    if(c %% 50 == 0) message(paste("Column ", c, " completed."))
  }
  if (exists("guest")) cbind(guest, df)
  df
}

# Precondition: data is the result of standardize
# standard deviation has been centered to be 1
domain_num <- function(data) {
  newSk <- c(rep(0, nrow(data)))
  # Checking if first column is "guest"
  if (class(data[, 1]) != "numeric") {
    guest <- data[, 1]
    data <- sapply(data[, -1], abs)
    if(nrow(data) == 1) {
      result <- mean(
        as.numeric(data[1, ], na.rm = T) + 1.28 * find_sd(as.numeric(data[1, ])))
      result <- data.frame(result)
    } else {
      result <-
        apply(data, 1, function(x)
          mean(as.numeric(x), na.rm = T) + 1.28 * find_sd(as.numeric(x))) %>%
        as.data.frame()
    }
    result <- data.frame(guest, result) %>%
      mutate(guest = as.character(guest))
    colnames(result)[2] <- "newSk"
  } else {
    data <- sapply(data[ , -1], abs)
    result <- apply(data, 1, function(x)
      mean(as.numeric(x), na.rm = T) + 1.28 * find_sd(as.numeric(x))) %>%
      as.data.frame()
    colnames(result)[1] <- "newSk"
  }
  max_ski <- apply(data, 1, max, na.rm = T)
  min_ski <- apply(data, 1, min, na.rm = T)
  result <- cbind(result, max_ski, min_ski)
  result %>%
    mutate(
      domain = ifelse(result$max_ski < 3, "inside", 
                      ifelse(result$min_ski > 3, "outside", 
                             ifelse(result$newSk > 3, "outside", "inside")))
    )
}

# Removes descriptors w/ very little variation (<= 2 unique values) in an 
# attempt to make applicability domain a little more useful
remove_binary <- function(data) {
  binary <- sapply(data, unique) %>% 
    sapply(., length) %>% 
    data.frame()
  binary_pred <- which(binary < 3)
  data[ , -binary.pred]
}