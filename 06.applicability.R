# Function based off of Roy 2015: Determining Applicability Domain of QSAR Models
# setwd("~/SREP LAB/qsar")
dir.create("./domain")

# Libraries and Packages --------------------------------------------------

library(caret)
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

# Testing on Some Data ----------------------------------------------------

df <- readRDS("./padel.pp.new.RDS") %>% select(., -DelG, -host, -data.source)
test <- df[ , 1:25]
sd.all <- lapply(test[ , -1], find.sd.desc) %>% unlist() %>% as.vector() # Works
stand.test <- initial.standardize(test)
stand.test2 <- standardize(test, sd.all)

# Analysis ----------------------------------------------------------------

# Test to make sure it works
df <- readRDS("./padel.pp.new.RDS") %>% select(., -DelG, -host, -data.source)
desc <- df[ , -1]
sdevs <- lapply(desc, find.sd.desc) %>% unlist() %>% as.vector() 
saveRDS(sdevs, "./domain/standard deviations.RDS")
# standardize(df, sdevs)
stand <- initial.standardize(df)
app.domain <- domain(stand)
saveRDS(app.domain, "applicability domain.RDS")

# Training Set (seed 25)
set.seed(25)
trn.ind <- sample(x = 1:nrow(df), size = round(0.7 * nrow(df)))
trn <- df[trn.ind, ]
tst <- df[-trn.ind, ]
trn.desc <- trn[ , -1]
trn.sdevs <- lapply(trn.desc, find.sd.desc) %>% unlist() %>% as.vector()
saveRDS(trn.sdevs, "./domain/trn.sdevs.RDS")
trn.stand <- standardize(trn, trn.sdevs)
trn.dmn <- domain(trn.stand) # All inside

# Test Set
tst.stand <- standardize(tst, trn.sdevs)
tst.dmn <- domain(tst.stand) # All inside

# Testing to see if "outside" even exists
samp <- c(rep(7, 712))
samp <- rbind(tst, samp)
samp.stand <- standardize(samp, trn.sdevs)
samp.dmn <- domain(samp.stand)
# OK, confirmed that "outside" the applicability domain can exist, which is reassuring