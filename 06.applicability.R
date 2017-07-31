# Function based off of Roy 2015: Determining Applicability Domain of QSAR Models

# Libraries and Packages --------------------------------------------------

library(caret)
library(stringr)
library(tidyverse)



# Functions ---------------------------------------------------------------

# Finds standard deviation for a single descriptor
# Requires a vector or single column; returns num
find.sd.desc <- function(data) {
  return(sum((data - mean(data)) ^ 2) / (length(data - 1)) %>%
           sqrt())
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

# Analysis ----------------------------------------------------------------


