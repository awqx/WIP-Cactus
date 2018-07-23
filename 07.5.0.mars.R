# MARS - Multivariate Area Regression 
dir.create("./models/mars")

# Libraries ---------------------------------------------------------------

if(!require("pacman"))
  install.packages("pacman")
library(pacman)
p_load(caret, tidyverse, earth)
source("./07.model.functions.R")

# Functions ---------------------------------------------------------------


# Test --------------------------------------------------------------------


