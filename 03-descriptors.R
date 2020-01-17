# Run this script after the SDFs have been run through PaDEL descriptor
# The .csv files should be saved under descriptors/ 

# install.packages("rcdk")
# install.packages("rJava")
# source("http://bioconductor.org/biocLite.R") # Sources the biocLite.R installation script. 
# biocLite("ChemmineR")
# library(rcdk)
if (!require(pacman)) install.packages("pacman")
require(pacman)
p_load(stringr, tidyverse)

# PaDEL-Descriptor --------------------------------------------------------

# Current settings:
#     1D, 2D
#     Fingerprints (PubChem)
#     Remove Salt
#     Detect aromaticity
#     Standardize nitro groups
#     Standardize tautomer
#     Convert to 3D, MM2, retain 3D coordinates
#     Use filename 
dg_data <- readRDS("cleaning/ri-suzuki.RDS")
dg_gamma <- readRDS("cleaning/singh.RDS")

# Total: 245/291 = 0.84192
alpha_dg <- inner_join(
  dg_data %>%
    filter(host == "alpha") %>%
    select(guest, dG),
  rename(
    read_csv("descriptors/alpha.csv"), 
    guest = Name),
  by = "guest")
# nrow(alpha_dg)
# nrow(filter(dg_data, host == "alpha"))

# Total: 372/407 = 0.91400 
beta_dg <- inner_join(
  dg_data %>%
    filter(host == "beta") %>%
    select(guest, dG),
  rename(
    read_csv("descriptors/beta.csv"), 
    guest = Name),
  by = "guest")

# 38/38 = 1
gamma_dg <- inner_join(
  dg_gamma %>%
    filter(host == "gamma") %>%
    select(guest, dG),
  rename(
    read_csv("descriptors/gamma.csv"), 
    guest = Name),
  by = "guest")

#     Saving Files --------------------------------------------------------

saveRDS(alpha_dg, "./descriptors/alpha-dg.RDS")
saveRDS(beta_dg, "./descriptors/beta-dg.RDS")
saveRDS(gamma_dg, "./descriptors/gamma-dg.RDS")
