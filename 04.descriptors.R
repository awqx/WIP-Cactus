# Run this script after the SDFs have been run through PaDEL descriptor
# The .csv files should be saved under descriptors/ 

# install.packages("rcdk")
# install.packages("rJava")
# source("http://bioconductor.org/biocLite.R") # Sources the biocLite.R installation script. 
# biocLite("ChemmineR")
# library(rcdk)
library(stringr)
library(tidyverse)

# Functions ---------------------------------------------------------------
# 
# # Write into a single .SDF for easy processing. Uses the filename as the
# # molecule name
# # mol.dir   a directory containing the SDF files that need to becompiled
# compile.sdf <- function(mol.dir) {
#   mol.files <- list.files(mol.dir, full.names = T)
#   mol.names <- list.files(mol.dir) %>%
#     str_remove(".SDF")
#   
#   # Removing empty files
#   mol.info <- file.info(mol.files) 
#   mol.files <- mol.files[!(mol.info$size == 0)]
#   
#   # Reading into a list
#   sdf.list <- lapply(mol.files, readRDS)
#   # Assigning the correct names
#   for(i in 1:length(sdf.list))
#     sdf.list[[i]][1] <- mol.names[i]
#   # Compiling into a single large data.frame
#   sdf <- do.call(rbind, sdf.list)
#   
#   return(sdf)
# }
# 
# Optimize energies -------------------------------------------------------
# 
# mol.files <- list.files("./molecules/alphaCD", full.names = T)
# mol.info <- file.info(mol.files)
# empties <- mol.info$size == 0
# mol.files <- mol.files[!empties]
# mol.names <- mol.files %>% str_remove("./molecules/alphaCD/") %>%
#   str_remove(".SDF")
# compiled.sdf <- lapply(mol.files, read.csv, 
#                        header = F, stringsAsFactors = F)
# for(n in 1:length(compiled.sdf)) {
#   temp <- compiled.sdf[[n]]
#   temp[1, 1] <- mol.names[[n]]
#   compiled.sdf[[n]] <- temp
# }
# compiled.sdf <- do.call(rbind, compiled.sdf)
# write.table(compiled.sdf, "./molecules/alpha.SDF", row.names = F, quote = F, 
#             col.names = F)
# # Run openBabel settings: 
# # Explicit hydrogens
# # continue past errors
# # add 3d coordinates
# 
# mol.files <- list.files("./molecules/betaCD", full.names = T)
# mol.info <- file.info(mol.files)
# empties <- mol.info$size == 0
# mol.files <- mol.files[!empties]
# mol.names <- mol.files %>% str_remove("./molecules/betaCD/") %>%
#   str_remove(".SDF")
# compiled.sdf <- lapply(mol.files, read.csv, 
#                        header = F, stringsAsFactors = F)
# for(n in 1:length(compiled.sdf)) {
#   temp <- compiled.sdf[[n]]
#   temp[1, 1] <- mol.names[[n]]
#   compiled.sdf[[n]] <- temp
# }
# compiled.sdf <- do.call(rbind, compiled.sdf)
# write.table(compiled.sdf, "./molecules/beta.SDF", row.names = F, quote = F, 
#             col.names = F)

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
dataset <- readRDS("cleaning/03.rename.RDS")

#     Alpha-CD ------------------------------------------------------------

# SDFs that fail to process:
#     Phenol?, p-cresol
alpha.dg <- dataset %>% filter(host == "alpha")
alpha.padel.raw <- read_csv("./descriptors/alpha.csv") %>%
  rename(guest = Name)
alpha.padel <- inner_join(alpha.dg, alpha.padel.raw, by = "guest")

# Total: 223/241 = 92.5%

#     Beta-CD -------------------------------------------------------------

# SDFs that fail to process:
#     Barbital, p-cresol, thianapthene, 4-hydroxyacetophenone
beta.dg <- dataset %>% filter(host == "beta")
beta.padel.raw <- read_csv("./descriptors/beta.csv") %>%
  rename(guest = Name)
beta.padel <- inner_join(beta.dg, beta.padel.raw, by = "guest")

# Total: 341/354 = 96.3% yield

#     Gamma-CD ------------------------------------------------------------

# No SDFs failed to process
gamma.dg <- dataset %>% filter(host == "gamma")
gamma.padel.raw <- read_csv("./descriptors/gamma.csv") %>%
  rename(guest = Name)
gamma.padel <- inner_join(gamma.dg, gamma.padel.raw, by = "guest")

# Total: 124/130 = 95.4% yield

# Right now, there isn't enough data for gamma-CD to create a reliable model
# so the descriptors won't be analyzed...yet

# # Suzuki Only -------------------------------------------------------------
# 
# suz <- readRDS("./dwnld/suzuki.only.RDS")
# suz.a <- suz %>% filter(host == "alpha")
# suz.b <- suz %>% filter(host == "beta")
# 
# suz.a.padel <- inner_join(suz.a, alpha.padel.raw)
# suz.b.padel <- inner_join(suz.b, beta.padel.raw)
# suz.padel <- rbind(suz.a.padel, suz.b.padel)

#     Saving Files --------------------------------------------------------

# all.padel <- rbind(alpha.padel, beta.padel, gamma.padel)
# saveRDS(all.padel, "./descriptors/all.padel.RDS")
# write.csv(all.padel, "./descriptors/all.padel.csv")
saveRDS(alpha.padel, "./descriptors/alpha.padel.RDS")
saveRDS(beta.padel, "./descriptors/beta.padel.RDS")
saveRDS(gamma.padel, "./descriptors/gamma.padel.RDS")

# saveRDS(suz.padel, "./descriptors/suz.padel.RDS")
# write.csv(suz.padel, "./descriptors/suz.padel.csv")

# Rcdk Descriptors --------------------------------------------------------

# mol <- load.molecules(c("./molecules/alphaCD/toluene.SDF", 
#                         "./molecules/alphaCD/quinoline.SDF"))
# view.molecule.2d(mol)
# mol2 <- get.murcko.fragments(mol) # works
# get.bonds(mol)
# get.bonds(mol2)
# get.exact.mass(mol)
# get.exact.mass(mol2)
# eval.atomic.desc(mol)
# eval.atomic.desc(mol2)
# do.aromaticity(mol)
# do.aromaticity(mol2)
# is.aromatic(mol)
# is.aromatic(mol2)
# 
# mol3 <- parse.smiles("C1=CC=C2C(=C1)C(=C(C(=O)O2)CC3=C(C4=CC=CC=C4OC3=O)O)O")
# get.fingerprint(mol, type = "standard")
# get.exact.mass(mol3)
# 
# is.aromatic(mol3)
# do.aromaticity(mol3)
# 
# desc.names <- unique(unlist(sapply(get.desc.categories(), get.desc.names)))
# temp <- eval.desc(mol, desc.names)
# fp <- get.fingerprint(mol3[[1]], type = "maccs")
# # ----
