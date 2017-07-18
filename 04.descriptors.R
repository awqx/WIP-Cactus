install.packages("rcdk")
library(rcdk)
library(tidyverse)

# PaDEL-Descriptor --------------------------------------------------------

# Current settings:
#     1D, 2D, and 3D descriptors
#     Fingerprints (PubChem)
#     Remove Salt
#     Detect aromaticity
#     Standardize tautomers
#     Convert to 3D, MMFF94, retain 3D coordinates
#     Use filename 

dataset <- readRDS("./bound/combined data.RDS")

#     Alpha-CD ------------------------------------------------------------

# SDFs that fail to process:
#     Phenol?, p-cresol
alpha.dg <- dataset %>% filter(host == "alpha")
alpha.padel.raw <- read_csv("./molecules/descriptors/2017-07-14 alpha predictors.csv") %>%
  rename(guest = Name)
alpha.padel <- inner_join(alpha.dg, alpha.padel.raw, by = "guest")
# Total: 180 guests

#     Beta-CD -------------------------------------------------------------

# SDFs that fail to process:
#     Barbital, p-cresol, thianapthene, 4-hydroxyacetophenone
beta.dg <- dataset %>% filter(host == "beta")
beta.padel.raw <- read_csv("./molecules/descriptors/2017-07-14 beta predictors.csv") %>%
  rename(guest = Name)
beta.padel <- inner_join(beta.dg, beta.padel.raw, by = "guest")
# Total: 15 guest molecules

#     Gamma-CD ------------------------------------------------------------

# No SDFs failed to process
gamma.dg <- dataset %>% filter(host == "gamma")
gamma.padel.raw <- read_csv("./molecules/descriptors/2017-07-14 gamma predictors.csv") %>%
  rename(guest = Name)
gamma.padel <- inner_join(gamma.dg, gamma.padel.raw, by = "guest")
# Total: 20 guest molecules

#     Saving Files --------------------------------------------------------

all.padel <- rbind(alpha.padel, beta.padel, gamma.padel)
saveRDS(all.padel, "./molecules/descriptors/04.all.padel.RDS")
saveRDS(alpha.padel, "./molecules/descriptors/04.alpha.padel.RDS")
saveRDS(beta.padel, "./molecules/descriptors/04.beta.padel.RDS")
saveRDS(gamma.padel, "./molecules/descriptors/04.gamma.padel.RDS")

write.csv(all.padel, "./molecules/descriptors/04.all.padel.csv")
write.csv(alpha.padel, "./molecules/descriptors/04.alpha.padel.csv")
write.csv(beta.padel, "./molecules/descriptors/04.beta.padel.csv")
write.csv(gamma.padel, "./molecules/descriptors/04.gamma.padel.csv")

# Rcdk Descriptors --------------------------------------------------------

mol <- load.molecules(c("C:/Users/Wei Xin/Documents/SREP LAB/Rekharsky and Inoue/Cactus/AlphaCD/2-butanol.SDF", 
                        "C:/Users/Wei Xin/Documents/SREP LAB/Rekharsky and Inoue/Cactus/AlphaCD/benzene.SDF"))
view.molecule.2d(mol)
mol2 <- get.murcko.fragments(mol) # works
get.bonds(mol)
get.bonds(mol2)
get.exact.mass(mol)
get.exact.mass(mol2)
eval.atomic.desc(mol)
eval.atomic.desc(mol2)
do.aromaticity(mol)
do.aromaticity(mol2)
is.aromatic(mol)
is.aromatic(mol2)

mol3 <- parse.smiles("C1=CC=C2C(=C1)C(=C(C(=O)O2)CC3=C(C4=CC=CC=C4OC3=O)O)O")
get.fingerprint(mol, type = "standard")
get.exact.mass(mol3)

is.aromatic(mol3)
do.aromaticity(mol3)

desc.names <- unique(unlist(sapply(get.desc.categories(), get.desc.names)))
temp <- eval.desc(mol, desc.names)
fp <- get.fingerprint(mol3[[1]], type = "maccs")
# ----
