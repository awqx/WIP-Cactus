install.packages("rcdk")
library(rcdk)
library(tidyverse)

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

dataset <- readRDS("./bound/combined data.RDS")

#     Alpha-CD ------------------------------------------------------------

# SDFs that fail to process:
#     Phenol?, p-cresol
alpha.dg <- dataset %>% filter(host == "alpha")
alpha.padel.raw <- read_csv("./molecules/descriptors/2017-07-19 alpha predictors.csv") %>%
  rename(guest = Name)
alpha.padel <- inner_join(alpha.dg, alpha.padel.raw, by = "guest")

# alpha.padel.2d <- read_csv("./molecules/descriptors/2017-07-19 alpha pred 2d.csv") %>%
#   rename(guest = Name)
# alpha.padel.2d <- inner_join(alpha.dg, alpha.padel.2d)
# Total: 209 guests

#     Beta-CD -------------------------------------------------------------

# SDFs that fail to process:
#     Barbital, p-cresol, thianapthene, 4-hydroxyacetophenone
beta.dg <- dataset %>% filter(host == "beta")
beta.padel.raw <- read_csv("./molecules/descriptors/2017-07-19 beta predictors.csv") %>%
  rename(guest = Name)
beta.padel <- inner_join(beta.dg, beta.padel.raw, by = "guest")

# beta.padel.2d <- read_csv("./molecules/descriptors/2017-07-19 beta pred 2d.csv") %>%
#   rename(guest = Name)
# beta.padel.2d <- inner_join(beta.dg, beta.padel.2d)
# Total: 321 guest molecules

#     Gamma-CD ------------------------------------------------------------

# No SDFs failed to process
gamma.dg <- dataset %>% filter(host == "gamma")
gamma.padel.raw <- read_csv("./molecules/descriptors/2017-07-19 gamma predictors.csv") %>%
  rename(guest = Name)
gamma.padel <- inner_join(gamma.dg, gamma.padel.raw, by = "guest")

# gamma.padel.2d <- read_csv("./molecules/descriptors/2017-07-19 gamma pred 2d.csv") %>%
#   rename(guest = Name)
# gamma.padel.2d <- inner_join(gamma.dg, gamma.padel.2d)
# Total: 17 guest molecules

#     Saving Files --------------------------------------------------------

all.padel <- rbind(alpha.padel, beta.padel, gamma.padel)
saveRDS(all.padel, "./molecules/descriptors/all.padel.RDS")
saveRDS(alpha.padel, "./molecules/descriptors/alpha.padel.RDS")
saveRDS(beta.padel, "./molecules/descriptors/beta.padel.RDS")
saveRDS(gamma.padel, "./molecules/descriptors/gamma.padel.RDS")

write.csv(all.padel, "./molecules/descriptors/all.padel.csv")
write.csv(alpha.padel, "./molecules/descriptors/alpha.padel.csv")
write.csv(beta.padel, "./molecules/descriptors/beta.padel.csv")
write.csv(gamma.padel, "./molecules/descriptors/gamma.padel.csv")

# all.padel.2d <- rbind(alpha.padel.2d, beta.padel.2d, gamma.padel.2d)
# saveRDS(all.padel.2d, "./molecules/descriptors/04.all.2d.RDS")
# saveRDS(alpha.padel.2d, "./molecules/descriptors/04.alpha.2d.RDS")
# saveRDS(beta.padel.2d, "./molecules/descriptors/04.beta.2d.RDS")
# saveRDS(gamma.padel.2d, "./molecules/descriptors/04.gamma.2d.RDS")

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
