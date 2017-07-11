install.packages("rcdk")
library(rcdk)
library(tidyverse)

# ------
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
