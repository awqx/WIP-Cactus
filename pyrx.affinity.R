install.packages("stringr")
install.packages("dplyr")
library(stringr)
library(dplyr)

source(file = paste0(folder, "cactus.functions.R"))
# Put all the PyRx results, as .csv files, here
# Feels like this could definitely be made into a cohesive function
affinity.path <- create.host.dir(folder, "Affinity")

gamma.aff.path <- paste0(affinity.path, "/Gamma")

aff.gamma.files <-
  list.files(path = gamma.aff.path, pattern = ".csv")

raw.gamma.aff <- read.csv(aff.gamma.files,
                          header = T)

gamma.aff <-
  do.call(rbind, lapply(
    paste0(gamma.aff.path, "/", aff.gamma.files),
    read.csv,
    header = T
  ))

sample.groups <- as.character(unique(raw.gamma.aff$Ligand))

min.aff.gamma <- raw.gamma.aff %>%
  group_by(Ligand) %>%
  summarize(min(Binding.Affinity))

min.aff.gamma$Ligand <-
  str_replace(min.aff.gamma$Ligand, "gamma-cd_", "")

min.aff.gamma$Ligand <-
  str_replace(min.aff.gamma$Ligand,
              "_uff_E=[[:digit:]]+\\.[[:digit:]]+",
              "")

colnames(min.aff.gamma)[colnames(min.aff.gamma) == "Ligand"] <-
  "formula"

min.gamma <- min.aff.gamma

test <- full_join(min.gamma, gamma.names.table)
