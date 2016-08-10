install.packages("stringr")
install.packages("dplyr")
install.packages("ggplot2")
library(stringr)
library(dplyr)
library(ggplot2)

source(file = paste0(folder, "cactus.functions.R"))
import.aff <- function(filename, path, host){
  file <- paste0(path, filename)
  raw.aff <- read.csv(file, header = T)
  chem.num <- length(raw.aff$Ligand)
  run.num <- str_extract(filename, "[[:digit:]]+")
  run.num <- rep(run.num, chem.num)
  host <- rep(host, chem.num)
  raw.aff$runs <- run.num
  raw.aff$Ligand <-
    str_replace(raw.aff$Ligand, paste0(host, "-cd_(uff|mmff94|ghemical)_E=[[:digit:]]+\\.[[:digit:]]+_"), "")
  raw.aff$Ligand <-
    str_replace(raw.aff$Ligand,
                "_(uff|mmff94|ghemical)_E=[[:digit:]]+\\.[[:digit:]]+",
                "")
  colnames(raw.aff)[colnames(raw.aff) == "Ligand"] <-
    "formula"
  return(raw.aff)
}
multi_join <- function(..., by){
  
  jointhis <- list(...) 
  joint_df <- Reduce(function(x,y) merge(x,y, by = by),jointhis)
  
  return(joint_df)}
# Put all the PyRx results, as .csv files, here
# Feels like this could definitely be made into a cohesive function
affinity.path <- create.host.dir(folder, "Affinity/")
sample.aff.path <- paste0(affinity.path, "Sample/")
# Function that reads .csv from PyRx
#Table containing numbers corresponding to PyRx variables
runs <- c(1:72)
forcefields <- c("uff", "mmff94", "ghemical")
min.alg <- c("conj", "steep")
vina <- c("default", "maximized")
vina <- rep(vina, 36)
host.ff <- rep(forcefields, each = 24)
host.alg <- rep(rep(min.alg, each = 12), 3)
ligand.ff <- rep(rep(forcefields, each = 4), 6)
ligand.alg <- rep(rep(min.alg, each = 2), 18)
runs.df <- data.frame(runs, host.ff, host.alg, ligand.ff, ligand.alg, vina)

sample.aff.files <- list.files(sample.aff.path, pattern = ".csv")
sample.aff <-
  do.call(rbind,
          lapply(sample.aff.files, import.aff, sample.aff.path, "alpha"))
sample.aff.all <- multi_join(sample.aff, runs.df, by = "runs")

plot1 <- ggplot(sample.aff.all, aes(x = Binding.Affinity, fill = vina)) +
  geom_histogram(
    data = subset(sample.aff.all, vina == "default"),
    binwidth = 0.2,
    alpha = 0.3, 
    center = -4.1
  ) +
  geom_histogram(
    data = subset(sample.aff.all, vina == "maximized"),
    binwidth = 0.2,
    alpha = 0.3, 
    center = -4.1
  ) +
  facet_grid(ligand.alg ~ ligand.ff) +
  scale_fill_manual(values=c("green", "blue")) +
  theme_bw()

plot2 <- ggplot(sample.aff.all, aes(x = Binding.Affinity, fill = ligand.ff)) +
  geom_histogram(
    data = subset(sample.aff.all, ligand.ff == "uff"), 
    binwidth = 0.2,
    alpha = 0.2, 
    center = -4.1
  ) + 
  geom_histogram(
    data = subset(sample.aff.all, ligand.ff == "mmff94"), 
    binwidth = 0.2,
    alpha = 0.2, 
    center = -4.1
  ) + 
  geom_histogram(
    data = subset(sample.aff.all, ligand.ff == "ghemical"), 
    binwidth = 0.2,
    alpha = 0.2, 
    center = -4.1
  ) + 
  facet_grid(ligand.alg ~ vina) +
  scale_fill_manual(values=c("green", "blue", "red")) +
  theme_bw()

plot3 <- ggplot(sample.aff.all, aes(x = Binding.Affinity, fill = ligand.alg)) +
  geom_histogram(
    data = subset(sample.aff.all, ligand.alg == "conj"), 
    binwidth = 0.2, 
    alpha = 0.2, 
    center = -4.1
  ) + 
  geom_histogram(
    data = subset(sample.aff.all, ligand.alg == "steep"), 
    binwidth = 0.2, 
    alpha = 0.2, 
    center = -4.1
  ) + 
  facet_grid(ligand.ff ~ vina) +
  scale_fill_manual(values=c("green", "blue")) +
  theme_bw()

min.aff.sample <- sample.aff.all %>%
  group_by(formula, ligand.ff, ligand.alg, vina) %>%
  summarize(min(Binding.Affinity))
colnames(min.aff.sample)[colnames(min.aff.sample) == "min(Binding.Affinity)"] <-
  "Binding.Affinity"

plot4 <- ggplot(min.aff.sample, aes(x = ligand.ff, y = Binding.Affinity, fill = formula)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(.~ligand.alg) +
  theme_bw()
