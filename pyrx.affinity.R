install.packages("stringr")
install.packages("dplyr")
install.packages("ggplot2")
library(stringr)
library(dplyr)
library(ggplot2)

source(file = paste0(folder, "cactus.functions.R"))
# Function that reads .csv from PyRx
import.aff <- function(filename, path, host){
  file <- paste0(path, filename)
  raw.aff <- read.csv(file, header = T, row.names = NULL)
  chem.num <- length(raw.aff$Ligand)
  run.num <- str_extract(filename, "[[:digit:]]+")
  run.num <- rep(run.num, chem.num)
  host <- rep(host, chem.num)
  raw.aff$runs <- run.num
  raw.aff$Ligand <-
    str_replace(raw.aff$Ligand, paste0(host, "-cd_(uff|mmff94|ghemical)_E=[[:digit:]]+\\.[[:digit:]]+_"), "")
  raw.aff$Ligand <-
    str_replace(raw.aff$Ligand,
                "_(uff|mmff94|ghemical)_E=//-?[[:digit:]]+\\.[[:digit:]]+",
                "")
  colnames(raw.aff)[colnames(raw.aff) == "Ligand"] <-
    "formula"
  return(raw.aff)
}
find.atoms <- function(chem.name, atom, symbol, atom.wt){
  pattern <- paste0(symbol, "[[:digit:]]*")
  atom <- str_extract(chem.name, pattern)
  if(is.na(atom)){
    atom <- 0
  }
  if(str_detect(atom, "[[:digit:]]")){
    atom.num <- str_extract(atom, "[[:digit:]]+") %>% as.integer
  } else if(atom == symbol){
    atom.num <- 1
  } else {
    atom.num <- 0
  }
  comb.wt <- atom.wt*atom.num
  return(comb.wt)
}
elements <- c("hydrogen", "carbon", "nitrogen", "oxygen", "fluorine", "sodium", "phosphorous", "sulfur", "chlorine", "bromine", "iodine", "platinum")
symbols <- c("H", "C", "N", "O", "F", "Na", "P", "S", "Cl", "Br", "I", "Pt")
atom.wts <- c(1.00794, 12.011, 14.0067, 15.9994, 18.9984, 22.98977, 30.97376, 32.066, 35.4527, 79.904, 126.9045, 195.08)

lapply(min.aff.sample$formula, find.atoms, "hydrogen", "H", 1.00794)
hydrogen <- min.aff.sample$formula %>%
  lapply(find.atoms, "hydrogen", "H", 1.00794) %>%
  unlist

carbon <- min.aff.sample$formula %>%
  lapply(find.atoms, "carbon", "C", 12.011) %>%
  unlist

mol.wt <- Map(sum, hydrogen, carbon) %>% unlist
mol.wt.sample <- data.frame(min.aff.sample$formula, hydrogen, carbon, mol.wt)

# Folder for results of chemical ID
chem.path <- create.host.dir(folder, "ChemIDs")

# Put all the PyRx results, as .csv files, in the following folders
affinity.path <- create.host.dir(folder, "Affinity/")
alpha.aff.path <- create.host.dir(affinity.path, "Alpha/")
beta.aff.path <- create.host.dir(affinity.path, "Beta/")
gamma.aff.path <- create.host.dir(affinity.path, "Gamma/")
sample.aff.path <- create.host.dir(affinity.path, "Sample/")
#Table containing numbers corresponding to PyRx variables
runs <- as.character(c(1:72))
forcefields <- c("uff", "mmff94", "ghemical")
min.alg <- c("conj", "steep")
vina <- c("default", "maximized")
vina <- rep(vina, 36)
host.ff <- rep(forcefields, each = 24)
host.alg <- rep(rep(min.alg, each = 12), 3)
ligand.ff <- rep(rep(forcefields, each = 4), 6)
ligand.alg <- rep(rep(min.alg, each = 2), 18)
runs.df <- data.frame(runs, host.ff, host.alg, ligand.ff, ligand.alg, vina)
#----------------------------------------------------------
#                       Gamma
#----------------------------------------------------------
gamma.aff.files <- list.files(gamma.aff.path, pattern = ".csv")
gamma.aff <-
  do.call(rbind,
          lapply(gamma.aff.files, import.aff, gamma.aff.path, "gamma"))
gamma.aff$Binding.Affinity <- as.numeric(gamma.aff$Binding.Affinity)
atom.num <- gamma.aff$formula %>%
  str_extract("^[[:alnum:]]+_") %>%
  str_extract_all("[[:digit:]]+") %>%
  lapply(as.numeric) %>%
  lapply(sum)
gamma.aff$atom.num <- atom.num
gamma.aff.all <- left_join(gamma.aff, runs.df)
# Plotting 
graph.path <-create.host.dir(folder, "Graphs")
gamma.graphs <- create.host.dir(graph.path, "/Gamma/")
plot0gamma <- ggplot(gamma.aff.all, aes(x = Binding.Affinity)) +
  geom_histogram(binwidth = 0.2, center = -4.1) +
  theme_bw() +
  xlab("Binding Affinity, kcal/mol") +
  ylab("Number of Instances")
ggsave(paste0(gamma.graphs, "0.png"), plot = plot0gamma)

plot1gamma <- ggplot(gamma.aff.all, aes(x = Binding.Affinity, fill = vina)) +
  geom_histogram(
    data = subset(gamma.aff.all, vina == "default"),
    binwidth = 0.2,
    alpha = 0.3, 
    center = -4.1
  ) +
  geom_histogram(
    data = subset(gamma.aff.all, vina == "maximized"),
    binwidth = 0.2,
    alpha = 0.3, 
    center = -4.1
  ) +
  facet_grid(ligand.alg ~ ligand.ff) +
  scale_fill_manual(values=c("green", "blue")) +
  theme_bw() +
  xlab("Binding Affinity, kcal/mol") +
  ylab("Number of Instances")
ggsave(paste0(gamma.graphs, "1.png"), plot = plot1gamma)

plot2gamma<- ggplot(gamma.aff.all, aes(x = Binding.Affinity, fill = ligand.ff)) +
  geom_histogram(
    data = subset(gamma.aff.all, ligand.ff == "uff"), 
    binwidth = 0.2,
    alpha = 0.2, 
    center = -4.1
  ) + 
  geom_histogram(
    data = subset(gamma.aff.all, ligand.ff == "mmff94"), 
    binwidth = 0.2,
    alpha = 0.2, 
    center = -4.1
  ) + 
  geom_histogram(
    data = subset(gamma.aff.all, ligand.ff == "ghemical"), 
    binwidth = 0.2,
    alpha = 0.2, 
    center = -4.1
  ) + 
  facet_grid(ligand.alg ~ vina) +
  scale_fill_manual(values=c("green", "blue", "magenta")) +
  theme_bw() +
  xlab("Binding Affinity, kcal/mol") +
  ylab("Number of Instances")

ggsave(paste0(gamma.graphs, "2.png"), plot = plot2gamma)

plot3gamma <- ggplot(gamma.aff.all, aes(x = Binding.Affinity, fill = ligand.alg)) +
  geom_histogram(
    data = subset(gamma.aff.all, ligand.alg == "conj"), 
    binwidth = 0.2, 
    alpha = 0.2, 
    center = -4.1
  ) + 
  geom_histogram(
    data = subset(gamma.aff.all, ligand.alg == "steep"), 
    binwidth = 0.2, 
    alpha = 0.2, 
    center = -4.1
  ) + 
  facet_grid(ligand.ff ~ vina) +
  scale_fill_manual(values=c("green", "blue")) +
  theme_bw() +
  xlab("Binding Affinity, kcal/mol") +
  ylab("Number of Instances")
ggsave(paste0(gamma.graphs, "3.png"), plot = plot3gamma)

plot4gamma <- ggplot(gamma.aff.all, aes(x = as.numeric(atom.num), y = Binding.Affinity, color = ligand.alg)) +
  geom_point(alpha = 0.1) +
  facet_grid(ligand.ff~vina) +
  scale_color_manual(values=c("seagreen1", "steelblue1")) +
  geom_smooth(method = lm) +
  xlab("Number of Atoms") +
  ylab("Binding Affinity, kcal/mol") +
  theme_bw()
ggsave(paste0(gamma.graphs, "4.png"), plot = plot4gamma)

plot5gamma <- ggplot(gamma.aff.all, aes(x = as.numeric(atom.num), y = Binding.Affinity, color = vina)) +
  geom_point(alpha = 0.1) +
  facet_grid(ligand.ff~ligand.alg) +
  scale_color_manual(values=c("seagreen1", "steelblue1")) +
  geom_smooth(method = lm) +
  xlab("Number of Atoms") +
  ylab("Binding Affinity, kcal/mol") +
  theme_bw()
ggsave(paste0(gamma.graphs, "5.png"), plot = plot5gamma)

plot6gamma <- ggplot(gamma.aff.all, aes(x = as.numeric(atom.num), y = Binding.Affinity, color = ligand.ff)) +
  geom_point(alpha = 0.1) +
  facet_grid(vina~ligand.alg) +
  scale_color_manual(values=c("seagreen1", "steelblue1", "orchid")) +
  geom_smooth(method = lm) +
  xlab("Number of Atoms") +
  ylab("Binding Affinity, kcal/mol") +
  theme_bw()
ggsave(paste0(gamma.graphs, "6.png"), plot = plot6gamma)

min.aff.gamma <- gamma.aff.all %>%
  group_by(formula, ligand.ff, ligand.alg, vina, runs, atom.num) %>%
  summarize(min(Binding.Affinity))
colnames(min.aff.gamma)[colnames(min.aff.gamma) == "min(Binding.Affinity)"] <-
  "Binding.Affinity"
min.aff.gamma$atom.num <- as.integer(min.aff.gamma$atom.num)
plot7gamma <- ggplot(min.aff.gamma, aes(x = atom.num, y = Binding.Affinity, color = vina)) +
  geom_point(alpha = 0.2) +
  scale_x_continuous() +
  facet_grid(ligand.ff~ligand.alg) +
  scale_color_manual(values=c("green", "blue")) +
  theme_bw()
plot8gamma <- ggplot(min.aff.gamma, aes(x = atom.num, y = Binding.Affinity, color = ligand.ff)) +
  geom_point(alpha = 0.2) +
  facet_grid(vina~ligand.alg) +
  scale_color_manual(values=c("green", "blue", "magenta")) +
  theme_bw()
plot9gamma <- ggplot(min.aff.gamma, aes(x = ligand.ff, y = Binding.Affinity, color = ligand.alg, shape = ligand.alg)) +
  geom_point(alpha = 0.3) +
  facet_grid(.~vina) +
  scale_color_manual(values=c("green", "blue")) +
  theme_bw()






# Sample
sample.aff.files <- list.files(sample.aff.path, pattern = ".csv")
sample.aff <-
  do.call(rbind,
          lapply(sample.aff.files, import.aff, sample.aff.path, "alpha"))

# plotting, looking for patterns
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

plot4 <- ggplot(min.aff.sample, aes(x = ligand.ff, y = Binding.Affinity, color = formula)) +
  geom_point() +
  facet_grid(vina~ligand.alg) +
  theme_bw()

plot5 <- ggplot(min.aff.sample, aes(x = vina, y = Binding.Affinity, color = formula)) +
  geom_point() +
  facet_grid(ligand.ff~ligand.alg) +
  theme_bw()

plot6 <- ggplot(min.aff.sample, aes(x = formula, y = Binding.Affinity, color = vina)) +
  geom_point() +
  facet_grid(ligand.ff~ligand.alg) +
  theme_bw()

plot7 <- ggplot(min.aff.sample, aes(x = formula, y = Binding.Affinity, color = ligand.alg)) +
  geom_point() +
  facet_grid(ligand.ff~vina) +
  theme_bw()

gamma.thing <- import.aff("1gamma.csv", paste0(affinity.path, "Gamma/"), "gamma")

# Cleaning up the RI.rds dataset 
dataset.clean <- dataset[dataset$T.K == 298,]
dataset.clean.2 <- dataset[dataset$T.K == 303,]
dataset.clean <- rbind(dataset.clean, dataset.clean.2)
dataset.clean <- dataset.clean[!str_detect(dataset.clean$ref, "b|c|g|i|j|m"),]
alpha.guest.clean <- unique(dataset.clean$guest[dataset.clean$host == "1\u03b1"])
beta.guest.clean  <- unique(dataset.clean$guest[dataset.clean$host == "1\u03b2"])
gamma.guest.clean <- unique(dataset.clean$guest[dataset.clean$host == "1γ"])
hosts <- c("1\u03b1", "1\u03b2", "1γ")
all.guests.clean <- unique(c(alpha.guest.clean, beta.guest.clean, gamma.guest.clean))
ds.cleaner <- dataset.clean[dataset.clean$guest %in% all.guests.clean, ]
ds.cleaner <- ds.cleaner[ds.cleaner$host %in% hosts, ]
convert.kj.kcal <- function(kJ){
  kcal <- kJ/4.184
  return(kcal)
}
binding.aff <- unlist(Map(convert.kj.kcal, ds.cleaner$DelG))
ds.cleaner$binding.affinity <- binding.aff
# Separating various parts of datasets into something graphable
gamma.part <- min.aff.gamma
gamma.part$formula <- str_replace(gamma.part$formula, "^[[:alnum:]]+_", "")
gamma.part$formula <- str_replace(gamma.part$formula, "_(ghemical|mmff94)_E=\\-[[:digit:]]+.[[:digit:]]+", "")
gamma.part$formula <- str_replace(gamma.part$formula, "_", " ")
gamma.ds <- ds.cleaner[ds.cleaner$host == "1γ",c("guest", "binding.affinity")]
gamma.ds$formula <- gamma.ds$guest
gamma.ds$guest <- NULL
predict.gamma <- left_join(gamma.part, gamma.ds, by = "formula")
perc.yield <- (predict.gamma$Binding.Affinity / predict.gamma$binding.affinity) * 100
predict.gamma$perc.yield <- perc.yield
gamma.graphs.affinity <- create.host.dir(gamma.graphs, "Affinity/")
plot1affinity <- ggplot(predict.gamma, aes(y = Binding.Affinity, x = binding.affinity)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "gold", size = 1) +
  geom_smooth(method = lm, se = T, color = "seagreen1") +
  xlab("Rekharsky and Inoue Experimental Affinity") +
  ylab("Predicted (Vina) Affinity")
ggsave(paste0(gamma.graphs.affinity, "1.png"), plot1affinity)

plot2affinity <- ggplot(predict.gamma, aes(y = Binding.Affinity, x = binding.affinity, color = vina)) +
  geom_point(alpha = 0.2) +
  geom_abline(slope = 1, intercept = 0, color = "gold", size = 1) +
  geom_smooth(method = lm, se = T) +
  xlab("Rekharsky and Inoue Experimental Affinity") +
  ylab("Predicted (Vina) Affinity") +
  scale_color_manual(values=c("seagreen1", "slateblue2"))
ggsave(paste0(gamma.graphs.affinity, "2.png"), plot2affinity)

plot3affinity <- ggplot(predict.gamma, aes(y = Binding.Affinity, x = binding.affinity, color = ligand.ff)) +
  geom_point(alpha = 0.2) +
  geom_abline(slope = 1, intercept = 0, color = "gold", size = 1) +
  geom_smooth(method = lm, se = T) +
  xlab("Rekharsky and Inoue Experimental Affinity") +
  ylab("Predicted (Vina) Affinity") +
  scale_color_manual(values=c("seagreen1", "slateblue2", "orchid"))
ggsave(paste0(gamma.graphs.affinity, "3.png"), plot3affinity)

plot3.5affinity <- ggplot(predict.gamma, aes(y = Binding.Affinity, x = binding.affinity, color = ligand.alg)) +
  geom_point(alpha = 0.2) +
  geom_abline(slope = 1, intercept = 0, color = "gold", size = 1) +
  geom_smooth(method = lm, se = T) +
  xlab("Rekharsky and Inoue Experimental Affinity") +
  ylab("Predicted (Vina) Affinity") +
  scale_color_manual(values=c("seagreen1", "slateblue2"))
ggsave(paste0(gamma.graphs.affinity, "3.5.png"), plot3.5affinity)

plot4affinity <- ggplot(predict.gamma, aes(x = atom.num, y = perc.yield, color = ligand.ff)) +
  geom_point(alpha = 0.2) +
  geom_abline(slope = 0, intercept = 100, color = "gold", size = 1) +
  geom_smooth(method = lm) +
  xlab("Atom Number") +
  ylab("% Yield") +
  scale_color_manual(values=c("seagreen1", "slateblue2", "orchid"))
ggsave(paste0(gamma.graphs.affinity, "4.png"), plot4affinity)

ggplot(predict.gamma, aes(y = Binding.Affinity, x = binding.affinity, color = ligand.ff)) +
  geom_point(alpha = 0.4) +
  geom_abline(slope = 0, intercept = 100, color = "gold", size = 1) +
  facet_grid(ligand.alg~vina) +
  geom_smooth(method = lm, se = T) +
  xlab("Expected Affinity") +
  ylab("Predicted (Vina) Affinity") +
  scale_color_manual(values=c("seagreen1", "slateblue2", "orchid"))

plot6affinity <- ggplot(predict.gamma, aes(x = MW, y = perc.yield, color = ligand.alg)) +
  geom_point(alpha = 0.4) +
  geom_abline(slope = 0, intercept = 100, color = "gold", size = 1) +
  facet_grid(ligand.ff~vina) +
  geom_smooth(method = lm) +
  xlab("Molecular Weight") +
  ylab("Percent Yield") +
  scale_color_manual(values=c("seagreen1", "slateblue2"))
ggsave(paste0(gamma.graphs.affinity, "6.png"), plot6affinity)

# look at all these chemical identifiers
gamma.chem <- read.csv(paste0(chem.path, "/GammaID.mmff94.csv"))
gamma.chem.mm2 <- read.csv(paste0(chem.path, "/GammaID.mm2.csv"))
gamma.mol.wt <- gamma.chem[ , c("Name", "MW")]
names(gamma.mol.wt)[names(gamma.mol.wt) == "Name"] <- "formula"
predict.gamma <- left_join(predict.gamma, gamma.mol.wt)
predict.gamma.mw <- predict.gamma[!is.na(predict.gamma$MW), ]
plot10gamma <- ggplot(predict.gamma.mw, aes(x = MW, y = Binding.Affinity, color = ligand.ff)) +
  geom_point(alpha = 0.1) +
  scale_color_manual(values=c("seagreen1", "steelblue1", "orchid")) +
  geom_smooth(method = lm) +
  facet_grid(ligand.alg~vina) +
  xlab("Molecular Weight") +
  ylab("Binding Affinity, kcal/mol") +
  theme_bw()
ggsave(paste0(gamma.graphs, "10.png"), plot = plot10gamma)


# ======================================================================
#                  Organizing solubility data
# ======================================================================
alpha.bindaff <- ds.cleaner[ds.cleaner$host == "1\u03b1", c("guest", "host")]
