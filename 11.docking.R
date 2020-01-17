# Importing data from PyRx and comparing with experimental results

# Libraries and Packages --------------------------------------------------

library(caret)
library(data.table)
library(stringr)
library(tidyverse)
# Chemical analysis
source("http://bioconductor.org/biocLite.R")
biocLite("ChemmineR")
library(ChemmineR)

# Functions ---------------------------------------------------------------

kcal.to.kj <- function(kcal)
{
  kj <- kcal / 0.23912003826
  return(kj)
}

# Data Cleaning -----------------------------------------------------------

# Experimental results

exp.df <- readRDS("./cleaning/03.rename.RDS") %>%
  select(., -data.source)
a.df <- exp.df %>% filter(host == "alpha") %>% 
  select(., -host)
b.df <- exp.df %>% filter(host == "beta") %>%
  select(., -host)

c.df <- exp.df %>% filter(host == "gamma") %>%
  select(., -host)
c.df2 <- readRDS("./cleaning/02.connors.RDS") %>%
  select(., -ka)
c.df <- rbind(c.df, c.df2)

# PyRx docking data
    # Alpha
    # Because PyRx crashes so much, the results are split in several files
a.1 <- read.csv("./data/docking/2017-12-22 alpha-cd docking affinity 1.csv")
a.2 <- read.csv("./data/docking/2017-12-22 alpha-cd docking affinity 2.csv")
a.docking <- rbind(a.1, a.2) %>%
  rename(DelG = Binding.Affinity, guest = Ligand) %>%
  select(guest:DelG) %>%
  mutate(DelG = kcal.to.kj(DelG)) %>%
  mutate(guest = str_replace(guest, "alpha-cd_uff_E=[[:digit:]]+\\.[[:digit:]]+_", "")) %>%
  mutate(guest = str_replace(guest, "_uff_E=[[:digit:]]+\\.[[:digit:]]+", "")) %>% 
  mutate(guest = str_replace_all(guest, "_", " ")) %>%
  data.table(., key = "guest")
a.docking <- a.docking[ , list(DelG = min(DelG)), by = guest]
    # Beta
b.1 <- read.csv("./data/docking/2017-12-22 beta-cd docking affinity 1.csv")
b.2 <- read.csv("./data/docking/2017-12-22 beta-cd docking affinity 2.csv")
b.3 <- read.csv("./data/docking/2017-12-22 beta-cd docking affinity 3.csv", header = F)
colnames(b.3) <- colnames(b.2)
b.docking <- rbind(b.1, b.2, b.3) %>%
  rename(DelG = Binding.Affinity, guest = Ligand) %>%
  select(guest:DelG) %>%
  mutate(DelG = kcal.to.kj(DelG)) %>% 
  mutate(guest = str_replace(guest, "beta-cyclodextrin_", "")) %>%
  mutate(guest = str_replace(guest, "_uff_E=[[:digit:]]+\\.[[:digit:]]+", "")) %>% 
  mutate(guest = str_replace_all(guest, "_", " ")) %>%
  data.table(., key = "guest")
b.docking <- b.docking[ , list(DelG = min(DelG)), by = guest]
    # Gamma
    # Only one file, so rbind not necessary
# Filename used to be 2017-12-21 gamma-cd docking affinity.csv
# Removing files that were messed up from commas
c.docking <- read.csv("./data/docking/gamma-cd docking affinity vina.csv") 
colnames(c.docking) <- c("guest", "DelG", "a", "b")
c.docking <- c.docking[!str_detect(c.docking$a, "[[:alpha:]]"), ]
c.docking <- c.docking[!str_detect(c.docking$DelG, "[[:alpha:]]"), ]
c.docking <- c.docking[!is.na(c.docking$DelG), ]
c.docking <- c.docking %>% 
  mutate(DelG = as.numeric(as.character(DelG))) %>% 
  select(guest:DelG) %>%
  mutate(DelG = kcal.to.kj(as.numeric(DelG))) %>%
  mutate(guest = str_replace(guest, "gamma-cd_uff_E=1768\\.99_", "")) %>%
  mutate(guest = str_replace(guest, "_uff_E=[[:digit:]]+\\.[[:digit:]]+", "")) %>%
  mutate(guest = str_replace_all(guest, "_", " ")) %>%
  mutate(guest = as.factor(guest)) %>%
  data.table(., key = "guest")
c.docking <- c.docking[, list(DelG = min(DelG)),
                       by = guest] 

# Compilation/Calculation -------------------------------------------------

# Using inner join because any independence data points wil be pointless
# Alpha
a.data <- inner_join(a.df, a.docking, by = "guest") %>% 
  rename(pred = DelG.y, obs = DelG.x) %>% 
  as.data.frame()
defaultSummary(a.data) 
# RMSE   Rsquared        MAE 
# 5.51380926 0.05633648 4.33843926 

# Beta
b.data <- inner_join(b.df, b.docking, by = "guest") %>%
  rename(pred = DelG.y, obs = DelG.x) %>%
  as.data.frame()
defaultSummary(b.data) 
# RMSE  Rsquared       MAE 
# 7.0252935 0.1641873 5.6096630 

# Gamma
c.data <- inner_join(c.df, c.docking, by = "guest") %>% 
  rename(pred = DelG.y, obs = DelG.x) %>%
  as.data.frame()
defaultSummary(c.data) 
# RMSE     Rsquared          MAE 
# 9.707894e+00 5.062259e-05 8.042973e+00 

# Compiled
a.temp <- a.data %>% mutate(host = "alpha")
b.temp <- b.data %>% mutate(host = "beta")
c.temp <- c.data %>% mutate(host = "gamma")
all.data <- rbind(a.temp, b.temp, c.temp) %>% group_by(host)

defaultSummary(as.data.frame(all.data)) 
# RMSE  Rsquared       MAE 
# 6.6293336 0.1686069 5.2244895 
# saveRDS(all.data, "./data/docking.RDS")

# Graphs ------------------------------------------------------------------

# dir.create("./graphs")

ggplot(a.data, aes(x = obs, y = pred)) +
  geom_point() + 
  geom_abline(slope = 1, intercept = 0) +
  theme_bw() + 
  coord_fixed() +
  labs(x = "Experimental DelG", y = "PyRx (Docking) DelG", 
       title = "Alpha-CD Docking Calculations")

ggplot(b.data, aes(x = obs, y = pred)) +
  geom_point() + 
  geom_abline(slope = 1, intercept = 0) +
  theme_bw() + 
  coord_fixed() +
  labs(x = "Experimental DelG", y = "PyRx (Docking) DelG", 
       title = "Beta-CD Docking Calculations")

ggplot(c.data, aes(x = obs, y = pred)) +
  geom_point() + 
  geom_abline(slope = 1, intercept = 0) +
  theme_bw() + 
  coord_fixed() + 
  labs(x = "Experimental DelG", y = "PyRx (Docking) DelG", 
       title = "Gamma-CD Docking Calculations")

ggplot(all.data, aes(x = obs, y = pred, color = host, shape = host)) +
  geom_point() + 
  geom_abline(slope = 1, intercept = 0) +
  theme_bw() + 
  coord_fixed() +
  scale_shape_manual(values=c(9, 16, 3)) +
  labs(x = "Experimental DelG", y = "PyRx (Docking) DelG", 
       title = "Cyclodextrin Docking Calculations", 
       color = "CD type", host = "CD type")
# ggsave("./graphs/cd docking.png")

# Analysis of molecules ---------------------------------------------------

docking <- all.data %>% select(., -guest.charge)
# Fairly simple percent error analysis
docking <- docking %>% mutate(perc.error = (pred-obs)/obs*100) %>%
  mutate(pyrx = ifelse(abs(perc.error) <= 15, "good",
                       ifelse(abs(perc.error) <= 60, "meh", "bad")))

# Analysis of alpha
alpha.guest <- docking %>% filter(host == "alpha")
alpha.desc <- readRDS("./pre-process/alpha/2/pp.RDS") %>% 
  select(., -DelG) %>% filter(guest %in% alpha.guest$guest)
alpha.pyrx <- inner_join(alpha.guest, alpha.desc, by = "guest") %>%
  as.data.frame()

# Analysis of beta
beta.guest <- docking %>% filter(host == "beta")
beta.desc <- readRDS("./pre-process/beta/2/pp.RDS") %>% 
  select(., -DelG) %>% filter(guest %in% beta.guest$guest)
beta.pyrx <- inner_join(beta.guest, beta.desc, by = "guest") %>%
  as.data.frame()

# Analysis of gamma
gamma.guest <- docking %>% filter(host == "gamma")
gamma.desc <- readRDS("./pre-process/gamma/2/pp.RDS") %>% 
  select(., -DelG) %>% filter(guest %in% gamma.guest$guest)
gamma.pyrx <- inner_join(gamma.guest, gamma.desc, by = "guest") %>%
  as.data.frame()

ggplot(docking %>% filter(pyrx == "good"), aes(x = obs, y = pred, color = host)) +
  geom_point() + 
  geom_abline(slope = 1, intercept = 0) +
  theme_bw() + 
  coord_fixed() +
  labs(x = "Experimental DelG", y = "PyRx (Docking) DelG", 
       title = "Cyclodextrin Docking Calculations", 
       color = "CD type")
ggplot(docking %>% filter(pyrx == "meh"), aes(x = obs, y = pred, color = host)) +
  geom_point() + 
  geom_abline(slope = 1, intercept = 0) +
  theme_bw() + 
  coord_fixed() +
  labs(x = "Experimental DelG", y = "PyRx (Docking) DelG", 
       title = "Cyclodextrin Docking Calculations", 
       color = "CD type")
ggplot(docking %>% filter(pyrx == "bad"), aes(x = obs, y = pred, color = host)) +
  geom_point() + 
  geom_abline(slope = 1, intercept = 0) +
  theme_bw() + 
  coord_fixed() +
  labs(x = "Experimental DelG", y = "PyRx (Docking) DelG", 
       title = "Cyclodextrin Docking Calculations", 
       color = "CD type")
ggplot(docking, aes(x = obs, y = pred, color = pyrx, shape = host)) +
  geom_point() + 
  geom_abline(slope = 1, intercept = 0) +
  theme_bw() + 
  coord_fixed() +
  labs(x = "Experimental DelG", y = "PyRx (Docking) DelG", 
       title = "Cyclodextrin Docking Calculations", 
       color = "PyRx prediction", 
       shape = "CD type") + 
  scale_shape_manual(values=c(9, 16, 3)) 

# Important features ------------------------------------------------------

ctrl <- rfeControl(functions = nbFuncs, 
                   method = "repeatedcv", 
                   repeats = 5, 
                   verbose = T)
subsets <- c(1, 5, 10, 15, 25, 50)


alpha.x <- alpha.pyrx %>% select(., -guest:-pyrx)
alpha.y <- alpha.pyrx$pyrx %>% as.factor()
alpha.rfe <- rfe(x = alpha.x, y = alpha.y, 
                 sizes = subsets, rfeControl = ctrl)
predictors(alpha.rfe)


beta.x <- beta.pyrx %>% select(., -guest:-pyrx)
beta.y <- beta.pyrx$pyrx %>% as.factor()
beta.rfe <- rfe(x = beta.x, y = beta.y, 
                 sizes = subsets, rfeControl = ctrl)
predictors(beta.rfe)

gamma.x <- gamma.pyrx %>% select(., -guest:-pyrx)
gamma.y <- gamma.pyrx$pyrx %>% as.factor()
gamma.rfe <- rfe(x = gamma.x, y = gamma.y, 
                sizes = subsets, rfeControl = ctrl)

# Chemical similarity -----------------------------------------------------

compile.sdf <- function(path, name) {
  sdf <- read.csv(paste0(path, "/", name, ".SDF"), header = F)
  sdf[ , 1] <- as.character(sdf[ , 1])
  sdf[1, 1] <- name
  return(sdf)
}

calc.similarity <- function(mol1, mol2, sdfset) {
  if (mol1 == mol2)
    score <- 1
  else
    score <- cmp.similarity(sdf2ap(sdfset[[mol1]]), sdf2ap(sdfset[[mol2]]))
  if (mol1 %% 10 == 0 && mol2 %% 10 == 0)
    message(paste0(mol1, " and ", mol2, " completed."))
  return(data.frame(
    mol1 = mol1, mol2 = mol2, similarity = score
  ))
}

#     Alpha ---------------------------------------------------------------

# For some reason, read.SDFset works best with a single .SDF, so 
# the "good" predictions on alpha must be compiled
alpha.guest.good <- alpha.guest %>% filter(pyrx == "good") %>% .$guest
alpha.good.sdf <- do.call(rbind, 
                          lapply(FUN = compile.sdf, 
                                 X = alpha.guest.good, 
                                 path = "./molecules/alphaCD"))
write.table(alpha.good.sdf, "./data/docking/alpha.good.SDF", quote = F, 
            row.names = F, col.names = F)
alpha.sdfset <- read.SDFset("./data/docking/alpha.good.SDF")

# Obtaining the atom pairs isn't working on the entire SDF list, 
# so "manual' comparison must be used instead
alpha.combos <- expand.grid(1:length(alpha.guest.good), 
                            1:length(alpha.guest.good)) %>%
  rename(mol1 = Var1, mol2 = Var2)
mol1.combos <- alpha.combos$mol1
mol2.combos <- alpha.combos$mol2
alpha.good.similarity <- mapply(FUN = calc.similarity, 
                                mol1 = mol1.combos, 
                                mol2 = mol2.combos, 
                                MoreArgs = list(sdfset = alpha.sdfset), 
                                SIMPLIFY = F)
alpha.good.similarity <- do.call(rbind, alpha.good.similarity)

# A graph of similarity in the dataset
ggplot(alpha.good.similarity, aes(x = mol1, y = mol2, fill = similarity)) + 
  geom_raster() + 
  coord_fixed() + 
  theme_bw() + 
  scale_fill_gradientn(colors = terrain.colors(10)) + 
  labs(x = "Molecule 1", y = "Molecule 2", 
       title = "Alpha-CD: Good PyRx predictions", 
       fill = 'Similarity score')
ggplot(alpha.good.similarity, aes(x = similarity)) + 
  geom_histogram() + 
  theme_bw() + 
  labs(x = "Similarity score", title = "Alpha-CD: Good PyRx predictions")

# Repeating with "meh"
alpha.guest.meh <- alpha.guest %>% filter(pyrx == "meh") %>% .$guest
alpha.meh.sdf <- do.call(rbind, 
                          lapply(FUN = compile.sdf, 
                                 X = alpha.guest.meh, 
                                 path = "./molecules/alphaCD"))
write.table(alpha.meh.sdf, "./data/docking/alpha.meh.SDF", quote = F, 
            row.names = F, col.names = F)
alpha.sdfset <- read.SDFset("./data/docking/alpha.meh.SDF")

alpha.combos <- expand.grid(1:length(alpha.guest.meh), 
                            1:length(alpha.guest.meh)) %>%
  rename(mol1 = Var1, mol2 = Var2)
mol1.combos <- alpha.combos$mol1
mol2.combos <- alpha.combos$mol2
alpha.meh.similarity <- mapply(FUN = calc.similarity, 
                                mol1 = mol1.combos, 
                                mol2 = mol2.combos, 
                                MoreArgs = list(sdfset = alpha.sdfset), 
                                SIMPLIFY = F)
alpha.meh.similarity <- do.call(rbind, alpha.meh.similarity)

ggplot(alpha.meh.similarity, aes(x = mol1, y = mol2, fill = similarity)) + 
  geom_raster() + 
  coord_fixed() + 
  theme_bw() + 
  scale_fill_gradientn(colors = terrain.colors(10)) + 
  labs(x = "Molecule 1", y = "Molecule 2", 
       title = "Alpha-CD: Meh PyRx predictions", 
       fill = 'Similarity score')
ggplot(alpha.meh.similarity, aes(x = similarity)) + 
  geom_histogram() + 
  theme_bw() + 
  labs(x = "Similarity score", title = "Alpha-CD: Meh PyRx predictions")

# Repeating with "bad"
alpha.guest.bad <- alpha.guest %>% filter(pyrx == "bad") %>% .$guest
alpha.bad.sdf <- do.call(rbind, 
                         lapply(FUN = compile.sdf, 
                                X = alpha.guest.bad, 
                                path = "./molecules/alphaCD"))
write.table(alpha.bad.sdf, "./data/docking/alpha.bad.SDF", quote = F, 
            row.names = F, col.names = F)
alpha.sdfset <- read.SDFset("./data/docking/alpha.bad.SDF")
valid <- validSDF(alpha.sdfset)
alpha.sdfset <- alpha.sdfset[valid]

n <- length(alpha.guest.bad) - sum(!valid)
alpha.combos <- expand.grid(1:n, 1:n) %>%
  rename(mol1 = Var1, mol2 = Var2)
mol1.combos <- alpha.combos$mol1
mol2.combos <- alpha.combos$mol2
alpha.bad.similarity <- mapply(FUN = calc.similarity, 
                               mol1 = mol1.combos, 
                               mol2 = mol2.combos, 
                               MoreArgs = list(sdfset = alpha.sdfset), 
                               SIMPLIFY = F)
alpha.bad.similarity <- do.call(rbind, alpha.bad.similarity)

ggplot(alpha.bad.similarity, aes(x = mol1, y = mol2, fill = similarity)) + 
  geom_raster() + 
  coord_fixed() + 
  theme_bw() + 
  scale_fill_gradientn(colors = terrain.colors(10)) + 
  labs(x = "Molecule 1", y = "Molecule 2", 
       title = "Alpha-CD: Bad PyRx predictions", 
       fill = 'Similarity score')
ggplot(alpha.bad.similarity, aes(x = similarity)) + 
  geom_histogram() + 
  theme_bw() + 
  labs(x = "Similarity score", title = "Alpha-CD: Bad PyRx predictions")

#     Beta ----------------------------------------------------------------

beta.guest.good <- beta.guest %>% filter(pyrx == "good") %>% .$guest
beta.good.sdf <- do.call(rbind, 
                          lapply(FUN = compile.sdf, 
                                 X = beta.guest.good, 
                                 path = "./molecules/betaCD"))
write.table(beta.good.sdf, "./data/docking/beta.good.SDF", quote = F, 
            row.names = F, col.names = F)
beta.sdfset <- read.SDFset("./data/docking/beta.good.SDF")

# Obtaining the atom pairs isn't working on the entire SDF list, 
# so "manual' comparison must be used instead
beta.combos <- expand.grid(1:92, 1:92) %>%
  rename(mol1 = Var1, mol2 = Var2)
mol1.combos <- beta.combos$mol1
mol2.combos <- beta.combos$mol2
beta.good.similarity <- mapply(FUN = calc.similarity, 
                                mol1 = mol1.combos, 
                                mol2 = mol2.combos, 
                                MoreArgs = list(sdfset = beta.sdfset), 
                                SIMPLIFY = F)
beta.good.similarity <- do.call(rbind, beta.good.similarity)

# A graph of similarity in the dataset
ggplot(beta.good.similarity, aes(x = mol1, y = mol2, fill = similarity)) + 
  geom_raster() + 
  coord_fixed() + 
  theme_bw() + 
  scale_fill_gradientn(colors = terrain.colors(10)) + 
  labs(x = "Molecule 1", y = "Molecule 2", 
       title = "beta-CD: Good PyRx predictions", 
       fill = 'Similarity score')
ggplot(beta.good.similarity, aes(x = similarity)) + 
  geom_histogram() + 
  theme_bw() + 
  labs(x = "Similarity score", title = "beta-CD: Good PyRx predictions")

# Repeating with "meh"
beta.guest.meh <- beta.guest %>% filter(pyrx == "meh") %>% .$guest
beta.meh.sdf <- do.call(rbind, 
                         lapply(FUN = compile.sdf, 
                                X = beta.guest.meh, 
                                path = "./molecules/betaCD"))
write.table(beta.meh.sdf, "./data/docking/beta.meh.SDF", quote = F, 
            row.names = F, col.names = F)
beta.sdfset <- read.SDFset("./data/docking/beta.meh.SDF")

beta.combos <- expand.grid(1:length(beta.guest.meh), 
                            1:length(beta.guest.meh)) %>%
  rename(mol1 = Var1, mol2 = Var2)
mol1.combos <- beta.combos$mol1
mol2.combos <- beta.combos$mol2
beta.meh.similarity <- mapply(FUN = calc.similarity, 
                               mol1 = mol1.combos, 
                               mol2 = mol2.combos, 
                               MoreArgs = list(sdfset = beta.sdfset), 
                               SIMPLIFY = F)
beta.meh.similarity <- do.call(rbind, beta.meh.similarity)

ggplot(beta.meh.similarity, aes(x = mol1, y = mol2, fill = similarity)) + 
  geom_raster() + 
  coord_fixed() + 
  theme_bw() + 
  scale_fill_gradientn(colors = terrain.colors(10)) + 
  labs(x = "Molecule 1", y = "Molecule 2", 
       title = "beta-CD: Meh PyRx predictions", 
       fill = 'Similarity score')
ggplot(beta.meh.similarity, aes(x = similarity)) + 
  geom_histogram() + 
  theme_bw() + 
  labs(x = "Similarity score", title = "beta-CD: Meh PyRx predictions")

# Repeating with "bad"
beta.guest.bad <- beta.guest %>% filter(pyrx == "bad") %>% .$guest
beta.bad.sdf <- do.call(rbind, 
                         lapply(FUN = compile.sdf, 
                                X = beta.guest.bad, 
                                path = "./molecules/betaCD"))
write.table(beta.bad.sdf, "./data/docking/beta.bad.SDF", quote = F, 
            row.names = F, col.names = F)
beta.sdfset <- read.SDFset("./data/docking/beta.bad.SDF")
valid <- validSDF(beta.sdfset)
beta.sdfset <- beta.sdfset[valid]

n <- length(beta.guest.bad) - sum(!valid)
beta.combos <- expand.grid(1:n, 1:n) %>%
  rename(mol1 = Var1, mol2 = Var2)
mol1.combos <- beta.combos$mol1
mol2.combos <- beta.combos$mol2
beta.bad.similarity <- mapply(FUN = calc.similarity, 
                               mol1 = mol1.combos, 
                               mol2 = mol2.combos, 
                               MoreArgs = list(sdfset = beta.sdfset), 
                               SIMPLIFY = F)
beta.bad.similarity <- do.call(rbind, beta.bad.similarity)

ggplot(beta.bad.similarity, aes(x = mol1, y = mol2, fill = similarity)) + 
  geom_raster() + 
  coord_fixed() + 
  theme_bw() + 
  scale_fill_gradientn(colors = terrain.colors(10)) + 
  labs(x = "Molecule 1", y = "Molecule 2", 
       title = "beta-CD: Bad PyRx predictions", 
       fill = 'Similarity score')
ggplot(beta.bad.similarity, aes(x = similarity)) + 
  geom_histogram() + 
  theme_bw() + 
  labs(x = "Similarity score", title = "beta-CD: Bad PyRx predictions")

# Descriptor vs. Perc Error -----------------------------------------------

desc.cols <- colnames(alpha.guest)
# alpha.amw <- alpha.pyrx %>% select(desc.cols, "AMW")
# ggplot(alpha.amw, aes(x = AMW, y = perc.error, color = pyrx)) + 
#   geom_point() + 
#   theme_bw()
alpha.xlogp <- alpha.pyrx %>% select(desc.cols, "XLogP")
# a slight pattern, especially with "bad" 
ggplot(alpha.xlogp, aes(x = XLogP, y = perc.error, color = pyrx, shape = pyrx)) + 
  geom_point() + 
  theme_bw() + 
  coord_cartesian(ylim = c(-100, 200))

alpha.nring <- alpha.pyrx %>% select(desc.cols, "nRing")
ggplot(alpha.nring, aes(x = nRing, y = perc.error, color = pyrx, shape = pyrx)) + 
  geom_jitter() + 
  theme_bw() + 
  coord_cartesian(ylim = c(-100, 200))

alpha.chain <- alpha.pyrx %>% select(desc.cols, "nAtomLC")
ggplot(alpha.chain, aes(x = nAtomLC, y = perc.error, color = pyrx, shape = pyrx)) + 
  geom_jitter() + 
  theme_bw() + 
  coord_cartesian(ylim = c(-100, 200))



beta.amw <- beta.pyrx %>% select(desc.cols, "AMW")
ggplot(beta.amw, aes(x = AMW, y = perc.error, color = pyrx)) +
  geom_point() +
  theme_bw()

beta.xlogp <- beta.pyrx %>% select(desc.cols, "XLogP")
ggplot(beta.xlogp, aes(x = XLogP, y = perc.error, color = pyrx, shape = pyrx)) + 
  geom_point() + 
  theme_bw() 

beta.chain <- beta.pyrx %>% select(desc.cols, "nAtomLC")
ggplot(beta.chain, aes(x = nAtomLC, y = perc.error, color = pyrx, shape = pyrx)) + 
  geom_jitter() + 
  theme_bw() + 
  coord_cartesian(ylim = c(-100, 200))
