# Importing data from PyRx and comparing with experimental results
# setwd("~/SREP LAB/docking")

# Libraries and Packages --------------------------------------------------

library(caret)
library(data.table)
library(stringr)
library(tidyverse)

# Functions ---------------------------------------------------------------

kcal.to.kj <- function(kcal)
{
  kj <- kcal / 0.23912003826
  return(kj)
}

# Data Cleaning -----------------------------------------------------------

# Experimental data
df <- readRDS("~/SREP LAB/qsar/padel.pp.new.RDS") %>%
  select(guest:DelG)
    # Sorting into alpha, beta, and gamma-cd
a.df <- df %>% filter(host == "alpha") %>% select(-host)
b.df <- df %>% filter(host == "beta") %>% select(-host)
c.df <- df %>% filter(host == "gamma") %>% select(-host)

# PyRx docking data
    # Alpha
    # Because PyRx crashes so much, the results are split in several files

a.1 <- read.csv("./results/2017-12-22 alpha-cd docking affinity 1.csv")
a.2 <- read.csv("./results/2017-12-22 alpha-cd docking affinity 2.csv")
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
b.1 <- read.csv("./results/2017-12-22 beta-cd docking affinity 1.csv")
b.2 <- read.csv("./results/2017-12-22 beta-cd docking affinity 2.csv")
b.3 <- read.csv("./results/2017-12-22 beta-cd docking affinity 3.csv", header = F)
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
c.docking <- read.csv("./results/2017-12-21 gamma-cd docking affinity.csv") %>%
  rename(DelG = Binding.Affinity, guest = Ligand) %>%
  select(guest:DelG) %>%
  mutate(DelG = kcal.to.kj(DelG)) %>%
  mutate(guest = str_replace(guest, "gamma-cd_uff_E=[[:digit:]]+\\.[[:digit:]]+_", "")) %>%
  mutate(guest = str_replace(guest, "_uff_E=[[:digit:]]+\\.[[:digit:]]+", "")) %>%
  mutate(guest = str_replace(guest, "_uff_E=294.49", "")) %>%
  mutate(guest = str_replace_all(guest, "_", " ")) %>%
  data.table(., key = "guest")
c.docking <- c.docking[, list(DelG = min(DelG)),
                       by = guest] 
# Compilation/Calculation -------------------------------------------------

# Using inner join because any independence data points wil be pointless
# Alpha
a.data <- inner_join(a.df, a.docking, by = "guest") %>% 
  rename(pred = DelG.y, obs = DelG.x)
defaultSummary(a.data) # R^2 0.0999244
# Beta
b.data <- inner_join(b.df, b.docking, by = "guest") %>%
  rename(pred = DelG.y, obs = DelG.x)
defaultSummary(b.data) # 0.1367831
# Gamma
c.data <- inner_join(c.df, c.docking, by = "guest") %>% 
  rename(pred = DelG.y, obs = DelG.x)
defaultSummary(c.data) # R^2 = 0.000194

# Compiled
a.temp <- a.data %>% mutate(host = "alpha")
b.temp <- b.data %>% mutate(host = "beta")
c.temp <- c.data %>% mutate(host = "gamma")
all.data <- rbind(a.temp, b.temp, c.temp) %>% group_by(host)
defaultSummary(as.data.frame(all.data)) #R2 = 0.1759711

# Graphs ------------------------------------------------------------------

# dir.create("./graphs")

ggplot(a.data, aes(x = obs, y = pred)) +
  geom_point() + 
  geom_abline(slope = 1, intercept = 0) +
  theme_bw() + 
  coord_fixed(xlim = c(-45,5), ylim = c(-45, 5)) +
  labs(x = "Experimental DelG", y = "PyRx (Docking) DelG", 
       title = "Alpha-CD Docking Calculations")
ggsave("./graphs/2017-12-22 alpha-cd docking.png")

ggplot(b.data, aes(x = obs, y = pred)) +
  geom_point() + 
  geom_abline(slope = 1, intercept = 0) +
  theme_bw() + 
  coord_fixed(xlim = c(-45,5), ylim = c(-45, 5)) +
  labs(x = "Experimental DelG", y = "PyRx (Docking) DelG", 
       title = "Beta-CD Docking Calculations")
ggsave("./graphs/2017-12-22 beta-cd docking.png")

ggplot(c.data, aes(x = obs, y = pred)) +
  geom_point() + 
  geom_abline(slope = 1, intercept = 0) +
  theme_bw() + 
  coord_fixed(xlim = c(-45,5), ylim = c(-45, 5)) + 
  labs(x = "Experimental DelG", y = "PyRx (Docking) DelG", 
       title = "Gamma-CD Docking Calculations")
ggsave("./graphs/2017-12-22 gamma-cd docking.png")

ggplot(all.data, aes(x = obs, y = pred, color = host)) +
  geom_point() + 
  geom_abline(slope = 1, intercept = 0) +
  theme_bw() + 
  coord_fixed(xlim = c(-45,5), ylim = c(-45, 5)) +
  labs(x = "Experimental DelG", y = "PyRx (Docking) DelG", 
       title = "Cyclodextrin Docking Calculations")
ggsave("./graphs/2017-12-22 cd docking.png")
