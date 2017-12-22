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
# Gamma
c.data <- inner_join(c.df, c.docking, by = "guest")
c.data %>% rename(pred = DelG.y, obs = DelG.x) %>% 
  defaultSummary(c.data) # R^2 = 0.000194

# Graphs ------------------------------------------------------------------
dir.create("./graphs")
ggplot(c.data, aes(x = DelG.x, y = DelG.y)) +
  geom_point() + 
  geom_abline(slope = 1, intercept = 0) +
  theme_bw() + 
  coord_fixed(ratio = 1) + 
  labs(x = "Experimental DelG", "PyRx (Docking) DelG", 
       title = "Gamma-CD Docking Calculations")
ggsave("./graphs/2017-12-22 gamma-cd docking.png")
