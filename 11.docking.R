# Importing data from PyRx and comparing with experimental results

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

# Experimental results
exp.df <- readRDS("./cleaning/02.combined.data.RDS") %>%
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
  rename(pred = DelG.y, obs = DelG.x)
defaultSummary(a.data) # R^2 0.087
# Beta
b.data <- inner_join(b.df, b.docking, by = "guest") %>%
  rename(pred = DelG.y, obs = DelG.x)
defaultSummary(b.data) # 0.161
# Gamma
c.data <- inner_join(c.df, c.docking, by = "guest") %>% 
  rename(pred = DelG.y, obs = DelG.x)
defaultSummary(c.data) # R^2 = 0.00018

# Compiled
a.temp <- a.data %>% mutate(host = "alpha")
b.temp <- b.data %>% mutate(host = "beta")
c.temp <- c.data %>% mutate(host = "gamma")
all.data <- rbind(a.temp, b.temp, c.temp) %>% group_by(host)
defaultSummary(as.data.frame(all.data)) #R2 = 0.078
saveRDS(all.data, "./data/docking.RDS")

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
