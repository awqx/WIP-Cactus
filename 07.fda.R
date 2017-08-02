# Libraries and Packages --------------------------------------------------

library(stringr)
library(tidyverse)

# Data --------------------------------------------------------------------

# Downloaded from http://zinc15.docking.org/substances/subsets/fda/, 
# analyzed via PaDEL-Descriptor
dir.create("./molecules/fdaSDFs")

download.cactus.results <- function(guest, path, chemical.format) {
  report <- tryCatch({
    destfile       <- paste0(path, "/", guest, ".SDF")
    # Chemical format must be parsed to match all the outputs from NCI cactus
    Guest.URL      <- unlist(lapply(guest, URLencode, reserved = T))
    URL            <- paste0(
      "https://cactus.nci.nih.gov/chemical/structure/",
      Guest.URL,
      "/",
      chemical.format
    )
    Map(download.file, url = URL, destfile = destfile)
    data.frame(
      guest = guest,
      downloaded = "yes",
      warning = "no",
      error = "no"
    )
  },
  warning = function(warn) {
    message("Warning: either URL error or already existing directory.")
    destfile       <- paste0(path, "/", guest, ".SDF")
    Guest.URL      <- unlist(lapply(guest, URLencode, reserved = T))
    URL            <- paste0(
      "https://cactus.nci.nih.gov/chemical/structure/",
      Guest.URL, "/", chemical.format
    )
    Map(download.file, url = URL, destfile = destfile)
    data.frame(
      guest = guest,
      downloaded = "yes",
      warning = "yes",
      error = "no"
    )
  },
  error = function(err) {
    message("An error occurred")
    data.frame(
      guest = guest,
      downloaded = "no",
      warning = "yes",
      error = "yes"
    )
  },
  finally = {
    message("Chemical name processed")
  })
  return(report)
}
fda <- read.csv("./molecules/fda/fda.csv")
fda.guest <- fda[ , 2]

guest <- fda.guest %>% as.vector()
guest1 <- fda.guest[1:100] %>% as.vector()

cactus1 <- 
  do.call(
    rbind,
    lapply(
      guest1,
      download.cactus.results,
      path = "./molecules/fdaSDFs",
      chemical.format = "SDF"
    )
  )

fda.results <- 
  do.call(
    rbind,
    lapply(
      guest,
      download.cactus.results,
      path = "./molecules/fdaSDFv2",
      chemical.format = "SDF"
    )
  )


# PaDEL-Descriptor --------------------------------------------------------

# Loading raw data
fda.padel1 <- read.csv("./molecules/fda/fda desc 0-3 kb.csv")
fda.padel2 <- read.csv("./molecules/fda/fda desc 4-5 kb.csv")
# fda.padel3 <- read.csv("./molecules/fda/fda desc 6-8 kb.csv")
# fda.padel4 <- read.csv("./molecules/fda/fda desc 9+ kb.csv")
fda.padel <- rbind(fda.padel1, fda.padel2)
fda.padel <- rbind(fda.padel1, fda.padel2, fda.padel3, fda.padel3)
fda.padel <- rbind(fda.padel, fda.padel, fda.padel)

# fda.alpha <- fda.padel %>% mutate(alpha = 1) %>%
#   mutate(beta = 0) %>% mutate(gamma = 0)
# fda.beta <- fda.padel %>% mutate(alpha = 0) %>%
#   mutate(beta = 1) %>% mutate(gamma = 0)
# fda.gamma <- fda.padel %>% mutate(alpha = 0) %>%
#   mutate(beta = 0) %>% mutate(gamma = 1)

pp.settings <- readRDS("./preprocess.settings.RDS")
zero.pred <- readRDS("./zero.pred.RDS") %>% str_replace(., "-", ".")
zero.pred2 <- readRDS("./zero.pred2.RDS")[1]
too.high <- readRDS("./high.cor.RDS") %>% str_replace(., "-", ".")

fda.desc <- fda.padel[ , !colnames(fda.padel) %in% zero.pred]

fda.desc <- fda.desc[ , -1] 
fda.desc <- fda.desc %>% mutate(alpha = c(rep(1, nrow(fda.desc)/3), rep(0, nrow(fda.desc)/3 * 2))) %>%
  mutate(beta = c(rep(0, nrow(fda.desc)/3), rep(1, nrow(fda.desc)/3), rep(0,nrow(fda.desc)/3))) %>%
  mutate(gamma = c(rep(0,nrow(fda.desc)/3 * 2), rep(1, nrow(fda.desc)/3)))

fda.pp <- predict(pp.settings, fda.desc)

fda.pp <- fda.pp[ , !colnames(fda.pp) %in% too.high]
fda.pp <- fda.pp[ , !colnames(fda.pp) %in% zero.pred2]

fda.data <- cbind(fda.padel[ , 1], fda.pp)
colnames(fda.data)[1] <- "guest"

# Models ------------------------------------------------------------------

#     SVM -----------------------------------------------------------------

svm.alpha <- readRDS("./models/svm/polysvm.alpha.RDS")
svm.beta <- readRDS("./models/svm/polysvm.beta.RDS")
svm.gamma <- readRDS("./models/svm/polysvm.gamma.RDS")

fda.alpha <- predict(svm.alpha, fda.pp)[1:830]
fda.beta <- predict(svm.beta, fda.pp)[831:1660]
fda.gamma <- predict(svm.gamma, fda.pp)[1661:2490]

fda.svm <- c(fda.alpha, fda.beta, fda.gamma)
guest <- fda.padel$Name %>% as.vector()

fda.svm <- data.frame(guest, fda.svm) %>%
  dplyr::rename(., pred = fda.svm) %>%
  mutate(cd.type = c(rep("alpha", 830), rep("beta", 830), rep("gamma", 830)))

ggplot(fda.svm, aes(x = guest, y = pred, color = cd.type, shape = cd.type)) +
  geom_point(alpha = 0.6) + theme_bw() + 
  scale_x_discrete(breaks = NULL) + 
  theme(
    axis.ticks.x=element_blank(),
    axis.line = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
  ) + 
  geom_hline(yintercept = -24, linetype = "dotted", color = "red", size = 1) +
  geom_hline(yintercept = -20, linetype = "dotted", color = "orange", size = 1) +
  labs(title = "Polynomial SVM - FDA Approved Drugs", 
       y = "Predicted DelG, kJ/mol", 
       color = "Cyclodextrin", shape = "Cyclodextrin")
