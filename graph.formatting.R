# Creating visually consistent graphs for posters, slides, etc.
# setwd("~/SREP LAB/qsar")

# Libraries and Packages --------------------------------------------------

library(extrafont)
library(tidyverse)

# Importing Data ----------------------------------------------------------

# Data from docking experiments
#   I'll fix the filepaths later
docking <- readRDS("~/SREP LAB/docking/docking data.RDS")

# Compiled a, b, c-CD data from QSARs
svm <- readRDS("./models/svm/polysvm.tst.results.RDS") %>% rename(host = cd.type)

# Graphs ------------------------------------------------------------------

dir.create("./graphs/2018 poster")
setwd("./graphs/2018 poster")
# Importing a preferred font
font_import(pattern = "[B/b]ahnschrift")
loadfonts(device = "win")

theme.2018 <- theme(
  plot.background = element_rect(fill = "#EFF0F5", color = NA), 
  panel.grid.major = element_line(color = "lightgray"),
  panel.background = element_rect(fill = "white", color = "lightgray"), 
  legend.background = element_rect(fill = "#EFF0F5"),
  legend.key = element_rect(fill = "#EFF0F5", color = NA),
  text = element_text(size = 16, family = "Bahnschrift")
)

plot.2018 <- function(data) {
  ggplot(data, aes(x = obs, y = pred, color = host)) +
    geom_point() + 
    geom_abline(slope = 1, intercept = 0, color = "maroon") +
    theme.2018 + 
    coord_fixed(xlim = c(-45,5), ylim = c(-45, 5)) +
    labs(x = "Experimental DelG", y = "PyRx (Docking) DelG", 
         title = "Docking Calculations",
         color = "CD Type")
}

plot.2018(docking)
ggsave("./2018-02-11 docking.png", scale = 0.85, dpi = 600)

plot.2018(svm)
ggsave("./2018-02-11 svm.png", scale = 0.85, dpi = 600)
