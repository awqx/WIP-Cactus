# Creates visually consistent graphs

# Libraries and Packages --------------------------------------------------

library(data.table)
library(extrafont)
library(tidyverse)

# Functions ---------------------------------------------------------------

font_import(pattern = "[B/b]ahnschrift")
font_import(pattern = "Average")
font_import(pattern = "ClearSans")
loadfonts(device = "win")

theme.2018 <- theme(
  plot.background = element_rect(fill = "#EFF0F5", color = NA), 
  panel.grid.major = element_line(color = "lightgray"),
  panel.background = element_rect(fill = "white", color = "lightgray"), 
  legend.background = element_rect(fill = "#EFF0F5", color = NA),
  legend.key = element_rect(fill = "#EFF0F5", color = NA),
  strip.background = element_rect(fill = "#EFF0F5", color = NA),
  text = element_text(size = 16, family = "Bahnschrift", color = "#404040")
)

theme.isef <- theme(
  plot.background = element_rect(fill = "#EFF0F5", color = NA), 
  panel.grid.major = element_line(color = "lightgray"),
  panel.background = element_rect(fill = "white", color = "lightgray"), 
  legend.background = element_rect(fill = "#EFF0F5", color = NA),
  legend.key = element_rect(fill = "#EFF0F5", color = NA),
  legend.key.size = unit(1.75, 'lines'),
  strip.background = element_rect(fill = "#EFF0F5", color = NA),
  text = element_text(size = 32, family = "Clear Sans Light", color = "#404040")
)

plot.2018 <- function(data) {
  ggplot(data, aes(x = obs, y = pred, color = cd.type)) +
    geom_point() + 
    geom_abline(slope = 1, intercept = 0, color = "maroon") +
    theme.2018 + 
    coord_fixed(xlim = c(-45,5), ylim = c(-45, 5)) +
    labs(x = "Experimental dG, kJ/mol", y = "Predicted dG, kJ/mol",
         color = "CD Type")
}

plot.isef <- function(data) {
  ggplot(data, aes(x = obs, y = pred, color = cd.type)) +
    geom_point() + 
    geom_abline(slope = 1, intercept = 0, color = "#404040") +
    theme.isef + 
    coord_fixed(xlim = c(-35,5), ylim = c(-35, 5)) +
    labs(x = "Experimental dG, kJ/mol", y = "Predicted dG, kJ/mol",
         color = "CD Type")
}

theme.paper.2018 <- theme(
  plot.background = element_rect(fill = "white", color = NA), 
  panel.grid.major = element_line(color = NA),
  panel.background = element_rect(fill = "NA", color = "black"), 
  panel.border = element_rect(fill = NA, color = "black"),
  legend.background = element_rect(fill = "white", color = NA),
  legend.key = element_rect(fill = "white", color = NA),
  strip.background = element_rect(fill = "white", color = NA),
  text = element_text(size = 16, family = "Clear Sans Light")
)

plot.paper.2018 <- function(data) {
  ggplot(data, aes(x = obs, y = pred, color = cd.type)) +
    geom_abline(slope = 1, intercept = 0, color = "black") +
    geom_point(size = 1) + 
    theme.paper.2018 + 
    coord_fixed(xlim = c(-45,5), ylim = c(-45, 5)) +
    labs(x = "Experimental dG, kJ/mol", y = "Predicted dG, kJ/mol",
         color = "CD Type")
}
