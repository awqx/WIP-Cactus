# Integrating DelG results from the QSARs with Fu's mathematical modeling of
# complexation dynamics. Meant to produce estimates of cumulative drug release
# over time
#
# Source: https://doi.org/10.1007/s10439-011-0336-z

# Libraries and Packages --------------------------------------------------

library(data.table)
library(deSolve)
library(dplyr)
library(tidyverse)

# Functions ---------------------------------------------------------------

# Well, if we know that DelG = RTlnKeq, then we should invert that
# In ths case, b/c drug + CD <--> Complex, Keq = [complex]/[cd][drug]
# And therefore, Keq = Ka

# assuming delg is in kj/mol
convert.delg.ka <- function(delg) {
  # exp(n) = e^n
  joules <- delg*1000
  return (exp(joules / (-8.314 * 298)))
}

convert.delg.kd <- function(delg) {
  return (1/convert.delg.ka(delg))
}

# Diffusivity can be related to molecular weight with Vulic's modified
# Stokes-Einstein-Sutherland equation
# 
# This is probably wrong in terms of how the units scale
get.diffusivity <- function(mol.wt) {
  # assuming the conditions are standardized to 298 K
  # and the release media is water
  diff <- 9.940 * 10^-15 * 298 / (0.8921 * mol.wt^(1/3))
  # assuming the output is in m^2*s^-1...
  # correcting units to cm^2*h^-1
  diff <- diff * 10000 * 3600
  return(diff)
}

accumulate.new.dg <- function(dg) {
  ka <- convert.delg.ka(dg)
  k1 <- k2*ka
  P1 <<- P1 <- k1/k2*c0
  out.temp <- ode(y = state, times = times, func = drug.model, parms = NULL)
  dmr.sum <- accumulate.dMr(out = out.temp) %>% .[-1, ] %>%
    mutate(dG = dg)
  return(dmr.sum)
} 
# Parameters --------------------------------------------------------------

# Known
ka <- convert.delg.ka(-18) # benzoic acid
ct <- 0.22 # at 1 g/4 mL (Fu)
D <- 0.1 # technically searchable, but this is arbitrary
# just using benzoic acid
D <- get.diffusivity(122.123)
width <- 0.05 # width of half of hydrogel
vol <- 0.0785 # cm^2, vol when swollen

# Reasonably derivable
k2 <- 35 # h^-1, according to Fu
k1 <- k2*ka # ka = k1/k2 

# Educated guesses
c0 <- 0.15 # total drug concentration
cc <- 0.08 # concentration of uncomplexed cd
clc.eq <- ct - cc
cl.eq <- c0 - clc.eq
# # Checking for approximate accuracy
# clc.eq
# c0 - c0/(1 + ka*cc)

# Dimensionless
P1 <- k1/k2*c0
P2 <- D/(k2*width^2)
P2 <- 0.8 # approximation from Fu
P3 <- ct/c0

# Dimensionless ODE Implementation ----------------------------------------

# Functions
drug.model <- function(t, state, parms) {
  drug <- state[1:n]
  comp <- state[(n+1):(2*n)]
  dcomp <- P1 * drug * (P3 - comp) - k2*comp
  
  ddrug <- rep(0, n)
  ddrug[1] <- P2 * (drug[2] - drug[1]) / (del^2) - dcomp[1]
  ddrug[n - 1] <-
    P2*(-2*drug[n-1]+drug[n-2])/(del^2) - dcomp[n - 1]
  for (i in 2:(n-2)) {
    ddrug[i] <- P2*(drug[i+1]-2*drug[i]+drug[i-1])/(del^2)-dcomp[i]
  }
  dt <- c(ddrug, dcomp)
  return(list(c(dt)))
}

# out is the result of ode with times 0 to t and layers 1 to n
get.dMr <- function(out) {
  times <- seq(1, t, 1)
  rr <- rep(0, (t-1))
  for(i in 1:t) {
    rr[i] <- -0.5 * P2 * (out[i, n] - out[i , (n - 1)])/del
  }
  return(data.frame(times, rr))
}

accumulate.dMr <- function(out) {
  times <- seq(1, t, 1)
  rr <- rep(0, (t-1))
  rr[1] <- -0.5 * P2 * (out[1, n] - out[1 , (n - 1)])/del
  for(i in 2:t) {
    rr[i] <- rr[i - 1] - 0.5 * P2 * (out[i, n] - out[i , (n - 1)])/del
  }
  return(data.frame(times, rr))
}

# Parameters

n <- 100 # number of layers
del <- 1/n # width of a layer
del <- 0.1 # I don't know why this value works 
state <- c(rep(cl.eq, n), 
           rep(clc.eq, n))
t <- 1000
times <- seq(0, t, 1) 
out <- ode(y = state, times = times, func = drug.model, parms = NULL)

drug.out <- out[ -1, c(1:n)] %>% as.data.frame()
drug.out <- data.table::melt(setDT(drug.out), id.vars = "time",
                             variable.name = "location")

ggplot(drug.out, aes(x = time, y = location, fill = value)) +
  geom_raster()
# 
# dmr.raw <- get.dMr(out = out) %>% .[-1, ]
# ggplot(dmr.raw, aes(x = times, y = rr)) +
#   geom_line()
dmr.sum <- accumulate.dMr(out = out) %>% .[-1, ]
ggplot(dmr.sum, aes(x = times, y = rr)) + 
  geom_line() + 
  theme_bw() + 
  labs(x = "Time, in hours", y = "Cumulative proportion of drug released")

# Compilation of several
drug.out1 <- drug.out
ggplot(drug.out, aes(x = time, y = location, fill = value)) +
  geom_raster() + 
  theme.paper.2018 + 
  theme(text = element_text(size = 12, family = "Clear Sans Light"), 
        panel.border = element_rect(fill = NA, color = "white"), 
        axis.text.y=element_blank(), axis.ticks.y = element_blank()) + 
  labs(fill = "[Free drug]",
       y = "Distance from center of hydrogel", 
       x = "Time")
ggsave("./graphs/presentation/diffusion.png", dpi = 600)


dmr.sum.all <- do.call(rbind, lapply(c(-5, -10, -15, -18, -20, -25), 
                                     accumulate.new.dg)) %>%
  mutate(dG = as.factor(dG))
ggplot(dmr.sum.all, aes(x = times, y = rr, color = dG)) + 
  geom_line() + 
  # theme.paper.2018 + 
  labs(x = "Time, in hours", y = "Cumulative proportion of drug released")
ggsave("./graphs/presentation/profiles.png", dpi = 600)


# Testing RMSE ranges -----------------------------------------------------

# Results from external validation
ev.results <- readRDS("./results/ensemble.RDS")
# (The following observations are limited to beta-CD)
# As an example, a strong binder would be protriptyline and a weak binder
# would be d-arabinose.
# Ensemble predictions: 
#     -23.143283 kJ/mol for protriptyline
#     -6.603071 kJ/mol for d-arabinose

# The RMSE of prediction for the beta-CD ensemble is 2.6439926 kJ/mol. 
# The RMSEof prediction for PyRx on beta-CD is 7.0252935 kJ/mol (weird, 
# that's higher than I remember. Will recalculate if necessary). 

# Adjusting parameters for the molecules
# protriptyline
set.strong <- function() { 
  # Known
  ka <<- convert.delg.ka(protriptyline) # benzoic acid
  # ct <<- 0.22 # at 1 g/4 mL (Fu)
  D <<- 0.1 # technically searchable, but this is arbitrary
  # just using benzoic acid
  D <<- get.diffusivity(263.384)
  # width <<- 0.05 # width of half of hydrogel
  # vol <<- 0.0785 # cm^2, vol when swollen
  
  # Reasonably derivable
  k2 <<- 35 # h^-1, according to Fu
  k1 <<- k2*ka # ka = k1/k2 
  
  # Educated guesses
  # c0 <<- 0.15 # total drug concentration
  # cc <<- 0.08 # concentration of uncomplexed cd
  # clc.eq <<- ct - cc
  # cl.eq <<- c0 - clc.eq
  # # Checking for approximate accuracy
  # clc.eq
  # c0 - c0/(1 + ka*cc)
  
  # Dimensionless
  P1 <<- k1/k2*c0
  # P2 <<- D/(k2*width^2)
  # P2  <<- 0.8 # approximation from Fu
  # P3 <<- ct/c0
}

set.weak <- function() { 
  # Known
  ka <<- convert.delg.ka(d_arabinose) 
  D <<- get.diffusivity(150.13)
  
  k2 <<- 35 # h^-1, according to Fu
  k1 <<- k2*ka # ka = k1/k2 
  
  # Dimensionless
  P1 <<- k1/k2*c0
  # P2 <<- D/(k2*width^2)
  # P2  <<- 0.8 # approximation from Fu
  # P3 <<- ct/c0
}

protriptyline <- -23.143
d_arabinose <- -6.603
qsar.sd <- 2.644
docking.sd <- 7.025
set.seed(1337)
strong.qsar.range <- rnorm(n = 15, mean = protriptyline, sd = qsar.sd)
weak.qsar.range <- rnorm(15, d_arabinose, qsar.sd)
strong.docking.range <- rnorm(15, protriptyline, docking.sd)
weak.docking.range <- rnorm(15, d_arabinose, docking.sd)

# Strong binding 
set.strong()
strong.qsar.release <- do.call(rbind, 
        lapply(strong.qsar.range, accumulate.new.dg)) %>% 
  mutate(dG = as.factor(dG)) %>% 
  mutate(method = "QSAR", strength = "strong")
strong.docking.release <- do.call(rbind, 
                              lapply(strong.docking.range, accumulate.new.dg)) %>% 
  mutate(dG = as.factor(dG)) %>% 
  mutate(method = "Docking", strength = "strong")
strong.release <- rbind(strong.qsar.release, strong.docking.release)

ggplot(strong.qsar.release, aes(x = times, y = rr, color = dG)) + 
  geom_line() + 
  theme_bw() + 
  # theme.paper.2018 + 
  labs(x = "Time, in hours", y = "Cumulative drug release",
       title = "QSAR RMSE range in release curves")
ggplot(strong.docking.release, aes(x = times, y = rr, color = dG)) + 
  geom_line() + 
  theme_bw() + 
  # theme.paper.2018 + 
  labs(x = "Time, in hours", y = "Cumulative drug release",
       title = "Docking RMSE range in relese curves")
ggplot(strong.release, aes(x = times, y = rr, color = method)) + 
  geom_line() + 
  theme_bw() + 
  # theme.paper.2018 + 
  labs(x = "Time", y = "Cumulative drug release", 
       title = "Release profiles for a strong binder (protriptyline)")
ggsave("./graphs/release profile, strong, block.png")

ggplot(strong.release, aes(x = times, y = rr, 
                           color = method, group = dG)) + 
  geom_line() + 
  theme_bw() + 
  # theme.paper.2018 + 
  labs(x = "Time", y = "Cumulative drug release", 
       title = "Release profiles for a strong binder (protriptyline)")
ggsave("./graphs/release profile, strong, line.png")

# Weak binding
set.weak() 
weak.qsar.release <- do.call(rbind, 
                               lapply(weak.qsar.range, accumulate.new.dg)) %>% 
  mutate(dG = as.factor(dG)) %>% 
  mutate(method = "QSAR", strength = "weak")
weak.docking.release <- do.call(rbind, 
                                  lapply(weak.docking.range, accumulate.new.dg)) %>% 
  mutate(dG = as.factor(dG)) %>% 
  mutate(method = "Docking", strength = "weak")
weak.release <- rbind(weak.qsar.release, weak.docking.release)

ggplot(weak.qsar.release, aes(x = times, y = rr, color = dG)) + 
  geom_line() + 
  theme_bw() + 
  # theme.paper.2018 + 
  labs(x = "Time, in hours", y = "Cumulative drug release",
       title = "QSAR RMSE range in release curves")
ggplot(weak.docking.release, aes(x = times, y = rr, color = dG)) + 
  geom_line() + 
  theme_bw() + 
  # theme.paper.2018 + 
  labs(x = "Time, in hours", y = "Cumulative drug release",
       title = "Docking RMSE range in relese curves")
ggplot(weak.release, aes(x = times, y = rr, color = method)) + 
  geom_line() + 
  theme_bw() + 
  # theme.paper.2018 + 
  labs(x = "Time", y = "Cumulative drug release", 
       title = "Release profiles for a weak binder")
ggsave("./graphs/release profile, weak, block.png")
ggplot(weak.release, aes(x = times, y = rr, color = method, group = dG)) + 
  geom_line() + 
  theme_bw() + 
  # theme.paper.2018 + 
  labs(x = "Time, in hours", y = "Cumulative drug release", 
       title = "Release profiles for a weak binder", 
       color = "Method")
ggsave("./graphs/release profile, weak, line.png")


# ODE with dimensions -----------------------------------------------------

# Testing values ----------------------------------------------------------

dir.create("./release")

# A function to reset the parameters to certain values
# because for some reason ODE works best with global parameters only
set.param <- function() {
  ka <<- convert.delg.ka(-15) # choose a value from the QSPR
  ct <<- 0.22 # at 1 g/4 mL (Fu)

  # Reasonably derivable
  k2 <<- 35 # h^-1, according to Fu
  k1 <<- k2*ka # ka = k1/k2 
  
  # Educated guesses
  c0 <<- 0.15 # total drug concentration
  cc <<- 0.08 # concentration of uncomplexed cd
  clc.eq <<- ct - cc
  # # Checking for approximate accuracy
  # clc.eq
  # c0 - c0/(1 + ka*cc)
  
  # Dimensionless
  P1 <<- k1/k2*c0
  P2 <<- 0.0002753496
  P3 <<- ct/c0
}

#     P1 ------------------------------------------------------------------

# Provides a cumulative release dataframe for a new value of P1
vary.P1 <- function(newP1) {
  P1 <<- newP1
  ode.out <- ode(y = state, times = times, func = drug.model, parms = NULL)
  result <- accumulate.dMr(out = ode.out) %>% .[-1, ] %>%
    mutate(param = as.factor(newP1))
  return(result)
}

P1.variance <- do.call(rbind, lapply(seq(20, 320, 50),
                                     FUN = vary.P1))

ggplot(P1.variance, aes(x = times, y = rr, color = param)) + 
  geom_line() + 
  theme_bw() 

# Simple modification of P1 that takes inputted DelG values
vary.delG <- function(delG) {
  P1 <<- convert.delg.ka(delG) * c0
  ode.out <- ode(y = state, times = times, func = drug.model, parms = NULL)
  result <- accumulate.dMr(out = ode.out) %>% .[-1, ] %>%
    mutate(param = as.factor(delG))
  return(result)
}
# P2 <- 0.01
delg.variance <- do.call(rbind, lapply(seq(-5, -40, -5),
                                       FUN = vary.delG)) 
ggplot(delg.variance, aes(x = times, y = rr, color = param)) + 
  geom_line() +
  labs(title = "Effect of delG on release profiles", x = "Time, hr", y = "Cumulative release", 
       color = "delG, kJ/mol") + 
  theme_bw() 
ggsave("./release/delg.variance.png")
saveRDS(delg.variance, "./release/delg.variance.RDS")

ggplot(delg.variance, aes(x = times, y = rr, color = param)) + 
  geom_line(size = 1) +
  labs( x = "Time, hr", y = "Cumulative release", 
       color = "dG, kJ/mol") +
  coord_fixed(ratio = 12500) +
  theme.isef
ggsave("./graphs/2018 isef/dg and release.png", dpi = 450)

#     P3 ------------------------------------------------------------------

# Varying P3 is a lot less simple than P1 or P2 because the constituent 
# variables, ct and c0, exist in a balance with clc.eq, cl.eq

set.param()

# Let's keep ct, the total amount of CD, constant
# (and ignore the fact that this set of equations contradicts basic
# chemical relationships in the model)
# c0 is the total drug concentration
vary.c0 <- function(newc0) {
  P3 <<- ct/newc0
  ode.out <- ode(y = state, times = times, func = drug.model, parms = NULL)
  result <- accumulate.dMr(out = ode.out) %>% .[-1, ] %>%
    mutate(param = as.factor(newc0))
  return(result)
}

c0.variance <- do.call(rbind, lapply(seq(0.01, 0.51, 0.1),
                                     FUN = vary.c0)) 
ggplot(c0.variance, aes(x = times, y = rr, color = param)) + 
  geom_line() +
  labs(title = "Effect of c0, total drug concentration, on release profiles", 
       x = "Time, hr", y = "Cumulative release", color = "c0") + 
  theme_bw() 
ggsave("./release/c0.variance.png")

set.param()
vary.ct <- function(newct) {
  P3 <<- newct/c0
  ode.out <- ode(y = state, times = times, func = drug.model, parms = NULL)
  result <- accumulate.dMr(out = ode.out) %>% .[-1, ] %>%
    mutate(param = as.factor(newct))
  return(result)
}

ct.variance <- do.call(rbind, lapply(seq(0.01, 0.51, 0.1),
                                     FUN = vary.ct)) 
ggplot(ct.variance, aes(x = times, y = rr, color = param)) + 
  geom_line() +
  labs(title = "Effect of ct, total CD concentration, on release profiles", 
       x = "Time, hr", y = "Cumulative release", color = "ct") + 
  theme_bw() 
ggsave("./release/ct.variance.png")

#     Diffusivity ---------------------------------------------------------

# First, finding relationships between diffusivity and kd
# 1. Find the data with descriptors
padel <- readRDS("./descriptors/all.padel.RDS") 
data.mw <- padel %>% select(., guest, MW)
#     i. convert molecular weight to diffusivity
data.diff <- data.mw %>% mutate(diffusivity = get.diffusivity(MW)) %>%
  select(., guest, diffusivity)

# 2. Find the data with DelG
binding <- readRDS("./dwnld/combined.data.RDS")
#     i. convert delG to Ka
data.ka <- binding %>% mutate(ka = convert.delg.ka(DelG)) %>%
  select(., guest, host, ka)
data.dg <- binding %>% select(., guest, host, DelG)

# 3. split the ka data by host and join to the diffusivity data
#     i. alpha-CD
ka.alpha <- data.ka %>% filter(host == "alpha")
ka.diff.alpha <- inner_join(ka.alpha, data.diff, by = "guest")
ka.diff.alpha <- ka.diff.alpha[!duplicated(ka.diff.alpha), ]
#     ii. beta-CD
ka.beta <- data.ka %>% filter(host == "beta")
ka.diff.beta <- inner_join(ka.beta, data.diff, by = "guest")
ka.diff.beta <- ka.diff.beta[!duplicated(ka.diff.beta), ]
#     iii. gamma-CD
ka.gamma <- data.ka %>% filter(host == "gamma")
ka.diff.gamma <- inner_join(ka.gamma, data.diff, by = "guest")
ka.diff.gamma <- ka.diff.gamma[!duplicated(ka.diff.gamma), ]
#     iv. join everything
data.ka.diff <- rbind(ka.diff.alpha, ka.diff.beta, ka.diff.gamma)
#     v. repeat for only delg
dg.alpha <- data.dg %>% filter(host == "alpha")
dg.diff.alpha <- inner_join(dg.alpha, data.diff, by = "guest")
dg.diff.alpha <- dg.diff.alpha[!duplicated(dg.diff.alpha), ]
#     ii. beta-CD
dg.beta <- data.dg %>% filter(host == "beta")
dg.diff.beta <- inner_join(dg.beta, data.diff, by = "guest")
dg.diff.beta <- dg.diff.beta[!duplicated(dg.diff.beta), ]
#     iii. gamma-CD
dg.gamma <- data.dg %>% filter(host == "gamma")
dg.diff.gamma <- inner_join(dg.gamma, data.diff, by = "guest")
dg.diff.gamma <- dg.diff.gamma[!duplicated(dg.diff.gamma), ]

data.dg.diff <- rbind(dg.diff.alpha, dg.diff.beta, dg.diff.gamma)


#         Graphing diffusivity vs. dG -------------------------------------

# On an ordinary scale, trends aren'ts very clear due to the exponential 
# character of Ka
ggplot(data.ka.diff, aes(x = ka, y = diffusivity,
  color = host, shape = host)) +
  geom_point() +
  theme_bw()
ggplot(data.ka.diff, aes(x = ka, y = diffusivity,
                         color = host, shape = host)) +
  geom_point() +
  theme_bw() + 
  coord_cartesian(xlim = c(0,700)) + 
  facet_grid(~host)

# plotting ka on an exponential scale
ggplot(data.ka.diff, aes(x = ka, y = diffusivity, color = host, shape = host)) + 
  geom_point() + 
  scale_x_log10() + 
  theme_bw()

# Graphing with labels
ggplot(data.ka.diff, aes(x = ka, y = diffusivity,
                         color = host, shape = host)) +
  geom_point() +
  theme_bw() + 
  scale_x_log10() + 
  # facet_grid(~host) + 
  geom_text(data = subset(data.ka.diff, ka > 10000 & diffusivity < (2 * 10^-5)), 
            aes(ka, diffusivity, label = guest))

# Graphing with labels
ggplot(data.dg.diff, aes(x = DelG, y = diffusivity,
                         color = host, shape = host)) +
  geom_point() +
  theme_bw() 
  # facet_grid(~host) + 


#     i. divide into quadrants
# search for some extreme values (not the most extreme, though, b/c reliability)
# quadrant 1: high binding, high diffusivity
subset(data.dg.diff, DelG < -18 & diffusivity > 2.35e-5)
# benzaldehyde (211) for alpha, cyclohexanecarboxylic acid (299) for beta
# nothing available for gamma
q1 <- data.dg.diff[rownames(data.dg.diff) %in% c("211", "299"), ] %>%
  mutate(quadrant = as.factor(1))

# quadrant 2: low binding, high diff
subset(data.dg.diff, DelG > -5 & diffusivity > 3e-5)
#     acetonitrile (200) for alpha, ethanol (339) for beta
subset(data.dg.diff, DelG > -8 & diffusivity > 2.5e-5) 
#     benzene (121) for gamma
q2 <- data.dg.diff[rownames(data.dg.diff) %in% c("200", "339", "121"), ] %>%
  mutate(quadrant = as.factor(2))

# quadrant 3: low binding, low diff
subset(data.dg.diff, DelG > -9 & diffusivity < 2.1e-5)
# Tyr-Ile-Gly-Ser-Arg (197) for alpha, griseofulvin (354) for beta
# 2-naphthalenesulfonate (2312) for gamma
q3 <- data.dg.diff[rownames(data.dg.diff) %in% c("197", "354", "2312"), ] %>%
  mutate(quadrant = as.factor(3))

# quadrant 4: high binding, low diff
subset(data.dg.diff, DelG < -20 & diffusivity < 2e-5)
# 3-[(4-hydroxyphenyl)azo]benzoate (83) for alpha, 
# 2-chloro-4-[(4-hydroxyphenyl)azo]benzoate (123) for gamma
subset(data.dg.diff, DelG < -25 & diffusivity < 1.8e-5)
# spironolactone (440) for beta
q4 <- data.dg.diff[rownames(data.dg.diff) %in% c("83", "123", "440"), ] %>%
  mutate(quadrant = as.factor(4))

dg.diff.subset <- rbind(q1, q2, q3, q4)
ggplot(data.dg.diff, aes(x = DelG, y = diffusivity, color = host)) +
  geom_point() +
  theme_bw() + 
  geom_hline(yintercept = mean(data.dg.diff$diffusivity)) + 
  geom_vline(xintercept = mean(data.dg.diff$DelG)) + 
  facet_grid(~host) + 
  geom_text(data = dg.diff.subset, aes(x = DelG, y = diffusivity, label = guest),
            color = "#161616", vjust = 1) + 
  geom_point(data = dg.diff.subset, aes(x = DelG, y = diffusivity)) + 
  # theme.2018 + 
  labs(y = "Diffusivity", x = "dG, kJ/mol", color = "CD Type")

#         Release curves --------------------------------------------------

# modifying the set.param function to be more usable

# set.param.p <- function(p1, p2, p3) {
#   P1 <<- p1
#   P2 <<- p2
#   P3 <<- p3
# }

set.param.dg.diff <- function(dg, diff) {
  P1 <<- convert.delg.ka(dg) * c0
  P2 <<- diff/(k2*width^2)
}

vary.dg.diff <- function(dg, diff) {
  set.param.dg.diff(dg, diff)
  ode.out <- ode(y = state, times = times, func = drug.model, parms = NULL)
  result <- accumulate.dMr(out = ode.out) %>%
    mutate(DelG = as.factor(dg)) %>%
    mutate(diffusivity = as.factor(diff))
  return(result)
}
q1.dd <- do.call(rbind, Map(f = vary.dg.diff, 
                            dg = q1$DelG, 
                            diff = q1$diffusivity)) %>%
  mutate(quadrant = "q1")
q2.dd <- do.call(rbind, Map(f = vary.dg.diff, 
                            dg = q2$DelG, 
                            diff = q2$diffusivity)) %>%
  mutate(quadrant = "q2")
q3.dd <- do.call(rbind, Map(f = vary.dg.diff, 
                            dg = q3$DelG, 
                            diff = q3$diffusivity)) %>%
  mutate(quadrant = "q3")
q4.dd <- do.call(rbind, Map(f = vary.dg.diff, 
                            dg = q4$DelG, 
                            diff = q4$diffusivity)) %>%
  mutate(quadrant = "q4")

allq.dd <- rbind(q1.dd, q2.dd, q3.dd, q4.dd) 
ggplot(allq.dd, aes(x = times, y = rr, color = quadrant, group = DelG)) + 
  geom_line() +
  theme_bw() 

# For the heck of it, why don't we just increase the diffusivity by 
# a factor of 10. Or something. 
set.param.dg.diff <- function(dg, diff) {
  P1 <<- convert.delg.ka(dg) * c0
  P2 <<- diff/(k2*width^2) * 10 
}

q1.dd <- do.call(rbind, Map(f = vary.dg.diff, 
                            dg = q1$DelG, 
                            diff = q1$diffusivity)) %>%
  mutate(quadrant = "q1")
q2.dd <- do.call(rbind, Map(f = vary.dg.diff, 
                            dg = q2$DelG, 
                            diff = q2$diffusivity)) %>%
  mutate(quadrant = "q2")
q3.dd <- do.call(rbind, Map(f = vary.dg.diff, 
                            dg = q3$DelG, 
                            diff = q3$diffusivity)) %>%
  mutate(quadrant = "q3")
q4.dd <- do.call(rbind, Map(f = vary.dg.diff, 
                            dg = q4$DelG, 
                            diff = q4$diffusivity)) %>%
  mutate(quadrant = "q4")

allq.dd <- rbind(q1.dd, q2.dd, q3.dd, q4.dd) 
ggplot(allq.dd, aes(x = times, y = rr, color = quadrant, group = DelG)) + 
  geom_line(size = 0.75) +
  # theme.2018 + 
  labs(x = "Time, hr", y = "Cumulative release", color = "Quadrant")
# ggsave("./release/presentation graph.png")

#     Diffusivity (~constant dG) ------------------------------------------

subset(data.dg.diff, DelG < -12.5 & DelG > -12.7)
# 346, 317, 67, and 2381 make a nice range
same.dg <- data.dg.diff[rownames(data.dg.diff) %in% c("346", "317", "67", "2381"), ]
same.dg.out <- do.call(rbind, Map(f = vary.dg.diff, 
                            dg = -12.6, 
                            diff = same.dg$diffusivity)) 
ggplot(same.dg.out, aes(x = times, y = rr, color = diffusivity)) + 
  geom_line() + 
  theme_bw()

# Checking an area with weak affinity
same.dg.low <- subset(data.dg.diff, DelG <= -3.5 & DelG > -3.7)
same.dg.low.out <- do.call(rbind, Map(f = vary.dg.diff, 
                                  dg = -3.6, 
                                  diff = same.dg.low$diffusivity)) 
ggplot(same.dg.low.out, aes(x = times, y = rr, color = diffusivity)) + 
  geom_line() + 
  theme_bw()

# checking an area with high affinity
same.dg.high <- subset(data.dg.diff, DelG < -24.5 & DelG > -25.5)
same.dg.high.out <- do.call(rbind, Map(f = vary.dg.diff, 
                                      dg = -25, 
                                      diff = same.dg.high$diffusivity)) 
ggplot(same.dg.high.out, aes(x = times, y = rr, color = diffusivity)) + 
  geom_line() + 
  theme_bw()

# Checking for conservation of mass ---------------------------------------

vary.dg.diff <- function(dg, diff) {
  set.param.dg.diff(dg, diff)
  ode.out <- ode(y = state, times = times, func = drug.model, parms = NULL)
  result <- accumulate.dMr(out = ode.out) %>%
    mutate(DelG = as.factor(dg)) %>%
    mutate(diffusivity = as.factor(diff))
  return(result)
}

temp <- ode(y = state, times = times, func = drug.model, parms = NULL)

length(times)
drug.total <- rep(0, length(times))
comp.total <- rep(0, length(times))

for(i in 1:length(times)) {
  drug.total[i] <- sum(temp[i, 2:102])
}
for(i in 1:length(times)) {
  comp.total[i] <- sum(temp[i, 102:201])
}
drug.comp.total <- drug.total + comp.total
drug.total.time <- data.frame(times, drug.total, comp.total, drug.comp.total)
ggplot(drug.total.time, aes(x = times, y = drug.total)) + 
  geom_line(color = "blue") + 
  geom_line(aes(y = comp.total), color = "red") + 
  geom_line(aes(y = drug.comp.total)) + 
  theme_bw()



# ODE on external validation ----------------------------------------------


