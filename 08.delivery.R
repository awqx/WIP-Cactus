# Integrating DelG results from the QSARs with Fu's mathematical modeling of
# complexation dynamics. Meant to produce estimates of cumulative drug release
# over time
#
# Source: https://doi.org/10.1007/s10439-011-0336-z

# Libraries and Packages --------------------------------------------------

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
# ODE ---------------------------------------------------------------------

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


# Parameters --------------------------------------------------------------

# Known
ka <- convert.delg.ka(-15) # choose a value from the QSPR
ct <- 0.22 # at 1 g/4 mL (Fu)
D <- 0.1 # technically searchable, but this is arbitrary
width <- 0.05 # width of half of hydrogel
vol <- 0.0785 # cm^2, vol when swollen

# Reasonably derivable
k2 <- 35 # h^-1, according to Fu
k1 <- k2*ka # ka = k1/k2 

# Educated guesses
c0 <- 0.15 # total drug concentration
cc <- 0.08 # concentration of uncomplexed cd
clc.eq <- ct - cc
# # Checking for approximate accuracy
# clc.eq
# c0 - c0/(1 + ka*cc)

# Dimensionless
P1 <- k1/k2*c0
P2 <- D/(k2*width^2)
P2 <- 0.8 # approximation from Fu
P3 <- ct/c0


# ODE Implementation ------------------------------------------------------

n <- 100 # number of layers
del <- 1/n # width of a layer
del <- 0.1 # I don't know why this value works 
state <- c(rep(cl.eq, n), 
           rep(clc.eq, n))
t <- 200
times <- seq(0, t, 1) 
out <- ode(y = state, times = times, func = drug.model, parms = NULL)

# drug.out <- out[ , c(1:n)] %>% as.data.frame()
# drug.out <- data.table::melt(setDT(drug.out), id.vars = "time", 
#                              variable.name = "location")
# 
# ggplot(drug.out, aes(x = time, y = location, fill = value)) + 
#   geom_tile()
# 
# dmr.raw <- get.dMr(out = out) %>% .[-1, ]
# ggplot(dmr.raw, aes(x = times, y = rr)) + 
#   geom_line()
dmr.sum <- accumulate.dMr(out = out) %>% .[-1, ]
ggplot(dmr.sum, aes(x = times, y = rr)) + 
  geom_line() + 
  theme_bw() + 
  coord_cartesian(ylim = c(0,2))

# Testing values ----------------------------------------------------------

#     P1 ------------------------------------------------------------------

# Provides a cumulative release dataframe for a new value of P1
vary.P1 <- function(newP1) {
  P1 <<- newP1
  ode.out <- ode(y = state, times = times, func = drug.model, parms = NULL)
  result <- accumulate.dMr(out = ode.out) %>% .[-1, ] %>%
    mutate(param = paste0("P1 = ", newP1))
  return(result)
}

P1.variance <- do.call(rbind, lapply(seq(20, 320, 50),
                                     FUN = vary.P1))

ggplot(P1.variance, aes(x = times, y = rr, color = param)) + 
  geom_line() + 
  theme_bw() 

