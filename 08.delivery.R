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

n <- 100 # number of layers
del <- 1
del <- 0.01 # I don't know why this value works 
del <- 1/n # width of a layer
# del <- 1
# width <- 1
width <- 0.05 # width of hydrogel, in cm
cl.eq <- 0.3 # conc of drug at equilibrium
clc.eq <- 0.5 # conc of complexes @ equilibrium
c0 <- cl.eq + clc.eq # total drug concentration
ct <- clc.eq + 0.6 # total cd conc
D <- 0.08 # diffusivity coef
k1 <- 36*20 # rate of binding
k2 <- 36 # rate of unbinding

P1 <- k1/k2*c0
P2 <- D/(k2*width^2)
P2 <- 0.9 # according to Fu
P3 <- ct/c0

# parameters <- c(P1 = p1, 
#                 P2 = p2, 
#                 P3 = p3)
state <- c(rep(cl.eq, n), 
           rep(clc.eq, n))

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

t <- 500
times <- seq(0, t, 1) 
out <- ode(y = state, times = times, func = drug.model, parms = NULL)

drug.out <- out[ , c(1:n)] %>% as.data.frame()
drug.out <- data.table::melt(setDT(drug.out), id.vars = "time", 
                             variable.name = "location")

ggplot(drug.out, aes(x = time, y = location, fill = value)) + 
  geom_tile()

dmr.raw <- get.dMr(out = out) %>% .[-1, ]
ggplot(dmr.raw, aes(x = times, y = rr)) + 
  geom_line()
dmr.sum <- accumulate.dMr(out = out) %>% .[-1, ]
ggplot(dmr.sum, aes(x = times, y = rr)) + 
  geom_line()

