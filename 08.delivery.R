# Integrating DelG results from the QSARs with Vulic's mathematical modeling of
# complexation dynamics Meant to produce estimates of cumulative drug release
# over time
#
# Source: https://doi.org/10.1016/j.jconrel.2014.10.032

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
 
# Functions (Fu, w/ dimensions) -------------------------------------------

# The ODE/PDE system is
# 
# dDRUG/dt = D*d2DRUG/d2z - Rb
# dCOMP/dt = Rb
# Rb = Kon*DRUG*CD-KoffCOMP = Kon*COMP*(TOTAL-COMP)-KoffCOMP
# 
# Boundary condition
# z = delta, DRUG = 0 
# 
# Initial conditions
# DRUG = 0; COMP = 0

k1 <- 0.1 
k2 <- 25
cd <- 1
D <- 10

get.Rb <- function(drug, comp, loc) {
  rb <- k1 * drug[loc] * (cd - comp[loc]) - k2*comp[loc]
}

ddrug.dt <- function(drug, comp) {
  ddrug <- c(rep(0, n)) # create a blank vector w/ correct number of layers
  # ddrug[0] (imaginary) = 0, ddrug[n] = 0 ... these are boundary conditions
  ddrug[1] <- D * (drug[2] - drug[1]) / (del^2) - get.Rb(drug = drug, 
                                                         comp = comp, 
                                                         loc = 1)
  ddrug[n - 1] <-
    D*(-2*drug[n-1]+drug[n-2])/(del^2) - get.Rb(drug = drug,
                                                comp = comp,
                                                loc = n - 1)
  for (i in 2:(n-2)) {
    ddrug[i] <- D*(drug[i+1]-2*drug[i]+drug[i-1])/(del^2)-get.Rb(drug = drug, 
                                                                 comp = comp,
                                                                 loc = i)
  }
  
  return(ddrug)
}

dcomp.dt <- function(drug, comp) {
  dcomp <- sapply(c(1:(n-1)), get.Rb, drug = drug, comp = comp)
  dcomp[n] <- 0
  return(dcomp)
}

dMr.dt <- function(drug, comp, time) {
  time.results <- list()
  initial <- list(drug, comp, 
                  ddrug.dt(drug = drug, comp = comp), 
                  dcomp.dt(drug = drug, comp = comp)) %>%
    setNames(., c("drug", "comp", "ddrug", "dcomp"))
  time.results[[1]] <- initial
  print(time.results)
  for (t in 2:time) 
  {
    # calculates new concentration of drug based on previous concentration
    # plus the calculated change
    drug.new <- time.results[[t-1]][[1]] + time.results[[t-1]][[3]]
    comp.new <- time.results[[t-1]][[2]] + time.results[[t-1]][[4]]
    ddrug.new <- ddrug.dt(drug.new, comp.new)
    dcomp.new <- dcomp.dt(drug.new, comp.new)
    list.new <- list(drug.new, comp.new, ddrug.new, dcomp.new) %>%
      setNames(., c("drug", "comp", "ddrug", "dcomp"))
    time.results[[t]] <- list.new
  }
  setNames(time.results, c(1:time))
  return(time.results)
}

dMr.dt2 <- function(time.results) {
  time <- length(time.results)
  change <- c(rep(0, time))
  
  for(i in 1:time)
    change[i] <- -D * 0.0005^2 * pi * (time.results[[i]]$drug[20] - time.results[[i]]$drug[19]) # 20 and 19 should be N and N-1

  return(change)
  }

temp2 <- dMr.dt(DRUG, COMP, 10)
dMr.dt2(temp2)
