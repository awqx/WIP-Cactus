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

k1 <- .1 
k2 <- .5
cd <- 1
D <- 1

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
    change[i] <- -D * 0.001^2 * pi * (time.results[[i]]$drug[20] - time.results[[i]]$drug[19]) # 20 and 19 should be N and N-1

  return(change)
  }

temp2 <- dMr.dt(DRUG, COMP, 10)
dMr.dt2(temp2)


# ODE ---------------------------------------------------------------------
library(deSolve)
# deSolve Paper
n <- 100
del <- 1
cl.eq <- 0.3
clc.eq <- 0.5
ct <- cl.eq + clc.eq
c0 <- clc.eq + 0.2
D <- 600
k1 <- 36*20
k2 <- 36

P1<- p1 <- k1/k2*c0
P2 <- p2 <- D/(k2*del^2)
P3 <- p3 <- ct/c0

parameters <- c(P1 = p1, 
                P2 = p2, 
                P3 = p3)
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

t <- 100
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

# ReacTran PDE ------------------------------------------------------------

# This probably isn't going to work
# ReacTran is difficult because it can't handle >1 variable at a time
# and only comes with examples of constant consumption

# Also, I'm pretty sure reactran doesn't give an output over time

#     Example: hyperbolic wave --------------------------------------------

library(ReacTran)
xgrid <- setup.grid.1D(0, 1, dx.1 = 0.01)
x <- xgrid$x.mid
N <- xgrid$N

drug <- rep(0.05, N) # initial concentration of drug
comp <- rep(0.75, N) # initial concentratino of complexes
c0 <- 0.05+0.75 # total drug loading
ct <- 0.85 # tital CD
yini <- c(drug, comp)
times <- seq(from = 0, to = 10, by = 1)

D <- 30 # Diffusion coef
D.grid<- setup.prop.1D(value = D, grid = grid)

k1 <- 36.9*40
k2 <- 36.9

drug.model <- function(t, y, parms) {
  u1 <- y[1:N]
  u2 <- y[-(1:N)]
  
  du2 <- k1/k2 * c0 * u1 * (ct/c0 - u2) - u2
  du1 <- tran.1D(C = u1, C.up = 0, C.down = 0, D = 30, dx = xgrid)$dC
  du1 <- du1 - du2
  return(list(c(du1, du2)))
}

# Warning message
# In rk(y, times, func, parms, method = "ode45", ...) :
# Number of time steps 22915 exceeded maxsteps at t = 1.89484
out <- ode.1D(func = drug.model, y = yini, times = times, parms = NULL, 
              nspec = 2, method = "ode45", dimens = N, names = c("drug", "comp"))

# release of drug
drug.model <- function(t = 0, drug, pars = NULL) {
  tran <- tran.1D(C = drug, D = D.grid, dx = grid)$dC
  reac <- -Rb
}

comp.model <- function(t = 0, comp, pars = NULL) {
  tran <- tran.1D(C)
}



#     O2 consumption in spherical aggregate -------------------------------

Aggregate.Model <- function(time, O2, pars) {
  tran <- tran.1D(C = O2, C.down = C.ow.O2, 
                  D = D.grid, A = A.grid, 
                  VF = por.grid, dx = grid)$dC
  reac <- -R.O2 * (O2/(Ks + O2)) * (O2 > 0)
  return(list(dCdt = tran + reac, consumption = -reac))
}

C.ow.O2 <- 0.25     # concentration O2 water [micromol cm-3]
por     <- 0.8      # porosity
D       <- 400      # diffusion coefficient O2 [cm2 yr-1]
v       <- 0        # advective velocity [cm yr-1]
R.O2    <- 1000000  # O2 consumption rate [micromol cm-3 yr-1]
Ks      <- 0.005    # O2 saturation constant [micromol cm-3]

# Grid definition
R <- 0.025           # radius of the agggregate [cm]
N <- 100             # number of grid layers
grid <- setup.grid.1D(x.up = 0, L = R, N = N)

# Volume fractions
por.grid <- setup.prop.1D(value = por, grid = grid)
D.grid   <- setup.prop.1D(value = D, grid = grid)
# Surfaces
A.mid <- 4*pi*grid$x.mid^2  # surface of sphere at middle of grid cells
A.int <- 4*pi*grid$x.int^2  # surface of sphere at interface
A.grid <- list(int = A.int, mid = A.mid)

O2.agg <- steady.1D(runif(N), func = Aggregate.Model, nspec = 1, 
                    atol = 1e-10, names = "O2")

par(mfrow = c(1,1))
plot(grid$x.mid, O2.agg$y, xlab = "distance", ylab = "mmol/m3", 
     main = "Diffusion-reaction")
