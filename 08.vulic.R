# Integrating DelG results from the QSARs with Vulic's mathematical modeling of
# complexation dynamics Meant to produce estimates of cumulative drug release
# over time
#
# Source: https://doi.org/10.1016/j.jconrel.2014.10.032
# setwd("~/SREP LAB/qsar")

# Libraries and Packages --------------------------------------------------

install.packages("ReacTran")
library(ReacTran)
library(tidyverse)

# Functions ---------------------------------------------------------------

# Well, if we know that DelG = RTlnKd, then we should invert that
convert.delg.kd <- function(delg) {
  # exp(n) = e^n
  return (exp(delg / (8.314 * 298)))
}

# These are dimensionless constants used by Vulic:
#   alpha = Cpro,0 / Cpep,T (guest / total host (peptide))
#   beta = L^2 * koff / D 
#   gamma = kon * Cpro,0 / koff

# Regime 1 calculations
#   for cases when gamma << 1 or ~ 1
regime1 <- function(gamma, alpha, L, D) {
  return(L^2*(1 + gamma/alpha)/D)
}

# Regime 2 calculations
#   for gamma >> 1; or, a large fraction of peptide initially bound
#   This one is more complicated b/c the timeframe/calculations actually change
#   as the concentration gradually drops, and the time frame calculation shifts
#   to regime 1

# Regime 3 (slow unbinding)
#   for cases where unbinding is slow compared to diffusion, or beta << 1 or 1/koff >> L^2/D

# Preliminary Trial -------------------------------------------------------

glm <- readRDS("./models/glmnet/glm.results.RDS")
glm.delg <- glm$pred
glm.kd <- map(glm.delg, convert.delg.kd)

# ReacTran Trial ----------------------------------------------------------

# Grid definition
r.N     <- 4   # number of cells in r-direction
theta.N <- 6   # number of cells in theta-direction
z.N     <- 3   # number of cells in z-direction

D       <- 100 # diffusion coefficient
# D can be calculated as = (9.940E-15*T)/
# (viscosity*cuberoot(molecular weight)) [Vulic]
# In this case, T = 298, viscosity = 8.90E-4 Pa 

r      <- seq(0,   8, len = r.N+1)       # cell size r-direction [cm]
theta  <- seq(0,2*pi, len = theta.N+1)   # theta-direction - theta: from 0, 2pi
phi    <- seq(0,2*pi, len = z.N+1)       # phi-direction (0,2pi)
z      <- seq(0,5, len = z.N+1)          # cell size z-direction [cm]

# Intial conditions 
C <- array(dim = c(r.N, theta.N, z.N), data = 0)

# Concentration boundary conditions
tran.cylindrical (C = C, D.r = D, D.theta = D, 
                  C.r.up = 1, C.r.down = 1,
                  C.theta.up = 1, C.theta.down = 1, 
                  C.z.up = 1, C.z.down = 1,
                  r = r, theta = theta, z = z )

# Flux boundary conditions
tran.cylindrical(C = C, D.r = D, r = r, theta = theta, z = z,
                 flux.r.up = 10, flux.r.down = 10,
                 flux.theta.up = 10, flux.theta.down = 10,
                 flux.z.up = 10, flux.z.down = 10)

# cyclic boundary conditions
tran.cylindrical(C = C, D.r = D, r = r, theta = theta, z = z,
                 cyclicBnd = 1:3)

# zero-gradient boundary conditions
tran.cylindrical(C = C, D.r = D, r = r, theta = theta, z = z)


# Model w/ Diffusion and First-Order Consumption
N     <- 50          # number of grid cells
XX    <- 4           # total size
rr    <- 0.005       # consumption rate
ini   <- 1           # initial value at x=0
D     <- 400
r     <- seq (2, 4, len = N+1)
theta   <- seq(0, 2*pi, len = N+1)
theta.m <- 0.5*(theta[-1]+theta[-(N+1)])
# The model equations
Diffpolar <- function (t, y, parms)  {
  CONC  <- matrix(nrow = N, ncol = N, data = y)
  tran  <- tran.polar(CONC, D.r = D, D.theta = D, r = r, theta = theta,
                      C.r.up = 0, C.r.down = 1*sin(5*theta.m),
                      cyclicBnd = 2, full.output=TRUE )
  dCONC <- tran$dC  - rr * CONC
  return (list(dCONC))
}
# solve to steady-state; cyclicBnd = 2, because of C.theta.up, C.theta.down
out <- steady.2D (y = rep(0, N*N), func = Diffpolar, parms = NULL,
                  dim = c(N, N), lrw = 1e6, cyclicBnd = 2)
image(out)
cart <- polar2cart(out, r = r, theta = theta,
                   x = seq(-4, 4, len = 100),
                   y = seq(-4, 4, len = 100))
image(cart)

# Messing Around
diff1 <- tran.cylindrical (C = C, D.r = D, D.theta = D, 
                  C.r.up = 1, C.r.down = 1,
                  C.theta.up = 1, C.theta.down = 1, 
                  C.z.up = 1, C.z.down = 1,
                  r = r, theta = theta, z = z )
diff1.1 <- diff1[1] %>% unlist()
diff1.1 <- array(dim = c(r.N, theta.N, z.N), data = diff1.1)
diff2 <- tran.cylindrical (C = diff1.1, D.r = D, D.theta = D, 
                           C.r.up = 1, C.r.down = 1,
                           C.theta.up = 1, C.theta.down = 1, 
                           C.z.up = 1, C.z.down = 1,
                           r = r, theta = theta, z = z )

# deSolve Trial -----------------------------------------------------------

# Following the Soetart documentation
# library(deSolve)
Aphid <- function(t, APHIDS, parameters) {
  deltax <- c(0.5, rep(1, numboxes -1), 0.5)
  Flux <- -D*diff(c(0,APHIDS, 0)) / deltax
  dAPHIDS <- -diff(Flux) / delx + APHIDS * r
  
  # the return value
  list(dAPHIDS)
}

D <- 0.3 # m2/day diffusion rate
r <- 0.01 # /day net growth rate
delx <- 1 # thickness of boxes
numboxes <- 60 
# distance of boxes on plant, m, 1 m intervals)
Distance <- seq(from = 0.5, by = delx, length.out = numboxes)

# Initial conditions # ind/m2
APHIDS <- rep(0, times = numboxes)
APHIDS[30:31] <- 1
state <- c(APHIDS = APHIDS) # initialize state variables

# the model is run for 200 days, producing output every day
times <- seq(0,200, by = 1)
out <- ode.1D(state, times, Aphid, parms = 0, nspec = 1, names = "Aphid")
image(out, method = "filled.contour", grid = Distance, 
      xlab = "Time, days", ylab = "Distance on plant, m", 
      main = "Aphid density")

data <- cbind(dist = c(0, 10, 20, 30, 40, 50, 60), 
              Aphid = c(0, 0.1, 0.25, 0.5, 0.25, 0.1, 0))
par(mfrow = c(1,2))
matplot.1D(out, grid = Distance, type = "l", mfrow = NULL, 
           subset = time %in% seq(0, 200, by = 10), 
           obs = data, obspar = list(pch = 18, cex = 2, col = "red"))
plot.1D(out, grid = Distance, type = "l", mfrow = NULL, subset = time == 100, 
        obs = data, obspar = list(pch = 18, cex = 2, col = "red"))

# deSolve and Fu, von Recum -----------------------------------------------

# Let's start by trying to model concentratino of free ligand in 
# the hydrogel
# total concentration of Beta-CD = Ct; free CD = Cc; bound CD = Clc
# Binding rate = Rb
# height of cyclodextrin 0<z<smalldelta
# Rb = Kon*Cl*Cc-KoffClc = Kon*Clc*(Ct-Clc)-KoffClc
# partialCl/partialt = D*secondpartialCl/secondpartialz^2-Rb
# partialClc/partialt = Rb

Hydrogel <- function(t, CONC, parameters){
  deltac <- c(0.5, rep(1, numboxes - 1), 0.5)
  Flux <- D*diff(c(0, CONC, 0)) / deltac
  Rb <- k1*CONC*(cd.total - COMP) - k2*COMP
  dCOMP <- Rb
  dCONC <- diff(Flux) / delx - Rb
  
  list(dCONC)
}
# Model parameters
D <- 0.2 # Diffusion rate
k1 <- 1.01218 # kon, or koff^-1
k2 <- 0.98796 # koff, or kon^-1
cd.total <- 0.4
delx <- 1 # thickness of boxes
numboxes <- 50
Distance <- seq(from = 0.5, by = delx, length.out = numboxes) # for image(out)

# Initial conditions
CONC <- rep(0.5, # concentration  of ligand/drug at equilibrium
            times = numboxes)
COMP <- rep(0.5, # concentration of complexes @ equilibrium
            times = numboxes)
state <- c(CONC = CONC)

times <- seq(0, 20, by = 1)
out <- ode.1D(state, times, Hydrogel, parms = 0, nspec = 1, 
              names = "Hydrogel")
image(out, grid = Distance)

hydrogel.mod <- function(t, state, parms, N, rr, ri, dr, dri) {
  with (as.list(parms), 
        {
          DRUG <- state[1:N] # ligand/drug concentration
          CD <- state[(N+1):(2*N)] # complexed CD
          
          # Fluxes due to diffusion
          # Only a flux for drug b/c CD remains stationary
          FluxDrug <- D * diff(c(DRUG[1], DRUG, DRUG[N])) / dri
          
          # Rb = net binding rate
          Rb <- k1 * DRUG * (cd.total - CD) - k2 * CD
          
          dDRUG <- diff(ri * FluxDrug)/rr/dr - Rb
          dCD <- Rb
          
          return(list(c(dDRUG, dCD)))
        })
}

R <- 20 # radius of surface
N <- 100 # number of layers
dr <- R/N # thickness of each layer
r <- seq(dr/2, by = dr, len = N) # distance of center to mid-layer
ri <- seq (0, by = dr, len = N+1) # distance to layer interface
dri <- dr # dispersion distances
parms <- c(D = 0.05, k1 = 0.5, k2 = 2, cd.total = 1)

# Distance <- seq(from = 0.5, by = delx, length.out = numboxes)

state <- rep(0, 2 * N)
state[1] <- state[N+1] <- 0.1

times <- seq(0, 10, by = 1)

out <- ode.1D(y = state, times = times, func = hydrogel.mod, parms = parms, nspec = 2, names = c("DRUG", "CD"), 
              N = N, rr = r, ri = ri, dr = dr, dri = dri)

DRUG  <- out[, 2:9]

filled.contour(x = times, y = c(1:8), DRUG, color = topo.colors,
               xlab = "time, days", ylab = "Distance",
               main = "Drug concentration")
