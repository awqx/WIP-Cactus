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
