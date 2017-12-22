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
