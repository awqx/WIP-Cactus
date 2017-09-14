# source("http://bioconductor.org/biocLite.R")
# biocLite("ChemmineR") 
# biocLite("fmcsR") 
install.packages("fmcsR")
library(ChemmineR)
library(data.table)
library(fmcsR)
library(stringr)
library(tidyverse)
# setwd("~/SREP LAB/qsar")
dir.create("./molecules/compiledSDFs")

# Functions ---------------------------------------------------------------

# Read SDF, replace header with name of file
# precondition: guest is an SDF file path
name.sdf <- function(guest) {
  name <- str_extract(guest, "/{1}[^/]+\\.SDF$") %>% 
    str_replace(., "\\.SDF", "") %>%
    str_replace(., "/", "")
  sdf <- read.csv(guest, header = F)
  sdf[ , 1] <- sdf[ , 1] %>% as.character()
  sdf[1, ] <- name
  return(sdf)
}

# Compiling SDFs ----------------------------------------------------------

# Compiling guests into single SDF

alpha.files <- list.files("./molecules/alphaCD", full.names = T)
sdf.all.alpha <- do.call(rbind, lapply(alpha.files, name.sdf))
write.table(sdf.all.alpha, "./molecules/compiledSDFs/alpha.SDF", 
            row.names = F, col.names = F, quote = F)

sdf.all.beta <- do.call(rbind, lapply(list.files("./molecules/betaCD", full.names = T), 
                                      name.sdf))
write.table(sdf.all.beta, "./molecules/compiledSDFs/beta.SDF",
  row.names = F, col.names = F, quote = F)

sdf.all.gamma <- do.call(rbind, lapply(list.files("./molecules/gammaCD", full.names = T), 
                                      name.sdf))
write.table(sdf.all.gamma, "./molecules/compiledSDFs/gamma.SDF",
            row.names = F, col.names = F, quote = F)

# Vignette ----------------------------------------------------------------

# Getting SDFset classes
alpha.sdf <- read.SDFset("./molecules/compiledSDFs/alpha.SDF")
valid <- validSDF(alpha.sdf)
alpha.sdf <- alpha.sdf[valid]
plot(alpha.sdf[1:4], regenCoords = T, print = F)
# sdf.visualize(fda.sdf[1:4])

beta.sdf <- read.SDFset("./molecules/compiledSDFs/beta.SDF")
# plot(beta.sdf[1:4], regenCoords = T, print = F)

gamma.sdf <- read.SDFset("./molecules/compiledSDFs/gamma.SDF")
# plot(gamma.sdf[1:4], regenCoords = T, print = F)

# Getting APset classes (atom pairs)
alpha.apset <- sdf2ap(alpha.sdf)
beta.apset <- sdf2ap(beta.sdf)
gamma.apset <- sdf2ap(gamma.sdf)


# Sample Comparisons ------------------------------------------------------

fmcs(alpha.sdf[4], alpha.sdf[10]) 
cmp.similarity(alpha.apset[1], alpha.apset[2])
alpha.clusters <- cmp.cluster(db = alpha.apset, cutoff = c(07, 0.8, 0.9, quiet = T))

# Structure Plot
data(fmcstest)
test <- fmcs(fmcstest[1], fmcstest[2], au = 2, bu = 1)
plotMCS(test)

test <- fmcs(alpha.sdf[8], alpha.sdf[64], au = 2, bu = 1)
plotMCS(test)
