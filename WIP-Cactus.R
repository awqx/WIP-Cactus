install.packages("stringr")
library("stringr")
install.packages("dplyr")
library("dplyr")
install.packages("Amelia")
library(Amelia)

dataset <-
  readRDS(file = "~/SREP LAB/Rekharsky and Inoue/Cactus/RI.rds")
# Fixing chemical names that affect URL at download-cactus stage
# Typos
dataset[48, "guest"] <- "biebrich scarlet"
dataset[1081, "guest"] <- "4-[(4-hydroxyphenyl)azo]benzoate"
dataset[23, "guest"] <- "3-(aminomethyl)-proxyl"
dataset[82, "guest"] <- "3-carbamoyl-proxyl"
# Cactus uses alternate identifying name for the following:
dataset[42, "guest"] <- "1-o-benzyl-rac-glycerol"
dataset[126, "guest"] <-
  "cis-diammine(1,1-cyclobutanedicarboxylato)platinum(II)"
dataset[c(142, 143, 727, 728), "guest"] <-
  "[(1R,2S)-1-hydroxy-1-phenylpropan-2-yl]-methylazanium"
dataset[c(144, 145, 729, 730), "guest"] <-
  "[(1S,2R)-1-hydroxy-1-phenylpropan-2-yl]-methylazanium"

source(file = "~/SREP LAB/Rekharsky and Inoue/Cactus/read_cactus.R")
folder  <- "~/SREP LAB/Rekharsky and Inoue/Cactus/"

# Unicode detects the characters alpha and beta
# Alpha == \u03b1, beta == \u03b2
alpha.guest <- unique(dataset$guest[dataset$host == "1\u03b1"])

beta.guest  <- unique(dataset$guest[dataset$host == "1\u03b2"])

gamma.guest <- unique(dataset$guest[dataset$host == "1γ"])

# Cactus is a chemical identifier resolver and can be used
# to find and download SD files used by PyRx
# Solution to function inspired by hfty's answer to
# http://stackoverflow.com/questions/33423725/
# using-trycatch-to-populate-a-data-frame-inside-a-loop-nicely

download.cactus.results <- function(guest, path, chemical.format) {
  report <- tryCatch({
    destfile       <- paste0(path, "/", guest, ".SDF")
    # Chemical format must be parsed to match all the outputs from NCI cactus
    Guest.URL      <- unlist(lapply(guest, URLencode, reserved = T))
    URL            <- paste0(
      "https://cactus.nci.nih.gov/chemical/structure/",
      Guest.URL,
      "/",
      chemical.format
    )
    Map(download.file, url = URL, destfile = destfile)
    data.frame(
      guest = guest,
      downloaded = "yes",
      warning = "no",
      error = "no"
    )
  },
  warning = function(warn) {
    message("Warning: either URL error or already existing directory.")
    destfile       <- paste0(path, "/", guest, ".SDF")
    Guest.URL      <- unlist(lapply(guest, URLencode, reserved = T))
    URL            <- paste0(
      "https://cactus.nci.nih.gov/chemical/structure/",
      Guest.URL, "/", chemical.format
    )
    Map(download.file, url = URL, destfile = destfile)
    data.frame(
      guest = guest,
      downloaded = "yes",
      warning = "yes",
      error = "no"
    )
  },
  error = function(err) {
    message("An error occurred")
    data.frame(
      guest = guest,
      downloaded = "no",
      warning = "yes",
      error = "yes"
    )
    
  },
  finally = {
    message("Chemical name processed")
  })
  return(report)
}

# Find files that are too small to be SDF
filter_filesize <- function(path, pattern, size, logical) {
  filenames <- list.files(path = path, pattern = pattern)
  location  <- dir(path = path,
                   full.names = T,
                   recursive = T)
  data      <- file.info(paste0(path, "/", filenames))[1]
  index     <- eval(call(logical, data, size))
  return(data.frame(
    filepath = location[index],
    molecule = filenames[index],
    stringsAsFactors = F
  ))
}

# Create a directory to store files

create.host.dir <- function(path, host) {
  host.directory <- paste0(path, host)
  dir.create(path = host.directory)
  return(host.directory)
}
# Read SDF files and assign a chemical name
fix.sdf <- function(guest, path){
  report <- tryCatch({
    file <- paste0(path, "/", guest, ".SDF")
    raw.sdf <- read.table(file = file, header = F, sep = "\t", stringsAsFactors = F)
    raw.sdf[1, ] <- paste0(raw.sdf[1,], " ", guest)
    return(raw.sdf)
  }, 
  warning = function(warn){
    message("A warning occurred:")
    message(last.warning)
    file <- paste0(path, "/", guest, ".SDF")
    raw.sdf <- read.table(file = file, header = F, sep = "\t", stringsAsFactors = F)
    raw.sdf[1, ] <- paste0(raw.sdf[1,], " ", guest)
    return(raw.sdf)
  }, 
  error = function(err){
    message("An error occurred")
  })
  return(report)
}

#--------------------------------------
#          Download Gamma
#--------------------------------------

gamma.dest <- create.host.dir(folder, "GammaCD")

results.gamma <-
  do.call(
    rbind,
    lapply(
      gamma.guest,
      download.cactus.results,
      path = gamma.dest,
      chemical.format = "SDF"
    )
  )


results.all <- create.host.dir(folder, "Download Results")

saveRDS(results.gamma, file = paste0(results.all, "/results.gamma.RDS"))

success.gamma <- results.gamma[results.gamma$downloaded == "yes", ]

fail.gamma <- results.gamma[results.gamma$downloaded == "no", ]

empty.gamma <- filter_filesize(
  path = gamma.dest,
  pattern = "SDF",
  size = 100,
  logical = "<"
)

saveRDS(empty.gamma, file = paste0(results.all, "/empty_gamma_sdfiles.RDS"))

file.remove(empty.gamma$filepath)

# This step adjusts names of molecules that failed to download
# Into a form more likely to be used by Cactus
empty.gamma$molecule <- str_replace(string      = empty.gamma$molecule,
                                    pattern     = "\\α",
                                    replacement = "alpha")

# Still don't know what to do with these
results.empty.gamma <-
  do.call(
    rbind,
    lapply(
      empty.gamma$molecule,
      download.cactus.results,
      path = gamma.dest,
      chemical.format = "SDF"
    )
  )

gamma.comb <- rbind(results.empty.gamma, results.gamma)

#--------------------------------------
#         Downloading Beta
#--------------------------------------
# Expect this step to take a while, depending on speed of internet and computer
beta.dest <- create.host.dir(folder, "BetaCD")

results.beta <-
  do.call(
    rbind,
    lapply(
      beta.guest,
      download.cactus.results,
      path = beta.dest,
      chemical.format = "SDF"
    )
  )

empty.beta <-
  filter_filesize(
    path = beta.dest,
    pattern = "SDF",
    size = 100,
    logical = "<"
  )
success.beta <- results.beta[results.beta$downloaded == "yes", ]
file.remove(empty.beta$filepath)

saveRDS(results.beta, file = paste0(beta.dest, "/results.beta.RDS"))

saveRDS(empty.beta, file = paste0(beta.dest, "/empty_beta_sdfiles.RDS"))

# Trying to modify failed downloads to work properly
empty.beta$molecule <-
  str_replace(string = empty.beta$molecule,
              pattern = "\\β",
              replacement = "beta")

empty.beta$molecule <-
  str_replace(string = empty.beta$molecule,
              pattern = "\\α",
              replacement = "alpha")

empty.beta$molecule <-
  str_replace(
    string = empty.beta$molecule,
    pattern =  "\\(\\-*±*[0-9A-Z]*,*
    [0-9A-Z]*\\)\\-\\(*\\+*\\-*\\)*\\-*|nor(?!t)|
    \\([a-z]*\\,*\\s*[a-z]*(I[0-9]\\-)*\\)"
    ,
    replacement = ""
    )

beta.dest.2 <- create.host.dir(paste0(beta.dest), "/2nd.Try")

results.beta.2 <- do.call(
  rbind,
  lapply(
    empty.beta$molecule,
    download.cactus.results,
    path = beta.dest.2,
    chemical.format = "SDF"
  )
)
# Nothing changed. Maybe remove this step?
#--------------------------------------
#          Download Alpha
#--------------------------------------
# Sometimes the code generates a table of just errors, even
# if the moleules download correctly
alpha.dest <- create.host.dir(folder, "AlphaCD")

results.alpha <-
  do.call(
    rbind,
    lapply(
      alpha.guest,
      download.cactus.results,
      path = alpha.dest,
      chemical.format = "SDF"
    )
  )

success.alpha <- results.alpha[results.alpha$downloaded == "yes", ]

fail.alpha <- results.alpha[results.alpha$downloaded == "no", ]

empty.alpha <- filter_filesize(
  path = alpha.dest,
  pattern = "SDF",
  size = 100,
  logical = "<"
)

file.remove(empty.alpha$filepath)

saveRDS(results.alpha, file = paste0(alpha.dest, "/results.alpha.RDS"))

saveRDS(empty.alpha, file = paste0(alpha.dest, "/empty_alpha_sdfiles.RDS"))

empty.alpha$molecule <-
  str_replace(string = empty.alpha$molecule,
              pattern =  "\\(\\-*\\±*[0-9A-Z]*\\,*[0-9A-Z]*\\)\\-\\(*\\+*\\-*\\)*\\-*|nor(?!t)|\\([a-z]*\\,*\\s*[a-z]*(I[0-9]\\-)*\\)",
              replacement = "")

empty.alpha$molecule <-
  str_replace(string = empty.alpha$molecule,
              pattern = "\\β",
              replacement = "beta")

empty.alpha$molecule <-
  str_replace(string = empty.alpha$molecule,
              pattern = "\\α",
              replacement = "alpha")

download_cactus(
  guest = empty.alpha$molecule,
  host = "AlphaCD",
  path = folder,
  chemical.format = "SDF"
)

alpha.dest.2 <- create.host.dir(paste0(alpha.dest), "/2nd.Try")

results.alpha.2 <- do.call(
  rbind,
  lapply(
    empty.alpha$molecule,
    download.cactus.results,
    path = alpha.dest.2,
    chemical.format = "SDF"
  )
)

# Combining all the SDF files into one file
# For Gamma
all.sdf <- create.host.dir(folder, "All.SDF")

gamma.success <- as.character(success.gamma$guest)
all.gamma <- do.call(bind_rows, lapply(gamma.success, fix.sdf, path = gamma.dest))
write.table(
  all.gamma,
  file = paste0(all.sdf, "/allgamma.SDF"),
  append = TRUE,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

# For Beta
beta.success <- as.character(success.beta$guest)
e.dansyl <- fix.sdf("e-dansyl-l-lysine", beta.dest)
d.valero <- fix.sdf("d-valerolactam", beta.dest)
all.beta <- do.call(bind_rows, lapply(beta.success, fix.sdf, path = beta.dest))
all.beta <- rbind(e.dansyl, d.valero, all.beta)
write.table(
  all.beta,
  file = paste0(all.sdf, "/allbeta.SDF"),
  append = TRUE,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)
# For Alpha
alpha.success <- as.character(success.alpha$guest)
a.methyl1 <- fix.sdf("(+)-a-methylbenzylamine", alpha.dest)
a.methyl2 <- fix.sdf("(-)-a-methylbenzylamine", alpha.dest)
d.valero.a <- fix.sdf("d-valerolactam", alpha.dest)
all.alpha <- do.call(bind_rows, lapply(alpha.success, fix.sdf, path = alpha.dest))
all.alpha <- rbind(a.methyl1, a.methyl2, d.valero.a, all.alpha)
write.table(
  all.alpha,
  file = paste0(all.sdf, "/allalpha.SDF"),
  append = TRUE,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

# Creating a small sample
sample.aff.dir <- create.host.dir(folder, "Affinity/Sample")

common.chem.names <- Reduce(intersect, list(alpha.names, beta.names, gamma.names))
set.seed(4)
sample.chem.list <- sample(common.chem.names, 3)
# Sample should include cinnarizine, (E)-stilbene, and 4-nitrophenol
sample.files <- c(paste0(gamma.dest, "/", sample.chem.list, ".SDF"))
sample.chem <- lapply(sample.files, read.table, header = FALSE, sep = "\t")
sample.chem <- bind_rows(sample.chem)
write.table(
  sample.chem,
  file = paste0(sample.aff.dir, "/sample.SDF"),
  append = TRUE,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)






#==============================================================================
#                              Missing Data
#==============================================================================
# The alpha, beta, and gamma guests from Rekharsky and Inoue will serve as a 
# "full" data set
all.guests <- unique(c(alpha.guest, beta.guest, gamma.guest))
blank.space <- rep(NA, 555)
thing <- data.frame(all.guests, blank.space)
thing2 <- data.frame(all.guests.clean.na, blank.space)
thing2$all.guests <- all.guests.clean.na
thing2 <- thing2[ , c(3, 2, 1)]
thing3 <- full_join(thing, thing2, by = "all.guests")
thing3 <- thing3[-c(556:574), -c(2,3)]
guest.success <- unique(c(alpha.success, beta.success, gamma.success))
guest.success.df <- data.frame(guest.success)
guest.success.df$all.guests <- guest.success
thing3 <- full_join(thing3, guest.success.df)
missmap(thing3, rank.order = F, col =  c("seagreen1", "slateblue"), main = "Data Loss per Step")
