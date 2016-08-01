install.packages("stringr")
library("stringr")
install.packages("dplyr")
library("dplyr")

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

file.remove(empty.beta$filepath)

saveRDS(results.beta, file = paste0(gamma.dest, "/results.beta.RDS"))

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
gamma.list <-
  list.files(path = gamma.dest, pattern = "SDF")

gamma.files <- c(paste0(gamma.dest, "/", gamma.list))
all.gamma <- lapply(gamma.files, read.table, header = FALSE, sep = "\t")
all.gamma <- bind_rows(all.gamma)
write.table(
  all.gamma,
  file = paste0(all.sdf, "/allgamma.SDF"),
  append = TRUE,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

# Extracting name because the SDFs don't ID by name
formula.gamma <- str_extract(all.gamma$V1, "C[[:digit:]]+[[:alnum:]]+")
formula.gamma <- formula.gamma[!is.na(formula.gamma)]

gamma.names <- str_replace(gamma.list, ".SDF", "")

gamma.names.table <- data.frame(formula.gamma, gamma.names)

# For Beta
beta.list <-
  list.files(path = beta.dest, pattern = ".SDF")

beta.files <- c(paste0(beta.dest, "/", beta.list))

all.beta <- lapply(beta.files, read.table, header = FALSE, sep = "\t")
all.beta <- bind_rows(all.beta)
write.table(
  all.beta,
  file = paste0(all.sdf, "/allbeta.SDF"),
  append = TRUE,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

formula.beta <- str_extract(all.beta$V1, "C[[:digit:]]+[[:alnum:]]+|CF3O3S|ClO4|F6P")
formula.beta <- formula.beta[!is.na(formula.beta)]
beta.names <- str_replace(beta.list, ".SDF", "")
beta.names.table <- data.frame(formula.beta, beta.names)



# For Alpha
alpha.list <-
  list.files(path = alpha.dest, pattern = ".SDF")

alpha.files <- c(paste0(alpha.dest, "/", alpha.list))
all.alpha <- lapply(alpha.files, read.table, header = FALSE, sep = "\t")
all.alpha <- bind_rows(all.alpha)
write.table(
  all.alpha,
  file = paste0(all.sdf, "/allalpha.SDF"),
  append = TRUE,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)
formula.alpha <- str_extract(all.alpha$V1, "C[[:digit:]]+[[:alnum:]]+|F6P|^[I]+$|ClO4|ClHO4|FO3P|CNS|CH2O2|CF3O3S")
formula.alpha <- formula.alpha[!is.na(formula.alpha)]
alpha.names <- str_replace(alpha.list, ".SDF", "")

alpha.names.table <- data.frame(formula.alpha, alpha.names)
