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
alpha.guest <- unique(dataset$guest[dataset$host == "1\u03b1"])

source(file = "~/SREP LAB/Rekharsky and Inoue/Cactus/read_cactus.R")
folder  <- "~/SREP LAB/Rekharsky and Inoue/Cactus/"

# Unicode detects the characters alpha and beta
# Alpha == \u03b1, beta == \u03b2
alpha.guest <- unique(dataset$guest[dataset$host == "1\u03b1"])

beta.guest  <- unique(dataset$guest[dataset$host == "1\u03b2"])

gamma.guest <- unique(dataset$guest[dataset$host == "1γ"])

# Cactus is a chemical identifier resolver and can be used
# to find and download SD files used by PyRx

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
    report         <- read_cactus(URL, destfile)
    files.written  <- list.files(path = destfile, pattern = ".SDF")
    lapply(files.written,
           str_extract,
           pattern = "\\.SDF",
           replacement = "")
    data.frame(
      guest = guest,
      downloaded = "yes",
      warning = "no",
      error = "no"
    )
  },
  warning = function(warn) {
    message("probably something irrelevant,
            like the directory already existing")
    Guest.URL      <- unlist(lapply(guest, URLencode, reserved = T))
    URL            <- paste0(
      "https://cactus.nci.nih.gov/chemical/structure/",
      Guest.URL,
      "/",
      chemical.format
    )
    Map(download.file, url = URL, destfile = destfile)
    report         <- read_cactus(URL, destfile)
    files.written  <- list.files(path = destfile, pattern = ".SDF")
    lapply(files.written,
           str_extract,
           pattern = "\\.SDF",
           replacement = "")
    data.frame(
      guest = guest,
      downloaded = "yes",
      warning = "yes",
      error = "no"
    )
  },
  error = function(err) {
    report         <- read_cactus(URL, destfile)
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


empty.gamma <- filter_filesize(
  path = gamma.dest,
  pattern = "SDF",
  size = 100,
  logical = "<"
)

file.remove(empty.gamma$filepath)

# This step adjusts names of molecules that failed to download
# Into a form more likely to be used by Cactus
empty.gamma$molecule <- str_replace(string      = empty.gamma$molecule,
                                    pattern     = "\\α",
                                    replacement = "alpha")

download_cactus(
  guest = empty.gamma$molecule,
  host = "GammaCD",
  path = folder,
  chemical.format = "SDF"
)


#--------------------------------------
#         Downloading Beta
#--------------------------------------

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
    size = 700,
    logical = "<"
  )

file.remove(empty.beta$filepath)

saveRDS(empty.beta, file = "Desktop/Postdoctoral Research/Rekharsky and Inoue/Beta.Guest/empty_beta_sdfiles.RDS")


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
    pattern =  "\\(\\-*\\±*[0-9A-Z]*\\,*
    [0-9A-Z]*\\)\\-\\(*\\+*\\-*\\)*\\-*|nor(?!t)|
    \\([a-z]*\\,*\\s*[a-z]*(I[0-9]\\-)*\\)"
    ,
    replacement = ""
    )


download_cactus(
  guest = empty.beta$molecule,
  host = "Beta.Guest",
  path = folder,
  chemical.format = "SDF"
)

#--------------------------------------
#          Download Alpha
#--------------------------------------
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

empty.alpha <- filter_filesize(
  path = "Desktop/Postdoctoral Research/Rekharsky and Inoue/AlphaCD",
  pattern = "SDF",
  size = 100,
  logical = "<"
)
file.remove(empty.alpha$filepath)

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

# Combining all the SDF files in Gamma CD into one file.
# So far only applies to the gamma CDs.
# Not sure if properly translates to PyRx
sdf.list <-
  list.files(path = gamma.dest, pattern = "SDF")

folder.gammacd <- paste0(folder, "GammaCD/")
sdf.files <- c(paste0(folder.gammacd, sdf.list))
all.sdf <- lapply(sdf.files, read.table, header = FALSE, sep = "\t")
all.sdf <- bind_rows(all.sdf)
write.table(
  all.sdf,
  file = paste0(folder.gammacd, "allgamma.SDF"),
  append = TRUE,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)
