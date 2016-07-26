#adding some packages
install.packages("stringr")
library("stringr")
install.packages("dplyr")
library("dplyr")

#Tilde on the file path
dataset <- readRDS(file = "~/SREP LAB/Rekharsky and Inoue/Cactus/RI.rds")
#There's a small typo in row 1081, column "guest"
#This will matter at the download_cactus stage
dataset[1081, "guest"] <- "4-[(4-hydroxyphenyl)azo]benzoate"
source(file = "~/SREP LAB/Rekharsky and Inoue/Cactus/read_cactus.R")
folder  <- "~/SREP LAB/Rekharsky and Inoue/Cactus/"

#fixed unicode to detect the characters alpha and beta
#alpha == \u03b1, beta == \u03b2
#changed values to all lowercase for style and ease of use
alpha.guest <- unique(dataset$guest[dataset$host == "1\u03b1"])

beta.guest  <- unique(dataset$guest[dataset$host == "1\u03b2"])

gamma.guest <- unique(dataset$guest[dataset$host == "1γ"])

download_cactus <- function(guest, host, path, chemical.format) {
  host.directory <- paste0(path, host)
  dir.create(path = host.directory)
  destfile       <- paste0(host.directory, "/", guest, ".SDF")
  # Chemical format must be parsed to match all the outputs from NCI cactus
  Guest.URL      <- unlist(lapply(guest, URLencode, reserved = T))
  URL            <- paste0(
      "https://cactus.nci.nih.gov/chemical/structure/",
      Guest.URL, "/", chemical.format
    )
  #Removed method = "curl" from Map, fixed nonzero exit status error when downloading
  Map(download.file, url = URL, destfile = destfile)
  report         <- read_cactus(URL, destfile)
  files.written  <- list.files(path = destfile, pattern = ".SDF")
  lapply(files.written, str_extract, pattern = "\\.SDF", replacement = "")
  return(report)
}

filter_filesize <- function(path, pattern, size, logical) {
  filenames <- list.files(path = path, pattern = pattern)
  location  <- dir(path = path, full.names = T, recursive = T)
  data      <- file.info(paste0(path, "/", filenames))[1]
  index     <- eval(call(logical, data, size))
  return(data.frame(
    filepath = location[index],
    molecule = filenames[index],
    stringsAsFactors = F
    )
  )
}

download_cactus(
  guest = gamma.guest,
  host = "GammaCD",
  path = folder,
  chemical.format = "SDF"
)
#Cactus can't ID the last 3 gamma.guests:
#cis-1,2,3,4-tetraphenylcyclobutane, trans-1,2,3,4-tetraphenylcyclobutane
#and a-(2,4,6-trimethoxyphenyl)benzyl tert-butyl nitroxide

empty.gamma <- filter_filesize(
  path = "~/SREP LAB/Rekharsky and Inoue/Cactus/GammaCD",
  pattern = "SDF",
  size = 100,
  logical = "<"
)

file.remove(empty.gamma$filepath)

#Removed redundancy in the code
#I don't know why this step even exists 
#b/c it replaces something that was removed last step
empty.gamma$molecule <- str_replace(
  string      = empty.gamma$molecule,
  pattern     = "\\α",
  replacement = "alpha"
)

#changed guest from empty.gamma$molecule to beta.guest
#b/c I thought that made more sense
#don't actually know if that was the goal
download_cactus(
  guest = beta.guest,
  host = "BetaCD",
  path = folder,
  chemical.format = "SDF"
)

empty.beta <-
  filter_filesize(
    path = "~SREP LAB/Cactus/BetaCD",
    pattern = "SDF",
    size = 700,
    logical = "<"
  )

file.remove(empty.beta$filepath)

saveRDS(empty.beta, file = "Desktop/Postdoctoral Research/Rekharsky and Inoue/Beta.Guest/empty_beta_sdfiles.RDS")


empty.beta$molecule <-
  str_replace(
    string = empty.beta$molecule,
    pattern = "\\β",
    replacement = "beta"
  )
empty.beta$molecule <-
  str_replace(
    string = empty.beta$molecule,
    pattern = "\\α",
    replacement = "alpha"
  )
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

download_cactus(
  guest = alpha.guest,
  host = "AlphaCD",
  path = folder,
  chemical.format = "SDF"
)
empty.alpha <- filter_filesize(
  path = "Desktop/Postdoctoral Research/Rekharsky and Inoue/AlphaCD",
  pattern = "SDF",
  size = 100,
  logical = "<"
)
file.remove(empty.alpha$filepath)

empty.alpha$molecule <-
  str_replace(
    string = empty.alpha$molecule,
    pattern =  "\\(\\-*\\±*[0-9A-Z]*\\,*[0-9A-Z]*\\)\\-\\(*\\+*\\-*\\)*\\-*|nor(?!t)|\\([a-z]*\\,*\\s*[a-z]*(I[0-9]\\-)*\\)",
    replacement = ""
  )
empty.alpha$molecule <-
  str_replace(
    string = empty.alpha$molecule,
    pattern = "\\β",
    replacement = "beta"
  )
empty.alpha$molecule <-
  str_replace(
    string = empty.alpha$molecule,
    pattern = "\\α",
    replacement = "alpha"
  )

download_cactus(
  guest = empty.alpha$molecule,
  host = "AlphaCD",
  path = folder,
  chemical.format = "SDF"
)

sdf.list <-
  list.files(path = paste0(folder, "GammaCD"), pattern = "SDF"
  ) 

#Combining all the SDF files in Gamma CD into one file
#So far only applies to the gamma CDs
#Not sure if properly translates to PyRx, but the .txt file looks
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