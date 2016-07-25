dataset <- readRDS(file = "~/SREP LAB/Rekharsky and Inoue/Cactus/RI.rds")
source(file = "~/SREP LAB/Rekharsky and Inoue/Cactus/read_cactus.R")
folder  <- "~/SREP LAB/Rekharsky and Inoue/Cactus/"

alpha.guest <- unique(Dataset$guest[Dataset$host == "1α"])

beta.guest  <- unique(Dataset$guest[Dataset$host == "1β"])

gamma.guest <- unique(Dataset$guest[Dataset$host == "1γ"])

download_cactus <- function(guest, host, path, chemical.format) {
  host.directory <- paste0(path , host)
  dir.create(path = host.directory)
  destfile       <- paste0(host.directory, "/", guest, ".SDF")
  # Chemical format must be parsed to match all the outputs from NCI cactus
  Guest.URL      <- unlist(lapply(guest, URLencode, reserved = T))
  URL            <-
    paste0(
      "https://cactus.nci.nih.gov/chemical/structure/"
      ,
      Guest.URL,
      "/",
      chemical.format
    )
  #Map(download.file, url = URL, destfile = destfile, method = "curl")
  report         <- read_cactus(URL, destfile)
  files.written  <- list.files(path = destfile, pattern = ".SDF",)
  lapply(files.written,
         str_extract,
         pattern = "\\.SDF",
         replacement = "")
  return(report)
}

filter_filesize <- function(path, pattern, size, logical) {
  filenames <- list.files(path = path, pattern = pattern, )
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

download_cactus(
  guest = Gamma.Guest,
  host = "GammaCD",
  path = folder,
  chemical.format = "SDF"
)

Empty.Gamma <- filter_filesize(
    path = "~SREP LAB/Cactus/Gamma.Guest/",
    pattern = "SDF",
    size = 100,
    logical = "<"
  )

file.remove(Empty.Gamma$filepath)

Empty.Gamma$molecule <- str_replace(
  string = Empty.Gamma$molecule,
  Empty.Gamma$molecule <-
    str_replace(
      string = Empty.Gamma$molecule,
      pattern = "\\α",
      replacement = "alpha"
    )
)
#Added a bracket that was messing up all the code after it
  download_cactus(
    guest = Empty.Gamma$molecule,
    host = "BetaCD",
    path = folder,
    chemical.format = "SDF"
  )
  
  Empty.Beta <-
    filter_filesize(
      path = "~SREP LAB/Cactus/Beta.Guest",
      pattern = "SDF",
      size = 700,
      logical = "<"
    )
  
  file.remove(Empty.Beta$filepath)
  saveRDS(Empty.Beta, file = "Desktop/Postdoctoral Research/Rekharsky and Inoue/Beta.Guest/empty_beta_sdfiles.RDS")
  
  
  Empty.Beta$molecule <-
    str_replace(
      string = Empty.Beta$molecule,
      pattern = "\\β",
      replacement = "beta"
    )
  Empty.Beta$molecule <-
    str_replace(
      string = Empty.Beta$molecule,
      pattern = "\\α",
      replacement = "alpha"
    )
  Empty.Beta$molecule <-
    str_replace(
      string = Empty.Beta$molecule,
      pattern =  "\\(\\-*\\±*[0-9A-Z]*\\,*
      [0-9A-Z]*\\)\\-\\(*\\+*\\-*\\)*\\-*|nor(?!t)|
      \\([a-z]*\\,*\\s*[a-z]*(I[0-9]\\-)*\\)"
      ,
      replacement = ""
)
  
  
  download_cactus(
    guest = Empty.Beta$molecule,
    host = "Beta.Guest",
    path = folder,
    chemical.format = "SDF"
  )
  
  download_cactus(
    guest = Alpha.Guest,
    host = "AlphaCD",
    path = folder,
    chemical.format = "SDF"
  )
  Empty.Alpha <- filter_filesize(
    path = "Desktop/Postdoctoral Research/Rekharsky and Inoue/AlphaCD",
    pattern = "SDF",
    size = 100,
    logical = "<"
  )
  file.remove(Empty.Alpha$filepath)
  
  Empty.Alpha$molecule <-
    str_replace(
      string = Empty.Alpha$molecule,
      pattern =  "\\(\\-*\\±*[0-9A-Z]*\\,*[0-9A-Z]*\\)\\-\\(*\\+*\\-*\\)*\\-*|nor(?!t)|\\([a-z]*\\,*\\s*[a-z]*(I[0-9]\\-)*\\)",
      replacement = ""
    )
  Empty.Alpha$molecule <-
    str_replace(
      string = Empty.Alpha$molecule,
      pattern = "\\β",
      replacement = "beta"
    )
  Empty.Alpha$molecule <-
    str_replace(
      string = Empty.Alpha$molecule,
      pattern = "\\α",
      replacement = "alpha"
    )
  
  download_cactus(
    guest = Empty.Alpha$molecule,
    host = "AlphaCD",
    path = folder,
    chemical.format = "SDF"
  )
  
  SDFlist <-
    list.files(path = paste0(folder, "Gamma.Guest"), pattern = "SDF"
    )