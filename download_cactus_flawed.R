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
  Map(download.file, url = URL, destfile = destfile)
  report         <- read_cactus(URL, destfile)
  files.written  <- list.files(path = destfile, pattern = ".SDF")
  lapply(files.written, str_extract, pattern = "\\.SDF", replacement = "")
  return(report)
}
#Attempt to fix download_cactus w/ tryCatch
download.cactus.stable <- function(guest, host, path, chemical.format) {
  report <- tryCatch(
    {
      host.directory <- paste0(path, host)
      dir.create(path = host.directory)
      destfile       <- paste0(host.directory, "/", guest, ".SDF")
      # Chemical format must be parsed to match all the outputs from NCI cactus
      Guest.URL      <- unlist(lapply(guest, URLencode, reserved = T))
      URL            <- paste0(
        "https://cactus.nci.nih.gov/chemical/structure/",
        Guest.URL, "/", chemical.format
      )
      Map(download.file, url = URL, destfile = destfile)
      report         <- read_cactus(URL, destfile)
      files.written  <- list.files(path = destfile, pattern = ".SDF")
      lapply(files.written, str_extract, pattern = "\\.SDF", replacement = "")
    }, 
    warning = function(warn) {
      message("probably something irrelevant, like the directory already existing")
      host.directory <- paste0(path, host)
      destfile       <- paste0(host.directory, "/", guest, ".SDF")
      # Chemical format must be parsed to match all the outputs from NCI cactus
      Guest.URL      <- unlist(lapply(guest, URLencode, reserved = T))
      URL            <- paste0(
        "https://cactus.nci.nih.gov/chemical/structure/",
        Guest.URL, "/", chemical.format
      )
      Map(download.file, url = URL, destfile = destfile)
      report         <- read_cactus(URL, destfile)
      files.written  <- list.files(path = destfile, pattern = ".SDF")
      lapply(files.written, str_extract, pattern = "\\.SDF", replacement = "")
    },
    error = function(err) {
      message("URL doesn't exist")
      return(NA)
    }, 
    finally ={
      message("Chemical processed")
    }
  )
  return(report)
}
