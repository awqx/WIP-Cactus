# Created a new function
create.host.dir <- function(path, host) {
  host.directory <- paste0(path, host)
  dir.create(path = host.directory)
  return(host.directory)
}

# This version of cactus is the best so far
# Never mind. 
#I'm gonna keep it anyway for reference
download.cactus.stable <- function(guest, path, chemical.format) {
  report <- tryCatch({
    destfile       <- paste0(path, "/", guest, ".SDF")
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
    message("probably something irrelevant, 
            like the directory already existing")
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
    report         <- read_cactus(URL, destfile)
  },
  finally = {
    message("Chemical name processed")
  })
  return(report)
}

#An attempt to note down results
note.results <- function(guest){
  good.col <- c(rep(NA, length(guest)))
  warning.col <- c(rep(NA, length(guest)))
  error.col <- c(rep(NA, length(guest)))
  data.frame(guest, good.col, warning.col, error.col)
}
guest.sample <- c(1:10)

results.sample <- note.results(guest.sample)

download.cactus.results <- function(guest, path, chemical.format, dataframe) {
  report <- tryCatch({
    destfile       <- paste0(path, "/", guest, ".SDF")
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
    location <- which(guest == guest)
    dataframe["good.col", location] <- "y"
    dataframe["warning.col", location] <- "n"
    dataframe["error.col", location] <- "n"
  },
  warning = function(warn) {
    message("probably something irrelevant, 
            like the directory already existing")
    Guest.URL      <- unlist(lapply(guest, URLencode, reserved = T))
    URL            <- paste0(
      "https://cactus.nci.nih.gov/chemical/structure/",
      Guest.URL, "/", chemical.format
    )
    Map(download.file, url = URL, destfile = destfile)
    report         <- read_cactus(URL, destfile)
    files.written  <- list.files(path = destfile, pattern = ".SDF")
    lapply(files.written, str_extract, pattern = "\\.SDF", replacement = "")
    location <- which(guest == guest)
    dataframe["good.col", location] <- "y"
    dataframe["warning.col", location] <- "y"
    dataframe["error.col", location] <- "n"
  },
  error = function(err) {
    report         <- read_cactus(URL, destfile)
    location <- which(guest == guest)
    dataframe["good.col", location] <- "n"
    dataframe["warning.col", location] <- "y"
    dataframe["error.col", location] <- "y"
  },
  finally = {
    message("Chemical name processed")
  })
  return(report)
}