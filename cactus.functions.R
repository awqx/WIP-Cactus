# Cactus is a chemical identifier resolver and can be used
# to find and download SD files used by PyRx
# But there's sometimes hiccups in downloading, so I used a tryCatch function
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
filter.filesize <- function(path, pattern, size, logical) {
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
# Read SDF files and assign a chemical name -- this makes things easier later
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