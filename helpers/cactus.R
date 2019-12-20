# This file is for loading the functions without having to go through the 
# trouble of accidentally sourcing 03.cactus.R and downloading 615
# molecules unnecessarily

# Libraries and Packages --------------------------------------------------

packages <- c("httr", "stringr", "tidyverse", "XML")
lapply(packages, require, character.only = T)

# Functions ---------------------------------------------------------------

# write.sdf <- function(name, loc) {
#   path <- paste0(loc, "/", name, ".SDF")
#   write.table(file = path)
# }

make_regex <- function(s) {
  str_replace_all(
    s, pattern = "\\(", 
    replacement = "\\\\("
    ) %>% 
    str_replace_all(
      pattern = "\\)", 
      replacement = "\\\\)"
      ) %>%
    str_replace_all(
      pattern = "\\-", 
      replacement = "\\\\-")
}

download_sdf <- function(guest, path, chemical_format) {
  report <- tryCatch({ 
    destfile   <- paste0(path, "/", guest, ".SDF")
    guest_url  <- unlist(lapply(guest, URLencode, reserved = T))
    cactus_url <- paste0(
      "https://cactus.nci.nih.gov/chemical/structure/",
      guest_url, "/", chemical_format
    )
    GET(cactus_url, write_disk(destfile, overwrite = T))
    data.frame(
      guest = guest,
      downloaded = "yes",
      warning = "no",
      error = "no"
    )
  },
  warning = function(warn) {
    message("Warning: either URL error or already existing directory.")
    destfile   <- paste0(path, "/", guest, ".SDF")
    guest_url  <- unlist(lapply(guest, URLencode, reserved = T))
    cactus_url <- paste0(
      "https://cactus.nci.nih.gov/chemical/structure/",
      guest_url, "/", chemical_format
    )
    GET(cactus_url, write_disk(destfile, overwrite = T))
    data.frame(
      guest = guest,
      downloaded = "yes",
      warning = "no",
      error = "no"
    )
  },
  error = function(err) {
    message("An error occurred:")
    message(err)
    data.frame(
      guest = guest,
      downloaded = "no",
      warning = "yes",
      error = "yes"
    )
    
  },
  finally = {
    message(paste0("Processed ", guest))
  })
  return(report)
}
