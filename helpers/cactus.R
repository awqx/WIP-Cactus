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

# escapes special regex characters in a string
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

# takes molecule name and downloads the SDF from NCI
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

# combines a directory of SDFs into one SDF w/ mol name = filename
combine_sdf <- function(mol_dir) {
  mol_files <- list.files(mol_dir, full.names = T)
  mol_names <- list.files(mol_dir) %>%
    str_remove(".SDF")

  # Removing empty files
  mol_info  <- lapply(mol_files, file.info) %>%
    lapply(function(x) x$size > 30) %>%
    unlist()
  mol_files <- mol_files[mol_info]
  mol_names <- mol_names[mol_info]

  # Reading into a list
  sdf_list <- lapply(
    mol_files, 
    read.csv, 
    header = F, 
    stringsAsFactors = F)
  # Assigning the correct names
  for(i in 1:length(sdf_list)) sdf_list[[i]][1, ] <- mol_names[i]

  # return combined list
  do.call(rbind, sdf_list)
}
