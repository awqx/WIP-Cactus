# Libraries ---------------------------------------------------------------
packages <- c("data.table", "RCurl", "XML")
mapply(install.packages, packages)
library(data.table)
library(RCurl) # Read webpages with special access requirements
library(XML)   # Import XML files

# Download from ACS -------------------------------------------------------

# Solution by spacedman from stackoverflow:
# http://stackoverflow.com/questions/3616406/ 

# Warning: the following code will not work if access is not granted
# Alternate code may be created later

# ---- Rekharsky and Inoue Data ----

url <- "http://pubs.acs.org/doi/full/10.1021/cr970015o"
# The success of this step depends on access to the file
# A paywall will cause this to crash. 
file <- getURL(url,
               .opts = curlOptions(followlocation = TRUE, 
                                   cookiefile = "nosuchfile"))
ri.table.list <-
  readHTMLTable(
    file,
    header = T,
    as.data.frame = T,
    stringAsFactors = F
  )

dir.create("./dwnld")
saveRDS(ri.table.list, file = "./dwnld/ri.table.list.RDS")

# ---- Suzuki Data ----

url   <- "http://pubs.acs.org/doi/full/10.1021/ci010295f"
file  <-
  getURL(url, 
         .opts = curlOptions(followlocation = TRUE, 
                             cookiefile = "nosuchfile"))
suzuki.list <-
  readHTMLTable(
    file,
    header = T,
    as.data.frame = T,
    stringAsFactors = F
  )

saveRDS(suzuki.list, "./dwnld/suzuki.list.RDS")