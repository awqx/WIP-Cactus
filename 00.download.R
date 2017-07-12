# Libraries ---------------------------------------------------------------
library(data.table)
library(RCurl) # Read webpages with special access requirements
library(XML)   # Import XML files

# Download from ACS -------------------------------------------------------

# Solution by spacedman from stackoverflow:
# http://stackoverflow.com/questions/3616406/ 
# object-moved-error-in-using-the-rcurl-geturl-function-in-order-to-access-an-as

# Warning: the following code will not work if access is not granted
# Alternate code may be created later

# ~~~~ Rekharsly and Inoue Data ~~~~

# Change the working directory location at your discretion
setwd("~/SREP LAB/qsar")
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
save(ri.table.list, file = "./dwnld/ri.table.list.RData")
saveRDS(ri.table.list, file = "./dwnld/ri.table.list.RDS")

# ~~~~ Suzuki Data ~~~~

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

save(suzuki.list, file = "./dwnld/suzuki.list.RData")
saveRDS(suzuki.list, "./dwnld/suzuki.list.RDS")