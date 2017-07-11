# Libraries ---------------------------------------------------------------
library(XML)   # Import XML files
library(RCurl) # Read webpages with special access requirements


# Download from ACS -------------------------------------------------------

# Solution by spacedman from stackoverflow
# http://stackoverflow.com/questions/3616406/
# object-moved-error-in-using-the-rcurl-geturl-function-in-order-to-access-an-as

url           <- "http://pubs.acs.org/doi/full/10.1021/cr970015o"
file          <-
  getURL(url, .opts = curlOptions(followlocation = TRUE, cookiefile = "nosuchfile"))
ri.table.list <-
  readHTMLTable(
    file,
    header = T,
    as.data.frame = T,
    stringAsFactors = F
  )
# date.accessed  <- date()
dwnld.dir <- "./downloaded/"
# if the previous URL does not work, try "~/downloaded/"
dir.create(dwnld.dir)
save(ri.table.list, file = paste0(dwnld.dir, "ri.table.list.RData"))
saveRDS(ri.table.list, file = paste0(dwnld.dir, "ri.table.list.RDS"))
