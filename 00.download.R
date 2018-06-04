# Libraries ---------------------------------------------------------------
packages <- c("data.table", "RCurl", "tabulizer", "XML")
mapply(install.packages, packages)

library(data.table)
library(RCurl) # Read webpages with special access requirements
library(tabulizer) # read tables from PDFs
library(tidyverse)
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

# Suzuki Data ----

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

#  Singh Data -------------------------------------------------------------

convert.ka.delg <- function(ka) {
  return(-8.314*298*log(ka)/1000)
}

# The file should be downloaded as a PDF at the filepath listed below
singh.raw <- extract_tables("./dwnld/singh.pdf") %>% as.data.frame()
singh1 <- singh.raw[[1]][-1, 1:4] %>% data.frame()
colnames(singh1) <- c("guest", "alpha1", "beta1", "gamma1")
singh1[ , 2:4] <- sapply(singh1[ , 2:4], as.character) %>% sapply(., as.numeric)
singh2 <- singh.raw[[2]][-1, 1:4] %>% data.frame()
colnames(singh2) <- c("guest", "alpha2", "beta2", "gamma2")
singh2[ , 2:4] <- sapply(singh2[ , 2:4], as.character) %>% sapply(., as.numeric)
singh <- inner_join(singh1, singh2) %>%
  mutate(alpha = alpha1/alpha2, beta = beta1/beta2, gamma = gamma1/gamma2) %>%
  select(guest, alpha, beta, gamma) 
singh[ , 2:4] <- sapply(singh[ , 2:4], convert.ka.delg)

saveRDS(singh, "./dwnld/singh.RDS")

