# downloads data on cyclodextrin binding from HTML tables

# Libraries ---------------------------------------------------------------

# If tabulizer doesn't load, make sure that the installed version of Java
# matches the version of R (32-bit or 64-bit)
packages <- c("data.table", "httr", "tabulizer", "tidyverse", "XML")
mapply(install.packages, packages)
lapply(packages, require, character.only = T)
if (!dir.exists("dwnld")) dir.create("dwnld")

# Rekharsky and Inoue -----------------------------------------------------

# the following code will not work if access is not granted
# a vpn granting academic access is suggested

ri_html <- GET("http://pubs.acs.org/doi/full/10.1021/cr970015o")
ri_list <- readHTMLTable(
  rawToChar(ri_html$content), 
  header = T,
  as.data.frame = T,
  stringAsFactors = F
)
# saveRDS(ri_list, file = "dwnld/ri-list.RDS")

# ---- Combining table ----------------------------------------------------

# entries in the list with necessary information
ri_index <- ri_list %>% 
  lapply(names) %>% 
  lapply(length) >= 10 %>% 
  as.vector()
# columns 1, 2, 6 correspond to host, guest, and deltaG
# columns 3, 4, and 10/11 (last) are for sorting later
# setNames ensures correct rbind
ri_df <- ri_list[ri_index] %>%
  lapply(
    function(x) 
      setNames(
        mutate_all(x, as.character)[, c(1, 2, 6, 3, 4, length(x))], 
        c("host", "guest", "dG", "solvent", "temp", "ref")
      )
    ) %>%
  bind_rows() %>%
  mutate(host = str_replace_all(host, "1\u03b1", "alpha")) %>%
  mutate(host = str_replace_all(host, "1\u03b2", "beta")) %>%
  mutate(host = str_replace_all(host, "1\u03b3", "gamma"))
saveRDS(ri_df, "dwnld/ri-df.RDS")

# Suzuki ------------------------------------------------------------------

suzuki_html <- GET("http://pubs.acs.org/doi/full/10.1021/ci010295f")
suzuki_list <-
  readHTMLTable(
    rawToChar(suzuki_html$content),
    header = T,
    as.data.frame = T,
    stringAsFactors = F
  )
# saveRDS(suzuki_list, "dwnld/suzuki-list.RDS")
# tables with the complete information
suzuki_index <- suzuki_list %>% 
  lapply(names) %>% 
  lapply(length) == 11 %>% 
  as.vector()
# the suzuki data is organized so that alpha comes first in the table
# and then beta comes second. 
# the table is wide, not long
suzuki_alpha <- suzuki_list[suzuki_index] %>%
  lapply(
    function(x)
      mutate_all(x, as.character)[-1, c(2, 4)]
  ) %>%
  bind_rows() %>%
  mutate(host = "alpha") %>%
  setNames(c("guest", "dG", "host")) %>%
  select(host, guest, dG)# style choice to reorder
suzuki_beta <- suzuki_list[suzuki_index] %>%
  lapply(
    function(x)
      mutate_all(x, as.character)[-1, c(2, 8)]
  ) %>%
  bind_rows() %>%
  mutate(host = "beta") %>%
  setNames(c("guest", "dG", "host")) %>%
  select(host, guest, dG)
suzuki_df <- rbind(suzuki_alpha, suzuki_beta)
saveRDS(suzuki_df, "dwnld/suzuki-df.RDS")


# Singh -------------------------------------------------------------------

# converts a Ka value into dG (kJ/mol)
convert_ka <- function(ka) (-8.314*298*log(ka)) /1000

# file provided as a pdf only; was only made available upon request
singh_raw <- extract_tables("dwnld/singh.pdf")# %>% as.data.frame()
# table 1 contains ka values
singh <- data.frame(singh_raw[[1]][-1, 1:4])
colnames(singh) <- c("guest", "alpha", "beta", "gamma")
singh[, 2:4] <- sapply(singh[, 2:4], as.character) %>%
  sapply(., as.numeric)
singh[, 2:4] <- sapply(singh[ , 2:4], convert_ka)
singh <- melt(
  data.table(singh),
  id.vars = "guest",
  measure.vars = c(2:4),
  variable.name = "cd",
  value.name = "dG"
)
saveRDS(singh, "dwnld/singh.RDS")