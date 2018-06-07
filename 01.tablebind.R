# Libraries and Packages --------------------------------------------------

library(tidyverse)
library(data.table)

# Version 1 (RI) ---------------------------------------------------------------
# NOTICE: Use Version 1 for files downloaded with 00.download.R

# Load the table from file created in 00.download.R
ri.table.list <- readRDS("./dwnld/ri.table.list.RDS")

# Find entries with 10 columns - missing DelCp column
index10 <- ri.table.list %>% 
  lapply(names) %>% 
  lapply(length) == 10 %>% 
  as.vector()
# Index11 have all the data
index11 <- ri.table.list %>% 
  lapply(names) %>% 
  lapply(length) == 11 %>% 
  as.vector()

ri.delcp    <- rbindlist(ri.table.list[index11]) %>%
  lapply(as.character) %>%
  as.data.frame(stringsAsFactors = F)
ri.delcp  <- ri.delcp[, c(1:8, 10, 11, 9)]

ri.delcp.none  <- rbindlist(ri.table.list[index10]) %>%
  lapply(as.character) %>%
  as.data.frame(stringsAsFactors = F)
ri.delcp.none[, 11] <- NA

ri.bound.df <- rbindlist(list(ri.delcp, ri.delcp.none))

saveRDS(ri.bound.df, "./dwnld/ri.bound.df.RDS")

# Version 1 (Suzuki) ------------------------------------------------------

suzuki.list <- readRDS("./dwnld/suzuki.list.RDS")

index11 <- suzuki.list %>% 
  lapply(names) %>% 
  lapply(length) == 11 %>% 
  as.vector()

suz.bound    <- rbindlist(suzuki.list[index11]) %>%
  lapply(as.character) %>%
  as.data.frame(stringsAsFactors = F)
suz.bound  <- suz.bound[, c(1:8, 10, 11, 9)]

saveRDS(suz.bound, "./dwnld/suz.bound.df.RDS")

# Version 2 (HTML Download) ---------------------------------------------------

# # NOTICE: Use this code in the event of downloading the HTML code manually, 
# # rather than with the 00.download.R code
# # Commenting it out for now, for clarity 
# 
# directory <- paste0("./dwnld/ri.table.list.rds")
# ri.table.list <- readRDS(directory)
# 
# # Cleaning HTML import prior to binding due to:
# # (1) Ref Column in list entry 26 is an integer when Ref columns in other 
# #     list entries are integer
# 
# ri.table.list[[26]][11] <- ri.table.list[[26]][11] %>% map(as.character)
# 
# # (2) Some list entries [14,17,21] are actually table entries that were imported 
# #     as empty table names. Some table list entries [15,18,22] had no title 
# #     and its colnames were actual results from the paper.
# 
# # ri.15 <- ri.table.list[[15]] %>% names()
# names(ri.table.list[[15]]) <- c("host",
#                                 "guest",
#                                 "solvent",
#                                 "T/K",
#                                 "logâK",
#                                 "ÎGÂ°/ kJâmol-1",
#                                 "ÎHÂ°/ kJâmol-1",
#                                 "TÎSÂ°/ kJâmol-1",
#                                 "methoda" ,
#                                 "ref")
# 
# # ri.18 <- ri.table.list[[18]] %>% names()
# names(ri.table.list[[18]]) <- c("host",
#                                 "guest",
#                                 "solvent",
#                                 "T/K",
#                                 "logâK",
#                                 "ÎGÂ°/ kJâmol-1",
#                                 "ÎHÂ°/ kJâmol-1",
#                                 "TÎSÂ°/ kJâmol-1",
#                                 "methoda" ,
#                                 "ref")
# 
# # ri.22 <- ri.table.list[[22]] %>% names()
# names(ri.table.list[[22]]) <- c("host",
#                                 "guest",
#                                 "solvent",
#                                 "T/K",
#                                 "logâK",
#                                 "ÎGÂ°/ kJâmol-1",
#                                 "ÎHÂ°/ kJâmol-1",
#                                 "TÎSÂ°/ kJâmol-1",
#                                 "methoda" ,
#                                 "ref")
# 
# # (2.3) List table 21 is missing the log K value
# 
# ri.table.list[[21]] <- ri.table.list[[21]] %>%
#   names() %>% 
#   data_frame(entries = .) %>%
#   mutate(key = c("host",
#                  "guest",
#                  "solvent",
#                  "T/K",
#                  "ÎGÂ°/ kJâmol-1",
#                  "ÎHÂ°/ kJâmol-1",
#                  "TÎSÂ°/ kJâmol-1",
#                  "methoda" ,
#                  "ref")) %>% 
#   spread(key = key, value = entries) %>%
#   mutate("T/K" = as.integer(`T/K`),
#          "logâK" = "NA")
# 
# 
# # (2.2) List table 17 is missing its reference (Likely to be reference  
# # 242 since all other 5-methoxyresorcinol are from reference 242)
# ri.table.list[[17]] <- ri.table.list[[17]] %>% 
#   names() %>% 
#   data_frame(entries = . ) %>%
#   mutate(key = c("host",
#                  "guest",
#                  "solvent",
#                  "T/K",
#                  "logâK",
#                  "ÎGÂ°/ kJâmol-1",
#                  "ÎHÂ°/ kJâmol-1",
#                  "TÎSÂ°/ kJâmol-1",
#                  "methoda")) %>%
#   spread(key = key, value = entries) %>%
#   mutate("T/K" = as.integer(`T/K`),
#          ref = "NA")
# 
# # (2.1) List table 14 is missing the log K value
# ri.table.list[[14]] <- ri.table.list[[14]] %>% 
#   names() %>% 
#   data_frame(entries = . ) %>%
#   mutate(key = c("host",
#                  "guest",
#                  "solvent",
#                  "T/K",
#                  "ÎGÂ°/ kJâmol-1",
#                  "ÎHÂ°/ kJâmol-1",
#                  "TÎSÂ°/ kJâmol-1",
#                  "methoda" ,
#                  "ref")) %>% 
#   spread(key = key, value = entries) %>%
#   mutate("T/K" = as.integer(`T/K`),
#          "logâK" = "NA")
# 
# # (3) Remove the summary output from the linear regression and save into its own
# # data frame
# 
# delThermPararm_delAlk <- ri.table.list[[28]] 
# 
# summary_table <- ri.table.list[[29]]
# 
# ri.table.list[[29]] <- NULL
# ri.table.list[[28]] <- NULL
# 
# ri.table.list %>%
#   tibble() %>%
#   unnest() -> ri_10and11col_allhost
# 
# # Saving efforts to file
# 
# save(ri_10and11col_allhost, file = "./Output Data/01-ri_boundDF.RData")
# saveRDS(ri_10and11col_allhost, file = "./Output Data/01-ri_boundDF.rds")
# save(delThermPararm_delAlk, file = "./Output Data/01.2-delThermodynamicParameter_delCcontent.RData")
# saveRDS(delThermPararm_delAlk, file = "./Output Data/01.2-delThermodynamicParameter_delCcontent.rds")
# save(summary_table, file = "./Output Data/01.3-DifferentCyclodextrins_SummaryTable.RData")
# save(summary_table, file = "./Output Data/01.3-DifferentCyclodextrins_SummaryTable.rds")
