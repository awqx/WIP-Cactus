require("stringr")
require("tidyverse")

source("03.1.molecule.renaming.R")
suzuki <- readRDS("./dwnld/02.suz.clean.RDS")
ri <- readRDS("./dwnld/02.ri.squeaky.clean.RDS")


# Cleaning names ----------------------------------------------------------

# Unconventional apostrophe \\â or â\u0080\u0098  
suzuki$guest <- suzuki$guest %>%
  str_replace("(\\â | â\u0080\u0098)", "\\'") 

# Requires manual replacement - 42
# Changed the previous method (vector, gather, etc.) to a .csv for easier 
# management
fixed.names <- read.csv("cleaning/fixed.names.csv", header = T) %>%
  mutate(pattern = as.character(pattern), 
         replacement = as.character(replacement))

# Converting the table to regex for str_replace
pattern.reg <- fixed.names %>%
  mutate(pattern = paste0(make.regex(fixed.names$pattern), "$")) %>%
  mutate(replacement = make.regex(fixed.names$replacement))

for(i in 1:nrow(pattern.reg)) {
  suzuki$guest <- str_replace(suzuki$guest, pattern = pattern.reg[i, 1], 
                               replacement = pattern.reg[i, 2])
}
# "-acid"has problems reading from the .csv
suzuki$guest <- str_replace(suzuki$guest, pattern = "-acid", 
                             replacement = " acid")

suzuki.alpha.guest <- suzuki %>% filter(host == "alpha") %>% .$guest
