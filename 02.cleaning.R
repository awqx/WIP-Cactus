library(tidyverse)
library(stringr)
# The following reading will vary based on user preference
ri.bound <- readRDS(file = paste0(bound.dir, "ri.bound.df.RDS"))


# Rename columns 
#     Note: for some reason, this doesn't work on Windows
#     Also, I'm pretty sure this is supposed to go at the end of the code

ri.bound <- ri.bound %>%
  rename(DelG = `ÎGÂ°/ kJâmol-1`, 
         DelH =`ÎHÂ°/ kJâmol-1`,
         TDelS = `TÎSÂ°/ kJâmol-1`,
         DelCp = `ÎCpÂ°/ Jâmol-1âK-1`,
         log.K = `logâK`,
         method = `methoda`)


# Splitting columns containing a variable value and its uncertainty
# Occasionally the naming of the columns will get messed up
# In this case, you would need colnames(ri.bound)[3], [5], [6], [7], and [11]

# ri.clean <- ri.bound %>% 
#   separate(solvent, c("solvent","solvent.specs"),
#            sep = "(?=\\s*\\()", extra = "merge", fill = "right") %>% 
#   separate(log.K, c("log.K", "log.K.Uncertainty"), 
#            sep = "\\s\\Â±\\s", extra = "merge", fill = "right") %>% 
#   separate(DelG, c("DelG", "DelG.Uncertainty"),
#            sep = "\\s\\Â±\\s",extra = "merge", fill = "right") %>%
#   separate(DelH, c("DelH", "DelH.Uncertainty"),
#            sep = "\\s\\Â±\\s", extra = "merge", fill = "right") %>%
#   separate(TDelS, c("TDelS", "TDelS.Uncertainty"), 
#            sep = "\\s\\Â±\\s",extra = "merge", fill = "right")  %>%
#   separate(DelCp, c("DelCp", "DelCp.Uncertainty"), 
#           sep = "\\s\\Â±\\s",extra = "merge", fill = "right")
# My version of the table
ri.clean <- ri.bound %>% 
  separate(solvent, c("solvent","solvent.specs"),
           sep = "(?=\\s*\\()", extra = "merge", fill = "right") %>% 
  separate(log.U.2009.K, c("log.K", "log.K.Uncertainty"), 
           sep = "\\s\\Â±\\s", extra = "merge", fill = "right")   %>% 
  separate(X.U.0394.G...kJ.U.2009.mol.1, c("DelG", "DelG.Uncertainty"),
           sep = "\\s\\Â±\\s",extra = "merge", fill = "right")    %>%
  separate(X.U.0394.H...kJ.U.2009.mol.1, c("DelH", "DelH.Uncertainty"),
           sep = "\\s\\Â±\\s", extra = "merge", fill = "right")   %>%
  separate(T.U.0394.S...kJ.U.2009.mol.1, c("TDelS", "TDelS.Uncertainty"), 
           sep = "\\s\\Â±\\s",extra = "merge", fill = "right")    %>%
  separate(X.U.0394.Cp...J.U.2009.mol.1.U.2009.K.1, c("DelCp", "DelCp.Uncertainty"), 
           sep = "\\s\\Â±\\s",extra = "merge", fill = "right")

# Cleaning strings that contain inconsistent patterns such as more than one space
# separation and unconventional separation symbols. Converting alpha, beta and
# gamma symbols to words for easier subset of tables 

ri.clean <- ri.clean %>%
  lapply(str_replace_all, pattern = "\\â", replacement = "-") %>%
  lapply(str_replace_all, pattern = "\\â+", replacement = " ") %>%
  lapply(str_replace_all, pattern = "\\s+", replacement = " ") %>%
  lapply(str_replace_all, pattern = "Â·", replacement = " ") %>% 
  as_tibble() %>%
  mutate(host = str_replace(host, pattern =  "\\Î±|\u03b1",replacement = "alpha" )) %>%
  mutate(host = str_replace(host, pattern = "\\Î²|\\B\\s+H\\+|\u03b2", replacement = "beta")) %>%
  mutate(host = str_replace(host, pattern = "\\Î³", replacement = "gamma"))






#=============================================================================== 
#                             pH Imputation                                    =
#===============================================================================
# Remove pH string and parenthesis to convert the pH Column to numerical. Assume 
# that columns with no value have a pH of 7.0. pH values with a range of 
# phValue1-phValue2 will be assigned an average pH value ((phValue1+phValue2)/2)
# ph's with an inequality (<,>) will be set at the value given. Guests with 
# the molarity of the acid given in the solvent composition will be set
# at the theoretical value of the solution

pattern.molarity <- gsub("\n", replacement = "", x = "(([0-9]*\\.*[0-9]*\\s+[M]
                         \\s+[A-Za-z]*[A-Za-z0-9]*\\s*[a-z]*)(\\;*\\,*\\+*\\s+[0-9]+\\.*[0-9]*\\s+[M]
                         \\s+[A-Za-z]*[A-Za-z0-9]*[a-z]*)*)")

pH.numeric <-  "[0-9]+\\.*[0-9]*"

ri.clean <- ri.clean %>% 
  mutate(ref.notes= str_extract(ref, pattern = "[:alpha:]"),
         ref = str_extract(ref, pattern = "\\d+(\\,\\s\\d+)*"),
         pH = str_extract(solvent.specs, 
                          pattern = "(pH\\s(\\<\\s)*[0-9]+(\\.[0-9]+)*)"),
         pH.range = str_extract(solvent.specs,
                                pattern = "pH\\s[0-9]+\\.[0-9]+\\-[0-9]+(\\.[0-9]+)*"),
         solvent.ratio = str_extract(solvent.specs,
                                     pattern = "[0-9]+(\\.[0-9]+)*\\:[0-9]+(\\.[0-9]+)*"),
         solvent.molarity = str_extract(solvent.specs, 
                                        pattern = pattern.molarity)) %>%
  mutate(pH = ifelse(!is.na(pH.range), NA, pH)) %>%
  mutate(pH = str_extract(pH, pattern = pH.numeric)%>% as.numeric())%>%
  separate(., pH.range, c("pH1", "pH2"), sep = "-") %>%
  mutate(pH1 = str_extract(pH1, pattern = pH.numeric) %>% as.numeric(),
         pH2 = as.numeric(pH2)) %>%
  mutate(pH = ifelse(!is.na(pH1), (pH1+pH2)/2, pH)) %>%
  mutate(pH = ifelse(is.na(pH), 7.0, pH))%>%
  select(-pH1, -pH2) 



# Setting pH to acidic value when the solution contained an acid

sulfuric.acid <- grep(pattern = "[M]+\\s+\\bH2SO4\\b", 
                      x = ri.clean$solvent.specs)

hcl.acid <- grep(pattern = "M+\\s+\\bHCl\\b", 
                 x = ri.clean$solvent.specs)

ri.clean$pH[sulfuric.acid] <- 1
ri.clean$pH[hcl.acid]  <-  1.21

# ------------------------------------------------------------------------------
#                          Reconsider this Section
# ------------------------------------------------------------------------------
# It would be nice to explain what each of these cleanup step is actually doing
# to the guest column and why is needed. I did a quick search in google 
# and was able to find many   of the structures without modifying the original 
# name. Maybe a guest.special.cases column could be more useful - ERD
# 
ri.clean <- ri.clean %>%
  #   mutate(clean.guest = str_replace_all(
  #     string = guest,
  #     pattern = "\\Â·+\\s*[A-Z]*[a-z]*",
  #     replacement = ""
  #   )) %>% 
  # mutate(
  #   clean.guest = str_replace_all(
  #     string = clean.guest,
  #     pattern = "[0-9]HCl",
  #     replacement = "hydrochloride"
  #   )
#   ) %>%
#   mutate(clean.guest = str_replace_all(
#     string = clean.guest,
#     pattern = "\\Î²",
#     replacement = "beta"
#   )) %>%
#   mutate(clean.guest = str_replace(
#     string = clean.guest,
#     pattern = "\\Î±",
#     replacement = "alpha"
#   )) %>%

mutate(guest = str_replace(
  string = guest,
  pattern = "HCl",
  replacement = "hydrochloride"
)) #  %>%
#   mutate(clean.guest = str_replace(
#     string = clean.guest,
#     pattern = "H2SO4",
#     replacement = "sulfonate"
#   )) %>% select(guest, clean.guest) %>% unique() %>% View()
#   tbl_df() %>%
#   bind_cols(ri.clean) %>%
#   
# 
# #fixed file path -AX
# clean_dir       <- "./Output Data/"
# saveRDS(ri_engineered, file = paste0(clean_dir, "02-ri_engineered.RDS"))
# save(ri_engineered, file = paste0(clean_dir, "ri_engineered.RData"))

# Useful Data -------------------------------------------------------------

ri.clean <- ri.clean[!str_detect(ri.clean$ref, "b|c|g|i|j|m"),]
single.complex <- str_detect(ri.clean$host, "1[[:alpha:]]")
ri.clean <- ri.clean[single.complex, ]

# Discriminating based on pH and T
ri.squeaky.clean <- ri.clean[ri.clean$pH > 6.9, ]
ri.squeaky.clean <- ri.squeaky.clean[ri.squeaky.clean$pH < 7.1, ]
ri.squeaky.clean <- ri.squeaky.clean[ri.squeaky.clean$T.K == 298,]


# ------------------------------------------------------------------------------
#                             Formatting Table           
# Reorganizing Table  into the following colunn families vector types
# Host <chr> | Guest <chr> | Solvent Specs <chr>| pH + Thermodynamic Values (dbl)
# | Method+References <chr>
ri.clean <- ri.clean %>%
  select(1:4, 20:21, pH, 5:13, DelCp, DelCp.Uncertainty, methoda, ref, ref.notes) %>%
  mutate_at(vars(7:18), as.numeric) 

ri.squeaky.clean <- ri.squeaky.clean %>%
  select(1:4, 20:21, pH, 5:13, DelCp, DelCp.Uncertainty, methoda, ref, ref.notes) %>%
  mutate_at(vars(7:18), as.numeric) 

# Nice to know:
# Total = 1227
# Alpha = 615
# Beta = 549
# Gamma = 63
# 564 unique

# Squeaky Clean:
# Total = 558
# Alpha = 266
# Beta = 244
# Gamma = 48
# 329 unique

# -----------------------------------------------------------------------------
#                     Save output of script 

saveRDS(ri.clean, file = paste0(bound.dir, "02.ri.clean.RDS"))
save(ri.clean, file = paste0(bound.dir, "02.ri.clean.RData"))
saveRDS(ri.squeaky.clean, file = paste0(bound.dir, "02.ri.squeaky.clean.RDS"))
save(ri.squeaky.clean, file = paste0(bound.dir, "02.ri.squeaky.clean.RData"))
