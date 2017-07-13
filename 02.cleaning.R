library(data.table)
library(tidyverse)
library(stringr)

# Rekharsky and Inoue -----------------------------------------------------

# The following reading will vary based on user preference
ri.bound <- readRDS("./bound/ri.bound.df.RDS")
# Rename columns 
#     Note: for some reason, this doesn't work on Windows
#     Keeping this, in case it's important for RStudio on iOS

# ri.bound <- ri.bound %>%
#   rename(DelG = `ÎGÂ°/ kJâmol-1`, 
#          DelH =`ÎHÂ°/ kJâmol-1`,
#          TDelS = `TÎSÂ°/ kJâmol-1`,
#          DelCp = `ÎCpÂ°/ Jâmol-1âK-1`,
#          log.K = `logâK`,
#          method = `methoda`)


# Splitting columns containing a variable value and its uncertainty
# Occasionally the naming of the columns will get messed up
# In this case, you would need colnames(ri.bound)[3], [5], [6], [7], and [11]
colnames(ri.bound) <- c("host", "guest", "solvent", 
                        "T.K", "log.K", 
                        "DelG", "DelH", "TDelS", 
                        "methoda", "ref", "DelCp")

ri.clean <- ri.bound %>%
  separate(solvent, c("solvent","solvent.specs"),
           sep = "(?=\\s*\\()", extra = "merge", fill = "right") %>%
  separate(log.K, c("log.K", "log.K.Uncertainty"),
           sep = "\\s\\±\\s", extra = "merge", fill = "right") %>% 
  separate(DelG, c("DelG", "DelG.Uncertainty"),
           sep = "\\s\\±\\s",extra = "merge", fill = "right") %>% 
  separate(DelH, c("DelH", "DelH.Uncertainty"),
           sep = "\\s\\±\\s", extra = "merge", fill = "right") %>%
  separate(TDelS, c("TDelS", "TDelS.Uncertainty"),
           sep = "\\s\\±\\s",extra = "merge", fill = "right")  %>%
  separate(DelCp, c("DelCp", "DelCp.Uncertainty"),
          sep = "\\s\\±\\s",extra = "merge", fill = "right")

# Cleaning strings that contain inconsistent patterns such as more than one space
# separation and unconventional separation symbols. Converting alpha, beta and
# gamma symbols to words for easier subset of tables 

ri.clean <- ri.clean %>%
  lapply(str_replace_all, pattern = "\\â", 
         replacement = "-") %>%
  lapply(str_replace_all, pattern = "\\â+", 
         replacement = " ") %>%
  lapply(str_replace_all, pattern = "\\s+", 
         replacement = " ") %>%
  lapply(str_replace_all, pattern = "Â·", 
         replacement = " ") %>% 
  as_tibble() %>%
  mutate(host = str_replace(host, pattern =  "\\Î±|\u03b1",
                            replacement = "alpha" )) %>%
  mutate(host = str_replace(host, pattern = "\\Î²|\\B\\s+H\\+|\u03b2", 
                            replacement = "beta")) %>%
  mutate(host = str_replace(host, pattern = "\\Î³|\u03b3", 
                            replacement = "gamma"))


#     pH Imputation -----------------------------------------------------------

# Remove pH string and parenthesis to convert the pH Column to numerical. Assume 
# that columns with no value have a pH of 7.0. pH values with a range of 
# phValue1-phValue2 will be assigned an average pH value ((phValue1+phValue2)/2)
# ph's with an inequality (<,>) will be set at the value given. Guests with 
# the molarity of the acid given in the solvent composition will be set
# at the theoretical value of the solution

# Creating regex for molarity detection
pattern.molarity <- gsub("\n", replacement = "", x = "(([0-9]*\\.*[0-9]*\\s+[M]
                         \\s+[A-Za-z]*[A-Za-z0-9]*\\s*[a-z]*)(\\;*\\,*\\+*\\s+[0-9]+\\.*[0-9]*\\s+[M]
                         \\s+[A-Za-z]*[A-Za-z0-9]*[a-z]*)*)")
# Regex for pH
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
  mutate(pH = str_extract(pH, pattern = pH.numeric) %>% as.numeric()) %>%
  separate(., pH.range, c("pH1", "pH2"), sep = "-") %>%
  mutate(pH1 = str_extract(pH1, pattern = pH.numeric) %>% as.numeric(),
         pH2 = as.numeric(pH2)) %>%
  mutate(pH = ifelse(!is.na(pH1), (pH1+pH2)/2, pH)) %>%
  mutate(pH = ifelse(is.na(pH), 7.0, pH)) %>%
  select(-pH1, -pH2) 

# Setting pH to acidic value when the solution contained an acid

sulfuric.acid <- grep(pattern = "[M]+\\s+\\bH2SO4\\b", 
                      x = ri.clean$solvent.specs)

hcl.acid <- grep(pattern = "M+\\s+\\bHCl\\b", 
                 x = ri.clean$solvent.specs)

ri.clean$pH[sulfuric.acid] <- 1
ri.clean$pH[hcl.acid]  <-  1.21


#     Questionable Section ------------------------------------------------


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

#     Useful Data ---------------------------------------------------------

# These annotations include data that is questionable or unreliable
ri.clean <- ri.clean[!str_detect(ri.clean$ref, "b|c|g|i|j|m"),]
# Filtering for 1:1 complexes
ri.clean <- ri.clean[str_detect(ri.clean$host, "^1[[:alpha:]]"), ]

# Discriminating based on pH and T
ri.squeaky.clean <- ri.clean[ri.clean$pH > 6.9, ]
ri.squeaky.clean <- ri.squeaky.clean[ri.squeaky.clean$pH < 7.1, ]
ri.squeaky.clean <- ri.squeaky.clean[ri.squeaky.clean$T.K == 298,]

#     Formatting  ---------------------------------------------------------

# Reorganizing Table  into the following colunn families vector types
# Host <chr> | Guest <chr> | Solvent Specs <chr>| pH + Thermodynamic Values (dbl)
# | Method+References <chr>
ri.clean <- ri.clean %>%
  mutate(DelG = str_replace(ri.clean$DelG, "(−|\u2212)", "-") %>% as.numeric()) %>%
  mutate(DelH = str_replace(ri.clean$DelH, "(−|\u2212)", "-") %>% as.numeric()) %>%
  mutate(TDelS = str_replace(ri.clean$TDelS, "(−|\u2212)", "-") %>% as.numeric()) %>%
  mutate(DelCp = str_replace(ri.clean$DelCp, "(−|\u2212)", "-") %>% as.numeric()) %>%
  select(1:4, 20:21, pH, 5:13, DelCp, DelCp.Uncertainty, methoda, ref, ref.notes) %>%
  mutate_at(vars(7:18), as.numeric) 

ri.squeaky.clean <- ri.squeaky.clean %>%
  mutate(DelG = str_replace(ri.squeaky.clean$DelG, "(−|\u2212)", "-") %>% as.numeric()) %>%
  mutate(DelH = str_replace(ri.squeaky.clean$DelH, "(−|\u2212)", "-") %>% as.numeric()) %>%
  mutate(TDelS = str_replace(ri.squeaky.clean$TDelS, "(−|\u2212)", "-") %>% as.numeric()) %>%
  mutate(DelCp = str_replace(ri.squeaky.clean$DelCp, "(−|\u2212)", "-") %>% as.numeric()) %>%
  select(1:4, 20:21, pH, 5:13, DelCp, DelCp.Uncertainty, methoda, ref, ref.notes) %>%
  mutate_at(vars(7:18), as.numeric) 

# Nice to know:
# Total = 1225
# Alpha = 615
# Beta = 547
# Gamma = 63
# 564 unique

# Squeaky Clean:
# Total = 558
# Alpha = 266
# Beta = 244
# Gamma = 48
# 329 unique

#     Save Output ---------------------------------------------------------
saveRDS(ri.clean, file = "./bound/02.ri.clean.RDS")

saveRDS(ri.squeaky.clean, file = "./bound/02.ri.squeaky.clean.RDS")
save(ri.squeaky.clean, file = "./bound/02.ri.squeaky.clean.RData")

#####
# Suzuki ------------------------------------------------------------------

suz.bound <- readRDS("./bound/suz.bound.df.RDS")
suz.clean <- suz.bound[-1, c(2, 4, 8)] # Selection of obsd data only
colnames(suz.clean) <- c("guest", "DelG.a", "DelG.b") # Easy renaming

# Sorting out characters that were read weirdly 
suz.clean <- suz.clean %>%
  mutate(DelG.a = str_replace(DelG.a, "â\u0088\u0092", "-")) %>%
  mutate(DelG.a = str_replace(DelG.a, "(â|Â)", "")) %>%
  mutate(DelG.b = str_replace(DelG.b, "â\u0088\u0092", "-")) %>%
  mutate(DelG.b = str_replace(DelG.b, "(â|Â)", "")) %>%
  mutate(guest = str_replace(guest, "â\u0080\u0089", "-")) 

# Detection of rows that ended up as column names
h.ind <- suz.clean$guest %>% str_detect("guest")
suz.clean <- suz.clean[!h.ind, ] 

# Converting the DelG columns to numeric
suz.clean[ , 2:3] <- sapply(suz.clean[, 2:3], as.numeric) 

# Separating the DelGs, creating a long df
suz.clean.dg.alpha <- suz.clean %>% 
  filter(!is.na(DelG.a)) %>% 
  select(guest, DelG.a) %>% 
  rename(DelG = DelG.a) %>% 
  mutate(host = "alpha") %>%
  mutate(data.source = "suzuki")
suz.clean.dg.beta <- suz.clean %>%  
  select(guest, DelG.b) %>%
  rename(DelG = DelG.b) %>% 
  mutate(host = "beta") %>%
  mutate(data.source = "suzuki")

suz.clean <- rbind(suz.clean.dg.alpha, suz.clean.dg.beta)

#     Save Output ---------------------------------------------------------
save(suz.clean, file = "./bound/02.suz.clean.RData")
saveRDS(suz.clean, file = "./bound/02.suz.clean.RDS")

#####
# Combining Datasets ------------------------------------------------------

ri <- ri.squeaky.clean %>%
  mutate(host = str_replace(ri.squeaky.clean$host, "1", "")) %>%
  select(host, guest, DelG) %>%
  mutate(data.source = "rekharsky.inoue")

# Splitting based on CD type
ri.a <- filter(ri, host == "alpha")
ri.b <- filter(ri, host == "beta")

suz.a <- filter(suz.clean, host == "alpha")
suz.b <- filter(suz.clean, host == "beta")

# Collapsing separate data points 
#     Adapted from GSee at stackoverflow.com/questions/12884695
#     As well as Sven Hohenstein at stackoverflow.com/questions/20854615
comb.a <- rbind(ri.a, suz.a) %>%
  data.table(., key = "guest")
comb.a <- comb.a[, list(host = host, DelG = mean(DelG),
                        data.source = paste(unlist(data.source),
                                            collapse = ", ")),
                 by = guest] 

comb.b <- rbind(ri.b, suz.b) %>%
  data.table(., key = "guest")
comb.b <- comb.b[ , list(host = host, DelG = mean(DelG), 
                         data.source = paste(unlist(data.source),
                                             collapse = ", ")), 
                  by = guest] 

comb.dg <- rbind(comb.a, comb.b)
saveRDS(comb.dg, "./bound/combined ri and suzuki.RDS")
