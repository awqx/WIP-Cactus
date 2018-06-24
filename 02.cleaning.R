# Libraries and Packages --------------------------------------------------

library(data.table)
library(tidyverse)
library(stringr)

# Rekharsky and Inoue -----------------------------------------------------

ri.bound <- readRDS("./dwnld/ri.bound.df.RDS")
# 1367 observations

#     Cosmetic cleaning ---------------------------------------------------

# 1. Rename columns 
#     Note: for some reason, this doesn't work on Windows
#     Keeping this, in case it's important for RStudio on iOS

# ri.bound <- ri.bound %>%
#   rename(DelG = `ÎGÂ°/ kJâmol-1`, 
#          DelH =`ÎHÂ°/ kJâmol-1`,
#          TDelS = `TÎSÂ°/ kJâmol-1`,
#          DelCp = `ÎCpÂ°/ Jâmol-1âK-1`,
#          log.K = `logâK`,
#          method = `methoda`)

colnames(ri.bound) <- c("host", "guest", "solvent", 
                        "T.K", "log.K", 
                        "DelG", "DelH", "TDelS", 
                        "methoda", "ref", "DelCp")

# 2. Splitting columns containing a variable value and its uncertainty
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

# 3. Cleaning strings that contain inconsistent patterns such as more than one space
#     separation and unconventional separation symbols. Converting alpha, beta and
#     gamma symbols to words for easier subset of tables 

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

# 1. Detecting molarity
molarity.regex <- paste0("[0-9]+\\.*[0-9]*[[:space:]]M",
                         "[[:space:]][A-Za-z0-9]*([[:space:]]", 
                         "*(;|\\,|\\+)[[:space:]][0-9]\\.*[0-9]*", 
                         "[[:space:]]M\\s[A-Za-z0-9]*)*")
# Regex for pH
pH.numeric <-  "[0-9]+\\.*[0-9]*"

ri.clean <- ri.clean %>% 
  mutate(# ref.notes= str_extract(ref, pattern = "[:alpha:]"),
         # ref = str_extract(ref, pattern = "\\d+(\\,\\s\\d+)*"),
         pH = str_extract(solvent.specs, 
                          pattern = "(pH\\s(\\<\\s)*[0-9]+(\\.[0-9]+)*)"),
         pH.range = str_extract(solvent.specs,
                                pattern = "pH\\s[0-9]+\\.[0-9]+\\-[0-9]+(\\.[0-9]+)*"),
         solvent.ratio = str_extract(solvent.specs,
                                     pattern = "[0-9]+(\\.[0-9]+)*\\:[0-9]+(\\.[0-9]+)*"),
         solvent.molarity = str_extract(solvent.specs, 
                                        pattern = molarity.regex)) %>%
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
ri.clean <- ri.clean[str_detect(ri.clean$solvent, "H2O$"), ]
ri.clean <- ri.clean[str_detect(ri.clean$solvent, "^H2O"), ]
# removing hydrochloride salts
ri.clean <- ri.clean[!str_detect(ri.clean$guest, pattern = "hydrochloride"), ]

# Discriminating based on pH and T
ri.squeaky.clean <- ri.clean[ri.clean$pH > 6.9, ]
ri.squeaky.clean <- ri.squeaky.clean[ri.squeaky.clean$pH < 7.1, ]
ri.squeaky.clean <- ri.squeaky.clean[ri.squeaky.clean$T.K == 298,]

#     Formatting  ---------------------------------------------------------

# Reorganizing Table  into the following column families vector types
# Host <chr> | Guest <chr> | Solvent Specs <chr>| pH + Thermodynamic Values (dbl)
# | Method+References <chr>
ri.clean <- ri.clean %>%
  mutate(DelG = str_replace(ri.clean$DelG, "(−|\u2212)", "-") %>% as.numeric()) %>%
  mutate(DelH = str_replace(ri.clean$DelH, "(−|\u2212)", "-") %>% as.numeric()) %>%
  mutate(TDelS = str_replace(ri.clean$TDelS, "(−|\u2212)", "-") %>% as.numeric()) %>%
  mutate(DelCp = str_replace(ri.clean$DelCp, "(−|\u2212)", "-") %>% as.numeric()) %>% 
  select(1:4, 19:20, pH, 5:13, DelCp, DelCp.Uncertainty, methoda, ref) %>%
  mutate_at(vars(7:18), as.numeric) 

ri.squeaky.clean <- ri.squeaky.clean %>%
  mutate(DelG = str_replace(ri.squeaky.clean$DelG, "(−|\u2212)", "-") %>% as.numeric()) %>%
  mutate(DelH = str_replace(ri.squeaky.clean$DelH, "(−|\u2212)", "-") %>% as.numeric()) %>%
  mutate(TDelS = str_replace(ri.squeaky.clean$TDelS, "(−|\u2212)", "-") %>% as.numeric()) %>%
  mutate(DelCp = str_replace(ri.squeaky.clean$DelCp, "(−|\u2212)", "-") %>% as.numeric()) %>%
  select(1:4, 19:20, pH, 5:13, DelCp, DelCp.Uncertainty, methoda, ref) %>%
  mutate_at(vars(7:18), as.numeric) 

# Nice to know:
# Total = 1069
# Alpha = 554
# Beta = 483
# Gamma = 32
# 504 unique

# Squeaky Clean:
# Total = 427
# Alpha = 221
# Beta = 185
# Gamma = 21
# 273 unique

#     Save Output ---------------------------------------------------------

# The difference between these files is that squeaky clean is more
# picky about the pH
saveRDS(ri.clean, file = "./dwnld/02.ri.clean.RDS")
saveRDS(ri.squeaky.clean, file = "./dwnld/02.ri.squeaky.clean.RDS")

#####
# Suzuki ------------------------------------------------------------------

suz.bound <- readRDS("./dwnld/suz.bound.df.RDS")
# 1. Selecting only observed data
suz.clean <- suz.bound[-1, c(2, 4, 8)] # Selection of obsd data only

# 2. Renaming columns
colnames(suz.clean) <- c("guest", "DelG.a", "DelG.b") 

# 3. Sorting out misread special characters
suz.clean <- suz.clean %>%
  mutate(DelG.a = str_replace(DelG.a, "â\u0088\u0092", "-")) %>%
  mutate(DelG.a = str_replace(DelG.a, "(â|Â)", "")) %>%
  mutate(DelG.b = str_replace(DelG.b, "â\u0088\u0092", "-")) %>%
  mutate(DelG.b = str_replace(DelG.b, "(â|Â)", "")) %>%
  mutate(guest = str_replace(guest, "â\u0080\u0089", "-")) 

# 4. Detection of rows that ended up as column names
h.ind <- suz.clean$guest %>% str_detect("guest")
suz.clean <- suz.clean[!h.ind, ] 

# 5. Converting the DelG columns to numeric
suz.clean[ , 2:3] <- sapply(suz.clean[, 2:3], as.numeric) 

# 6. Creation of long df
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

# Nice to know:
# Total 320 complexes
# Alpha: 102
# Beta: 218
# Unique guests: 218

#     Save Output ---------------------------------------------------------

saveRDS(suz.clean, file = "./dwnld/02.suz.clean.RDS")

#####
# Combining Datasets ------------------------------------------------------

# renaming for convenience
ri <- ri.squeaky.clean %>%
  mutate(host = str_replace(ri.squeaky.clean$host, "1", "")) %>%
  select(host, guest, DelG) %>%
  mutate(data.source = "rekharsky.inoue")

# 1. Splitting based on CD type
ri.a <- filter(ri, host == "alpha")
ri.b <- filter(ri, host == "beta")
ri.c <- filter(ri, host == "gamma") # Gamma is RI-exclusive

suz.a <- filter(suz.clean, host == "alpha")
suz.b <- filter(suz.clean, host == "beta")

# 2. Collapsing separate data points 
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

comb.dg <- rbind(comb.a, comb.b, ri.c)
comb.dg.nodup <- comb.dg[!duplicated(comb.dg), ] # filtering for unique

# The full combined data contains everything, where replicated points
# are all included
saveRDS(comb.dg, "./dwnld/02.full.combined.ri.suzuki.RDS")
saveRDS(comb.dg.nodup, "./dwnld/02.combined.data.RDS")


# Information about Data
# Total: 615
# Alpha: 241
# Beta: 354
# Gamma: 20
# Unique: 457

#####
# Duplicate Cases ---------------------------------------------------------

# Just a quick visualization of differences between duplicated cases
# Not necessary for model-building

ri.clean.sub <- ri.clean %>% select(., guest, DelG, host, pH, T.K) %>%
  mutate(data.source = "rekharsky.inoue") %>%
  mutate(host = str_replace(host, "1", ""))
suz.clean.add <- suz.clean %>% mutate(pH = 7) %>%
  mutate(T.K = 298)

ri.clean.sub.a <- ri.clean.sub %>% filter(host == "alpha")
suz.clean.add.a <- suz.clean.add %>% filter(host == "alpha")
all.a <- rbind(suz.clean.add.a, ri.clean.sub.a)
dup.a.guest <- all.a$guest[duplicated(all.a$guest)] %>% unique()
dup.a <- all.a[all.a$guest %in% dup.a.guest, ] %>% group_by(guest)

dup.a.tk <- dup.a %>% select(-pH)
dup.a.guest2 <- dup.a.tk[duplicated(dup.a.tk$guest), ]$guest %>% unique()
dup.a.tk <- dup.a.tk[dup.a.tk$guest %in% dup.a.guest2, ] %>% group_by(guest)
ggplot(dup.a.tk, aes(x = T.K, y = DelG, # color = pH, 
                  shape = data.source, 
                  group = guest)) + 
  geom_point() + 
  geom_line() 
  # scale_colour_gradient2(low = "orangered2", mid = "seagreen2",
  #                      high = "slateblue", midpoint = 7)

ggplot(dup.a, aes(x = pH, y = DelG, # color = pH, 
                  shape = data.source, 
                  group = guest)) + 
  geom_point() + 
  geom_line() 


# Graphing ----------------------------------------------------------------

# Some analysis of different conditions
# ri.clean before questionable section
ri.clean$DelG <- as.numeric(ri.clean$DelG)
ggplot(ri.clean, aes(x = DelG)) + 
  geom_histogram(stat = "bin") + 
  facet_wrap(~solvent)

ri.dmf <- ri.clean %>% 
  filter(solvent == "DMF") %>% 
  select(., host, guest, solvent, DelG) 
ri.h2o <- ri.clean %>%
  filter(solvent == "H2O")  %>% 
  select(., host, guest, solvent, DelG) 
ri.dmfh2o.guests <- inner_join(ri.dmf, ri.h20, by = c("host", "guest")) %>%
  .$guest
ri.dmfh2o <- rbind(
  ri.dmf %>% filter(guest %in% ri.dmfh2o.guests), 
  ri.h2o %>% filter(guest %in% ri.dmfh2o.guests)
)
ggplot(ri.dmfh2o, aes(x = DelG)) + 
  geom_histogram(stat = "bin") + 
  facet_wrap(~solvent)

ri.dmfh2o.wide <- inner_join(ri.dmf, ri.h20, by = c("host", "guest")) %>%
  as.data.frame()
ggplot(ri.dmfh2o.wide, aes(x = DelG, y = dG.h20)) + 
  geom_point() + 
  labs(x = "DMF", y = "H2O") + 
  geom_abline(intercept = 0, slope = 1) + 
  coord_fixed() + 
  theme_bw()
