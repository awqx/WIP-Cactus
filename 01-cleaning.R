# Libraries and Packages --------------------------------------------------

packages <- c("data.table", "stringr", "tidyverse")
lapply(packages, require, character.only = T)

if(!dir.exists("cleaning")) dir.create("cleaning")

# Rekharsky and Inoue -----------------------------------------------------

ri_df <- readRDS("dwnld/ri-df.RDS") # 1367 observations

# split columns with error
ri_clean <- ri_df %>%
  separate(
    solvent, c("solvent","solvent_specs"),
    sep = "(?=\\s*\\()", extra = "merge", fill = "right") %>%
  separate(
    dG, c("dG", "dG_error"),
    sep = "\\s\\±\\s",extra = "merge", fill = "right") %>%
  mutate(temp = as.numeric(temp))

# cleaning strings with unusual separators or double spaces
ri_clean <- ri_clean %>%
  lapply(str_replace_all, pattern = "\\â", 
         replacement = "-") %>%
  lapply(str_replace_all, pattern = "\\â+", 
         replacement = " ") %>%
  lapply(str_replace_all, pattern = "\\s+", 
         replacement = " ") %>%
  lapply(str_replace_all, pattern = "Â·", 
         replacement = " ") %>% 
  as_tibble() 

# pH
# Remove pH string and parenthesis to convert the pH Column to numerical. 
# pH values averaged across ranges 
# ph's with an inequality (<,>) will be set at the value given. 
molarity_regex <- paste0("[0-9]+\\.*[0-9]*[[:space:]]M",
                         "[[:space:]][A-Za-z0-9]*([[:space:]]", 
                         "*(;|\\,|\\+)[[:space:]][0-9]\\.*[0-9]*", 
                         "[[:space:]]M\\s[A-Za-z0-9]*)*")
pH_regex <-  "[0-9]+\\.*[0-9]*"

ri_clean <- ri_clean %>% 
  mutate(
    pH = str_extract(
      solvent_specs, 
      pattern = "(pH\\s(\\<\\s)*[0-9]+(\\.[0-9]+)*)"),
    pH_range = str_extract(
      solvent_specs,
      pattern = "pH\\s[0-9]+\\.[0-9]+\\-[0-9]+(\\.[0-9]+)*"),
    solvent.ratio = str_extract(
      solvent_specs,
      pattern = "[0-9]+(\\.[0-9]+)*\\:[0-9]+(\\.[0-9]+)*"),
    molarity = str_extract(
      solvent_specs, 
      pattern = molarity_regex)) %>%
  mutate(pH = ifelse(!is.na(pH_range), NA, pH)) %>%
  mutate(pH = as.numeric(str_extract(pH, pattern = pH_regex))) %>%
  separate(., pH_range, c("pH1", "pH2"), sep = "-") %>%
  mutate(
    pH1 = as.numeric(str_extract(pH1, pattern = pH_regex)),
    pH2 = as.numeric(pH2)) %>%
  mutate(pH = ifelse(!is.na(pH1), (pH1 + pH2)/2, pH)) %>%
  mutate(pH = ifelse(is.na(pH), 7.0, pH)) %>%
  select(-pH1, -pH2) 

# Setting pH to acidic value when the solution contained an acid
sulfuric_acid <- grep(pattern = "[M]+\\s+\\bH2SO4\\b", 
                      x = ri_clean$solvent_specs)
hcl_acid <- grep(pattern = "M+\\s+\\bHCl\\b", 
                 x = ri_clean$solvent_specs)
ri_clean$pH[sulfuric_acid] <- 1
ri_clean$pH[hcl_acid]  <-  1.21

# filter by H20 solvent, pH between 6.8-7.5, temperature
ri_clean <- ri_clean %>%
  filter(
    is.na(molarity), 
    solvent == "H2O", 
    temp == 298, 
    pH > 6.8 & pH < 7.5, 
    !str_detect(ref, "b|c|g|i|j|m"), # unreliable references
    !str_detect(guest, "HCl") # salts
    ) %>%
  select(host, guest, dG) %>%
  data.frame()

# changing dG to numeric
# the dash in the table is unusual and is not being properly cast 
# so it must be replaced
ri_clean <- ri_clean %>%
  mutate(dG = str_replace(dG, "\u2212", "-")) %>%
  mutate(dG = as.numeric(dG))
saveRDS(ri_clean, "cleaning/ri-clean.RDS")


# Suzuki ------------------------------------------------------------------

suz_clean <- readRDS("dwnld/suzuki-df.RDS") %>%
  mutate(dG = str_replace(dG, "\u2212", "-")) %>%
  mutate(dG = as.numeric(dG)) %>%
  mutate(guest = str_replace(guest, "\u2009", " ")) %>%
  filter(!is.na(dG))
saveRDS(suz_clean, "cleaning/suzuki-clean.RDS")

# Combining Datasets ------------------------------------------------------

ri   <- mutate(ri_clean, paper = "ri")
ri_a <- filter(ri, host == "alpha")
ri_b <- filter(ri, host == "beta")
ri_c <- filter(ri, host == "gamma")
suz   <- mutate(suz_clean, paper = "suzuki")
suz_a <- filter(suz, host == "alpha")
suz_b <- filter(suz, host == "beta")

# Collapsing separate data points 
#     Adapted from GSee at stackoverflow.com/questions/12884695
#     As well as Sven Hohenstein at stackoverflow.com/questions/20854615
combined_a <- rbind(ri_a, suz_a) %>%
  data.table(., key = "guest") %>%
  .[, list(
    host = host, 
    dG = mean(dG), 
    paper = paste(unlist(paper), collapse = ", ")), 
    by = guest
    ] %>%
  filter(!duplicated(guest))

combined_b <- rbind(ri_b, suz_b) %>%
  data.table(., key = "guest") %>%
  .[, list(
    host = host, 
    dG = mean(dG), 
    paper = paste(unlist(paper), collapse = ", ")), 
    by = guest
    ] %>%
  filter(!duplicated(guest))

combined <- rbind(combined_a, combined_b, ri_c) 
saveRDS(combined, "cleaning/ri-suzuki.RDS")

# Information about Data
# Total: 727
# Alpha: 294
# Beta: 413
# Gamma: 20
# Unique: 528

# Connors (gamma-CD) ------------------------------------------------------

# UPDATE: modeling on this data did not yield significant results
# Check older versions of GitHub for reference
# Connors, K.A. Feb 24, 1995, School of Pharmacy, University of Wisconsin.
# Population characteristics of cyclodextrin complex stabilities in aqueous
# solution
# File found from University of Wisconsin archives
# Thanks to Debra King and Joni Mitchell 


# Singh (gamma-CD) --------------------------------------------------------

singh <- readRDS("dwnld/singh.RDS") %>%
  select(guest, gamma) %>%
  mutate(host = "gamma", guest = tolower(guest)) %>%
  rename(dG = gamma) %>%
  select(guest, host, dG)
saveRDS(singh, "cleaning/singh.RDS")