# download guest molecules from NCI/CADD
# https://cactus.nci.nih.gov/chemical/structure
source("helpers/cactus.R")

# Cleaning data -----------------------------------------------------------

if(!dir.exists("molecules")) {
  dir.create("molecules")
  dir.create("molecules/alphaCD")
  dir.create("molecules/betaCD")
  dir.create("molecules/gammaCD")  
}

ri_suzuki <-readRDS("cleaning/ri-suzuki.RDS") %>% select(-paper)
singh <- readRDS("cleaning/singh.RDS") 
dataset <- rbind(ri_suzuki, singh) %>%
  filter(!is.na(dG))

# clean anionic species
dataset <- dataset %>%
  rowwise() %>%
  mutate(
    guest_charge = str_extract(guest, pattern = "anion|monoanion|dianion"),
    guest = gsub(
      x = guest,
      pattern = "\\(anion\\)|\\(monoanion\\)|\\(dianion\\)",
      replacement = ""
      )
    ) %>% 
  # fix unconventional apostrophes
  mutate(guest = str_replace(guest, "(\\â | â\u0080\u0098)", "\\'"))

# 42 molecules require manual replacement of names
fixed_names <- read.csv("cleaning/fixed-names.csv", header = T) %>%
  mutate(
    pattern = as.character(pattern), 
    replacement = as.character(replacement)
    )

# Converting the table to regex for str_replace
pattern_reg <- fixed_names %>%
  mutate(pattern = paste0(make_regex(fixed_names$pattern), "$")) %>%
  mutate(replacement = make_regex(fixed_names$replacement))
for(i in 1:nrow(pattern_reg)) {
  dataset$guest <- str_replace(
    dataset$guest, 
    pattern = pattern_reg[i, 1], 
    replacement = pattern_reg[i, 2])
}
# "-acid"has problems reading from the .csv
dataset$guest <- str_replace(
  dataset$guest, pattern = "-acid", 
  replacement = " acid")

# Problem with beta replacement
dataset <- dataset %>%
  mutate(guest = str_replace(guest, "\u03b2", "beta")) %>%
  mutate(
    guest = str_replace_all(
      guest,
      c(
        "4-nitrophenyl-beta-d-glucoside" = 
          "4-nitrophenyl beta-d-glucoside",
        "4-nitrophenyl-beta-d-xyloside" = 
          "(2S,3R,4S,5R)-2-(4-nitrophenoxy)oxane-3,4,5-triol",
        "4-nitrophenyl-beta-d-galactoside" = 
          "4-Nitrophenylgalactoside",
        "4-nitrophenyl-beta-d-glucosamide" =
          paste0("N-[(2R,3R,4R,5S,6R)-4,5-dihydroxy-6-(hydroxymethyl)-",
          "2-[(4-nitrophenyl)methoxy]oxan-3-yl]acetamide")
      )
    )
  )

# Reading dataset for guest molecules specific to host

beta.guest <- dataset %>% filter(host == "beta") %>% .$guest
gamma.guest <- dataset %>% filter(host == "gamma") %>% 
  .$guest %>% unique()

# Cactus Download ---------------------------------------------------------

# AlphaCD
#     Creates a dataframe of results; SDFs downloaded to disk
results_a <-
  do.call(
    rbind,
    lapply(
      filter(dataset, host == "alpha")$guest,
      download_sdf,
      path = "./molecules/alphaCD",
      chemical_format = "SDF"
    )
  ) %>% 
  mutate(host = "alpha")

# Beta-CD
results_b <-
  do.call(
    rbind,
    lapply(
      filter(dataset, host == "beta")$guest,
      download_sdf,
      path = "./molecules/betaCD",
      chemical_format = "SDF"
    )
  ) %>% 
  mutate(host = "beta")

# Gamma-CD
results_c <-
  do.call(
    rbind,
    lapply(
      filter(dataset, host == "gamma")$guest,
      download_sdf,
      path = "./molecules/gammaCD",
      chemical_format = "SDF"
    )
  ) %>% 
  mutate(host = "gamma")

results <- rbind(results_a, results_b, results_c)
saveRDS(results, "cleaning/cactus-results.RDS")

# failed_sdfs <- results %>% 
#   filter(downloaded == "no")