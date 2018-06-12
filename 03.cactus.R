# Libraries and Packages --------------------------------------------------

library(tidyverse)
library(stringr)
library(XML)
library(RCurl)

# Consider using WebChem?

# Functions ---------------------------------------------------------------

# write.sdf <- function(name, loc) {
#   path <- paste0(loc, "/", name, ".SDF")
#   write.table(file = path)
# }

make.regex <- function(string) {
  new.string <- str_replace(string, pattern = "\\(", 
                            replacement = "\\\\(")
  new.string <- str_replace(new.string, pattern = "\\)", 
                            replacement = "\\\\)")
  # new.string <- str_replace(new.string, pattern = "\\Î²|\\B\\s+H\\+|\u03b2", # alternatives: \\Î²|\\B\\s+H\\+
  #                           replacement = "\\(Î²|\u03b2\\)")
  return(new.string)
}

download.cactus.results <- function(guest, path, chemical.format) {
  report <- tryCatch({
    destfile       <- paste0(path, "/", guest, ".SDF")
    # Chemical format must be parsed to match all the outputs from NCI cactus
    Guest.URL      <- unlist(lapply(guest, URLencode, reserved = T))
    URL            <- paste0(
      "https://cactus.nci.nih.gov/chemical/structure/",
      Guest.URL,
      "/",
      chemical.format
    )
    Map(download.file, url = URL, destfile = destfile)
    data.frame(
      guest = guest,
      downloaded = "yes",
      warning = "no",
      error = "no"
    )
  },
  warning = function(warn) {
    message("Warning: either URL error or already existing directory.")
    destfile       <- paste0(path, "/", guest, ".SDF")
    Guest.URL      <- unlist(lapply(guest, URLencode, reserved = T))
    URL            <- paste0(
      "https://cactus.nci.nih.gov/chemical/structure/",
      Guest.URL, "/", chemical.format
    )
    Map(download.file, url = URL, destfile = destfile)
    data.frame(
      guest = guest,
      downloaded = "yes",
      warning = "yes",
      error = "no"
    )
  },
  error = function(err) {
    message("An error occurred")
    data.frame(
      guest = guest,
      downloaded = "no",
      warning = "yes",
      error = "yes"
    )
    
  },
  finally = {
    message("Chemical name processed")
  })
  return(report)
}


# Cactus Download (1) -----------------------------------------------------

dir.create("./molecules")
dir.create("./molecules/alphaCD")
dir.create("./molecules/betaCD")
dir.create("./molecules/gammaCD")
dataset <-readRDS("./dwnld/02.combined.data.RDS")

# Reading dataset for guest molecules specific to host
alpha.guest <- dataset %>% filter(host == "alpha") %>% .$guest
beta.guest <- dataset %>% filter(host == "beta") %>% .$guest
gamma.guest <- dataset %>% filter(host == "gamma") %>% .$guest

# AlphaCD
#     Creates a dataframe of results; SDFs downloaded to disk
results.alpha <-
  do.call(
    rbind,
    lapply(
      alpha.guest,
      download.cactus.results,
      path = "./molecules/alphaCD",
      chemical.format = "SDF"
    )
  ) %>% mutate(host = "alpha")

# Beta-CD
results.beta <-
  do.call(
    rbind,
    lapply(
      beta.guest,
      download.cactus.results,
      path = "./molecules/betaCD",
      chemical.format = "SDF"
    )
  ) %>% mutate(host = "beta")

# Gamma-CD
results.gamma <-
  do.call(
    rbind,
    lapply(
      gamma.guest,
      download.cactus.results,
      path = "./molecules/gammaCD",
      chemical.format = "SDF"
    )
  ) %>% mutate(host = "gamma")

results.all <- rbind(results.alpha, results.beta, results.gamma)
saveRDS(results.all, "./molecules/cactus.dwnld.results.RDS")

# Failed SDFs -------------------------------------------------------------

# 32 undownloaded alpha guests
alpha.fail <- results.alpha %>% 
  filter(downloaded == "no") %>%
  select(guest)

# 34 undownloaded beta guests
beta.fail <- results.beta %>% 
  filter(downloaded == "no") %>%
  select(guest)

# 3 undownloaded gamma guests
gamma.fail <- results.gamma %>% 
  filter(downloaded == "no") %>%
  select(guest)

problem.sdf <- results.all %>% 
  filter(downloaded == "no") %>%
  select(guest, host)

# Cleaning Issues ---------------------------------------------------------

# 1. Presence of anionic species - 18
wip.sdf <- problem.sdf %>%
  rowwise() %>%
  mutate(guest.charge = str_extract(guest, pattern = "anion|monoanion|dianion"),
         guest = gsub(x = guest,
                      pattern = "\\(anion\\)|\\(monoanion\\)|\\(dianion\\)",
                      replacement = "")) # %>%

# 2. Unconventional apostrophe \\â or â\u0080\u0098  - 2 # 
wip.sdf$guest <- wip.sdf$guest %>%
  str_replace("(\\â | â\u0080\u0098)", "\\'") 

# 3. Requires manual replacement - 42

#     Manual Replacement of Names -----------------------------------------

# Reassign name
#     Commented sections are compounds that may appear in another cleaning
#     sample, but are not present with the current configuration 
fixed.sdf <- data.frame(
  pattern = "biebricht scarlet",
  replacement = "biebrich scarlet",
  # pattern = "(S)-1,1â-binaphthyl-2,2â-dicarboxylic acid",
  # replacement = "1,1â-binaphthyl-2,2â-dicarboxylic acid",
  # pattern = "4-([(4-hydroxyphenyl)azo]benzoate",
  # replacement = "4-[(4-hydroxyphenyl)azo]benzoate",
  pattern = "3-(aminomethyl)proxyl",
  replacement = "3-(aminomethyl)-proxyl",
  pattern = "3-carbamoylproxyl",
  replacement = "3-carbamoyl-proxyl",
  # pattern = "l-Î±-O-benzylglycerol",
  # replacement = "1-o-benzyl-rac-glycerol",
  pattern = "4-nitrophenyl-Î²-d-glucoside",
  replacement = "4-nitrophenyl-beta-D-glucopyranoside",
  # pattern = "bromodiphenhydramine hydrochloride",  # Consider permanent removal
  # replacement = "Bromdiphenhydramine hydrochloride",
  pattern = "diammine(1,1-cyclobutane- dicarboxylato)platinum(II)" ,
  replacement = "cis-diammine(1,1-cyclobutanedicarboxylato)platinum(II)",
  # pattern = "(1R,2S)-(-)-ephedrine",
  # replacement = "L-Ephedrine",
  # pattern = "(1S,2R)-(+)-ephedrine",
  # replacement = "D-Ephedrine",
  pattern = "sulfoisomidine",
  replacement = "sulfisomidine",
  # pattern = "diphenyhydramine hydrochloride", # Consider permanent removal
  # replacement = "diphenhydramine hydrochloride",
  pattern = "trans,trans-2,4- hexadienedioic (muconic) acid",
  replacement = "muconic acid",
  # pattern = "hexyl-Î²-d-glucopyranoside",
  # replacement = "hexyl beta-d-glucopyranoside",
  # pattern = "2-thiophenobarbital",
  # replacement = "5-Ethyl-5-phenyl-2-thioxodihydro-4,6(1H,5H)-pyrimidinedione",
  # pattern = "4-(hydroxyphenethyl)ammonium" ,
  # replacement = "2-(4-hydroxyphenyl)ethylazanium",
  # pattern = "2hydrochloride", # Consider permanent removal
  # replacement = "dihydrochloride",
  # pattern = 'd-mandelate', replacement = '(R)-Mandelate',
  # pattern = 'l-mandelate', replacement = '(S)-Mandelate',
  # pattern = 'methapyriline hydrochloride', replacement = 'Methoxylene', # Consider permanent removal
  # pattern = 'methyl red (cation, protonated)', replacement = 'methyl red',
  pattern = '4-nitrophenyl-ß-d-galactoside', 
  replacement = '4-Nitrophenylgalactoside',
  pattern = '4-nitrophenyl-ß-d-glucosamide', 
  replacement = 'N-[(2R,3R,4R,5S,6R)-4,5-dihydroxy-6-(hydroxymethyl)-2-[(4-nitrophenyl)methoxy]oxan-3-yl]acetamide',
  pattern = '4-nitrophenyl-ß-d-xyloside', 
  replacement = '(2S,3R,4S,5R)-2-(4-nitrophenoxy)oxane-3,4,5-triol',
  # pattern = '(Â±)-norphenylephrine', 
  # replacement = '3-(2-amino-1-hydroxyethyl)phenol',
  pattern = 'pentakis(ethyleneglycol) monohexyl ether', 
  replacement = '2-hexoxyethanol',
  # pattern = '(R)-(-)-phenylephrine', 
  # replacement = '[(2R)-2-hydroxy-2-(3-hydroxyphenyl)ethyl]-methylazanium',
  # pattern = "phenyl-Î²-d-glucopyranoside",
  # replacement = "phenyl-beta-d-glucopyranoside", 
  # pattern = "(Â±)-octopamine",
  # replacement = "4-(2-amino-1-hydroxyethyl)phenol", 
  # pattern = 'd-phenyltrifluoroethanol', 
  # replacement = ' 2,2,2-trifluoro-1-phenylethanol',
  pattern = 'pyrilammonium maleate',
  replacement = 'pyrilamine maleate',
  # pattern = 'terfenadine hydrochloride', # Consider permanent removal
  # replacement = '1-(4-tert-butylphenyl)-4-[4-[hydroxy(diphenyl)methyl]piperidin-1-yl]butan-1-ol hydrochloride',
  # pattern = 'triiodide (I3-)', replacement = 'triiodide', # Consider permanent removal
  pattern = 'l-tryprophan', replacement = 'l-tryptophan',
  pattern = 'Tyr-Gly-Gly-Phe-Leu', replacement = 'Leucine enkephalin',
  # pattern = '1-adamantaneammonium', replacement = '1-Adamantanaminium',
  # pattern = '(Â±)-anisodamine', replacement = 'anisodamine',
  # pattern = '(-)-anisodamine', replacement = 'anisodamine',
  # pattern = '(-)-anisodamine HBr', 
  # replacement = '[(3S,6S)-6-hydroxy-8-methyl-8-azabicyclo[3.2.1]octan-3-yl] (2R)-3-hydroxy-2-phenylpropanoate',
  # pattern = '(-)-anisodine HBr', replacement = 'Anisodine hydrobromide',
  # pattern = '(Â±)-atropine H2SO4', replacement = 'HOBWAPHTEJGALG-UHFFFAOYSA-N',
  # pattern = '(4Z,15Z)-bilirubin-IXa', 
  # replacement = 'BPYKTIZUTYGOLE-KDUUSRDASA-N',
  pattern = 'dansyl-l-hydroxyplorine', 
  replacement = 'dansyl-l-hydroxyproline',
  # pattern = "methapyriline hydrochloride", # Consider permanent removal
  # replacement = "methapyrilene hydrochloride", 
  # pattern = 'dicumarol', 
  # replacement = '4-hydroxy-3-[(4-hydroxy-2-oxochromen-3-yl)methyl]chromen-2-one',
  # pattern = 'ethylbis(coumacetate)', 
  # replacement = 'Ethyl biscoumacetate',
  # pattern = 'ethylthiobarbituic acid', 
  # replacement = '5-Ethyl-2-thioxodihydropyrimidine-4,6(1H,5H)-dione',
  # pattern = '4-[(4-hydroxy-1-naphthyl)azo]- naphthalene-1-sulfonate', 
  # replacement = '4-[(4-hydroxy-1-naphthyl)azo]- naphthalene-1-sulfonic acid',
  pattern = 'maprotilin', replacement = 'maprotiline',
  # pattern = 'naproxenate', replacement = 'naproxen',
  pattern = 'nortriptylin', replacement = ' Aventyl',
  pattern = 'protriptylin', replacement = 'protriptyline',
  # pattern = '(-)-scopolamine HBr', 
  # replacement = '(-)-scopolamine hydrobromide',
  pattern = 'sulfasnilamide', 
  replacement = '4-aminobenzenesulfonamide',
  pattern = 'sulfathidole', 
  replacement = 'sulfaethidole',
  # pattern = 'thiophenobarbital', 
  # replacement = '5-(1,1,2,2,2-pentadeuterioethyl)-5-phenyl-2-sulfanylidene-1,3-diazinane-4,6-dione',
  # pattern = '6-(p-toluidinyl)-2-naphthalenesufonate', 
  # replacement = '6-(p-toluidinyl)-2-naphthalenesulfonic acid',
  # pattern = '4-([(4-hydroxyphenyl)azo]benzoate', 
  # replacement = '(E)-4-[(4-Hydroxyphenyl)azo]benzoic acid',
  pattern = "modant yellow 7", 
  replacement = "disodium 3-methyl-5-((4-sulphonatophenyl)azo)salicylate", 
  pattern = "-acid", 
  replacement = " acid"
) 
pattern.replacement <- fixed.sdf %>%
  gather("argument", "chemical.name") %>% 
  separate(argument, into = c("argument", "attribute"),sep = "\\.", fill ="right") %>% 
  spread(key = argument, value = "chemical.name") %>%
  dplyr::select(., -attribute) 

# Converting the table to regex for str_replace
pattern.reg <- pattern.replacement %>%
  mutate(pattern = make.regex(pattern)) %>%
  mutate(replacement = make.regex(replacement))

for(i in 1:nrow(pattern.replacement)) {
  wip.sdf$guest <- str_replace(wip.sdf$guest, pattern = pattern.reg[i, "pattern"], 
              replacement = pattern.reg[i, "replacement"])
}
wip.sdf$guest <- str_replace(wip.sdf$guest, "\u03b2", "beta")
# Problem with beta replacement
wip.sdf$guest <- str_replace(wip.sdf$guest, "4-nitrophenyl-beta-d-glucoside", 
                             "4-nitrophenyl beta-d-glucoside")
wip.sdf$guest <- str_replace(wip.sdf$guest, "4-nitrophenyl-beta-d-xyloside", 
                             "(2S,3R,4S,5R)-2-(4-nitrophenoxy)oxane-3,4,5-triol")
wip.sdf$guest <- str_replace(wip.sdf$guest, pattern = '4-nitrophenyl-beta-d-galactoside', 
                             '4-Nitrophenylgalactoside')
wip.sdf$guest <- str_replace(wip.sdf$guest, pattern = '4-nitrophenyl-beta-d-glucosamide', 
                             'N-[(2R,3R,4R,5S,6R)-4,5-dihydroxy-6-(hydroxymethyl)-2-[(4-nitrophenyl)methoxy]oxan-3-yl]acetamide')

# Cactus Download (2) -----------------------------------------------------

alpha.guest2 <- wip.sdf %>% filter(host == "alpha") %>% .$guest
beta.guest2 <- wip.sdf %>% filter(host == "beta") %>% .$guest
gamma.guest2 <- wip.sdf %>% filter(host == "gamma") %>% .$guest

# Alpha-CD
results2.alpha <-
  do.call(
    rbind,
    lapply(
      alpha.guest2,
      download.cactus.results,
      path = "./molecules/alphaCD",
      chemical.format = "SDF"
    )
  ) %>% mutate(host = "alpha")

# Beta-CD
results2.beta <-
  do.call(
    rbind,
    lapply(
      beta.guest2,
      download.cactus.results,
      path = "./molecules/betaCD",
      chemical.format = "SDF"
    )
  ) %>% mutate(host = "beta")

# Gamma-CD
results2.gamma <-
  do.call(
    rbind,
    lapply(
      gamma.guest2,
      download.cactus.results,
      path = "./molecules/gammaCD",
      chemical.format = "SDF"
    )
  ) %>% mutate(host = "gamma")

results2.all <- rbind(results2.alpha, results2.beta, results2.gamma)
results2.fail <- results2.all %>% filter(downloaded == "no")
saveRDS(results2.all, "./molecules/cactus.dwnld.2.RDS")

# 4 alpha, 13 beta, 3 gamma did not download 

#####
# Old Functions that may be Useful ----------------------------------------

# This takes around 3-4 minutes on 8GB RAM
# guest.sdf <- dataset %>% 
#   dplyr::select(guest) %>%
#   unique() %>% 
#   rowwise() %>%
#   mutate(encoded = URLencode(guest, reserved = T)) %>% 
#   mutate(url = paste0("https://cactus.nci.nih.gov/chemical/structure/", 
#                       encoded, "/sdf")) %>%
#   mutate(sdf = try(getURL(url = url))) 
# 
# guest.sdf %>% write.table("./molecules/02.guestSDF.csv", 
#                           quote = F, row.names = F)
# saveRDS(guest.sdf, "./molecules/02.guestSDF.RDS")
# guest.sdf.success <- guest.sdf %>% 
#   filter(!str_detect(string = sdf, pattern = "Page"))

# for(i in 1:nrow(guest.sdf)) {
#   filename <- paste("./molecules/", guest.sdf[i, 1], ".SDF", sep = "")
#   write.table(guest.sdf[i, 4], filename, 
#               col.names = F, row.names = F, quote = F)
# }

# Cactus resolved 395 out of 457 molecules for an 86% success rate 
# 62 failures 

# # 62 undownloaded molecules
# guest.sdf %>% filter(str_detect(string = sdf, pattern = "Page")) %>% nrow()
# 
# problem.sdf <- guest.sdf %>%
#   filter(str_detect(sdf, "Page")) %>%
#   select(guest) 
# 
# # This creates 40 entries
# problem.sdf.filter <- problem.sdf %>% 
#   filter(!str_detect(guest, pattern = "anion")) %>%
#   filter(!str_detect(guest, pattern = "carboxylate")) %>%
#   filter(!str_detect(guest, pattern = "[0-9][Hh]ydrochloride")) %>%
#   filter(!str_detect(guest, pattern = "\\â")) %>%
#   filter(!str_detect(guest, pattern = "ferrocen")) 
# Don't know why this is here, but I'm keeping it just in case
# replace.sdf <- problem.sdf %>%
#   filter(!str_detect(guest, pattern = "anion")) %>%
#   filter(!str_detect(guest, pattern = "carboxylate")) %>%
#   filter(!str_detect(guest, pattern = "[0-9][Hh]ydrochloride")) %>%
#   filter(!str_detect(guest, pattern = "\\â")) %>%
#   filter(!str_detect(guest, pattern = "ferrocen")) %>%
#   full_join(., fixed.replace.sdf, by = c("guest" = "pattern")) %>%
#   mutate(copythis = paste0("pattern = '", guest, "', ", "replacement = '", replacement, "',")) 
# 
# replace.sdf %>% dplyr::select(copythis) %>%
#   write.table("./molecules/02.1-problemSDF.txt",quote = F, row.names = F, col.names = F)
#
# problem.sdf.dwnld <- wip.sdf %>% 
#   select(guest) %>%
#   rowwise() %>%
#   mutate(encoded = URLencode(guest, reserved = T)) %>% 
#   mutate(url = paste0("https://cactus.nci.nih.gov/chemical/structure/", encoded, "/sdf")) %>%
#   mutate(sdf = try(getURL(url = url))) 
# 
# still.problem <- problem.sdf.dwnld %>%
#   filter(str_detect(string = sdf, pattern = "Page"))
# 
# sdf.success <- problem.sdf.dwnld %>%
#   filter(!str_detect(string = sdf, pattern = "Page")) %>%
#   rbind(guest.sdf.success)

# #     Saving SDFs ---------------------------------------------------------
# 
# #     Alpha - 240 molecules
# guest.a <- dataset[dataset$host == "alpha", ]
# guest.a <- guest.a$guest
# 
# sdf.a <- sdf.success %>%
#   filter(guest %in% guest.a)
# 
# for(i in 1:nrow(sdf.a)) {
#   filename <- paste("./molecules/alpha/", sdf.a[i, "guest"], ".SDF", sep = "")
#   write.table(guest.sdf[i, "sdf"], filename,
#               col.names = F, row.names = F, quote = F)
# }
# 
# #     Beta - 326 molecules
# guest.b <- dataset[dataset$host == "beta", ]
# guest.b <- guest.b$guest
# 
# sdf.b <- sdf.success %>%
#   filter(guest %in% guest.b)
# 
# for(i in 1:nrow(sdf.b)) {
#   filename <- paste("./molecules/beta/", sdf.b[i, "guest"], ".SDF", sep = "")
#   write.table(guest.sdf[i, "sdf"], filename,
#               col.names = F, row.names = F, quote = F)
# }
# 
# #     Gamma - 17 molecules
# guest.c <- dataset[dataset$host == "gamma", ]
# guest.c <- guest.c$guest
# 
# sdf.c <- sdf.success %>%
#   filter(guest %in% guest.c)
# 
# for(i in 1:nrow(sdf.c)) {
#   filename <- paste("./molecules/gamma/", sdf.c[i, "guest"], ".SDF", sep = "")
#   write.table(guest.sdf[i, "sdf"], filename,
#               col.names = F, row.names = F, quote = F)
# }
