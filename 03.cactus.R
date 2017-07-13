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


# Cactus Download (1) -----------------------------------------------------

dir.create("./molecules")
dataset <-readRDS("./bound/combined ri and suzuki.RDS")

# This takes around 3-4 minutes on 8GB RAM
guest.sdf <- dataset %>% 
  dplyr::select(guest) %>%
  unique() %>% 
  rowwise() %>%
  mutate(encoded = URLencode(guest, reserved = T)) %>% 
  mutate(url = paste0("https://cactus.nci.nih.gov/chemical/structure/", 
                      encoded, "/sdf")) %>%
  mutate(sdf = try(getURL(url = url))) 

guest.sdf %>% write.table("./molecules/02.guestSDF.csv", 
                          quote = F, row.names = F)
guest.sdf.success <- guest.sdf %>% 
  filter(!str_detect(string = sdf, pattern = "Page"))

for(i in 1:nrow(guest.sdf)) {
  filename <- paste("./molecules/", guest.sdf[i, 1], ".SDF", sep = "")
  write.table(guest.sdf[i, 4], filename, 
              col.names = F, row.names = F, quote = F)
}

# Cactus resolved 420 out of 511 molecules for an 82% success rate 
# 91 failures 

# Cactus Cleaning ---------------------------------------------------------

# 91 undownloaded molecules
guest.sdf %>% filter(str_detect(string = sdf, pattern = "Page")) %>% nrow()

problem.sdf <- guest.sdf %>%
  filter(str_detect(sdf, "Page")) %>%
  select(guest) 

# This creates 73 entries
problem.sdf %>% 
  filter(!str_detect(guest, pattern = "anion")) %>%
  filter(!str_detect(guest, pattern = "carboxylate")) %>%
  filter(!str_detect(guest, pattern = "[0-9][Hh]ydrochloride")) %>%
  filter(!str_detect(guest, pattern = "\\â")) %>%
  filter(!str_detect(guest, pattern = "ferrocen")) %>% write.csv("./clean.data/02.1.problem.sdf.csv")

write.csv(problem.sdf, "./clean.data/02.problem.sdf.all.csv")

# Cleaning issues with download :

# 1. Presence of anionic species - 26
# 29 cases - AX
problem.sdf <- problem.sdf %>%
  rowwise() %>%
  mutate(guest.charge = str_extract(guest, pattern = "anion|monoanion|dianion"),
         guest = gsub(x = guest,
                      pattern = "\\(anion\\)|\\(monoanion\\)|\\(dianion\\)",
                      replacement = "")) # %>%

# 3. Carboxylic acid should be carboxylate - 2
# Doesn't make a difference if they're converted. Dang. 

# 4. 2HCl should be dihydrochloride - 4 # None for me - AX
problem.sdf$guest <- problem.sdf$guest %>% 
  str_replace(., pattern = "2HCl", replacement = "dihydrochloride")

# 5. Unconventional apostrophe \\â  - 7 # None for me - AX
problem.sdf %>% filter(str_detect(guest, "\\â")) 

# 6. Non covalent structures(ferrocenes), using inchikey instead- 10
problem.sdf$guest <- problem.sdf$guest  %>%  
  str_replace(., pattern = "(S)-1-ferrocenylethanol|(R)-1-ferrocenylethanol",
              replacement = "1-(Ferrocenyl)ethanol")

# 7. Name not compatible by any search engine (Cactus, Google, Chemspider). Assumed 
#    to be close to the one replaced - 1


# 8. viologen compounds - 3

problem.sdf %>% filter(str_detect(guest, "viologen"))

# 8. Misspellings by hand? -  8
#-------------
# Option 1 detect and replace
# Reassign name
fixed.sdf <- data.frame(
  pattern = "biebricht scarlet",
  replacement = "biebrich scarlet",
  pattern = "(S)-1,1â-binaphthyl-2,2â-dicarboxylic acid",
  replacement = "1,1â-binaphthyl-2,2â-dicarboxylic acid",
  pattern = "4-([(4-hydroxyphenyl)azo]benzoate",
  replacement = "4-[(4-hydroxyphenyl)azo]benzoate",
  pattern = "3-(aminomethyl)proxyl",
  replacement = "3-(aminomethyl)-proxyl",
  pattern = "3-carbamoylproxyl",
  replacement = "3-carbamoyl-proxyl",
  pattern = "l-Î±-O-benzylglycerol",
  replacement = "1-o-benzyl-rac-glycerol",
  pattern = "4-nitrophenyl-Î²-d-glucoside",
  replacement = "4-nitrophenyl-beta-D-glucopyranoside",
  pattern = "3-(aminomethyl)proxyl",
  replacement = "3-aminomethyl-proxyl",
  pattern = "bromodiphenhydramine hydrochloride",
  replacement = "Bromdiphenhydramine hydrochloride",
  pattern = "diammine(1,1-cyclobutane- dicarboxylato)platinum(II)" ,
  replacement = "cis-diammine(1,1-cyclobutanedicarboxylato)platinum(II)",
  pattern = "(1R,2S)-(-)-ephedrine",
  replacement = "L-Ephedrine",
  pattern = "(1S,2R)-(+)-ephedrine",
  replacement = "D-Ephedrine",
  pattern = "sulfoisomidine",
  replacement = "sulfisomidine",
  pattern = "diphenyhydramine hydrochloride",
  replacement = "diphenhydramine hydrochloride",
  pattern = "trans,trans-2,4- hexadienedioic (muconic) acid",
  replacement = "muconic acid",
  pattern = "hexyl-Î²-d-glucopyranoside",
  replacement = "hexyl beta-d-glucopyranoside",
  pattern = "pentakis(ethyleneglycol) monohexyl ether)",
  replacement = "2-(Hexyloxy)ethanol",
  pattern = "2-thiophenobarbital",
  replacement = "5-Ethyl-5-phenyl-2-thioxodihydro-4,6(1H,5H)-pyrimidinedione",
  pattern = "4-(hydroxyphenethyl)ammonium" ,
  replacement = "2-(4-hydroxyphenyl)ethylazanium",
  pattern = "2hydrochloride",
  replacement = "dihydrochloride",
  pattern = '(1R,2S)-(-)-ephedrine', replacement = 'L-Ephedrine',
  pattern = '(1S,2R)-(+)-ephedrine', replacement = 'D-Ephedrine',
  pattern = 'd-mandelate', replacement = '(R)-Mandelate',
  pattern = 'l-mandelate', replacement = '(S)-Mandelate',
  pattern = 'methapyriline hydrochloride', replacement = 'Methoxylene',
  pattern = 'methyl red (cation, protonated)', replacement = 'methyl red',
  pattern = '4-nitrophenyl-Î²-d-galactoside', 
  replacement = ' 4-Nitrophenylgalactoside',
  pattern = '4-nitrophenyl-Î²-d-glucosamide', 
  replacement = 'N-[(2R,3R,4R,5S,6R)-4,5-dihydroxy-6-(hydroxymethyl)-2-[(4-nitrophenyl)methoxy]oxan-3-yl]acetamide',
  pattern = '4-nitrophenyl-Î²-d-xyloside', 
  replacement = '(2S,3R,4S,5R)-2-(4-nitrophenoxy)oxane-3,4,5-triol',
  pattern = '(Â±)-norphenylephrine', 
  replacement = '3-(2-amino-1-hydroxyethyl)phenol',
  pattern = 'pentakis(ethyleneglycol) monohexyl ether', 
  replacement = '2-hexoxyethanol',
  pattern = '(R)-(-)-phenylephrine', 
  replacement = '[(2R)-2-hydroxy-2-(3-hydroxyphenyl)ethyl]-methylazanium',
  pattern = "phenyl-Î²-d-glucopyranoside",
  replacement = "phenyl-beta-d-glucopyranoside", 
  pattern = "(Â±)-octopamine",
  replacement = "4-(2-amino-1-hydroxyethyl)phenol", 
  pattern = 'd-phenyltrifluoroethanol', 
  replacement = ' 2,2,2-trifluoro-1-phenylethanol',
  pattern = 'pyrilammonium maleate',
  replacement = 'pyrilamine maleate',
  pattern = 'terfenadine hydrochloride', 
  replacement = '1-(4-tert-butylphenyl)-4-[4-[hydroxy(diphenyl)methyl]piperidin-1-yl]butan-1-ol hydrochloride',
  pattern = 'triiodide (I3-)', replacement = 'triiodide',
  pattern = 'l-tryprophan', replacement = 'l-tryptophan',
  pattern = 'Tyr-Gly-Gly-Phe-Leu', replacement = 'Leucine enkephalin',
  pattern = '1-adamantaneammonium', replacement = '1-Adamantanaminium',
  pattern = '(Â±)-anisodamine', replacement = 'anisodamine',
  pattern = '(-)-anisodamine', replacement = 'anisodamine',
  pattern = '(-)-anisodamine HBr', 
  replacement = '[(3S,6S)-6-hydroxy-8-methyl-8-azabicyclo[3.2.1]octan-3-yl] (2R)-3-hydroxy-2-phenylpropanoate',
  pattern = '(-)-anisodine HBr', replacement = 'Anisodine hydrobromide',
  pattern = '(Â±)-atropine H2SO4', replacement = 'HOBWAPHTEJGALG-UHFFFAOYSA-N',
  pattern = '(4Z,15Z)-bilirubin-IXa', 
  replacement = 'BPYKTIZUTYGOLE-KDUUSRDASA-N',
  pattern = '(4Z,15Z)-bilirubin-IXa + cyclooctanol', 
  replacement = 'BPYKTIZUTYGOLE-KDUUSRDASA-N',
  pattern = 'dansyl-l-hydroxyplorine', 
  replacement = 'dansyl-l-hydroxyproline',
  pattern = "methapyriline hydrochloride", 
  replacement = "methapyrilene hydrochloride", 
  pattern = 'dicumarol', 
  replacement = '4-hydroxy-3-[(4-hydroxy-2-oxochromen-3-yl)methyl]chromen-2-one',
  pattern = 'ethylbis(coumacetate)', 
  replacement = 'Ethyl biscoumacetate',
  pattern = 'ethylthiobarbituic acid', 
  replacement = '5-Ethyl-2-thioxodihydropyrimidine-4,6(1H,5H)-dione',
  pattern = '4-[(4-hydroxy-1-naphthyl)azo]- naphthalene-1-sulfonate', 
  replacement = '4-[(4-hydroxy-1-naphthyl)azo]- naphthalene-1-sulfonic acid',
  pattern = 'maprotilin', replacement = 'maprotiline',
  pattern = 'naproxenate', replacement = 'naproxen',
  pattern = 'nortriptylin', replacement = ' Aventyl',
  pattern = 'protriptylin', replacement = 'protriptyline',
  pattern = '(-)-scopolamine HBr', 
  replacement = '(-)-scopolamine hydrobromide',
  pattern = 'sulfasnilamide', 
  replacement = '4-aminobenzenesulfonamide',
  pattern = 'sulfathidole', 
  replacement = 'sulfaethidole',
  pattern = 'thiophenobarbital', 
  replacement = '5-(1,1,2,2,2-pentadeuterioethyl)-5-phenyl-2-sulfanylidene-1,3-diazinane-4,6-dione',
  pattern = '6-(p-toluidinyl)-2-naphthalenesufonate', 
  replacement = '6-(p-toluidinyl)-2-naphthalenesulfonic acid',
  pattern = '4-([(4-hydroxyphenyl)azo]benzoate', 
  replacement = '(E)-4-[(4-Hydroxyphenyl)azo]benzoic acid',
  pattern = "modant yellow 7", 
  replacement = "disodium 3-methyl-5-((4-sulphonatophenyl)azo)salicylate"
) 
pattern.replacement <- fixed.sdf %>%
  gather("argument", "chemical.name") %>% 
  separate(argument, into = c("argument", "attribute"),sep = "\\.", fill ="right") %>% 
  spread(key = argument, value = "chemical.name") %>%
  dplyr::select(., -attribute) 

fixed.replace.sdf <- fixed.sdf %>%
  gather("argument", "chemical.name") %>% 
  separate(argument, into = c("argument", "attribute"),sep = "\\.", fill ="right") %>% 
  spread(key = argument, value = "chemical.name") 
# Str-replace ------
#  spread(key = argument, value = chemical.name)
# replace.1 <- guest.replace %>% filter(argument == "pattern")
# replace.2 <- guest.replace %>% filter(argument == "replacement")
# figure this out later
# replace.2[17, 2] <- ""
# chem.pattern <- replace.1$chemical.name %>% as.vector()
# chem.replacement <- replace.2$chemical.name %>% as.vector()
# 
# problem.sdf$guest <- problem.sdf$guest %>% str_replace(., chem.pattern, chem.replacement)
problem.sdf$guest <- problem.sdf$guest %>%
  str_replace(pattern = "biebricht scarlet", replacement = "biebrich scarlet") %>%
  str_replace(pattern = "\\(S\\)-1,1â-binaphthyl-2,2â-dicarboxylic acid", 
              replacement = "1,1â-binaphthyl-2,2â-dicarboxylic acid") %>%
  str_replace(pattern = "4-\\(\\[\\(4-hydroxyphenyl\\)azo\\]benzoate", 
              replacement = "4-[(4-hydroxyphenyl)azo]benzoate") %>%
  str_replace(pattern = "3-\\(aminomethyl\\)proxyl", 
              replacement = "3-(aminomethyl)-proxyl") %>%
  str_replace(pattern = "3-carbamoylproxyl", 
              replacement = "3-carbamoyl-proxyl") %>%
  str_replace(pattern = "l-.-O-benzylglycerol", 
              replacement = "1-o-benzyl-rac-glycerol") %>%
  str_replace(pattern = "4-nitrophenyl-.-d-glucoside", 
              replacement = "4-Nitrophenyl-beta-D-glucopyranoside") %>%
  str_replace(pattern = "3-\\(aminomethyl\\)proxyl", 
              replacement = "3-aminomethyl-proxyl") %>%
  str_replace(pattern = "bromodiphenhydramine hydrochloride", 
              replacement = "Bromdiphenhydramine hydrochloride") %>%
  str_replace(pattern = "diammine\\(1,1-cyclobutane- dicarboxylato\\)platinum\\(II\\)" ,
              replacement = "cis-diammine(1,1-cyclobutanedicarboxylato)platinum(II)") %>%
  str_replace(pattern = "\\(1R,2S\\)-\\(\\.\\)-ephedrine", replacement = "L-Ephedrine") %>%
  str_replace(pattern = "\\(1S,2R\\)-\\(\\+\\)-ephedrine", replacement = "D-Ephedrine") %>%
  str_replace(pattern = "sulfoisomidine", replacement = "sulfisomidine") %>%
  str_replace(pattern = "diphenyhydramine hydrochloride", 
              replacement = "diphenhydramine hydrochloride") %>%
  str_replace(pattern = "trans,trans-2,4- hexadienedioic \\(muconic\\) acid", 
              replacement = "muconic acid") %>%
  str_replace(pattern = "hexyl-.-d-glucopyranoside", 
              replacement = "hexyl beta-d-glucopyranoside") %>%
  str_replace(pattern = "4-\\(hydroxyphenethyl\\)ammonium", replacement = "") %>%
  str_replace(pattern = "pentakis\\(ethyleneglycol\\) monohexyl ether", 
              replacement = "2-(Hexyloxy)ethanol") %>%
  str_replace(pattern = "2-thiophenobarbital",
              replacement = "5-Ethyl-5-phenyl-2-thioxodihydro-4,6(1H,5H)-pyrimidinedione") %>%
  str_replace(pattern = "4-\\(hydroxyphenethyl\\)ammonium", 
              replacement = "2-(4-hydroxyphenyl)ethylazanium") %>%
  str_replace(pattern = "2hydrochloride",
              replacement = "dihydrochloride") %>%
  str_replace(pattern = "\\(.\\)-octopamine",
              replacement = "(S)-octopamine") %>%
  str_replace(pattern = "methapyriline hydrochloride", 
              replacement = "methapyrilene hydrochloride") %>%
  str_replace(pattern = "modant yellow 7", 
              replacement = "disodium 3-methyl-5-((4-sulphonatophenyl)azo)salicylate") %>%
  str_replace(pattern = "\\(R\\)-\\(\\-\\)-phenylephrine", 
              replacement = "l-phenylephrine")
# -----------
replace.sdf <- problem.sdf %>%
  filter(!str_detect(guest, pattern = "anion")) %>%
  filter(!str_detect(guest, pattern = "carboxylate")) %>%
  filter(!str_detect(guest, pattern = "[0-9][Hh]ydrochloride")) %>%
  filter(!str_detect(guest, pattern = "\\â")) %>%
  filter(!str_detect(guest, pattern = "ferrocen")) %>%
  full_join(., fixed.replace.sdf, by = c("guest" = "pattern")) %>%
  mutate(copythis = paste0("pattern = '", guest, "', ", "replacement = '", replacement, "',")) 

replace.sdf %>% dplyr::select(copythis) %>%
  write.table("./bound/02.1-problemSDF.txt",quote = F, row.names = F, col.names = F)

read.table("./bound/02.1-problemSDF.txt", sep = "\n") %>% View()

# Cactus Download (2) -----------------------------------------------------

problem.sdf.dwnld <- problem.sdf %>% 
  select(guest) %>%
  rowwise() %>%
  mutate(encoded = URLencode(guest, reserved = T)) %>% 
  mutate(url = paste0("https://cactus.nci.nih.gov/chemical/structure/", encoded, "/sdf")) %>%
  mutate(sdf = try(getURL(url = url))) 

problem.sdf %>% 
  filter(!str_detect(guest, pattern = "anion")) %>%
  filter(!str_detect(guest, pattern = "carboxylate")) %>%
  filter(!str_detect(guest, pattern = "[0-9][Hh]ydrochloride")) %>%
  filter(!str_detect(guest, pattern = "\\â")) %>%
  filter(!str_detect(guest, pattern = "ferrocen")) %>%  
  mutate(copythis = paste0("pattern = ", guest, ", ", "replacement = `` ")) %>% View()

