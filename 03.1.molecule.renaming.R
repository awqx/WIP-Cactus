make.regex <- function(string) {
  new.string <- str_replace_all(string, pattern = "\\(", 
                                replacement = "\\\\(") %>% 
    str_replace_all(pattern = "\\)", replacement = "\\\\)") %>%
    str_replace_all(pattern = "\\-", replacement = "\\\\-")
  # new.string <- str_replace(new.string, pattern = "\\ÃÂ²|\\B\\s+H\\+|\u03b2", # alternatives: \\ÃÂ²|\\B\\s+H\\+
  #                           replacement = "\\(ÃÂ²|\u03b2\\)")
  return(new.string)
}

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

# for(i in 1:nrow(pattern.replacement)) {
#   wip.sdf$guest <- str_replace(wip.sdf$guest, pattern = pattern.reg[i, "pattern"], 
#                                replacement = pattern.reg[i, "replacement"])
# }
# wip.sdf$guest <- str_replace(wip.sdf$guest, "\u03b2", "beta")
# # Problem with beta replacement
# wip.sdf$guest <- str_replace(wip.sdf$guest, "4-nitrophenyl-beta-d-glucoside", 
#                              "4-nitrophenyl beta-d-glucoside")
# wip.sdf$guest <- str_replace(wip.sdf$guest, "4-nitrophenyl-beta-d-xyloside", 
#                              "(2S,3R,4S,5R)-2-(4-nitrophenoxy)oxane-3,4,5-triol")
# wip.sdf$guest <- str_replace(wip.sdf$guest, pattern = '4-nitrophenyl-beta-d-galactoside', 
#                              '4-Nitrophenylgalactoside')
# wip.sdf$guest <- str_replace(wip.sdf$guest, pattern = '4-nitrophenyl-beta-d-glucosamide', 
#                              'N-[(2R,3R,4R,5S,6R)-4,5-dihydroxy-6-(hydroxymethyl)-2-[(4-nitrophenyl)methoxy]oxan-3-yl]acetamide')
