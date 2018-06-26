source("03.cactus.functions.R")

# Cleaning data -----------------------------------------------------------

dir.create("./molecules")
dir.create("./molecules/alphaCD")
dir.create("./molecules/betaCD")
dir.create("./molecules/gammaCD")
dataset <-readRDS("./cleaning/02.combined.data.RDS")
connors <- readRDS("./cleaning/02.connors.RDS") %>%
  mutate(host = "gamma", data.source = "connors") %>%
  select(-ka)
dataset <- rbind(dataset, connors)


# 1. Presence of anionic species 
dataset <- dataset %>%
  rowwise() %>%
  mutate(guest.charge = str_extract(guest, pattern = "anion|monoanion|dianion"),
         guest = gsub(x = guest,
                      pattern = "\\(anion\\)|\\(monoanion\\)|\\(dianion\\)",
                      replacement = "")) # %>%

# 2. Unconventional apostrophe \\â or â\u0080\u0098  
dataset$guest <- dataset$guest %>%
  str_replace("(\\â | â\u0080\u0098)", "\\'") 

# 3. Requires manual replacement - 42
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
  dataset$guest <- str_replace(dataset$guest, pattern = pattern.reg[i, 1], 
                               replacement = pattern.reg[i, 2])
}
# "-acid"has problems reading from the .csv
dataset$guest <- str_replace(dataset$guest, pattern = "-acid", 
                             replacement = " acid")

# Problem with beta replacement
dataset$guest <- str_replace(dataset$guest, "\u03b2", "beta")
dataset$guest <- str_replace(dataset$guest, "4-nitrophenyl-beta-d-glucoside", 
                             "4-nitrophenyl beta-d-glucoside")
dataset$guest <- str_replace(dataset$guest, "4-nitrophenyl-beta-d-xyloside", 
                             "(2S,3R,4S,5R)-2-(4-nitrophenoxy)oxane-3,4,5-triol")
dataset$guest <- str_replace(dataset$guest, pattern = '4-nitrophenyl-beta-d-galactoside', 
                             '4-Nitrophenylgalactoside')
dataset$guest <- str_replace(dataset$guest, pattern = '4-nitrophenyl-beta-d-glucosamide', 
                             'N-[(2R,3R,4R,5S,6R)-4,5-dihydroxy-6-(hydroxymethyl)-2-[(4-nitrophenyl)methoxy]oxan-3-yl]acetamide')

# Reading dataset for guest molecules specific to host
alpha.guest <- dataset %>% filter(host == "alpha") %>% .$guest
beta.guest <- dataset %>% filter(host == "beta") %>% .$guest
gamma.guest <- dataset %>% filter(host == "gamma") %>% 
  .$guest %>% unique()

# Cactus Download ---------------------------------------------------------

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
saveRDS(results.all, "./cleaning/cactus.dwnld.results.RDS")

# Failed SDFs -------------------------------------------------------------

# 3 undownloaded alpha guests
alpha.fail <- results.alpha %>% 
  filter(downloaded == "no") %>%
  select(guest)

# 23 undownloaded beta guests
beta.fail <- results.beta %>% 
  filter(downloaded == "no") %>%
  select(guest)

# 6 undownloaded gamma guests
gamma.fail <- results.gamma %>% 
  filter(downloaded == "no") %>%
  select(guest)

# Total 21 undownloaded
problem.sdf <- results.all %>% 
  filter(downloaded == "no") %>%
  select(guest, host)

# Saving the dataset with the renamed molecules
saveRDS(dataset, "cleaning/03.rename.RDS")

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
