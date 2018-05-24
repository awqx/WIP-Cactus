library(caret)
library(tidyverse)

# source("./03.cactus.functions.R")
dir.create("./predictions")
dir.create("./predictions/euroasia")
dir.create("./predictions/euroasia/molecules")

# Running molecules from ERD's data
# euroasia.mol <- molecules$molecule
# euroasia.mol <- euroasia.mol %>% str_replace(., "[[:space:]]$", "") %>%
#   str_replace(., "^[[:space:]]", "")
# saveRDS(euroasia.mol, "./predictions/euroasia/mol.names.RDS")

euroasia.mol <- readRDS("./predictions/euroasia/mol.names.RDS")

cactus.ea <-
  do.call(
    rbind,
    lapply(
      euroasia.mol,
      download.cactus.results,
      path = "./predictions/euroasia/molecules",
      chemical.format = "SDF"
    )
  )
# View(cactus.ea)
ea.fail <- cactus.ea %>% filter(downloaded == "no") %>% .$guest %>% as.character()
saveRDS(cactus.ea, "./molecules/euroasia/results.RDS")

# After running PaDEL-Descriptor
ea.padel <- read.csv("./molecules/euroasia.csv")
ea.desc <- ea.padel[complete.cases(ea.padel), ]
ea.desc <- ea.desc[ , !colnames(ea.desc) %in% zero.pred] 
ea.guest <- ea.desc[ , 1]

ea.desc <- ea.desc[ , -1] 
ea.desc <- rbind(ea.desc %>% mutate(alpha = 1, beta = 0, gamma = 0), 
                  ea.desc %>% mutate(alpha = 0, beta = 1, gamma = 0), 
                  ea.desc %>% mutate(alpha = 0, beta = 0, gamma = 1))
ea.pp <- predict(pp.settings, ea.desc, na.remove = T)

ea.pp <- ea.pp[ , !colnames(ea.pp) %in% too.high]
ea.pp <- ea.pp[ , !colnames(ea.pp) %in% zero.pred2]

ea.data <- cbind(ea.guest, ea.pp)
colnames(ea.data)[1] <- "guest"

# Applicability domain ----------------------------------------------------

ea.nobin <- remove.binary(ea.data[1:(nrow(ea.data)/3), ])
ea.nobin <- sapply(ea.nobin, as.numeric)
ea.ad <- domain.num(ea.nobin) %>% mutate(guest = as.character(ea.guest))

ea.inside <- ea.ad %>% filter(domain == "inside") %>% .$guest

ea <- ea.data %>% filter(guest %in% ea.inside)
saveRDS(ea, "./data/euroasia.RDS")

# SVM ---------------------------------------------------------------------

ea.svm.a <- predict(svm.a, ea[ , -1] %>% filter(alpha > 0)) %>%
  as.data.frame() %>% mutate(cd.type = "alpha")
ea.svm.b <- predict(svm.b, ea[ , -1] %>% filter(beta > 0)) %>%
  as.data.frame() %>% mutate(cd.type = "beta")
ea.svm.c <- predict(svm.c, ea[ , -1] %>% filter(gamma > 0)) %>%
  as.data.frame() %>% mutate(cd.type = "gamma")

ea.svm <-
  data.frame(guest <-
               ea[, 1], rbind(ea.svm.a, ea.svm.b, ea.svm.c)) %>%
  rename(., pred = `.`)

ggplot(ea.svm, aes(x = guest, y = pred, color = cd.type)) + 
  geom_point() + 
  theme_bw() + 
  facet_grid(~cd.type) + 
  theme(axis.text.x = NULL)

# GLMNet ------------------------------------------------------------------


