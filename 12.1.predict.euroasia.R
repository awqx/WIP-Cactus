library(caret)
library(tidyverse)

source("./master-process.R")
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
saveRDS(cactus.ea, "./predictions/euroasia/cactus-results.RDS")

# After running PaDEL-Descriptor
# Settings used:
  # Remove salts
  # Detect aromaticity
  # Use filename as molecule name
  # Maximum time of 200000 (200 seconds)
  # 1D/2D descriptors, fingerprints

ea.padel <- read.csv("./predictions/euroasia/descriptors.csv")
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

saveRDS(ea.data, "./predictions/euroasia/ea.all.pp.RDS")

# Applicability domain ----------------------------------------------------

ea.nobin <- remove.binary(ea.data[1:(nrow(ea.data)/3), ])
ea.nobin <- sapply(ea.nobin, as.numeric)
ea.ad <- domain.num(ea.nobin) %>% mutate(guest = as.character(ea.guest))

ea.inside <- ea.ad %>% filter(domain == "inside") %>% .$guest

ea <- ea.data %>% filter(guest %in% ea.inside)
saveRDS(ea, "./predictions/euroasia/ea.RDS")

# Predictions -------------------------------------------------------------

# Sorting data
ea.a <- ea[ , -1] %>% filter(alpha > 0)
ea.b <- ea[ , -1] %>% filter(beta > 0)
ea.c <- ea[ , -1] %>% filter(gamma > 0)

#     SVM -----------------------------------------------------------------

ea.svm.a <- predict(svm.a, ea.a) %>%
  as.data.frame() %>% mutate(cd.type = "alpha")
ea.svm.b <- predict(svm.b, ea.b) %>%
  as.data.frame() %>% mutate(cd.type = "beta")
ea.svm.c <- predict(svm.c, ea.c) %>%
  as.data.frame() %>% mutate(cd.type = "gamma")

ea.svm <-
  data.frame(guest <-
               ea[, 1], rbind(ea.svm.a, ea.svm.b, ea.svm.c)) %>%
  rename(., pred = `.`)
colnames(ea.svm)[1] <- "guest"

ggplot(ea.svm, aes(x = guest, y = pred, color = cd.type)) + 
  geom_point() + 
  theme_bw() + 
  facet_grid(~cd.type) + 
  theme(axis.text.x = NULL)

#     GLMNet --------------------------------------------------------------

ea.glm.a <- predict.glmnet(glm.a, as.matrix(ea.a),
                           s = tail(glm.a$lambda, n = 1)) %>%
  as.data.frame() %>% mutate(cd.type = "alpha")
ea.glm.b <- predict.glmnet(glm.b, as.matrix(ea.b),
                           s = tail(glm.b$lambda, n = 1)) %>%
  as.data.frame() %>% mutate(cd.type = "beta")
ea.glm.c <- predict.glmnet(glm.c, as.matrix(ea.c),
                           s = tail(glm.c$lambda, n = 1)) %>%
  as.data.frame() %>% mutate(cd.type = "gamma")

ea.glm <- data.frame(ea[ , 1], rbind(ea.glm.a, ea.glm.b, ea.glm.c)) 
colnames(ea.glm)[1] <- "guest"
colnames(ea.glm)[2] <- "pred"  

ggplot(ea.glm, aes(x = guest, y = pred, color = cd.type)) + 
  geom_point() + 
  theme_bw() + 
  facet_grid(~cd.type) + 
  theme(axis.text.x = NULL)

# Compiled ----------------------------------------------------------------

ea.ensemble <- inner_join(ea.glm %>% rename(glm = pred),
                          ea.svm %>% rename(svm = pred),
                          by = c("guest", "cd.type")) %>%
  mutate(pred = (glm + svm)/2)

ea.ensemble.long <- rbind(ea.glm %>% mutate(QSPR = "GLMNet"), 
                          ea.svm %>% mutate(QSPR = "SVM"), 
                          ea.ensemble %>% select(., guest, pred, cd.type) %>% mutate(QSPR = "Ensemble"))

ggplot(ea.ensemble, aes(x = guest, y = pred, color = cd.type)) + 
  geom_point() + 
  theme_bw() + 
  facet_grid(~cd.type) + 
  theme(axis.text.x = NULL)
ggsave("./predictions/euroasia/ensemble-graph.png")

ggplot(ea.ensemble.long, aes(x = guest, y = pred, shape = QSPR, color = cd.type)) + 
  geom_point() + 
  theme_bw() + 
  facet_grid(~cd.type)
ggsave("./predictions/euroasia/all-predictions.png")

saveRDS(ea.ensemble, "./predictions/euroasia/predictions.RDS")
write.csv(ea.ensemble[c(1, 3, 2, 4, 5)], "./predictions/euroasia/predictions.csv")
