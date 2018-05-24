# Libraries and Packages --------------------------------------------------

library(data.table)
library(e1071)
library(Matrix)

source("./06.1.ad.functions.R")

# Data --------------------------------------------------------------------

pp.settings <- readRDS("./pre-process/pp.settings.RDS")
zero.pred <- readRDS("./pre-process/zero.pred.RDS") %>% str_replace(., "-", ".")
zero.pred2 <- readRDS("./pre-process/zero.pred2.RDS")[1]
too.high <- readRDS("./pre-process/high.cor.RDS") %>% str_replace(., "-", ".")

fda.desc <- fda.padel[ , !colnames(fda.padel) %in% zero.pred]

fda.desc <- fda.desc[ , -1] 
fda.desc <- rbind(fda.desc %>% mutate(alpha = 1, beta = 0, gamma = 0), 
                  fda.desc %>% mutate(alpha = 0, beta = 1, gamma = 0), 
                  fda.desc %>% mutate(alpha = 0, beta = 0, gamma = 1))
fda.pp <- predict(pp.settings, fda.desc, na.remove = T)

fda.pp <- fda.pp[ , !colnames(fda.pp) %in% too.high]
fda.pp <- fda.pp[ , !colnames(fda.pp) %in% zero.pred2]

fda.data <- cbind(fda.padel[ , 1], fda.pp)
colnames(fda.data)[1] <- "guest"

# Applicability domain ----------------------------------------------------



# Models ------------------------------------------------------------------

#     SVM -----------------------------------------------------------------

svm.a <- readRDS("./models/svm/polysvm.alpha.RDS")
svm.b <- readRDS("./models/svm/polysvm.beta.RDS")
svm.c <- readRDS("./models/svm/polysvm.gamma.RDS")

fda.svm.a <- predict(svm.a, fda.pp %>% filter(alpha > 0)) %>%
  as.data.frame() %>% mutate(cd.type = "alpha")
fda.svm.b <- predict(svm.b, fda.pp %>% filter(beta > 0)) %>%
  as.data.frame() %>% mutate(cd.type = "beta")
fda.svm.c <- predict(svm.c, fda.pp %>% filter(gamma > 0)) %>%
  as.data.frame() %>% mutate(cd.type = "gamma")

guest <- fda.padel$Name %>% as.vector() %>% rep(., 3)
fda.svm <- data.frame(guest, rbind(fda.svm.a, fda.svm.b, fda.svm.c)) %>%
  rename(., pred = `.`)

# Removing outliers
fda.svm <- fda.svm[-which.max(fda.svm$pred), ]
fda.svm <- fda.svm[-which.max(fda.svm$pred), ]
fda.svm <- fda.svm[-which.min(fda.svm$pred), ]
fda.svm <- fda.svm[-which.min(fda.svm$pred), ]

ggplot(fda.svm, aes(x = guest, y = pred, color = cd.type)) +
  geom_point(alpha = 0.6) + theme_bw() + 
  scale_x_discrete(breaks = NULL) + 
  theme.2018 +
  # geom_hline(yintercept = -24, linetype = "dotted", color = "red", size = 1) +
  # geom_hline(yintercept = -20, linetype = "dotted", color = "orange", size = 1) +
  labs(title = "Polynomial SVM - FDA Approved Drugs", 
       y = "Predicted DelG, kJ/mol", 
       color = "CD Type", shape = "Cyclodextrin", 
       x = NULL)


#     Cubist --------------------------------------------------------------

cube.a <- readRDS("./models/cubist/cubist.alpha.RDS")
cube.b <- readRDS("./models/cubist/cubist.beta.RDS")
cube.c <- readRDS("./models/cubist/cubist.gamma.RDS")

fda.cube.a <- predict(cube.a, fda.pp %>% filter(alpha > 0)) %>%
  as.data.frame() %>% mutate(cd.type = "alpha")
fda.cube.b <- predict(cube.b, fda.pp %>% filter(beta > 0)) %>%
  as.data.frame() %>% mutate(cd.type = "beta")
fda.cube.c <- predict(cube.c, fda.pp %>% filter(gamma > 0)) %>%
  as.data.frame() %>% mutate(cd.type = "gamma")

fda.cube <- data.frame(guest, rbind(fda.cube.a, fda.cube.b, fda.cube.c)) %>%
  rename(., pred = `.`)

ggplot(fda.cube, aes(x = guest, y = pred, color = cd.type)) +
  geom_point(alpha = 0.6) + theme_bw() + 
  scale_x_discrete(breaks = NULL) + 
  theme.2018 +
  # geom_hline(yintercept = -24, linetype = "dotted", color = "red", size = 1) +
  # geom_hline(yintercept = -20, linetype = "dotted", color = "orange", size = 1) +
  labs(title = "Cubist - FDA Approved Drugs", 
       y = "Predicted DelG, kJ/mol", 
       color = "CD Type", shape = "Cyclodextrin", 
       x = NULL)


#     GLMnet --------------------------------------------------------------

glm.a <- readRDS("./models/glmnet/glm.alpha.RDS")
glm.b <- readRDS("./models/glmnet/glm.beta.RDS")
glm.c <- readRDS("./models/glmnet/glm.gamma.RDS")

fda.a <- fda.pp %>% filter(alpha > 0) %>% as.matrix()
fda.b <- fda.pp %>% filter(beta > 0) %>% as.matrix()
fda.c <- fda.pp %>% filter(gamma > 0) %>% as.matrix()
fda.glm.a <- predict.glmnet(glm.a, fda.a, s = tail(glm.a$lambda, n = 1)) %>%
  as.data.frame() %>% mutate(cd.type = "alpha")
fda.glm.b <- predict.glmnet(glm.b, fda.b, s = tail(glm.b$lambda, n = 1)) %>%
  as.data.frame() %>% mutate(cd.type = "beta")
fda.glm.c <- predict.glmnet(glm.c, fda.c, s = tail(glm.c$lambda, n = 1)) %>%
  as.data.frame() %>% mutate(cd.type = "gamma")

fda.glm <- data.frame(guest, rbind(fda.glm.a, fda.glm.b, fda.glm.c)) %>%
  rename(., pred = X1)

ggplot(fda.glm, aes(x = guest, y = pred, color = cd.type)) +
  geom_point(alpha = 0.6) + theme_bw() + 
  scale_x_discrete(breaks = NULL) + 
  theme.2018 +
  # geom_hline(yintercept = -24, linetype = "dotted", color = "red", size = 1) +
  # geom_hline(yintercept = -20, linetype = "dotted", color = "orange", size = 1) +
  labs(title = "GLMNet FDA Approved Drugs", 
       y = "Predicted DelG, kJ/mol", 
       color = "CD Type", shape = "Cyclodextrin", 
       x = NULL)

#     PLS -----------------------------------------------------------------

pls.a <- readRDS("./models/pls/pls.alpha.RDS")
pls.b <- readRDS("./models/pls/pls.beta.RDS")
pls.c <- readRDS("./models/pls/pls.gamma.RDS")

fda.pls.a <- predict(pls.a, fda.pp %>% filter(alpha > 0), ncomp = 4) %>%
  as.data.frame() %>% mutate(cd.type = "alpha") %>%
  rename(pred = `DelG.4 comps`)
fda.pls.b <- predict(pls.b, fda.pp %>% filter(beta > 0), ncomp = 4) %>%
  as.data.frame() %>% mutate(cd.type = "beta") %>%
  rename(pred = `DelG.4 comps`)
fda.pls.c <- predict(pls.c, fda.pp %>% filter(gamma > 0), ncomp = 5) %>%
  as.data.frame() %>% mutate(cd.type = "gamma") %>%
  rename(pred = `DelG.5 comps`)

fda.pls <- data.frame(guest, rbind(fda.pls.a, fda.pls.b, fda.pls.c))

ggplot(fda.pls, aes(x = guest, y = pred, color = cd.type)) +
  geom_point(alpha = 0.6) + theme_bw() + 
  scale_x_discrete(breaks = NULL) + 
  theme.2018 +
  # geom_hline(yintercept = -24, linetype = "dotted", color = "red", size = 1) +
  # geom_hline(yintercept = -20, linetype = "dotted", color = "orange", size = 1) +
  labs(title = "Cubist - FDA Approved Drugs", 
       y = "Predicted DelG, kJ/mol", 
       color = "CD Type", shape = "Cyclodextrin", 
       x = NULL)

#     Glmnet and SVM ------------------------------------------------------

glm.svm.fda <- rbind(fda.svm %>% mutate(Model = "SVM"), 
                     fda.glm %>% mutate(Model = "GLMnet"))
glm.svm.nogamma <- glm.svm.fda[glm.svm.fda$cd.type != "gamma", ]

ggplot(glm.svm.fda, aes(x = guest, y = pred, color = cd.type, shape = cd.type)) +
  geom_point(alpha = 0.6) + theme_bw() + 
  scale_x_discrete(breaks = NULL) + 
  theme(
    axis.ticks.x=element_blank(),
    axis.line = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
  ) + 
  facet_wrap(~Model) + 
  geom_hline(yintercept = -24, linetype = "dotted", color = "red", size = 1) +
  geom_hline(yintercept = -20, linetype = "dotted", color = "orange", size = 1) +
  labs(title = "GLMnet and SVM - FDA Approved Drugs", 
       y = "Predicted DelG, kJ/mol", 
       color = "Cyclodextrin", shape = "Cyclodextrin")

ggplot(glm.svm.nogamma, aes(x = guest, y = pred, color = cd.type, shape = cd.type)) +
  geom_point(alpha = 0.6) + theme_bw() + 
  scale_x_discrete(breaks = NULL) + 
  theme(
    axis.ticks.x=element_blank(),
    axis.line = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
  ) + 
  facet_wrap(~Model) + 
  geom_hline(yintercept = -24, linetype = "dotted", color = "red", size = 1) +
  geom_hline(yintercept = -20, linetype = "dotted", color = "orange", size = 1) +
  labs(title = "GLMnet and SVM - FDA Approved Drugs, No Gamma-CD", 
       y = "Predicted DelG, kJ/mol", 
       color = "Cyclodextrin", shape = "Cyclodextrin")

# Compiling ---------------------------------------------------------------

fda.comb <- rbind(fda.cube %>% mutate(model = "Cubist"), 
                  fda.glm %>% mutate(model = "GLMNet"), 
                  fda.svm %>% mutate(model = "SVM"))
ggplot(fda.comb, aes(x = guest, y = pred, color = cd.type)) +
  geom_point() + 
  theme.2018 + 
  facet_wrap(~model) + 
  labs(title = "QSPR predictions on FDA-approved drugs", 
       x = NULL, 
       y = "Predicted dG, kJ/mol", 
       color = "CD Type") + 
  ylim(-60, 25)  +
  scale_x_discrete(breaks = NULL) + 
  theme(
    axis.ticks.x=element_blank(),
    axis.line = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
  ) 
ggsave("./graphs/fda.png", dpi = 600)

# SDF splitting -----------------------------------------------------------

db.comp <- read.csv("./drugbank_approved_structures/structures.SDF", 
                    header = F, sep = "\n")

# index where SDF encoding ends, 
# aka where the new SDF starts in the next line
end.index <- which(db.comp$V1 == "$$$$") 
# index where names should go
name.index <- c(1, end.index[1:(length(end.index) - 1)] + 1) %>% as.integer()

# index where generic names can be found
generic.index <- which(db.comp$V1 == "> <GENERIC_NAME>") + 1
generic.names <- db.comp$V1[generic.index] %>% as.character()

# simple enough replacement at the correct indices
db.comp$V1[name.index] <- generic.names

# location where useful encoding ends
mend.index <- which(db.comp$V1 == "M  END")

# Just a function to split SDFs
# name is generic names, mend is the list of "M END" indices
subset.sdf <- function(name, mend, sdfs) {
  result <- sdfs[name:mend, 1] 
  header <- rep(" ", 3) # this is necessary b/c the header section syntax 
  header.end <- result %>% as.character() %>% str_detect(., "V2000") %>% which(.)
  header.temp <- result[1:(header.end - 1)]
  result <- result[header.end:length(result)]
  
  for(i in 1:length(header.temp))
    header[i] <- header.temp[i]
  
  result <- c(header, result, "$$$$")
  result <- data.frame(result)
  return(result)
} 

db.char <- db.comp %>% mutate(V1 = as.character(V1))
db.all <-
  mapply(
    FUN = subset.sdf,
    name = name.index,
    mend = mend.index,
    MoreArgs = list(sdfs = db.char),
    SIMPLIFY = F
  ) %>% rbindlist()

# write.table(temp, file = "./drugbank/fixed-sample.SDF", quote = F, row.names = F, col.names = F)
write.table(db.all, file = "./drugbank/all.SDF", quote = F, row.names = F, col.names = F)
write.table(db.all, file = "./fda/all.fda.SDF", quote = F, row.names = F, col.names = F)


# PaDEL Descriptor results ------------------------------------------------

fda.padel <- read.csv("./drugbank/all.csv")

# Pre-processing
pp.settings <- readRDS("./pre-process/pp.settings.RDS")
zero.pred <- readRDS("./pre-process/zero.pred.RDS") %>% str_replace(., "-", ".")
zero.pred2 <- readRDS("./pre-process/zero.pred2.RDS")[1]
too.high <- readRDS("./pre-process/high.cor.RDS") %>% str_replace(., "-", ".")

fda.desc <- fda.padel[ , !colnames(fda.padel) %in% zero.pred] 

# removes guest molecules that are missing values: 2327 --> 2037 guests
fda.desc <- fda.desc[complete.cases(fda.desc), ]

fda.guest <- fda.desc[ , 1]
fda.desc <- fda.desc[ , -1] 
fda.desc <- rbind(fda.desc %>% mutate(alpha = 1, beta = 0, gamma = 0), 
                  fda.desc %>% mutate(alpha = 0, beta = 1, gamma = 0), 
                  fda.desc %>% mutate(alpha = 0, beta = 0, gamma = 1))

fda.pp <- predict(pp.settings, fda.desc, na.ignore = T)

fda.pp <- fda.pp[ , !colnames(fda.pp) %in% too.high]
fda.pp <- fda.pp[ , !colnames(fda.pp) %in% zero.pred2]

saveRDS(fda.pp, "./fda/fda.pp.RDS")

fda.data <- cbind(fda.guest, fda.pp)
colnames(fda.data)[1] <- "guest"

saveRDS(fda.data, "./fda/fda.data.RDS")

# Applicability domain ----------------------------------------------------

# do to a weird conflict of pre-processing, the fda.pp dataframe actually 
# contains alpha, beta, and gamma data that doesn't need to exist
# so remove.binary is being applied to only 1/3 of fda.pp
fda.nobin <- remove.binary(fda.pp[1:(nrow(fda.pp)/3) , ])
fda.ad <- domain.num(fda.nobin)
fda.ad <- fda.ad %>% mutate(guest = as.character(fda.guest))

# There are 224 drugs that are outside of the model's applicability domain
# making up around 11% of the database
# fda.ad %>% filter(newSk > 3) %>% nrow() / 2037

# there are several outliers in fda.ad (newSk > 50)
fda.outliers <- fda.ad %>% filter(newSk >= 50) %>% .$guest

fda.ad.nooutliers <- fda.ad %>% filter(!(guest %in% fda.outliers))
ggplot(fda.ad.nooutliers, aes(x = guest, y = newSk, color = domain, shape = domain)) + 
  geom_point() + 
  theme.isef + 
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        panel.grid.major = element_blank()) + 
  geom_hline(yintercept = 3, size = 1, color = "#404040") + 
  labs(x = NULL, y = "Similarity score", color = "Domain", shape = "Domain") + 
  coord_fixed(ratio = 72)
ggsave("./fda/app.domain.png", scale = 1.35, dpi = 450)


# Sorting applicability domain drugs --------------------------------------

guest.in.ad <- fda.ad %>% filter(domain == "inside") %>% .$guest
fda <- fda.data %>% filter(guest %in% guest.in.ad)

# SVM predictions on DrugBank ---------------------------------------------

svm.a <- readRDS("./models/svm/polysvm.alpha.RDS")
svm.b <- readRDS("./models/svm/polysvm.beta.RDS")
svm.c <- readRDS("./models/svm/polysvm.gamma.RDS")

fda.svm.a <- predict(svm.a, fda[ , -1] %>% filter(alpha > 0)) %>%
  as.data.frame() %>% mutate(cd.type = "alpha")
fda.svm.b <- predict(svm.b, fda[ , -1] %>% filter(beta > 0)) %>%
  as.data.frame() %>% mutate(cd.type = "beta")
fda.svm.c <- predict(svm.c, fda[ , -1] %>% filter(gamma > 0)) %>%
  as.data.frame() %>% mutate(cd.type = "gamma")


fda.svm <- data.frame(guest <- fda[ , 1], rbind(fda.svm.a, fda.svm.b, fda.svm.c)) %>%
  rename(., pred = `.`)

ggplot(fda.svm, aes(x = guest, y = pred, color = cd.type)) + 
  geom_point() + 
  theme.isef + 
  facet_grid(~cd.type) + 
  theme(axis.text.x = )

# GLMNet ------------------------------------------------------------------

glm.a <- readRDS("./models/glmnet/glm.alpha.RDS")
glm.b <- readRDS("./models/glmnet/glm.beta.RDS")
glm.c <- readRDS("./models/glmnet/glm.gamma.RDS")

fda.a <- fda[ , -1] %>% filter(alpha > 0) %>% as.matrix()
fda.b <- fda[ , -1] %>% filter(beta > 0) %>% as.matrix()
fda.c <- fda[ , -1] %>% filter(gamma > 0) %>% as.matrix()
fda.glm.a <- predict.glmnet(glm.a, fda.a, s = tail(glm.a$lambda, n = 1)) %>%
  as.data.frame() %>% mutate(cd.type = "alpha")
fda.glm.b <- predict.glmnet(glm.b, fda.b, s = tail(glm.b$lambda, n = 1)) %>%
  as.data.frame() %>% mutate(cd.type = "beta")
fda.glm.c <- predict.glmnet(glm.c, fda.c, s = tail(glm.c$lambda, n = 1)) %>%
  as.data.frame() %>% mutate(cd.type = "gamma")

fda.glm <- data.frame(guest <- fda[ , 1], rbind(fda.glm.a, fda.glm.b, fda.glm.c)) %>%
  rename(., pred = X1)

ggplot(fda.glm, aes(x = guest, y = pred, color = cd.type)) +
  geom_point() + 
  theme.isef + 
  scale_x_discrete(breaks = NULL) + 
  labs(title = "GLMNet FDA Approved Drugs", 
       x = NULL, 
       y = "Predicted DelG, kJ/mol", 
       color = "CD Type", shape = "Cyclodextrin") + 
  facet_wrap(~cd.type)

ggsave("./fda/glmnet.png", dpi = 450, scale = 1.5)

# Compiled FDA predictions ------------------------------------------------

fda.ensemble <- inner_join(fda.svm, fda.glm, by = c("guest....fda...1.", "cd.type")) %>%
  mutate(pred  = (pred.x + pred.y)/2) %>%
  select(-pred.x) %>%
  select(-pred.y)
fda.all <- rbind(fda.svm %>% mutate(QSPR = "SVM"), 
                 fda.glm %>% mutate(QSPR = "GLMNet"), 
                 fda.ensemble %>% mutate(QSPR = "Ensemble"))
colnames(fda.all)[1] <- "guest"
ggplot(fda.all, aes(x = guest, y = pred, color = cd.type)) + 
  geom_point() + 
  theme.paper.2018 + 
  scale_x_discrete(breaks = NULL) + 
  labs(x = NULL, 
       y = "Predicted delF, kJ/mol", 
       color = "CD") + 
  facet_grid(QSPR ~ cd.type)
ggsave("./fda/qsars and ensemble.png")
ggsave("./graphs/2018 paper/qsar and ensemble.png")

# Promising alpha-CD drugs
# docosanol, oleic acid, cetyl alcohol
fda.a.drugs <- c("Docosanol", "Oleic Acid")
fda.rows <- fda.all %>% mutate(row = rep(1:1813, 6))
fda.a.drugs <- fda.rows %>% filter(guest %in% fda.a.drugs) %>%
  filter(cd.type == "alpha")

# Promising beta-CD drugs
# Vitamin E, cholesterol, testosterone undecanoate
fda.b.drugs <- c("Vitamin E", "Cholesterol")
fda.b.drugs <- fda.rows %>% filter(guest %in% fda.b.drugs) %>%
  filter(cd.type == "beta")

# Promising gamma-CD drugs
# balsalazide, montelukast
fda.c.drugs <- c("Montelukast")
fda.c.drugs <- fda.rows %>% filter(guest %in% fda.c.drugs) %>%
  filter(cd.type == "gamma")

fda.drugs <- rbind(fda.a.drugs, fda.b.drugs, fda.c.drugs)
# don't show up well on graph
# fda.drugs <- fda.drugs %>% filter(!(guest %in% c("Testosterone undecanoate", "Cetyl alcohol")))
ggplot(fda.rows, aes(x = row, y = pred)) + 
  geom_point(aes(color = cd.type)) +
  scale_x_discrete(breaks = NULL) + 
  labs(x = NULL, 
       y = "Predicted dG, kJ/mol", 
       color = "CD") + 
  facet_grid(QSPR ~ cd.type) + 
  geom_point(data = fda.drugs, aes(x = row, y = pred)) + 
  theme.isef + 
  geom_label(data = fda.drugs, aes(x = row, y = pred, label = guest, family = "Clear Sans Light"), 
             hjust = -0.05, vjust = 1) + 
  theme(text = element_text(size = 24), 
        legend.key.size = unit(1.25, 'lines'))
ggsave("./fda/predicted drugs.png", scale = 1.25, dpi = 450)  

# Applicability domain ----------------------------------------------------


source("./06.1.ad.functions.R")
fda.pp <- readRDS("./fda/fda.pp.RDS")
fda.guest <- read.csv("./drugbank/all.csv") %>% select(., Name)
fda.nobin <- fda.pp %>% remove.binary()
# fda.small <- fda.nobin[1:50, ]
# domain.num(fda.small)

fda.ad <- domain.num(fda.nobin)
fda.ad <- data.frame(fda.guest, fda.ad)
ggplot(fda.ad, aes(x = fda.guest, y = newSk)) + 
  geom_point()
