# Creating visually consistent graphs for posters, slides, etc.
# setwd("~/SREP LAB/qsar")
source("./graph.functions.R")

# Importing Data ----------------------------------------------------------

# Data from docking experiments
#   I'll fix the filepaths later
docking <- readRDS("~/SREP LAB/docking/docking data.RDS") %>% rename(cd.type = host)

# Compiled a, b, c-CD data from QSARs
cubist <- readRDS("./models/cubist/cubist.results.RDS") %>% 
  mutate(qsar.type = "Cubist") 
glmnet <- readRDS("./models/glmnet/glmnet.tst.results.RDS") %>%
  mutate(qsar.type = "GLMNet") %>%
  select(-resid)
pls <- readRDS("./models/pls/pls.results.RDS") %>%
  mutate(qsar.type = "PLS")
rforest <- readRDS("./models/rforest/rf.results.RDS") %>%
  rename(residual = tst.resid) %>%
  mutate(qsar.type = "Random Forest")
svm <- readRDS("./models/svm/polysvm.tst.results.RDS") %>%
  rename(residual = resid) %>%
  mutate(qsar.type = "SVM")
qsar.results <- rbind(cubist, glmnet, pls, rforest, svm)
# comparing qsar and docking
dock.temp <- docking %>% select(obs, pred, cd.type) %>% 
  mutate(type = "PyRx Docking") %>% data.table()
qsar.temp <- qsar.results %>% select(obs, pred, cd.type, qsar.type) %>% 
  rename (type = qsar.type) %>% data.table()
all.results <- rbind(dock.temp, qsar.temp)
all <- transform(all.results, type = factor(
  type,
  levels = c("PyRx Docking", "GLMNet", "PLS", "SVM", "Random Forest", "Cubist")
))

# Graphs ------------------------------------------------------------------

# For the poster
dir.create("./graphs/2018 poster")

plot.2018(docking) + 
  labs(title = "Docking (PyRx) predictions")
ggsave("./graphs/2018 poster/2018-02-11 docking.png", dpi = 600)

plot.2018(svm) + 
  labs(title = "SVM predictions")
ggsave("./graphs/2018 poster/2018-02-11 svm.png", dpi = 600)

plot.2018(rforest) + 
  labs(title = "Random Forest predictions")
ggsave("./graphs/2018 poster/2018-02-11 rforest.png", dpi = 600)

plot.2018(cubist) + 
  labs(title = "Cubist predictions")
ggsave("./graphs/2018 poster/2018-02-11 cubist.png", dpi = 600)

plot.2018(glmnet) + 
  labs(title = "GLMNet predictions")
ggsave("./graphs/2018 poster/2018-02-11 glm.png", dpi = 600)

plot.2018(pls) + 
  labs(title = "PLS predictions")
ggsave("./graphs/2018 poster/2018-02-11 pls.png", dpi = 600)

plot.2018(qsar.results) + 
  facet_wrap(~qsar.type) + 
  labs(title = "QSAR predictions")
ggsave("./graphs/2018 poster/2018-02-11 qsar.png", dpi = 600)

plot.2018(all) +
  facet_wrap(~type) + 
  labs(title = "Docking vs. Various QSAR Predictions") + 
  geom_point(size = 0.5)
ggsave("./graphs/2018 poster/2018-02-12 docking and qsar.png", dpi = 600, scale = 1.5)

# For the manuscript, 2018
dir.create("./graphs/2018 paper")

plot.paper.2018(docking) 
ggsave("./graphs/2018 paper/docking results.png")

plot.paper.2018(qsar.results) +
  facet_wrap(~qsar.type)
ggsave("./graphs/2018 paper/qsar results.png", scale = 1.5)

plot.paper.2018(all) + 
  facet_wrap(~type) 
ggsave("./graphs/2018 paper/docking and qsar.png", scale = 1.5)

temp.all <- all %>% mutate(Data = "Test")
ext.val <- readRDS("./ext val results.RDS") %>% mutate(Data = "External validation") %>%
  rename(type = model)
tst.ev.all <- rbind(temp.all, ext.val) %>%
  filter(type != "PyRx Docking")

# plot.paper.2018(tst.ev.all) + 
#   facet_wrap(~type) + 
#   geom_point(aes(shape = Data)) + 

shapes <- c("16" = "Test", "0" = "External validation")  
ggplot(tst.ev.all, aes(x = obs, y = pred, color = cd.type, shape = Data)) +
  geom_point() + 
  theme.paper.2018 + 
  facet_wrap(~type) + 
  coord_fixed(xlim = c(-45,5), ylim = c(-45, 5)) +
  geom_abline(intercept = 0, slope = 1) + 
  labs(x = "Experimental dG, kJ/mol", y = "Predicted dG, kJ/mol",
       color = "CD Type") + 
  scale_shape_manual(values = c(0,16)) 
ggsave("./graphs/2018 paper/test and extval.png")

qsar.type <- c("Cubist", "GLMNet", "PLS", "Random Forest", "SVM")
qsar.time <- c(445.91, 4.314, 176.01, 635.8, 9.59)
times <- data.frame(temp.types, qsar.times)

ggplot(times, aes(x = qsar.type, y = qsar.time)) + 
  geom_bar(stat = "identity") + 
  theme.paper.2018 + 
  labs(x = "QSAR type", y = "Time used, seconds")
ggsave("./graphs/2018 paper/time.png", scale = 0.75)
