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
  levels = c("PyRx Docking", "GLMnet", "PLS", "SVM", "Random Forest", "Cubist")
))

# Graphs ------------------------------------------------------------------

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
ggsave("./graphs/2018 poster/2018-02-12 docking and qsar.png", dpi = 600, scale = 1.05)
