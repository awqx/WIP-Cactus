# Libraries and Packages --------------------------------------------------
library(data.table)
library(e1071)
library(RCurl)
library(stringr)
library(tidyverse)
library(XML)

# Loading Suzuki ----------------------------------------------------------
url   <- "http://pubs.acs.org/doi/full/10.1021/ci010295f"
file  <-
  getURL(url, .opts = curlOptions(followlocation = TRUE, cookiefile = "nosuchfile"))
suzuki.list <-
  readHTMLTable(
    file,
    header = T,
    as.data.frame = T,
    stringAsFactors = F
  )
index11 <- suzuki.list %>% lapply(names) %>% lapply(length) == 11 %>% as.vector()
suzuki.raw <- rbindlist(suzuki.list[index11], fill = T) %>% 
  lapply(as.character) %>% 
  as.data.frame(stringsAsFactors = F)
# Column names got shifted to row one, only relevant columns are guest, and DelG obsd
suzuki <- suzuki.raw[-1, c(2, 4, 8)]
colnames(suzuki) <- c("guest", "DelG.a", "DelG.b")

suzuki <- suzuki %>%
  mutate(DelG.a = str_replace(DelG.a, "Ã¢\u0088\u0092", "-")) %>%
  mutate(DelG.a = str_replace(DelG.a, "(Ã¢|Ã)", "")) %>%
  mutate(DelG.b = str_replace(DelG.b, "Ã¢\u0088\u0092", "-")) %>%
  mutate(DelG.b = str_replace(DelG.b, "(Ã¢|Ã)", "")) %>%
  mutate(guest = str_replace(guest, "Ã¢\u0080\u0089", "-")) 
h.ind <- suzuki$guest %>% str_detect("guest")
suzuki <- suzuki[!h.ind, ] 
suzuki[ , 2:3] <- sapply(suzuki[, 2:3], as.numeric) 
suzuki.dg.alpha <- suzuki %>% 
  filter(!is.na(DelG.a)) %>% 
  select(guest, DelG.a) %>% 
  rename(DelG = DelG.a)
suzuki.dg.beta <- suzuki %>%  
  select(guest, DelG.b) %>%
  rename(DelG = DelG.b)


suzuki.beta <- read_csv("C:/Users/Wei Xin/Desktop/2017 06 16 Suzuki beta.csv") %>% 
  mutate(alpha = c(rep(0, 185))) %>%
  mutate(beta = c(rep(1, 185))) %>%
  mutate(gamma = c(rep(0, 185))) %>%
  rename(guest = Name) %>%
  inner_join(suzuki.dg.beta, by = "guest", .)
suzuki.alpha <- read_csv("C:/Users/Wei Xin/Desktop/2017 06 16 Suzuki alpha comp.csv") %>%
  mutate(alpha = c(rep(1, 55))) %>%
  mutate(beta = c(rep(0, 55))) %>%
  mutate(gamma = c(rep(0, 55))) %>%
  rename(guest = Name) %>%
  inner_join(suzuki.dg.alpha, by = "guest", .)
suzuki.padel <- rbind(suzuki.alpha, suzuki.beta)
saveRDS(suzuki.padel, "./suzuki.padel.RDS")
colnames(suzuki.padel) <- colnames(suzuki.padel) %>% str_replace("-", "\\.")

sp.no.zero <- suzuki.padel[, !colnames(suzuki.padel) %in% zero.pred.names]

st <- preProcess(sp.no.zero[ , -1:-2], na.remove = T, 
                 method = c("knnImpute", "center", "scale")) %>%
  predict(., sp.no.zero[ , -1:-2])

st <- st[ , !colnames(st) %in% zero.pred2.names]
st.nocor <- st[ , !colnames(st) %in% corr.pred]
suzuki.transformed <- cbind(sp.no.zero[ , 1:2], st.nocor)
saveRDS(suzuki.transformed, "./suzuki.clean.RDS")

suzuki.x <- suzuki.transformed[ , -1:-2]
suzuki.y <- suzuki.transformed[ , 2]
suz.sprse <- sparse.model.matrix(~., suzuki.transformed[ , -1])

suz.sprse.a <- suz.sprse[suz.sprse[ , 752] > 0, ]
suz.sprse.b <- suz.sprse[suz.sprse[ , 753] > 0, ]

suz.a.trn <- predict(svm.a, suz.sprse.a[ , -1:-2]) %>% 
  cbind(suz.sprse.a[ , 2]) %>%
  data.frame() %>%
  rename(., pred = `.`, obs = V2)
ggplot(suz.a.trn, aes(x = obs, y = pred)) + 
  geom_point() + 
  coord_fixed() + 
  theme_bw() + 
  geom_abline(slope = 1, intercept = 0)
defaultSummary(suz.a.trn)

suz.b.trn <- predict(svm.b, suz.sprse.b[ , -1:-2]) %>%
  cbind(., suz.sprse.b[ , 2]) %>%
  data.frame() %>%
  rename(., pred = `.`, obs = V2)
ggplot(suz.b.trn, aes(x = obs, y = pred)) + 
  geom_point() + 
  coord_fixed() + 
  theme_bw() + 
  geom_abline(intercept = 0, slope = 1) 
defaultSummary(suz.b.trn)
