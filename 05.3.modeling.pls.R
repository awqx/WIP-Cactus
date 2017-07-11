install.packages("pls")
library(pls)
library(tidyverse)


pls.data <- readRDS("./rpt.RDS") %>% data.frame()
pls.data <- pls.data[!is.na(pls.data$log.K), c(7, 18:758)]
pls.data <- rfImpute(x = pls.data[ , -1], 
                     y = pls.data[ , 1], 
                     ntree = 400, 
                     na.action = na.omit)
colnames(pls.data)[1] <- "logk"

set.seed(11)
trn.ind <- sample(x = 1:nrow(pls.data), 
                  size = round(0.8 * nrow(pls.data)))
pls.trn <- pls.data[trn.ind, ]
pls.tst <- pls.data[-trn.ind, ]

# Creating a basic PLS model
pls.mod <- plsr(logk ~ ., 
                ncomp = 10, 
                data = pls.trn, 
                validation = "LOO", 
                method = "oscorespls")
summary(pls.mod)
# Plotting number of components vs. RMSEP
plot(RMSEP(pls.mod), legendpos = "topright") # 3 comps best
# Plotting one of the validation sets
plot(pls.mod, ncomp = 3, asp = 1, line = T)
# Test data prediction plot
predict(pls.mod, ncomp = 3, newdata = pls.tst) %>%
  cbind(pls.tst[ , 1]) %>%
  data.frame() %>%
  ggplot(., aes(x = ., y = V2)) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0) + 
  coord_fixed()
# Testing R-squared
pls.results <- predict(pls.mod, ncomp = 3, newdata = pls.tst) %>%
  cbind(pls.tst[ , 1]) %>%
  data.frame()
colnames(pls.results)[1] <- "pred"
colnames(pls.results)[2] <- "obs"
defaultSummary(pls.results) # R-squared = 0.553

pls.results <- predict(pls.mod, ncomp = 3, newdata = pls.trn) %>%
  cbind(pls.trn[ , 1]) %>%
  data.frame()
colnames(pls.results)[1] <- "pred"
colnames(pls.results)[2] <- "obs"
defaultSummary(pls.results) # R-squared = 0.560

# Creating a basic mvr model
mvr.mod <- mvr(logk ~ ., 
               ncomp = 20, 
               data = pls.vip.trn, 
               validation = "LOO")
summary(mvr.mod)
# Plotting number of components vs. RMSEP
plot(RMSEP(mvr.mod), legendpos = "topright") # 3 comps best
# Plotting one of the validation sets
plot(mvr.mod, ncomp = 3, asp = 1, line = T)
# Test data prediction plot
predict(mvr.mod, ncomp = 3, newdata = pls.tst) %>%
  cbind(pls.tst[ , 1]) %>%
  data.frame() %>%
  ggplot(., aes(x = ., y = V2)) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0) + 
  coord_fixed()
# Testing R-squared
mvr.results <- predict(mvr.mod, ncomp = 3, newdata = pls.tst) %>%
  cbind(pls.tst[ , 1]) %>%
  data.frame()
colnames(mvr.results)[1] <- "pred"
colnames(mvr.results)[2] <- "obs"
defaultSummary(mvr.results) # R-squared = 0.553

mvr.results <- predict(mvr.mod, ncomp = 3, newdata = pls.trn) %>%
  cbind(pls.trn[ , 1]) %>%
  data.frame()
colnames(mvr.results)[1] <- "pred"
colnames(mvr.results)[2] <- "obs"
defaultSummary(mvr.results) # R-squared = 0.560

# VIP Analysis ------------------------------------------------------------
# Function sourced from http://mevik.net/work/software/pls.html
VIP <- function(object) {
  if (object$method != "oscorespls")
    stop("Only implemented for orthogonal scores algorithm.  Refit with 'method = \"oscorespls\"'")
  if (nrow(object$Yloadings) > 1)
    stop("Only implemented for single-response models")
  
  SS <- c(object$Yloadings)^2 * colSums(object$scores^2)
  Wnorm2 <- colSums(object$loading.weights^2)
  SSW <- sweep(object$loading.weights^2, 2, SS / Wnorm2, "*")
  sqrt(nrow(SSW) * apply(SSW, 1, cumsum) / cumsum(SS))
}

sort.vip <- function(vip.df, row.num) {
  long <- vip.df[row.num, ] %>% gather()
  return(long[long[ , 2] > 1.2, ])
}
pls.vip <- VIP(pls.mod) %>% data.frame()
pls.vip1 <- sort.vip(pls.vip, 1)
pls.vip2 <- sort.vip(pls.vip, 2)
pls.vip3 <- sort.vip(pls.vip, 3)
pls.vip4 <- sort.vip(pls.vip, 4)
pls.vip5 <- sort.vip(pls.vip, 5)
pls.vip6 <- sort.vip(pls.vip, 6)
pls.vip7 <- sort.vip(pls.vip, 7)
pls.vip8 <- sort.vip(pls.vip, 8)
pls.vip9 <- sort.vip(pls.vip, 9)
pls.vip10 <- sort.vip(pls.vip, 10)

highest.vip <- inner_join(pls.vip1, pls.vip2, by = "key") %>%
  inner_join(., pls.vip3, by = "key") %>%
  inner_join(., pls.vip4, by = "key") %>%
  inner_join(., pls.vip5, by = "key") %>%
  inner_join(., pls.vip6, by = "key") %>%
  inner_join(., pls.vip7, by = "key") %>%
  inner_join(., pls.vip8, by = "key") %>%
  inner_join(., pls.vip9, by = "key") %>%
  inner_join(., pls.vip10, by = "key") 
highest.vip.desc <- highest.vip$key

# VIP Model ---------------------------------------------------------------

pls.vip.data <- pls.data[ , colnames(pls.data) %in% highest.vip.desc] %>% cbind(pls.data[ , 1], .)
colnames(pls.vip.data)[1] <- "logk"
set.seed(201)
trn.ind <- sample(x = 1:nrow(pls.vip.data), 
                  size = round(0.8 * nrow(pls.vip.data)))
pls.vip.trn <- pls.vip.data[trn.ind, ]
pls.vip.tst <- pls.vip.data[-trn.ind, ]

pls.vip.mod <- plsr(logk ~ ., 
                    ncomp = 10, 
                    data = pls.vip.trn, 
                    validation = "LOO", 
                    method = "oscorespls")
summary(pls.vip.mod)
# Plotting number of components vs. RMSEP
plot(RMSEP(pls.vip.mod), legendpos = "topright") # 2 comps best
# Plotting one of the validation sets
plot(pls.vip.mod, ncomp = 2, asp = 1, line = T)
# Test data prediction plot
predict(pls.vip.mod, ncomp = 2, newdata = pls.vip.tst) %>%
  cbind(pls.vip.tst[ , 1]) %>%
  data.frame() %>%
  ggplot(., aes(x = ., y = V2)) + 
  geom_point() + 
  theme_bw() +
  geom_abline(slope = 1, intercept = 0) + 
  coord_fixed()
# Testing R-squared
pls.vip.results <- predict(pls.vip.mod, ncomp = 3, newdata = pls.vip.tst) %>%
  cbind(pls.vip.tst[ , 1]) %>%
  data.frame()
colnames(pls.vip.results)[1] <- "pred"
colnames(pls.vip.results)[2] <- "obs"
defaultSummary(pls.vip.results) # 0.582

VIP(pls.vip.mod)
sort.vip2 <- function(vip.df, row.num) {
  long <- vip.df[row.num, ] %>% gather()
  return(long[long[ , 2] > 0.9, ])
}

pls.vip <- VIP(pls.vip.mod) %>% data.frame()
pls.vip1 <- sort.vip2(pls.vip, 1)
pls.vip2 <- sort.vip2(pls.vip, 2)
pls.vip3 <- sort.vip2(pls.vip, 3)
pls.vip4 <- sort.vip2(pls.vip, 4)
pls.vip5 <- sort.vip2(pls.vip, 5)
pls.vip6 <- sort.vip2(pls.vip, 6)
pls.vip7 <- sort.vip2(pls.vip, 7)
pls.vip8 <- sort.vip2(pls.vip, 8)
pls.vip9 <- sort.vip2(pls.vip, 9)
pls.vip10 <- sort.vip2(pls.vip, 10)

highest.next.vip <- inner_join(pls.vip1, pls.vip2, by = "key") %>%
  inner_join(., pls.vip3, by = "key") %>%
  inner_join(., pls.vip4, by = "key") %>%
  inner_join(., pls.vip5, by = "key") %>%
  inner_join(., pls.vip6, by = "key") %>%
  inner_join(., pls.vip7, by = "key") %>%
  inner_join(., pls.vip8, by = "key") %>%
  inner_join(., pls.vip9, by = "key") %>%
  inner_join(., pls.vip10, by = "key") 
highest.next.vip.desc <- highest.next.vip$key

