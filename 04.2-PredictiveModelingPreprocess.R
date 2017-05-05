# =============================================

#                   PreProcess

# =============================================

# Checks for predictors that don't vary much
# and then removes them
# note: rp stands for ri.apdel
# and rpt stands for ri.padel.trans 
set.seed(999)
zero.pred <- nearZeroVar(ri.padel)
rp.no.zero <- ri.padel[, -zero.pred] # View(rp.no.zero)

rpt <- preProcess(rp.no.zero, na.remove = T, 
                  method = c("center", "scale"
                             # , "BoxCox" throws an error: Na/NaN/Inf in foreign function call
                  )) %>%
  predict(., ri.padel) 
rpt.bp <-  rpt %>%
  filter(host == "Beta-CD") %>% 
  filter(!is.na(log.K)) %>%
  dplyr::select(., -X1:-Solvent.composition, -log.K.Uncertainty:-`bind.aff, kcal/mol`)

rpt.btherm <-  rpt %>%
  filter(host == "Beta-CD") %>%
  dplyr::select(X1:Solvent.composition, log.K.Uncertainty:`bind.aff, kcal/mol`)
rpt.sprse <- sparse.model.matrix(~., rpt.bp)
rpt.trn.ind <- sample(x = 1:nrow(rpt.sprse), size = round(0.75*nrow(rpt.sprse)))

rpt.trn <- rpt.sprse[rpt.trn.ind, ]
rpt.tst <- rpt.sprse[-rpt.trn.ind,]

rpt.glm.2 <- cv.glmnet(x = rpt.trn[,-1:-6], y = rpt.trn[, 6])

predict(rpt.glm.2, rpt.tst[,-1:-6]) %>% 
  cbind(rpt.tst[,6]) %>% 
  data.frame() %>% 
  rename(predicted = X1, experimental = V2) %>%
  ggplot(., aes(x = experimental, y = predicted))+
  geom_point()+
  geom_abline(slope = 1, intercept = 0)
