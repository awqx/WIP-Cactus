library(caret)
library(keras)
install_tensorflow()

# Data Organization -------------------------------------------------------

df <- readRDS("./rpt.RDS") %>% 
  select(., -X1:-log.K.Uncertainty) %>%
  select(., -DelG.Uncertainty:-ref) %>%
  select(., -`bind.aff, kcal/mol`)
mat <- as.matrix(df)

set.seed(12)
trn.ind <- sample(x = 1:nrow(df), size = round(0.7 * nrow(df)))
trn.x <- mat[trn.ind, -1] 
trn.y <- mat[trn.ind, 1]
tst.x <- mat[-trn.ind, -1]
tst.y <- mat[-trn.ind, 1]


# Initial Keras model -----------------------------------------------------

model <- keras_model_sequential()
model %>%
  layer_dense(units = 128, activation = "relu", input_shape = c(752)) %>%
  layer_dense(units = 1, activation = "linear")
model %>% compile(
  loss = "mean_squared_error", 
  optimizer = "adam", 
  metrics = "accuracy"
)
model %>% fit(
  trn.x, 
  trn.y, 
  epochs = 400, 
  batch_size = 16, 
  validation_split = 0.2
)

result <- model %>%
  predict_proba(tst.x, batch_size = 128)
keras.df <- data.frame(tst.y, result) %>%
  rename(pred = result, obs = tst.y)
defaultSummary(keras.df)
ggplot(keras.df, aes(x = obs, y = pred)) +
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) + 
  coord_fixed() + 
  theme_bw()

