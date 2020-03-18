## @knitr load

library(tidymodels)
library(glmnet)

load("datasetDiff.RData")
load("glmn_tune.RData")

show_best(glmn_tune)

lasso_glmn <- show_best(glmn_tune)[4,1:2]
best_glmn <- select_best(glmn_tune)
mix <- show_best(glmn_tune)[3,1:2]

## @knitr Finalize_model

wfl_final <- 
  wfl %>%
  finalize_workflow(best_glmn) %>%
  fit(data = trainSetDiff)

## @knitr trainset_prediction

train_probs <- 
  predict(wfl_final, type = "prob", new_data = trainSetDiff) %>%
  bind_cols(obs = trainSetDiff$PIK3CA_T) %>%
  bind_cols(predict(wfl_final, new_data = trainSetDiff))

confusion_matrix <- conf_mat(train_probs, obs, .pred_class)

roc_curve_train <- autoplot(roc_curve(train_probs, obs, .pred_Mutant))

roc_auc_train <- roc_auc(train_probs, obs, .pred_Mutant)

## @knitr testset_prediction

test_probs <- 
  predict(wfl_final, type = "prob", new_data = testSetDiff) %>%
  bind_cols(obs = testSetDiff$PIK3CA_T) %>%
  bind_cols(predict(wfl_final, new_data = testSetDiff))

conf_mat(test_probs, obs, .pred_class)
autoplot(roc_curve(test_probs, obs, .pred_Mutant))

roc_auc(test_probs, obs, .pred_Mutant)
