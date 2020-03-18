## @knitr all_histology_diffset

load('datasetHistology.RData')
library(purrr)
library(tidyverse)
library(ggplot2)
library(tidymodels)

init_split <- initial_split(dataset, strata = PIK3CA_T)
trainset_ah <- training(init_split)
testset_ah <- testing(init_split)
  
variance <- trainset_ah %>%
  select(-PIK3CA_T, -ID, -HISTOLOGICAL_DIAGNOSIS) %>%
  map_dbl(var)

sum(is.na(variance))
summary(variance)
summary(variance[variance > median(variance)])

trainset_ahLargeVariance <- trainset_ah[,-1:-3][,variance > median(variance)]
trainset_ahLargeVariance$PIK3CA_T <- trainset_ah$PIK3CA_T
#trainset_ahLargeVariance$HISTOLOGICAL_DIAGNOSIS <- trainset_ah$HISTOLOGICAL_DIAGNOSIS

testset_ahLargeVariance <- testset_ah[,-1:-3][,variance > median(variance)]
testset_ahLargeVariance$PIK3CA_T <- testset_ah$PIK3CA_T
testset_ahLargeVariance$HISTOLOGICAL_DIAGNOSIS <- testset_ah$HISTOLOGICAL_DIAGNOSIS

stasticsWilcox <- vector(mode = "list", length = dim(trainset_ahLargeVariance)[2]-1)
pvalueWilcox <- vector(mode = "list", length = dim(trainset_ahLargeVariance)[2]-1)

for (i in 1:(dim(trainset_ahLargeVariance)[2]-1)) {
  a <- wilcox.test(trainset_ahLargeVariance[,i] ~ trainset_ahLargeVariance$PIK3CA_T)
  stasticsWilcox[[i]] <- a$statistic
  pvalueWilcox[[i]] <- a$p.value
}

colnames(trainset_ahLargeVariance)[log(unlist(pvalueWilcox), base = 10) < -8]

trainset_ahDiff <- trainset_ahLargeVariance[,log(unlist(pvalueWilcox), base = 10) < -8]

trainset_ahDiff$PIK3CA_T <- trainset_ah$PIK3CA_T

trainset_ahDiff$HISTOLOGICAL_DIAGNOSIS <- trainset_ah$HISTOLOGICAL_DIAGNOSIS

testset_ahDiff <- testset_ahLargeVariance[,log(unlist(pvalueWilcox), base = 10) < -8]

testset_ahDiff$PIK3CA_T<- testset_ah$PIK3CA_T
testset_ahDiff$HISTOLOGICAL_DIAGNOSIS <- testset_ah$HISTOLOGICAL_DIAGNOSIS

testset_lobular <- testset_ahDiff %>%
  filter(HISTOLOGICAL_DIAGNOSIS == "Infiltrating Lobular Carcinoma")# %>%
  #select(-HISTOLOGICAL_DIAGNOSIS)


testset_ductal <- testset_ahDiff %>%
  filter(HISTOLOGICAL_DIAGNOSIS == "Infiltrating Ductal Carcinoma")# %>%
  #select(-HISTOLOGICAL_DIAGNOSIS)

  
summary(trainset_ah$PIK3CA_T)
summary(testset_ah$PIK3CA_T)
summary(testset_ductal$PIK3CA_T)
summary(testset_lobular$PIK3CA_T)

save(trainset_ahDiff, testset_ahDiff, testset_ductal, testset_lobular,  file = "dataset_all_test.RData")

## @knitr all_histology_modeling


load("dataset_all_test.RData")


library(tidymodels)
library(glmnet)
library(skimr)

#skim(trainset_ahDiff)

set.seed(930093)

cv_splits <- rsample::vfold_cv(trainset_ahDiff, strata = PIK3CA_T, repeats = 5)

mod <- logistic_reg(penalty = tune(),
                    mixture = tune()) %>%
  set_engine("glmnet")

rec <- recipe(PIK3CA_T ~ ., data = trainset_ahDiff) %>%
  step_BoxCox(all_numeric()) %>%
  step_dummy(HISTOLOGICAL_DIAGNOSIS) %>%
  step_center(all_numeric()) %>%
  step_scale(all_numeric()) %>%
  step_downsample(PIK3CA_T)

wfl <- 
  workflow() %>%
  add_recipe(rec) %>%
  add_model(mod)

glmn_set <- parameters(penalty(range = c(-5,3), trans = log10_trans()),
                       mixture())

glmn_grid <- 
  grid_regular(glmn_set, levels = c(9, 5))

ctrl <- control_grid(save_pred = TRUE, verbose = TRUE)

## @knitr fitting

glmn_tune <- 
  tune_grid(wfl,
            resamples = cv_splits,
            grid = glmn_grid,
            metrics = metric_set(roc_auc),
            control = ctrl)


show_best(glmn_tune)

best_glmn <- select_best(glmn_tune)

## @knitr Finalize_model

wfl_final <- 
  wfl %>%
  finalize_workflow(best_glmn) %>%
  fit(data = trainset_ahDiff)

model <- extract_model(wfl_final)
model[2]

tidy(extract_model(wfl_final)) %>%
  filter(lambda > 0.98 & lambda < 1.01)

## @knitr trainset_prediction
train_predict <- stats::predict(wfl_final, type = "prob", new_data = trainset_ahDiff)
train_predict$
train_probs <- 
  predict(wfl_final, type = "prob", new_data = trainset_ahDiff) %>%
  bind_cols(obs = trainset_ahDiff$PIK3CA_T) %>%
  bind_cols(predict(wfl_final, new_data = trainset_ahDiff))

confusion_matrix <- conf_mat(train_probs, obs, .pred_class)

roc_curve_train <- autoplot(roc_curve(train_probs, obs, .pred_Mutant))

roc_auc_train <- roc_auc(train_probs, obs, .pred_Mutant)

confusion_matrix

roc_curve_train

roc_auc_train

tidypredict::tidypredict_fit(wfl_final$fit$fit)

## @knitr testset_prediction

ductal_probs <- 
  predict(wfl_final, type = "prob", new_data = testset_ductal) %>%
  bind_cols(obs = testset_ductal$PIK3CA_T) %>%
  bind_cols(predict(wfl_final, new_data = testset_ductal))

conf_mat(ductal_probs, obs, .pred_class)
autoplot(roc_curve(ductal_probs, obs, .pred_Mutant))

roc_auc(ductal_probs, obs, .pred_Mutant)



lobular_probs <- 
  predict(wfl_final, type = "prob", new_data = testset_lobular) %>%
  bind_cols(obs = testset_lobular$PIK3CA_T) %>%
  bind_cols(predict(wfl_final, new_data = testset_lobular))

conf_mat(lobular_probs, obs, .pred_class)
autoplot(roc_curve(lobular_probs, obs, .pred_Mutant))

roc_auc(lobular_probs, obs, .pred_Mutant)
