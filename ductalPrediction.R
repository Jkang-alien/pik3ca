## @knitr ductal_diffset

load('ductalset.RData')
library(purrr)
library(tidyverse)
library(ggplot2)
library(tidymodels)

init_split <- initial_split(ductalSet, strata = PIK3CA_T)
trainset_ah <- training(init_split)
testset_ah <- testing(init_split)

variance <- trainset_ah %>%
  select(-PIK3CA_T) %>%
  map_dbl(var)

sum(is.na(variance))
summary(variance)
summary(variance[variance > median(variance)])

trainset_ahLargeVariance <- trainset_ah[,-1][,variance > median(variance)]
trainset_ahLargeVariance$PIK3CA_T <- trainset_ah$PIK3CA_T

testset_ahLargeVariance <- testset_ah[,-1][,variance > median(variance)]
testset_ahLargeVariance$PIK3CA_T <- testset_ah$PIK3CA_T

lobularsetLargeVariance <- lobularSet[,-1][,variance > median(variance)]
lobularsetLargeVariance$PIK3CA_T <- lobularSet$PIK3CA_T

stasticsWilcox <- vector(mode = "list", length = dim(trainset_ahLargeVariance)[2])
pvalueWilcox <- vector(mode = "list", length = dim(trainset_ahLargeVariance)[2])

for (i in 1:(dim(trainset_ahLargeVariance)[2]-1)) {
  a <- wilcox.test(trainset_ahLargeVariance[,i] ~ trainset_ahLargeVariance$PIK3CA_T)
  stasticsWilcox[[i]] <- a$statistic
  pvalueWilcox[[i]] <- a$p.value
}

colnames(trainset_ahLargeVariance)[log(unlist(pvalueWilcox), base = 10) < -8]

trainset_ahDiff <- trainset_ahLargeVariance[,log(unlist(pvalueWilcox), base = 10) < -8]

trainset_ahDiff$PIK3CA_T<- trainset_ah$PIK3CA_T

testset_ahDiff <- testset_ahLargeVariance[,log(unlist(pvalueWilcox), base = 10) < -8]

testset_ahDiff$PIK3CA_T<- testset_ah$PIK3CA_T

lobularset_Diff <- lobularsetLargeVariance[,log(unlist(pvalueWilcox), base = 10) < -8]

lobularset_Diff$PIK3CA_T<- lobularSet$PIK3CA_T

summary(trainset_ah$PIK3CA_T)
summary(testset_ah$PIK3CA_T)
summary(lobularSet$PIK3CA_T)

save(trainset_ahDiff, testset_ahDiff, lobularset_Diff, file = "ductalset_ah_Diff.RData")

## @knitr all_histology_modeling


load("ductalset_ah_Diff.RData")


library(tidymodels)
library(glmnet)
library(skimr)

# skim(trainset_ahDiff)

set.seed(930093)

cv_splits <- rsample::vfold_cv(trainset_ahDiff, strata = PIK3CA_T, repeats = 5)

mod <- logistic_reg(penalty = tune(),
                    mixture = tune()) %>%
  set_engine("glmnet")

rec <- recipe(PIK3CA_T ~ ., data = trainset_ahDiff) %>%
  step_BoxCox(all_predictors()) %>%
  step_center(all_predictors()) %>%
  step_scale(all_predictors()) %>%
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

lasso_glmn <- show_best(glmn_tune)[4,1:2]
best_glmn <- select_best(glmn_tune)
mix <- show_best(glmn_tune)[3,1:2]

## @knitr Finalize_model

wfl_final <- 
  wfl %>%
  finalize_workflow(best_glmn) %>%
  fit(data = trainset_ahDiff)

## @knitr trainset_prediction

train_probs <- 
  predict(wfl_final, type = "prob", new_data = trainset_ahDiff) %>%
  bind_cols(obs = trainset_ahDiff$PIK3CA_T) %>%
  bind_cols(predict(wfl_final, new_data = trainset_ahDiff))

confusion_matrix <- conf_mat(train_probs, obs, .pred_class)

roc_curve_train <- autoplot(roc_curve(train_probs, obs, .pred_Mutant))

roc_auc_train <- roc_auc(train_probs, obs, .pred_Mutant)

## @knitr testset_prediction

test_probs <- 
  predict(wfl_final, type = "prob", new_data = testset_ahDiff) %>%
  bind_cols(obs = testset_ahDiff$PIK3CA_T) %>%
  bind_cols(predict(wfl_final, new_data = testset_ahDiff))

conf_mat(test_probs, obs, .pred_class)
autoplot(roc_curve(test_probs, obs, .pred_Mutant))

roc_auc(test_probs, obs, .pred_Mutant)

## @knitr lobular_prediction

lobular_probs <- 
  predict(wfl_final, type = "prob", new_data = lobularset_Diff) %>%
  bind_cols(obs = lobularset_Diff$PIK3CA_T) %>%
  bind_cols(predict(wfl_final, new_data = lobularset_Diff))

conf_mat(lobular_probs, obs, .pred_class)
autoplot(roc_curve(lobular_probs, obs, .pred_Mutant))

roc_auc(lobular_probs, obs, .pred_Mutant)

