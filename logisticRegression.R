## @knitr modeling

load("datasetDiff.RData")
load("glmn_tune.RData")

library(tidymodels)
library(glmnet)
library(skimr)

# skim(trainSetDiff)

set.seed(930093)

cv_splits <- rsample::vfold_cv(trainSetDiff, strata = PIK3CA_T, repeats = 5)

mod <- logistic_reg(penalty = tune(),
                    mixture = tune()) %>%
  set_engine("glmnet")

rec <- recipe(PIK3CA_T ~ ., data = trainSetDiff) %>%
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


