library(tidymodels)
set.seed(930093)

dataset <- iris %>% 
  filter (Species != "virginica") %>%
  mutate(Species = factor(Species))

initSplit <- initial_split(dataset)
trainset <- training(initSplit)
testset <- testing(initSplit)

cv_splits <- rsample::vfold_cv(trainset)

mod <- logistic_reg(penalty = tune(),
                    mixture = tune()) %>%
  set_engine("glmnet")

rec <- recipe(Species ~ ., data = trainset) %>%
  step_center(all_numeric()) %>%
  step_scale(all_numeric())


wfl <- 
  workflow() %>%
  add_recipe(rec) %>%
  add_model(mod)

glmn_set <- parameters(penalty(range = c(-5,2), trans = log10_trans()),
                       mixture())

glmn_grid <- 
  grid_regular(glmn_set, levels = c(8, 5))

ctrl <- control_grid(save_pred = FALSE, verbose = TRUE)

## @knitr fitting

glmn_tune <- 
  tune_grid(wfl,
            resamples = cv_splits,
            grid = glmn_grid,
            metrics = metric_set(roc_auc),
            control = ctrl)

## finalize model

best_glmn <- select_best(glmn_tune, metric = "roc_auc")

wfl_final <- 
  wfl %>%
  finalize_workflow(best_glmn) %>%
  fit(data = trainset)

## @knitr save

save(wfl_final, file="wfl_final.RData")

## @knitr load

load("wfl_final.RData")

tidy(extract_model(wfl_final))
