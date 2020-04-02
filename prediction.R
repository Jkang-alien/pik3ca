## @knitr exploreCVResult
library(tidyverse)
library(tidymodels)
load(file = "wflFinal.RData")
load(file = "glmn_tune_pancancer.RData")


glmn_tune_tidy <- glmn_tune %>%
  select(-splits)

show_best(glmn_tune, metric = "roc_auc")
show_best(glmn_tune, metric = "pr_auc")

lasso_glmn <- show_best(glmn_tune)[4,1:2]
best_glmn <- select_best(glmn_tune, metric = "roc_auc")
mix <- show_best(glmn_tune)[3,1:2]

## @knitr Finalize_model

wfl_final <- 
  wfl %>%
  finalize_workflow(best_glmn) %>%
  fit(data = trainset)

roc_vals <-
  collect_metrics(glmn_tune) %>%
  filter(.metric == "roc_auc")

roc_vals %>%
  mutate(mixture = format(mixture)) %>%
  ggplot(aes(x = penalty, y = mean, col = mixture)) +
  geom_line() +
  geom_point() +
  scale_x_log10()

## @knitr coefficiency

coeff <- tidy(extract_model(wfl_final)) %>%
  filter(lambda > 0.09 & lambda < 0.104) %>%
  arrange(-abs(estimate)) %>%
  filter(!grepl("type_", term))

head(coeff, n=30) %>%
  filter(term != "(Intercept)") %>%
  ggplot(aes(x=reorder(term, -estimate), y=estimate)) +
    geom_bar(stat = "identity") +
    coord_flip()

tidy(extract_model(wfl_final)) %>%
  filter(lambda > 0.09 & lambda < 0.104) %>%
  arrange(-abs(estimate)) %>%
  filter(grepl("type_", term)) %>%
  ggplot(aes(x=reorder(term, -estimate), y=estimate)) +
  geom_bar(stat = "identity") +
  coord_flip()

## @knitr trainset_prediction

train_probs <- 
  predict(wfl_final, type = "prob", new_data = trainset) %>%
  bind_cols(obs = trainset$variant) %>%
  bind_cols(predict(wfl_final, new_data = trainset)) %>%
  bind_cols(type = trainset$type)

confusion_matrix <- train_probs %>%
  #group_by(type) %>%
  conf_mat(obs, .pred_class)

confusion_matrix

confusion_matrix_type <- train_probs %>%
  group_by(type) %>%
  conf_mat(obs, .pred_class)

confusion_matrix[1,2][[1]]

roc_curve_train <- autoplot(roc_curve(train_probs, obs, .pred_TRUE))

roc_auc_train <- train_probs %>%
  roc_auc(obs, .pred_TRUE)

roc_auc_train$type <- "ALL"

roc_auc_train_type <- train_probs %>%
  group_by(type) %>%
  roc_auc(obs, .pred_TRUE) %>%
  mutate(type = as.character(type))

roc_plot_type <- roc_auc_train_type %>%
  bind_rows(roc_auc_train) %>%
  ggplot(aes(x=reorder(type, -.estimate), y=.estimate)) +
  geom_bar(stat="identity") + 
  coord_flip()

roc_auc_train <- train_probs %>%
  roc_auc(obs, .pred_TRUE)

roc_auc_train

## @knitr correlation

mutationRate <- trainset %>%
  select(type, variant) %>%
  group_by(type) %>%
  count(variant == TRUE) %>%
  rename(variant = 'variant == TRUE') %>%
  spread(variant, n) %>%
  mutate(rate = `TRUE`/(`TRUE`+`FALSE`)) %>%
  bind_cols(roc_auc_train_type)

mutationRate %>%
  ggplot(aes(x=rate, y=.estimate)) +
  geom_point() +
  geom_smooth(method='lm', formula= y~x) +
  ggrepel::geom_text_repel(data = mutationRate, aes(label = type))

## @knitr testset_prediction

test_probs <- 
  predict(wfl_final, type = "prob", new_data = testset) %>%
  bind_cols(obs = testset$variant) %>%
  bind_cols(predict(wfl_final, new_data = testset)) %>%
  bind_cols(type = testset$type)

confusion_matrix_test <- test_probs %>%
  #group_by(type) %>%
  conf_mat(obs, .pred_class)

confusion_matrix_test

confusion_matrix_test_type <- test_probs %>%
  group_by(type) %>%
  conf_mat(obs, .pred_class)

confusion_matrix_test_type[1,2][[1]]

roc_curve_test <- autoplot(roc_curve(test_probs, obs, .pred_TRUE))

roc_curve_test

roc_auc_test <- test_probs %>%
  roc_auc(obs, .pred_TRUE)

roc_auc_test$type <- "ALL"

roc_auc_test_type <- test_probs %>%
  group_by(type) %>%
  roc_auc(obs, .pred_TRUE) %>%
  mutate(type = as.character(type))

roc_plot_test_type <- roc_auc_test_type %>%
  bind_rows(roc_auc_test) %>%
  ggplot(aes(x=reorder(type, -.estimate), y=.estimate)) +
  geom_bar(stat="identity") + 
  coord_flip()

plot(roc_auc_train_type$.estimate, roc_auc_test_type$.estimate)

roc_plot_test_type

roc_auc_test <- test_probs %>%
  roc_auc(obs, .pred_TRUE)

roc_auc_test

## @knitr correlation testset

mutationRate <- testset %>%
  select(type, variant) %>%
  group_by(type) %>%
  count(variant == TRUE) %>%
  rename(variant = 'variant == TRUE') %>%
  spread(variant, n) %>%
  mutate(rate = `TRUE`/(`TRUE`+`FALSE`)) %>%
  bind_cols(roc_auc_test_type)

mutationRate %>%
  ggplot(aes(x=rate, y=.estimate)) +
  geom_point() +
  geom_smooth(method='lm', formula= y~x) +
  ggrepel::geom_text_repel(data = mutationRate, aes(label = type))
