## @knitr exploreCVResult
library(tidyverse)
library(tidymodels)

load(file = "glmn_tune_pancancer.RData")


glmn_tune_tidy <- glmn_tune %>%
  select(-splits)

best_roc <- show_best(glmn_tune, metric = "roc_auc")
best_pr <- show_best(glmn_tune, metric = "pr_auc")

best_roc_down <- show_best(glmn_tune_down, metric = "roc_auc")
best_pr_down <- show_best(glmn_tune_down, metric = "pr_auc")

lasso_glmn <- show_best(glmn_tune)[4,1:2]
best_glmn <- select_best(glmn_tune, metric = "roc_auc")
mix <- show_best(glmn_tune)[3,1:2]

best_glmn_down <- select_best(glmn_tune_down, metric = "roc_auc")
## @knitr Finalize_model

wfl_final <- 
  wfl %>%
  finalize_workflow(best_glmn) %>%
  fit(data = trainset)

wfl_down_final <- 
  wfl_down %>%
  finalize_workflow(best_glmn_down) %>%
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

coeffPlot <- head(coeff, n=30) %>%
  filter(term != "(Intercept)") %>%
  ggplot(aes(x=reorder(term, -estimate), y=estimate)) +
    geom_bar(stat = "identity") +
    coord_flip()

typeCoeff <- tidy(extract_model(wfl_final)) %>%
  filter(lambda > 0.09 & lambda < 0.104) %>%
  arrange(-abs(estimate)) %>%
  filter(grepl("type_", term))
typeCoeff$term <- gsub("type_", "", typeCoeff$term)

typePlot <- typeCoeff %>%
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
  roc_auc(obs, .pred_FALSE)

roc_auc_train$type <- "ALL"

roc_auc_train_type <- train_probs %>%
  group_by(type) %>%
  roc_auc(obs, .pred_Mutant) %>%
  mutate(type = as.character(type))

roc_plot_type <- roc_auc_train_type %>%
  filter(is.na(.estimate) == FALSE) %>%
  bind_rows(roc_auc_train) %>%
  ggplot(aes(x=reorder(type, -.estimate), y=.estimate)) +
  geom_bar(stat="identity") + 
  coord_flip()

#options(yardstick.event_first = TRUE)

pr_auc_train <- train_probs %>%
  pr_auc(obs, .pred_Mutant)

pr_auc_train

pr_curve_train <- autoplot(pr_curve(train_probs, obs, .pred_Mutant))


pr_auc_train$type <- "ALL"

pr_auc_train_type <- train_probs %>%
  group_by(type) %>%
  pr_auc(obs, .pred_Mutant) %>%
  mutate(type = as.character(type))

pr_plot_type <- pr_auc_train_type %>%
  filter(is.na(.estimate) == FALSE) %>%
  bind_rows(pr_auc_train) %>%
  ggplot(aes(x=reorder(type, -.estimate), y=.estimate)) +
  geom_bar(stat="identity") + 
  coord_flip()

## @knitr correlation

mutationRate <- trainset %>%
  select(type, variant) %>%
  group_by(type) %>%
  count(variant == "Mutant") %>%
  rename(variant = 'variant == "Mutant"') %>%
  spread(variant, n) %>%
  mutate(rate = `TRUE`/(`TRUE`+`FALSE`)) %>%
  bind_cols(roc_auc_train_type)

mutatioPlot <- mutationRate %>%
  ggplot(aes(x=rate, y=.estimate)) +
  geom_point() +
  geom_smooth(method='lm', formula= y~x) +
  ggrepel::geom_text_repel(data = mutationRate, aes(label = type))


mutationRatePR <- trainset %>%
  select(type, variant) %>%
  group_by(type) %>%
  count(variant == "Mutant") %>%
  rename(variant = 'variant == "Mutant"') %>%
  spread(variant, n) %>%
  mutate(rate = `TRUE`/(`TRUE`+`FALSE`)) %>%
  bind_cols(pr_auc_train_type)

mutatioPlotPR <- mutationRatePR %>%
  ggplot(aes(x=rate, y=.estimate)) +
  geom_point() +
  geom_smooth(method='lm', formula= y~x) +
  ggrepel::geom_text_repel(data = mutationRatePR, aes(label = type))

## @knitr testset_prediction

test_probs <- 
  predict(wfl_final, type = "prob", new_data = testset) %>%
  bind_cols(obs = testset$variant) %>%
  bind_cols(predict(wfl_final, new_data = testset)) %>%
  bind_cols(type = testset$type)

test_probs_down <- 
  predict(wfl_down_final, type = "prob", new_data = testset) %>%
  bind_cols(obs = testset$variant) %>%
  bind_cols(predict(wfl_down_final, new_data = testset)) %>%
  bind_cols(type = testset$type)

confusion_matrix_test <- test_probs_down %>%
  #group_by(type) %>%
  conf_mat(obs, .pred_class)

confusion_matrix_test

confusion_matrix_test_type <- test_probs %>%
  group_by(type) %>%
  conf_mat(obs, .pred_class)

confusion_matrix_test_type[1,2][[1]]

roc_curve_test <- autoplot(roc_curve(test_probs, obs, .pred_Mutant))

roc_curve_test

roc_auc_test <- test_probs %>%
  roc_auc(obs, .pred_Mutant)

roc_auc_test

pr_auc_test <- test_probs %>%
  pr_auc(obs, .pred_Mutant)

pr_auc_test

pr_curve_test <- autoplot(pr_curve(test_probs, obs, .pred_Mutant))

pr_auc_test_down <- test_probs_down %>%
  pr_auc(obs, .pred_Mutant)

pr_auc_test

pr_curve_test_down <- autoplot(pr_curve(test_probs_down, obs, .pred_Mutant))

pr_curve_test

pr_auc_test$type <- "ALL"

pr_auc_test_type <- test_probs %>%
  group_by(type) %>%
  pr_auc(obs, .pred_Mutant) %>%
  mutate(type = as.character(type))

pr_plot_test_type <- pr_auc_test_type %>%
  bind_rows(pr_auc_test) %>%
  filter(is.na(.estimate) == FALSE) %>%
  ggplot(aes(x=reorder(type, -.estimate), y=.estimate)) +
  geom_bar(stat="identity") + 
  coord_flip()

pr_train_test <- merge(pr_auc_train_type,
                       pr_auc_test_type,
                       by = "type")

pr_test_trainPlot <- pr_train_test %>%
  ggplot(aes(x=.estimate.x, y=.estimate.y)) +
  geom_point() +
  geom_smooth(method='lm', formula= y~x) +
  ggrepel::geom_text_repel(data = pr_train_test,
                           aes(label = type))

pr_test_trainPlot

roc_train_test <- merge(roc_auc_train_type,
                       roc_auc_test_type,
                       by = "type")

roc_test_trainPlot <- roc_train_test %>%
  ggplot(aes(x=.estimate.x, y=.estimate.y)) +
  geom_point() +
  geom_smooth(method='lm', formula= y~x) +
  ggrepel::geom_text_repel(data = roc_train_test,
                           aes(label = type))

roc_test_trainPlot



