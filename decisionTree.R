load("trainset.RData")
load("cart_tune.RData")

library(tidymodels)
library(rpart.plot)
library(partykit)

set.seed(930093)

cv_splits <- rsample::vfold_cv(trainSetDiff, strata = PIK3CA_T)

mod <- decision_tree(mode = "classification",
                     cost_complexity = tune(),
                     tree_depth = tune(),
                     min_n = tune()) %>%
    set_engine("rpart")

rec <- recipe(PIK3CA_T ~ ., data = trainSetDiff)

wfl <- 
  workflow() %>%
  add_recipe(rec) %>%
  add_model(mod)

cart_set <- parameters(cost_complexity(range = c(-10, -1)),
                       tree_depth(range = c(1L, 5L)),
                       min_n(range = c(5L, 20L)))

cart_grid <- 
  grid_regular(cart_set, levels = c(5, 5, 10))


ctrl <- control_grid(save_pred = TRUE, verbose = TRUE)

cart_tune <- 
  tune_grid(wfl,
            resamples = cv_splits,
            grid = cart_grid,
            metrics = metric_set(roc_auc),
            control = ctrl)

show_best(cart_tune)

bestTune <- select_best(cart_tune)

cart_rec_final <- prep(rec)

cart_mod_final <- finalize_model(mod, bestTune)

cart_fit <- 
  cart_mod_final %>% 
  fit(PIK3CA_T ~ ., data = juice(cart_rec_final))

cart_fit$fit$call

cart_fit$fit$call <- rlang::call2("rpart", PIK3CA_T ~ .,
                                  data = expr(trainSetDiff),
                                  cp = ~1e-10, maxdepth = ~5L, 
                                  minsplit = ~15L)

as.party()
plot( as.party(cart_fit$fit) )
rpart.plot(cart_fit$fit)

save(cart_tune, file = "cart_tuen.RData")
