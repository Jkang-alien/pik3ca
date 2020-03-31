library(GEOquery)
library(tidyverse)
library(tidymodels)

## https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE14860

gse <- getGEO("GSE60788", GSEMatrix = TRUE)
show(gse)

save(gse, file = "gse.RData")
load("gse.RData")

summary(gse)

dim(pData(gse[[1]]))
pData(gse[[1]])

# df1 <- getGSEDataTables("GSE60788")

gunzip("GSE60788_UCSC_Human_hg19_knownGenes_20120910_addInfo.txt.gz")

gds <- getGEO(filename=system.file("extdata/GSE22035.soft.gz",package="GEOquery"))

gseRNA <- read_delim("GSE60788_rnaseq_gex_normalized.txt", delim = '\t')

testsetGSE <- data.table::data.table(t(gseRNA[,-1]))
colnames(testsetGSE) <- gseRNA$`Gene Symbol`

## @knitr testsetGEO

glmn_tune_tidy <- glmn_tune %>%
  select(-splits)

show_best(glmn_tune)

lasso_glmn <- show_best(glmn_tune)[4,1:2]
best_glmn <- select_best(glmn_tune)
mix <- show_best(glmn_tune)[3,1:2]

best_glmn$mixture <- 1
best_glmn$penalty <- .1
## @knitr Finalize_model

wfl_final <- 
  wfl %>%
  finalize_workflow(best_glmn) %>%
  fit(data = trainset)

predict(wfl_final, type = "prob", new_data = testsetGSE)

sum(colnames(trainset) %in% colnames(testsetGSE))

test_probs <- 
  predict(wfl_final, type = "prob", new_data = testsetGSE) %>%
  bind_cols(obs = testset$variant) %>%
  bind_cols(predict(wfl_final, new_data = testset)) %>%
  bind_cols(type = testset$type)

test_probs %>%
  filter(type == "BLCA") %>%
  conf_mat(obs, .pred_class)

autoplot(roc_curve(test_probs, obs, .pred_TRUE))

roc_auc(test_probs, obs, .pred_TRUE)