library(tidyverse)
library(stringr)

rna <- read_delim("D:/rproject/XCIpancancer/EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2-v2.geneExp.tsv",
                  delim = '\t')

gene <- stringr::str_split_fixed(rna$gene_id, "\\|", n =2)[,1]
sample_ID <- colnames(rna)[-1]

data_rna <- data.table::data.table(t(rna[,-1]))
colnames(data_rna) <- gene
data_rna$ID <- sample_ID
data_rna[1:5,1:5]
data_rna$ID
data_rna$ID <- gsub("[A-Z]{1}-[0-9]{2}[A-Z]{1}-[0-9A-Z]{4}-[0-9]{2}", '', data_rna$ID)

data_rna_unique <- data_rna[!duplicated(data_rna$ID), c(-1:-29)]
data_rna_unique <- data_rna_unique[, c(-16273, -16274)]
data_rna_unique$ID  <- gsub("-[0-9]{2}$", "", data_rna_unique$ID)
sum(colnames(data_rna_unique) == `SLC35E2`)

data_rna_unique <- data_rna_unique %>%
  select_(as.name("SLC35E2"))

grep(as.name("SLC35E2"), colnames(data_rna_unique))

save(data_rna_unique, file = "panrna.RData")
## Mutation

mut <- read_delim("mutations.txt", delim = '\t')
summary(factor(mut$PIK3CA))

## https://www.accessdata.fda.gov/scripts/cdrh/cfdocs/cfpma/pma.cfm?id=P190001

variants <- c('C420R', 'E542K', 'E545A', 'E545D', 'E545G', 'E545K', 'Q546E', 'Q546R', 'H1047L', 'H1047R', 'H1047Y')

mut %>%
  filter(PIK3CA == variants)
mut$PIK3CA %in% variants

mut$PIK3CA

variant_dual <- stringr::str_split_fixed(mut$PIK3CA, " ", n =2)
variant_first <- mut$PIK3CA %in% variant_dual[,1]
variant_second <- mut$PIK3CA %in% variant_dual[,2]

mut$variant <- variant_first | variant_second
mut <- mut %>%
  mutate(ID = gsub("-[0-9]{2}$", "", SAMPLE_ID)) %>%
  select(ID, variant)

data <- inner_join(mut, data_rna_unique, by = "ID")
save(data, file = "panrna.RData")

load('panrna.RData')

library(purrr)
library(tidyverse)
library(ggplot2)

valueMad <- data %>%
  select(-variant, -ID) %>%
  map_dbl(mad, na.rm = TRUE)

summary(valueMad)

dataLV <- data[,-1:-2][,valueMad > median(valueMad)]

dataLV$variant <- data$variant
dataLV$ID <- data$ID


clinical <- read_delim("clinical.txt", delim = '\t')

## @knitr modeling

library(tidymodels)
library(glmnet)
library(skimr)



skim(clinical)
summary(factor(clinical$type))
head(clinical$bcr_patient_barcode)

hist <- clinical %>%
  select(bcr_patient_barcode, type) %>%
  rename(ID = "bcr_patient_barcode")

dataset <- inner_join(hist, dataLV, by = "ID" )

set.seed(930093)

initSplit <- initial_split(dataset, strata = variant)
trainset <- training(initSplit)
testset <- testing(initSplit)

cv_splits <- rsample::vfold_cv(trainset, strata = variant , repeats = 5)

mod <- logistic_reg(penalty = tune(),
                    mixture = tune()) %>%
  set_engine("glmnet")

rec <- recipe(variant ~ ., data = trainset) %>%
  update_role(ID, new_role = "id variable") %>%
  step_BoxCox(all_numeric()) %>%
  step_center(all_numeric()) %>%
  step_scale(all_numeric()) %>%
  step_dummy(type)%>%
  step_downsample(variant)

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

