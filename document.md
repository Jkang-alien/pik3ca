PIK3CA mutation prediction
================

### Purpose

PIK3CA mutation prediction by gene expression data

### Introduction

Breast cancer with PIK3CA mutation has been approved to use PIK3CA
inhibitor in hormone receptor positive HER2 negative subtype. Prediction
of PIK3CA mutation was done by gene expression data of TCGA.

### Dataset

PIK3CA mutation data was get using cgdsr rpackage. Gene expression data
was get from GDAC firehose using RTCGAToolbox R package. ER immunostain
postive HER2 immunostain negative and/or SISH negative breast cancers
are included. Data of invasive ductal carcinomas were used for training
set and data of invasive lobular carcinoma were used for test set.
Number of observations were 530 in training set and 188 in test set.

Table 1. Frequency of PIK3CA mutation in training set and test set

| TrainSet | TestSet |
| -------: | ------: |
|      331 |     117 |
|      199 |      71 |

### Selecting variable for modeling

To narrow down potential predictors, Genes with a large variance (more
than median) and differently expressed by PIK3CA mutation status (wilcox
test, *P*
![\<10^{-8}](https://latex.codecogs.com/png.latex?%3C10%5E%7B-8%7D
"\<10^{-8}")) were selected. 111 out of 20502 genes were included in the
modeling process.

### Preprocessing

BoxCox transformation was done to correct skewness. Centering and
scaling were done. All preprocessing was done using recipe r package.

### Modeling

Penalized logistic regression was applied to prediction modeling.
10-fold cross-validation with targe variable stratification was done
over the hyperparameter grid:
![lambda](https://latex.codecogs.com/png.latex?lambda "lambda")
{![10^{-5}](https://latex.codecogs.com/png.latex?10%5E%7B-5%7D
"10^{-5}"),
![10^{-4}](https://latex.codecogs.com/png.latex?10%5E%7B-4%7D
"10^{-4}"),![10^{-3}](https://latex.codecogs.com/png.latex?10%5E%7B-3%7D
"10^{-3}"),![10^{-2}](https://latex.codecogs.com/png.latex?10%5E%7B-2%7D
"10^{-2}"),![10^{-1}](https://latex.codecogs.com/png.latex?10%5E%7B-1%7D
"10^{-1}"), ![10^{0}](https://latex.codecogs.com/png.latex?10%5E%7B0%7D
"10^{0}")}, ![alpha](https://latex.codecogs.com/png.latex?alpha "alpha")
{0.0, 0.25, 0.5, 0.75}. Model performance was evealuated with the
receiver operating characteristic (AUROC).

The model showed best performance at lambda = 0.01 and alpha = 1.0
(Ridge regression).

|        | Wild | Mutant |
| ------ | ---: | -----: |
| Wild   |  289 |     56 |
| Mutant |   42 |    143 |
