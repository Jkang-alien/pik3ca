load('dataset.RData')
library(purrr)
library(tidyverse)
library(ggplot2)

variance <- trainSet %>%
  select(-PIK3CA_T) %>%
  map_dbl(var)

sum(is.na(variance))
summary(variance)
summary(variance[variance > median(variance)])

trainsetLargeVariance <- trainSet[,-1][,variance > median(variance)]
trainsetLargeVariance$PIK3CA_T <- trainSet$PIK3CA_T

testsetLargeVariance <- testSet[,-1][,variance > median(variance)]
testsetLargeVariance$PIK3CA_T <- testSet$PIK3CA_T

stasticsWilcox <- vector(mode = "list", length = dim(trainsetLargeVariance)[2])
pvalueWilcox <- vector(mode = "list", length = dim(trainsetLargeVariance)[2])

for (i in 1:(dim(trainsetLargeVariance)[2]-1)) {
  a <- wilcox.test(trainsetLargeVariance[,i] ~ trainsetLargeVariance$PIK3CA_T)
  stasticsWilcox[[i]] <- a$statistic
  pvalueWilcox[[i]] <- a$p.value
}

colnames(trainsetLargeVariance)[log(unlist(pvalueWilcox), base = 10) < -8]

trainSetDiff <- trainsetLargeVariance[,log(unlist(pvalueWilcox), base = 10) < -8]

trainSetDiff$PIK3CA_T<- trainSet$PIK3CA_T

testSetDiff <- testsetLargeVariance[,log(unlist(pvalueWilcox), base = 10) < -8]

testSetDiff$PIK3CA_T<- testSet$PIK3CA_T

heatmap(cor(trainSetDiff))


summary(trainSet$PIK3CA_T)
summary(testSet$PIK3CA_T)

save(trainSetDiff, testSetDiff,  file = "datasetDiff.RData")
