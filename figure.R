load('result.RData')

library(Cairo)
library(tidyverse)
library(tidymodels)
library(gridExtra)

source("theme_publication.R")

### ROC PR Plots
ROCPlot <- bind_rows("Trainset" = train_probs,
          "Testset" = test_probs,
          .id = "dataset") %>%
  group_by(dataset) %>%
  roc_curve(obs, .pred_Mutant) %>%
    ggplot(aes(x=1-specificity, y=sensitivity, color = dataset)) +
  geom_line() +
  geom_text(aes(x=.6, y=.6, label = "AUROC = 0.93")) + 
  geom_text(aes(x=.6, y=.4, label = "AUROC = 0.84"), color = "#386cb0") +
  labs(title = "ROC curve") +
  scale_colour_Publication() +
  theme_Publication()
  


prPlot <-  bind_rows("Trainset" = train_probs,
                     "Testset" = test_probs,
                     .id = "dataset") %>%
  group_by(dataset) %>%
  pr_curve(obs, .pred_Mutant) %>%
  ggplot(aes(x=recall, y=precision, color = dataset)) +
  geom_line() +
  geom_hline(yintercept = 0.11, color = 3)+
  geom_text(aes(x=.7, y=1.0, label = "AUPR = 0.66")) + 
  geom_text(aes(x=.7, y=.8, label = "AUPR = 0.39"), color = "#386cb0") +
  labs(title = "PR curve") +
  scale_colour_Publication() +
  theme_Publication()

CairoPDF(file = 'Figure1',
         width = 6, height = 3)
grid.arrange(ROCPlot, prPlot, nrow = 1)
dev.off()

typePlot +
  labs(title = "Coefficiency of cancer type",
       x = "Coefficiency",
       y = "Cancer type") +
  scale_colour_Publication() +
  theme_Publication()

### Train_test ROC 

Train_test_ROC <- test_trainPlot +
  labs(title = "AUROC of trainset and testset",
       x = "Trainset",
       y = "Testset") +
  scale_colour_Publication() +
  theme_Publication()
  
bind_rows("Trainset" = train_probs,
          "Testset" = test_probs,
          .id = "dataset") %>%
  skimr::skim()
  

