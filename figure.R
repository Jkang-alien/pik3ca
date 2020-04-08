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

cor_ROC <- roc_test_trainPlot +
  labs(title = "AUROC",
       x = "Trainset",
       y = "Testset") +
  scale_colour_Publication() +
  theme_Publication()

cor_PR <- pr_test_trainPlot +
  labs(title = "AUPR",
       x = "Trainset",
       y = "Testset") +
  scale_colour_Publication() +
  theme_Publication()

mutationROC <- mutatioPlotTest +
  labs(title = "AUROC",
       x = "Mutation rate",
       y = "ROC") +
  scale_colour_Publication() +
  theme_Publication()

mutationPR <- mutatioPlotPRTest +
  labs(title = "AUPR",
       x = "Mutation rate",
       y = "PR") +
  scale_x_continuous(limits = c(0, 0.6)) +
  scale_y_continuous(limits = c(0, 0.6)) +
  scale_colour_Publication() +
  geom_abline(slope = 1, intercept = 0, color = 3) +
  theme_Publication()

geneCoeff <- coeffPlot +
  labs(title = "Coefficiency of genes",
       x = "Coefficiency",
       y = "Gene") +
  scale_colour_Publication() +
  theme_Publication()

typeCoeff <- typePlot +
  labs(title = "Coefficiency of cancer type",
       x = "Coefficiency",
       y = "Cancer type") +
  scale_colour_Publication() +
  theme_Publication()


CairoPDF(file = 'Figure1',
         width = 8, height = 12, pointsize=10)
grid.arrange(ROCPlot, prPlot,
             cor_ROC, cor_PR,
             mutationROC, mutationPR,
             nrow = 3)
dev.off()

CairoPDF(file = 'Figure2',
         width = 8, height = 12)
grid.arrange(geneCoeff, typeCoeff,
             nrow = 1)
dev.off()



typePlot +
  labs(title = "Coefficiency of cancer type",
       x = "Coefficiency",
       y = "Cancer type") +
  scale_colour_Publication() +
  theme_Publication()