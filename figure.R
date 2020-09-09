load('result.RData')
load('dataset.RData')
library(Cairo)
library(tidyverse)
library(tidymodels)
library(gridExtra)

source("theme_publication.R")
### PIK3CA prevalence across cancer type

CairoPDF(file = 'Figure1',
         width = 4, height = 6, pointsize=8)
dataset %>%
  select(variant, type) %>%
  group_by(type) %>%
  count(variant) %>%
  spread(variant, n) %>%
  replace_na(list(Mutant=0)) %>%
  mutate(prevalence = Mutant/(Mutant + Wild)) %>%
  ggplot(aes(x=reorder(type, -prevalence), y=prevalence)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_Publication() + 
  labs(title = "PIK3CA prevalence",
       y = "Prevalence rate",
       x = "Cancer type") 
dev.off()
## Prevalence of PIK3CA 


dfprevalence <- dataset %>%
  select(variant, type) %>%
  group_by(type) %>%
  count(variant) %>%
  spread(variant, n) %>%
  replace_na(list(Mutant=0)) %>%
  mutate(prevalence = Mutant/(Mutant + Wild))

summary(dfprevalence$prevalence)
summary(dataset$variant)
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
  labs(title = "Coefficient of genes",
       y = "Coefficient",
       x = "Gene") +
  scale_colour_Publication() +
  theme_Publication()

typeCoeff <- typePlot +
  labs(title = "Coefficient of cancer type",
       y = "Coefficient",
       x = "Cancer type") +
  scale_colour_Publication() +
  theme_Publication()


CairoPDF(file = 'Figure2',
         width = 8, height = 12, pointsize=10)
grid.arrange(ROCPlot, prPlot,
             cor_ROC, cor_PR,
             mutationROC, mutationPR,
             nrow = 3)
dev.off()

CairoPDF(file = 'Figure3',
         width = 8, height = 12)
grid.arrange(geneCoeff, typeCoeff,
             nrow = 1)
dev.off()



typePlot +
  labs(title = "Coefficient of cancer type",
       x = "Coefficient",
       y = "Cancer type") +
  scale_colour_Publication() +
  theme_Publication()
library(tidyverse)
library(tidymodels)
library(latex2exp)
library(Cairo)

load(here::here("glmn_tune_pancancer.RData"))

roc_vals <- 
  collect_metrics(glmn_tune) %>% 
  filter(.metric == "roc_auc")

CairoPDF(file = 'S1Figure.pdf',
         width = 5, height = 5, pointsize=10)
roc_vals %>% 
  mutate(mixture = format(mixture)) %>% 
  ggplot(aes(x = penalty, y = mean, col = mixture)) + 
  geom_line() + 
  geom_point() + 
  geom_errorbar(aes(ymin = mean - std_err, ymax = mean + std_err, width = 0.05)) +
  scale_x_log10() +
  xlab(TeX("Penalty ($\\lambda$)")) +
  ylab("Area under the receiver operating characteristic (AUROC)") +
  labs(color = TeX("Mixture ($\\alpha$)"))
dev.off()
