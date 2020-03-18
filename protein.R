library(readr)
library(tidyverse)

protein <- read_delim("D:/Studies/TheHumanProteinAtlas/pathology.tsv/pathology.tsv", delim = "\t")

protein %>% 
  filter(Cancer == "breast cancer" & `Gene name` == " WBP1L")
