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

ID <- gsub("[A-Z]{1}-[0-9]{2}[A-Z]{1}-[0-9A-Z]{4}-[0-9]{2}", '',
           data_rna$ID)
data_rna_unique <- data_rna[duplicated(ID) == FALSE,]

duplicated(ID)
