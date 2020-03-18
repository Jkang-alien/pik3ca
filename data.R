## @knitr cbioportal

library(cgdsr)
library(tidyverse)
library(RTCGAToolbox)
library(data.table)
library(readr)
# Create CGDS object
mycgds = CGDS("https://www.cbioportal.org/")

test(mycgds)

# Get list of cancer studies at server
cancerStudy <- getCancerStudies(mycgds)

breastPanCancer <- cancerStudy[48,1]
# Get available case lists (collection of samples) for a given cancer study
caselist = getCaseLists(mycgds,breastPanCancer)[1,1]

# Get available genetic profiles
geneticprofile = getGeneticProfiles(mycgds,breastPanCancer)

mutation <- geneticprofile[5,1]

# Get data slices for a specified list of genes, genetic profile and case list
pik3ca <- getProfileData(mycgds,'PIK3CA', mutation, caselist)
pik3caDataTable <- data.table(pik3ca)
pik3caDataTable$ID <- gsub('\\.', '-',rownames(pik3ca))
pik3caDataTable$PIK3CA_T <- is.na(pik3caDataTable$PIK3CA)
pik3caDataTable$PIK3CA_T <- factor(pik3caDataTable$PIK3CA_T,
                                   levels = c(TRUE, FALSE),
                                   labels = c("Wild", "Mutant"))
pik3caDataTable <- pik3caDataTable[,-1]

# Get clinical data for the case list
clinicaldata = getClinicalData(mycgds,caselist)

clinicaldata <- clinicaldata %>%
  mutate(ID = rownames(clinicaldata))



breastFirehoseCancer <- cancerStudy[46,1]
# Get available case lists (collection of samples) for a given cancer study
caselistFirehose = getCaseLists(mycgds,breastFirehoseCancer)[1,1]

# Get clinical data for the case list

clinicaldataFirehose = getClinicalData(mycgds,caselistFirehose)
clinicaldataFirehose <-  clinicaldataFirehose %>%
  mutate(ID = rownames(clinicaldataFirehose))

clinicaldata <- inner_join(clinicaldata, clinicaldataFirehose, by = "ID")


clinicalIncluded <- clinicaldata %>%
  filter(ER_STATUS_BY_IHC == "Positive") %>%
  mutate(HER2_STATUS = (IHC_HER2 == "Negative" | HER2_FISH_STATUS == "Negative"))

## @knitr Firehose

# Firehose 
datasetBRCA <- getFirehoseDatasets()[3]

stddata <- getFirehoseRunningDates()

gisticDate <- getFirehoseAnalyzeDates(last=3)
gisticDate

rnaseqData <- getFirehoseData(dataset=datasetBRCA, runDate=stddata[1],
                            forceDownload=TRUE, clinical=FALSE, Mutation=FALSE,
                            RNASeq2GeneNorm = TRUE)

# save(rnaseqData, file = "rnaseq.RData" )

load("rnaseq.RData")

# pancancer 
rnaseqMatrix <- getData(rnaseqData, "RNASeq2GeneNorm")
row.names(rnaseqDataTable)[1:4]
colnames(rnaseqDataTable)[1:4]
rnaseqDataTable <- data.table(t(rnaseqMatrix))

colnames(rnaseqDataTable) <- row.names(rnaseqMatrix)
rnaseqDataTable$ID <- gsub("[A-Z]{1}-[0-9]{2}[A-Z]{1}-[0-9A-Z]{4}-[0-9]{2}", '',
                           colnames(rnaseqMatrix))

rnaseqDataTable$ID %in% gsub("\\.", "-", clinicalIncluded$ID)

rnaIncluded <- rnaseqDataTable[rnaseqDataTable$ID %in% gsub("\\.", "-", clinicalIncluded$ID),]

dataset <- inner_join(pik3caDataTable, rnaIncluded, by='ID')

dataset[1:5,1:5]

colnames(clinicalIncluded)

summary(factor(clinicalIncluded$HISTOLOGICAL_DIAGNOSIS))

histology <- clinicalIncluded %>%
  select(ID, HISTOLOGICAL_DIAGNOSIS) 

histology$ID <- gsub('\\.','-', histology$ID)

dataset <- inner_join(histology, dataset, by = "ID")

dataset <- dataset %>%
  filter(HISTOLOGICAL_DIAGNOSIS == "Infiltrating Ductal Carcinoma" |
           HISTOLOGICAL_DIAGNOSIS == "Infiltrating Lobular Carcinoma")

ductalSet <- dataset %>%
  filter(HISTOLOGICAL_DIAGNOSIS == "Infiltrating Ductal Carcinoma") %>%
  select(-ID, -HISTOLOGICAL_DIAGNOSIS)

lobularSet <- dataset %>%
  filter(HISTOLOGICAL_DIAGNOSIS == "Infiltrating Lobular Carcinoma") %>%
  select(-ID, -HISTOLOGICAL_DIAGNOSIS)

save(dataset, ductalSet, lobularSet, file = "datsetHistology.RData")
