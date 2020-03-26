library(GEOquery)


## https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE14860

gse <- getGEO("GSE14860", GSEMatrix = TRUE)
show(gse)

save(gse, file = "gse.RData")
