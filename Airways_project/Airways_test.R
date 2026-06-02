## Airways test ################################################################

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("airway")
BiocManager::install("DESeq2")

library(airway)
library(DESeq2)

data("airway")

## basic info ##################################################################

colData(airway) # sample information

dim(assay(airway)) # dimentions of the count matrix

head(assay(airway)) # look at the first few genes and samples

ncol(airway) # how many samples are there?

nrow(airway) # how many genes are there?

table(colData(airway)$dex) # look at the dex column from the metadata (treatment status), then count how many times this appears

head(rownames(airway)) # what does each row represent?

colnames(airway) # what does each column represent?


