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

## plotting ####################################################################

counts <- assay(airway) # extracts the number of genes matched to each gene

library_sizes <- colSums(counts) # sum of reads for each sample (sequencing depth)

barplot(
  library_sizes,
  las = 2,
  main = "Reads per sample",
  ylab = "Total counts"
)

library_sizes <- colSums(assay(airway)) # see if read numbers correlate to treatment groups

data.frame(
  Sample = colnames(airway),
  Treatment = colData(airway)$dex,
  Reads = library_sizes
) # no correlation

summary(library_sizes) # very large difference between the largest and smallest (requires normalisation)
library_sizes

## DESeq2 for normalisation ####################################################

dds <- DESeqDataSet(airway, design = ~ cell + dex) # Creates a DESeq2 object with the airway data set and specifying the study's design. 
dds <- estimateSizeFactors(dds) # how large the library is on average (>1 = above average or smaller library than normal). Use by DESeq2 to normalise counts.

sizeFactors(dds) # used for technical scaling corrections

# DESeq2 uses 'normalized count = raw count / size factor' internally with these numbers.

# Do samples now cluster by treatment?

## PCA plotting ################################################################

vsd <- vst(dds)
plotPCA(vsd, intgroup = "dex") # shows treated vs untreated samples. they separate clearly on PC1, so treatment is the strongest driver of variation.

## finding the genes that cause the separation #################################

dds <- DESeq(dds)
res <- results(dds)








