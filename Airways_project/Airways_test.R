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

## finding the genes that cause the separation after dex treatment #############

dds <- DESeq(dds) # differential expression analysis. Fits a statistical model for every gene
res <- results(dds) # extracting the results

res
summary(res)
head(res)

sig <- res[which(res$padj < 0.05), ] # filter significant genes

nrow(sig) # number of significant genes

res_sorted <- res[order(res$padj), ] # sort genes by significance
head(res_sorted) # look at most significant genes. The genes that responded most strongly to treatment. baseMean is the average normalised expression across all samples.

write.csv(
  as.data.frame(res_sorted),
  "Airways_project/deseq2_results.csv"
) # created an output file to be shared with collaborators

## Volcano plots ###############################################################

plot(
  res$log2FoldChange,
  -log10(res$padj),
  pch = 20,
  xlab = "Log2 Fold Change",
  ylab = "-log10 Adjusted P-value",
  main = "Differential Expression"
) # without colour

sig <- !is.na(res$padj) & res$padj < 0.05

plot(
  res$log2FoldChange,
  -log10(res$padj),
  col = ifelse(sig, "red", "black"),
  pch = 20,
  xlab = "Log2 Fold Change",
  ylab = "-log10 Adjusted P-value"
) # with colour. The treatment is upregulating more genes than it is repressing, most genes a not largely affected by treatment 

rownames(res_sorted)[1] # number one most significant gene

plotCounts(
  dds,
  gene = rownames(res_sorted)[1],
  intgroup = "dex"
) # normalised counts of treated vs untreated for a single gene

head(rowData(dds)) # gives information about the top few genes

## Which biological processes are changing? ####################################

head(rownames(res_sorted), 20) # identify the top 20 genes

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("org.Hs.eg.db")































