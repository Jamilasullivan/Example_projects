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

BiocManager::install("org.Hs.eg.db") # getting gene annotation information

library(org.Hs.eg.db)
library(AnnotationDbi)

top20 <- head(rownames(res_sorted), 20) # separate the top 20 genes

gene_info <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = top20,
  keytype = "ENSEMBL",
  columns = c("SYMBOL", "GENENAME")
) # change the first 20 gene names

gene_info # names of the first 20 genes

## Pathways enrichment analysis with DEGs ######################################

# What biological pathways are overrepresented in my data?

# clusterProfiler, ReactomePA or fgsea (GSEA-style ranking)

# test GO biological processes, KEGG pathways, Reactome pathways

sig_genes <- rownames(res)[res$padj < 0.05] # separating all significant genes

converted <- bitr(sig_genes,
                  fromType="ENSEMBL",
                  toType="ENTREZID",
                  OrgDb=org.Hs.eg.db) # converting all the gene IDs

ego <- enrichGO(converted$ENTREZID,
                OrgDb=org.Hs.eg.db,
                keyType="ENTREZID",
                ont="BP") # running enrichment analysis on all genes

library(enrichplot)

dotplot(ego, showCategory = 10)
barplot(ego, showCategory = 10)
cnetplot(ego)

## GSEA analysis ###############################################################

library(DESeq2)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

# res = your DESeq2 results object
res_df <- as.data.frame(res)

# remove NA logFC values
res_df <- res_df %>%
  filter(!is.na(log2FoldChange))

# rank the gene list
gene_list <- res_df$log2FoldChange
names(gene_list) <- rownames(res_df)

# sort descending (required for GSEA)
gene_list <- sort(gene_list, decreasing = TRUE)

# clean ENSEMLE IDs 
names(gene_list) <- gsub("\\..*", "", names(gene_list))

# map genes 
mapped <- bitr(
  names(gene_list),
  fromType = "ENSEMBL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

# keep only mapped genes
gene_list <- gene_list[mapped$ENSEMBL]

# assign ENTREZ IDs as names
names(gene_list) <- mapped$ENTREZID

# remove duplicates 
gene_list <- gene_list[!duplicated(names(gene_list))]

# final safety checks
cat("Number of genes:", length(gene_list), "\n")
cat("Duplicates:", sum(duplicated(names(gene_list))), "\n")
cat("Any NA values:", sum(is.na(gene_list)), "\n")

# run GSEA 
gse <- gseGO(
  geneList     = gene_list,
  OrgDb        = org.Hs.eg.db,
  keyType      = "ENTREZID",
  ont          = "BP",
  minGSSize    = 10,
  maxGSSize    = 500,
  pvalueCutoff = 0.05,
  verbose      = FALSE
)

# inspect results 
df <- as.data.frame(gse)

cat("GSEA terms found:", nrow(df), "\n")

head(df)

# plot results 
# main overview plot
dotplot(gse, showCategory = 15)

# ridgeplot (nice for presentation)
ridgeplot(gse)

# single pathway enrichment curve 
if (nrow(df) > 0) {
  gseaplot2(gse, geneSetID = 1)
}
