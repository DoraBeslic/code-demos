# runDESeq2.R
# script to perform differential gene expression analysis using DESeq2 package

# Goal: Visualize transcriptional changes in human airway smooth muscle cells 
# after treatment with dexamathesone

library(DESeq2)
library(tidyverse)
library(airway)
library(ggplot2)
library(ggrepel)
library(org.Hs.eg.db)

# read in count data and sample info
counts_data <- read.csv('~/Desktop/code_demos/DESeq2_demo/Data/counts_data.csv')
colData <- read.csv('~/Desktop/code_demos/DESeq2_demo/Data/sample_info.csv')

# make sure the row names in colData match column names in counts_data
all(colnames(counts_data) %in% rownames(colData))

# are they in the same order?
all(colnames(counts_data) == rownames(colData))

# construct a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = counts_data,
                              colData = colData, 
                              design = ~ dexamethasone)

# remove rows with low gene counts (less than 10 reads)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# set the factor level
dds$dexamethasone <- relevel(dds$dexamethasone, ref = 'untreated')

# run DESeq 
dds <- DESeq(dds)
res <- results(dds)

# explore results
summary(res)
res0.01 <- results(dds, alpha = 0.01)
summary(res0.01)

# convert res to dataframe for plotting
res_df = as.data.frame(res)

# color-code data points for volcano plot
res_df$color <- ifelse(res_df$padj < 0.05 & res_df$log2FoldChange > 0, "red", 
                       ifelse(res_df$padj < 0.05 & res_df$log2FoldChange < 0, 
                              "blue", "black"))

# map ensembl gene ids to gene names
# keytypes(org.Hs.eg.db)
# columns(org.Hs.eg.db)
gene_names_org <- mapIds(org.Hs.eg.db, 
                         keys = rownames(res_df),
                         keytype = "ENSEMBL",
                         column = "SYMBOL")

res_df$gene_name <- gene_names_org

# extract top 10 differentially expressed genes to label on volcano plot
top_genes <- head(res_df[order(-abs(res_df$log2FoldChange)), ], 10)

# create volcano plot
res_df %>%
  ggplot(., aes(x = log2FoldChange, y = -log10(padj), color = color)) +
  geom_point(alpha = 0.7) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray") +
  scale_color_identity() +
  geom_text_repel(
    data = top_genes,
    aes(label = gene_name, color = "black"),
    size = 2,
  ) +
  labs(title = "DGE of human airway smooth muscle cells 
       treated with 1 uM dex for 18 hours",
       x = "Log2 Fold Change",
       y = "-log10(Adjusted P-value)")


