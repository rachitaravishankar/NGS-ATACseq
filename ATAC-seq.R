#### Set Working Directory ####
setwd("./ATAC-seq")

#### QC ####
# Work w/ aligned reads directly

##### Install Packages ####
BiocManager::install("ATACseqQC")
BiocManager::install("soGGi")

##### Load Packages ####
library(ggplot2)
library(tidyverse)
library(BiocManager)
library(ATACseqQC)
library(soGGi) #visualization package for ATAC-seq


##### Data #####
# Two conditions - omni and standard 
# BAM - aligned seq
# BAM.bai - table of contents to navigate bam 

##### Set up Parse #####
OMNI_bam_path <- "./BAM_files/OMNI.bam"
OMNI_bai_path <- "./BAM_files/OMNI.bam.bai"

STD_bam_path <- "./BAM_files/STD.bam"
STD_bai_path <- "./BAM_files/STD.bam.bai"

###### Inspect OMNI BAM files #####
  # Sub-setting Data 
bamQC(OMNI_bam_path, index = OMNI_bai_path, outPath = NULL)
# not mapped reads - taken out chrm y for sub-setting and ignored all of the other chrms, no mapped reads other than chrm y

  # Visualize data - Fragment Size Distribution Plot
# Good quality ATC-seq data should display a characteristic pattern
fragSize_OMNI <- fragSizeDist(OMNI_bam_path, "OMNI", index = OMNI_bai_path)
# not perfect but typical banding pattern of nucleosome 
# Multi QC and R plot differences? - Here in coding using reads from chrm y not all so any plot or values from QC not rep of reads of entire dataset

###### Inspect STD BAM files ####
  # Sub-setting Data 
bamQC(STD_bam_path, index = STD_bai_path, outPath = NULL)

  # Visualize data - Fragment Size Distribution Plot
fragSize_STD <- fragSizeDist(STD_bam_path, "STD", index = STD_bai_path)

##### Load Annotations #####
  # Load Human Genome Object
BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

  # Extract gene locations
genesLocations <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)

  # Extract Transcription start sites TSS locations
tssLocations <- resize(genesLocations, fix="start", width = 1)

# BAM is only mapped to y chrm - but TSS across all chrm - subset TSS to include them in chrm Y

  # Subset TSS to include chrm Y 
tssLocations <- tssLocations[seqnames(tssLocations) %in% "chrY"]
# TSS are transcriptionally active = compared to other regions expect ATAC-seq mapped there = check for quality

##### TSS enrichment plot #####
OMNI_plot <- regionPlot(bamFile = OMNI_bam_path, testRanges = tssLocations, paired = TRUE)
plotRegion(OMNI_plot) + theme_bw()

STD_plot <- regionPlot(bamFile = STD_bam_path, testRanges = tssLocations, paired = TRUE)
plotRegion(STD_plot) + theme_bw()

#### Differential Gene Analysis ####
# Call reads to find where a large number of genes map 
# Peak calling algo - max2
# Can find individually or together 
# good quality data - specific regions - higher FRiP score 

##### Install Packages ####
BiocManager::install("DESeq2")
BiocManager::install("apeglm")
install.packages("data.table")

##### Load Packages #####
library(data.table)
library(dplyr)
library(DESeq2)
library(ggplot2)
library(apeglm)
library(tidyverse)

##### Output Directory ####
output_path <- './output_DESeq2/'
dir.create(output_path, recursive = TRUE)

##### Peak Files #####
# bed files - annotated and non-annotated
# each peak = 'Interval'

peaks_bed_path <- "./Peak_files/consensus_peaks_annotated.bed"
peaks_bed <- as.data.frame(read.delim2("~/Peak_files/consensus_peaks_annotated.bed"))

# manipulate df to find how many peaks in data
# each row = peak = 104659
nrow(peaks_bed)

# identify smallest and biggest peak sizes - (bp long)
# use start and end co-ords

df <- peaks_bed %>%
  mutate(peak_size = peaks_bed$End - peaks_bed$Start + 1) %>%
  arrange(peak_size)
head(df)
tail(df)

# Annotated peak file
# The data provided already contains annotations - peak_bed
# Nearest Promoter ID = when looking for enhancers - use ATAC-seq to identify which regions are open and find if these regions are enhancers that are controlling expression - closest gene to the peak - to make predictions 

# Count Matrix - accessibility differences b/w sample A and sample B
# do not have the count matrix data set

# OMNI vs STD w/ DESeq2
# pick peak IDs that are most differnet from each condition
# everything below this # is in case of having the matrix file 

##### Counts Matrix ####
peak_count_matrix

counts <- peak_count_matrix[,7:12]
dim(counts)


  # clean the data 
# remove endings of colnames
colnames(counts) <- gsub("\\..*", "",  colnames(counts))

# add peak names to the first column of the count df
counts <- counts %>% mutate(peak_ID = peak_count_matrix$Geneid)

# reorder columns
counts <- counts %>% select(peak_ID, GM12878_STD_REP1, GM12878_STD_REP2, GM12878_FAST_REP1,
                            GM12878_FAST_REP2, GM12878_OMNI_REP1, GM12878_OMNI_REP2)

# make peak IDs the rownames
counts2 <- counts[,-1]
rownames(counts2) <- counts[,1]
counts <- counts2

# check cleaned df
print(head(counts))

  # initialize DESeq2 object 
# set up metadata table
col_data <- data.frame(colnames(counts))
colnames(col_data) <- "sample_name"
col_data <- col_data %>% 
  mutate(condition = sapply(strsplit(as.character(col_data$sample_name), "_"), `[`, 2)) %>%
  mutate(replicate = sapply(strsplit(as.character(col_data$sample_name), "_"), `[`, 3))

print(col_data)

##### Run DESeq2 ####

deseq <- DESeqDataSetFromMatrix(counts, design = ~ condition, colData = col_data)

  # design 
deseq$condition <- relevel(deseq$condition, ref = "STD")

  # run
deseq <- DESeq(deseq)

  # normalize the data
rld <- rlog(deseq)

  # run PCA
plotPCA(rld, intgroup = "condition")

  # differential accessibility 
res <- results(deseq)
res

resultsNames(deseq)

res <- results(deseq, name=c("condition_FAST_vs_STD"))
res

  # set reference level 
deseq$condition <- relevel(deseq$condition, ref = "OMNI")
deseq <- DESeq(deseq)

resultsNames(deseq)

res <- results(deseq, name=c("condition_FAST_vs_OMNI"))
res

  # calculate more accurate log2FC
res <- lfcShrink(deseq, coef=c("condition_FAST_vs_OMNI"), type = 'apeglm')
res


  # identify peak with the highest log2FC
deseq$condition <- relevel(deseq$condition, ref = "FAST")
deseq <- DESeq(deseq)

resultsNames(deseq)

res <- results(deseq, name=c("condition_STD_vs_FAST"))
head(res[order(res$log2FoldChange, decreasing = TRUE),])

  # summary table with adjusted p-values 
summary(res, alpha = 0.05)

  # save the results to output 
saveRDS(deseq, paste0(output_path, 'deseq_output.RDS'))

#### Visualizations ####
# need to read in the data saved in the output directory 

##### Load Libraries ####
library(pheatmap)
library(ggplot2)
library(DESeq2)
library(dplyr)

  # set output path if needed 

##### Read DESeq2 output #####
deseq <- readRDS(paste0(output_path, 'deseq_output.RDS'))
res <- lfcShrink(deseq, coef=c("condition_STD_vs_FAST"), type = 'apeglm')

##### Volcano Plot ####
  #normalize the data 
rld <- rlog(deseq)

  # Format data for visualization
volcano_data <- as.data.frame(res)
volcano_data <- mutate(volcano_data, '-log10(padj)' = -log10(padj))

ggplot(volcano_data, aes(log2FoldChange, `-log10(padj)`)) +
  geom_point() +
  theme_classic()

# small p value = more significance = dots on the higher points of the plot have smaller p values thus more significant as -log10 adjusted.

  # can set thresholds on the data 
volcano_data <- volcano_data %>%
  filter(!is.na(padj)) %>%
  mutate(sig = case_when((padj < 0.001) == 'TRUE' ~ 'significant',
                         (padj >= 0.001) == 'TRUE' ~ 'not sigificant')) %>%
  arrange(abs(padj))

  # color peaks according to significance
ggplot(volcano_data, aes(log2FoldChange, `-log10(padj)`)) +
  geom_point(aes(colour = sig, fill = sig), size = 0.5) +
  scale_color_manual(breaks = c("not sigificant", "significant"),
                     values = alpha(c("#c1c1c1", "#48d1cc"))) +
  geom_hline(yintercept=-log10(0.001), linetype="dashed", color = "black") +
  theme_classic() +
  theme(legend.position = "none")


  #  p < 0.001 , log2FC > 1.5 = blue = unregulated peaks 
  #  p < 0.001 , log2FC < -1.5 = red = down regulated peaks 
  #  p > 0.001 , log3FC > abs(1.5) = grey = non significant peaks
volcano_data <- volcano_data %>%
  filter(!is.na(padj)) %>%
  mutate(sig = case_when((padj < 0.001 & log2FoldChange > 1.5) == 'TRUE' ~ 'upregulated',
                         (padj < 0.001 & log2FoldChange < -1.5) == 'TRUE' ~ 'downregulated',
                         (padj >= 0.001 | abs(log2FoldChange) <= 1.5) == 'TRUE' ~ 'not sig')) %>%
  arrange(abs(padj))

ggplot(volcano_data, aes(log2FoldChange, `-log10(padj)`)) +
  geom_point(aes(colour = sig, fill = sig), size = 0.5) +
  scale_color_manual(breaks = c("not sig", "downregulated", "upregulated"),
                     values = alpha(c("#c1c1c1", "#f55f20", "#48d1cc"))) +
  geom_hline(yintercept=-log10(0.001), linetype="dashed", color = "black") +
  geom_vline(xintercept=c(-1.5, 1.5), linetype="dashed", color = "black") +
  theme_classic() +
  theme(legend.position = "none")


##### Heat map of top differentially accessible peaks ####
res_sub <- filter(as.data.frame(res), padj < 0.001 & abs(log2FoldChange) > 2.5)
res_sub <- res_sub[order(-res_sub$log2FoldChange),]
res_sub[1:10,]

pheatmap(assay(rld)[rownames(res_sub),], color = colorRampPalette(c("#191d73", "white", "#ed7901"))(n = 100), cluster_rows=T, show_rownames=FALSE,
         show_colnames = F, cluster_cols=T, annotation_col=as.data.frame(colData(deseq)["condition"]), scale = "row", treeheight_row = 0, treeheight_col = 25,
         border_color = NA)









