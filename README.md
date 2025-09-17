# ATAC-seq Pipeline in R

This repository contains an ATAC-seq analysis workflow implemented in R.  
It was written as part of the **NGS: ATAC-seq course from King’s Health Partners**.

## Features
- Quality control with `ATACseqQC` and `soGGi`
- Fragment size distribution and TSS enrichment plots
- Peak annotation and inspection from BED files
- Differential accessibility analysis using `DESeq2` and `apeglm`
- Visualization: PCA, Volcano plots, Heatmaps

## Notes
- The BAM files used in this course were chromosome Y–restricted subsets.  
- A count matrix was not available; DESeq2 steps are included as a template.  

## Requirements
- R (≥ 4.0)
- Bioconductor packages: `ATACseqQC`, `soGGi`, `DESeq2`, `apeglm`, `TxDb.Hsapiens.UCSC.hg19.knownGene`
- Other CRAN packages: `tidyverse`, `ggplot2`, `pheatmap`, `data.table`, `dplyr`

## Acknowledgements
Based on materials from the **King’s Health Partners NGS: ATAC-seq course**.

