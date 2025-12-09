# DIFFERENTIAL-GENE-EXPRESSION-ANALYSIS
Differential gene expression is a tool used in bioinformatics to identify the genes which are expressed differently in different conditions.
This pipeline:

✔ Loads a raw count matrix
✔ Automatically identifies sample replicates using prefix matching
✔ Runs DESeq2 for each comparison
✔ Extracts both full and significant DEGs
✔ Generates:
    • Volcano plots
    • Heatmaps (top 30 DEGs)
✔ Saves all results into a single Excel workbook

This script is highly useful for experiments with multiple time points, treatments, or recombinant constructs (e.g., EtEC81 vs GFP, ZmEIP vs GFP, etc.).

├── data/
│   └── combined.csv
│
├── src/
│   └── deseq2_pipeline.R
│
├── results/
│   ├── volcano_<comparison>.png
│   ├── heatmap_<comparison>.png
│   └── All_DEG_Comparisons.xlsx
│
└── README.md

Input Requirements
1. Count Matrix (combined.csv)

*First column = gene IDs

*Other columns = raw integer counts

*Column names must begin with group prefixes 

2) Install Packages
   install.packages(c("tidyverse", "ggplot2", "openxlsx", "pheatmap"))

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")

3) pipeline
   *Place your count file in the data/ directory
   * Add ( change ) the comparison groups
   * Run the script
     
4) Output
   * All_DEG_Comparisons.xlsx containing:

      Comparison_full: all genes with DESeq2 stats
      Comparison_sig: significant DEGs (padj < 0.05 & |log2FC| > 1)
  * Volcano Plot
      Showing
          x-axis: log2 fold change
          y-axis: –log10 adjusted p-value
          red = significant DEGs
    * Heatmap - variance-stabilized counts (VST) for top 30 DEGs.
      
5) 🔍 Pipeline Logic Summary

    Load count data
    Clean column names (replace - with _)
    Set first column as rownames
    For each comparison:
    Detect samples based on prefix
    Subset count matrix
    Build metadata (colData)
    Create DESeq2 dataset
    Run differential expression
    Filter significant DEGs
    Save results + plots
    Compile all outputs into a single Excel file
    Print completion message
