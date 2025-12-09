# Load Libraries
library(DESeq2)
library(tidyverse)
library(openxlsx)
library(pheatmap)
library(ggplot2)

# Load read count data from csv file
counts_data <- read_csv("C:/Users/rambe/Downloads/combined.csv", locale = locale(encoding = "UTF-8"))  # Change path accordingly
counts_data <- as.data.frame(counts_data)
rownames(counts_data) <- counts_data[[1]]
counts_data <- counts_data[, -1]

# Fix sample names: replace - with _
colnames(counts_data) <- gsub("-", "_", colnames(counts_data))

# Define comparison groups manually (prefixes of sample names)
#  Add this in this format -> “ label “ = c(“Coloumn_1”,”Coloumn_2”)
comparisons <- list(
  "EtEC81_10d_vs_GFP_10d" = c("EtEC81_10d", "GFP_10d"),
  "EtEC81_14d_vs_GFP_14d" = c("EtEC81_14d", "GFP_14d"),
  "EtEC81_21d_vs_GFP_21d" = c("EtEC81_21d", "GFP_21d"),
  "ZmEIP_10d_vs_GFP_10d"  = c("ZmEIP_10d", "GFP_10d"),
  "ZmEIP_14d_vs_GFP_14d"  = c("ZmEIP_14d", "GFP_14d"),
  "ZmEIP_21d_vs_GFP_21d"  = c("ZmEIP_21d", "GFP_21d"),
  "EtEC81_10d_vs_ZmEIP_10d" = c("EtEC81_10d", "ZmEIP_10d"),
  "EtEC81_14d_vs_ZmEIP_14d" = c("EtEC81_14d", "ZmEIP_14d"),
  "EtEC81_21d_vs_ZmEIP_21d" = c("EtEC81_21d", "ZmEIP_21d")
)

# Output workbook
wb <- createWorkbook()

# Loop through each comparison
for (comp in names(comparisons)) {
  cat("Processing:", comp, "\n")
  
  group1 <- comparisons[[comp]][1]
  group2 <- comparisons[[comp]][2]
  
  # Get matching samples by prefix
  samples_g1 <- grep(paste0("^", group1), colnames(counts_data), value = TRUE)
  samples_g2 <- grep(paste0("^", group2), colnames(counts_data), value = TRUE)
  selected_samples <- c(samples_g1, samples_g2)
  
  # Subset counts
  sub_counts <- counts_data[, selected_samples]
  
  # Create metadata manually
  group_labels <- c(rep(group1, length(samples_g1)), rep(group2, length(samples_g2)))
  coldata <- data.frame(row.names = selected_samples, group = factor(group_labels))
  
  # Run DESeq2
  dds <- DESeqDataSetFromMatrix(countData = round(as.matrix(sub_counts)),
                                colData = coldata,
                                design = ~ group)
  dds$group <- relevel(dds$group, ref = group2)  # Set group2 as reference
  dds <- DESeq(dds)
  
  res <- results(dds)
  res <- res[order(res$padj), ]
  res_sig <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
  
  # Save to Excel (add full + sig as two sheets per comparison)
  addWorksheet(wb, paste0(comp, "_full"))
  writeData(wb, paste0(comp, "_full"), as.data.frame(res), rowNames = TRUE)
  
  addWorksheet(wb, paste0(comp, "_sig"))
  writeData(wb, paste0(comp, "_sig"), as.data.frame(res_sig), rowNames = TRUE)
  
  # Volcano plot
  res_df <- as.data.frame(res)
  res_df$geneid <- rownames(res_df)
  res_df$threshold <- as.factor(abs(res_df$log2FoldChange) > 1 & res_df$padj < 0.05)
  g <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = threshold)) +
    geom_point(alpha = 0.7, size = 1.5) +
    scale_color_manual(values = c("grey", "red")) +
    theme_minimal() +
    labs(title = paste("Volcano Plot:", comp), x = "Log2 Fold Change", y = "-log10(padj)")
  
  ggsave(paste0("C:/Users/rambe/Downloads/volcano_", comp, ".png"), g, width = 6, height = 5)
  
  # Heatmap of top 30 genes
  if (nrow(res_sig) >= 2) {
    top_genes <- head(rownames(res_sig), 30)
    vsd <- vst(dds, blind = TRUE)
    mat <- assay(vsd)[top_genes, ]
    pheatmap(mat, cluster_rows = TRUE, cluster_cols = TRUE,
             main = paste("Top 30 DEGs:", comp),
             filename = paste0("C:/Users/rambe/Downloads/heatmap_", comp, ".png"))
  }
}

# Save Excel file with all DEGs
saveWorkbook(wb, "C:/Users/rambe/Downloads/All_DEG_Comparisons_test.xlsx", overwrite = TRUE)   #Add the desired path and file name

cat("✅ All comparisons completed and saved.\n")
