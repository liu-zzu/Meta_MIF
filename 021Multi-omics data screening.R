# Required libraries
library(ggplot2)
library(reshape2)
library(ggsci)
library(dplyr)

# Target genes list
target_genes <- c("ADIPOQ", "AFP", "ALB", "ALPL", "B2M", "BIRC5", "CALCA", "CCL2", 
                  "CD274", "CD34", "CD44", "CDH1", "CEACAM5", "CGA", "CGB", "CHGA", 
                  "CHI3L1", "COL18A1", "CP", "CRP", "CSF2", "CSF3", "CXCL10", "CXCL8", 
                  "CXCR4", "EGF", "EGFR", "ENG", "ENO2", "ENOX2", "EPCAM", "F3", 
                  "FGA", "FGB", "FGG", "FGF2", "FTH1", "FTL", "FUT3", "GDF15", 
                  "HAVCR2", "HBA1", "HBA2", "HBB", "HGF", "HIF1A", "HMGB1", "HP", 
                  "HSP90AA1", "ICAM1", "IFNG", "IGF1", "IGF2", "IGFBP2", "IGFBP3", 
                  "IGHA1", "IGHA2", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHM", 
                  "IL10", "IL12A", "IL12B", "IL17A", "IL18", "IL1B", "IL2", "IL2RA", 
                  "IL4", "IL6", "IL7", "KRT18", "KRT19", "KRT8", "LCN2", "LDHA", 
                  "LDHB", "LEP", "LGALS3", "MIF", "MIR21", "MKI67", "MMP2", "MMP7", 
                  "MMP9", "MUC1", "MUC16", "PDCD1", "PLAT", "PTHLH", "S100A9", 
                  "SAA1", "SAA2", "SERPINA1", "SERPINE1", "SOD1", "SOD2", "SOD3", 
                  "SPP1", "TERT", "TF", "TGFB1", "TGFB2", "TGFB3", "TIMP1", "TIMP2", 
                  "TNF", "TP53", "TTR", "VEGFA", "VEGFC", "VIM")

# Data processing
df <- read.table("final_cate_all.txt", header=TRUE, sep="\t", check.names=FALSE)

df_filtered <- df[df$gene %in% target_genes, ] %>%
  arrange(Tumor)

df_melted <- melt(df_filtered, id.vars = "gene")

# Create visualization
p <- ggplot(df_melted, aes(
  x = factor(gene, levels = df_filtered$gene),
  y = ifelse(variable == "Tumor", value, -value),
  fill = variable
)) +
  geom_bar(stat = 'identity') +
  geom_text(
    aes(
      label = abs(value),
      hjust = ifelse(variable == "Tumor", -0.4, 1.1)
    ),
    size = 3
  ) +
  scale_y_continuous(
    labels = abs,
    expand = expansion(mult = c(0.1, 0.1))
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    plot.title = element_text(size = 18, face = "bold")
  ) +
  labs(
    title = "TCGA Differential Expression Analysis",
    x = NULL, 
    y = "Number of Cancer Types"
  ) +
  scale_fill_manual(values = c(Tumor = "#D73027", Normal = "#4575B4"))

# Save output
ggsave("differential_expression.pdf", p, width = 14, height = 4)





# Required libraries
library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)
library(stringr)

# Data processing and median calculations
data <- fread("filtered_tumor_gene_data.txt")

gene_medians <- data %>%
  pivot_longer(cols = -c(sample, tissue, type, type2), 
               names_to = "Gene name", 
               values_to = "expression") %>%
  group_by(tissue, `Gene name`) %>%
  summarise(median_expression = median(expression, na.rm = TRUE)) %>%
  ungroup()

# Calculate overall medians for sorting
gene_overall_medians <- gene_medians %>%
  group_by(`Gene name`) %>%
  summarise(overall_median = median(median_expression, na.rm = TRUE)) %>%
  arrange(overall_median)

# Handle missing data
missing_data_genes <- gene_overall_medians$`Gene name`[is.na(gene_overall_medians$overall_median)]
gene_order <- c(missing_data_genes, 
                gene_overall_medians$`Gene name`[!is.na(gene_overall_medians$overall_median)])

gene_medians$`Gene name` <- factor(gene_medians$`Gene name`, levels = gene_order)

# Color gradient definition
blue_gradient_colors <- colorRampPalette(
  c("#BFD3E6", "#9ECAE1", "#6BAED6", "#4292C6", "#2171B5", "#08519C")
)(length(gene_order))

# Create visualization
p <- ggplot(gene_medians, aes(x = `Gene name`, y = median_expression, fill = `Gene name`)) +
  geom_boxplot() +
  scale_fill_manual(values = blue_gradient_colors) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, 
                               size = 10, face = "bold"),
    axis.title.x = element_blank(),
    legend.position = "none",
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 14, face = "bold")
  ) +
  labs(
    y = "Gene Expression (TPM)", 
    title = "Gene Expression Distribution Across TCGA Cancer Cohorts",
    subtitle = " "
  ) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 15))

# Add missing data annotations if needed
if (length(missing_data_genes) > 0) {
  p <- p + geom_text(
    data = data.frame(
      `Gene name` = missing_data_genes,
      median_expression = rep(max(gene_medians$median_expression, na.rm = TRUE), 
                              length(missing_data_genes))
    ),
    aes(x = `Gene name`, y = median_expression, label = "Missing Data"),
    vjust = 0.5, 
    angle = 90, 
    size = 3, 
    color = "red"
  )
}

# Save output
ggsave("gene_expression_distribution.pdf", 
       p, 
       width = 14, 
       height = 4, 
       limitsize = FALSE, 
       device = cairo_pdf, 
       dpi = 300)






# Required libraries
library(readxl)
library(ggplot2)
library(dplyr)
library(viridis)
library(tidyr)

# Data processing
data <- read_excel("pathology_with_NES.xlsx")

# Identify genes with missing data
missing_data_genes <- data %>%
  group_by(`Gene_name`) %>%
  summarize(all_missing = all(is.na(normalized_score))) %>%
  filter(all_missing) %>%
  pull(`Gene_name`)

# Calculate gene medians and sort
gene_medians <- data %>%
  filter(!`Gene_name` %in% missing_data_genes) %>%
  group_by(`Gene_name`) %>%
  summarize(median_score = median(normalized_score, na.rm = TRUE)) %>%
  arrange(median_score)

# Set gene order and convert to factor
gene_order <- c(missing_data_genes, gene_medians$`Gene_name`)
data$`Gene_name` <- factor(data$`Gene_name`, levels = gene_order)

# Create color gradient
orange_gradient_colors <- colorRampPalette(
  c("#FFF5EB", "#FDD0A2", "#FD8D3C", "#E65100", "#CB4335", "#922B21")
)(length(gene_order))

# Create visualization
p <- ggplot(data, aes(x = `Gene_name`, y = normalized_score, fill = `Gene_name`)) +
  geom_boxplot() +
  scale_fill_manual(values = orange_gradient_colors) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10, face = "bold"),
    axis.title.x = element_blank(),
    legend.position = "none",
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 14, face = "bold")
  ) +
  labs(
    y = "Normalized Score", 
    title = "Protein Expression Distribution in Cancer Tissues",
    subtitle = "Based on Human Protein Atlas IHC Scores"
  ) +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 15))

# Add missing data annotations
if (length(missing_data_genes) > 0) {
  p <- p + geom_text(
    data = data.frame(
      `Gene_name` = missing_data_genes, 
      y = rep(1, length(missing_data_genes))
    ),
    aes(x = `Gene_name`, y = y, label = " "),
    vjust = 0.5, 
    angle = 90, 
    size = 3, 
    color = "red"
  )
}

# Save output
ggsave("protein_expression_distribution.pdf", 
       p, 
       width = 14, 
       height = 4, 
       limitsize = FALSE, 
       device = cairo_pdf, 
       dpi = 300)








# Required libraries
library(readxl)
library(tidyr)
library(ggplot2)
library(dplyr)

rm(list=ls())

# Data processing
data <- read_excel("all_statistical_results.xlsx")

df <- data %>%
  select(File, 
         tumor = Tumor_median_gt_Normal_median, 
         normal = Tumor_median_lt_Normal_median) %>%
  arrange(tumor) %>%
  pivot_longer(cols = c(tumor, normal), 
               names_to = "variable", 
               values_to = "value")

# Create visualization
p <- ggplot(df, aes(
  x = factor(File, levels = unique(df$File[df$variable == "tumor"])),
  y = ifelse(variable == "tumor", value, -value),
  fill = variable
)) +
  geom_bar(stat = 'identity') +
  geom_text(
    aes(
      label = abs(value),
      hjust = ifelse(variable == "tumor", -0.4, 1.1)
    ),
    size = 3
  ) +
  scale_y_continuous(
    labels = abs,
    expand = expansion(mult = c(0.1, 0.1))
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    plot.title = element_text(size = 18, face = "bold")
  ) +
  labs(
    title = "Differential Protein Expression in CPTAC Cancer Cohorts",
    x = NULL, 
    y = "Number of Significant Proteins"
  ) +
  scale_fill_manual(
    values = c(tumor = "#D73027", normal = "#4575B4"),
    labels = c("Normal", "Tumor")
  )

# Save output
ggsave("protein_expression_distribution.pdf", 
       p, 
       width = 14, 
       height = 4, 
       units = "in")







