# Required packages
library(ggplot2)
library(tidyr)
library(dplyr)
library(gridExtra)
library(ggsci)
library(scales)

rm(list=ls())

# Data import
boxplot_data <- read.delim("boxplot_data.txt", header=TRUE, sep="\t")
sample_size <- read.delim("sample_size.txt", header=TRUE, sep="\t", check.names=FALSE)

# Create statistical results dataframe with significance markers
stat_results <- data.frame(
  Organ = unique(boxplot_data$Organ),
  p.value = rep(0.001, length(unique(boxplot_data$Organ))),
  y.position = tapply(boxplot_data$MAX, boxplot_data$Organ, max) * 1.0
)

# Create sample size labels
sample_labels <- sample_size_long <- pivot_longer(sample_size, 
                                                  cols = -ID, 
                                                  names_to = "Organ", 
                                                  values_to = "Size") %>%
  group_by(Organ) %>%
  summarize(label = paste0("n=", paste(Size, collapse="/")))

# Set y-axis limits
y_min <- min(boxplot_data$MIN)
y_max <- 17500

# Create base plot
p <- ggplot(boxplot_data, aes(x=Organ, y=Median, fill=Types)) +
  geom_boxplot(
    aes(
      ymin = MIN,
      lower = Q1,
      middle = Median,
      upper = Q3,
      ymax = MAX
    ),
    stat = "identity",
    position = position_dodge(0.8),
    width = 0.7
  ) +
  scale_fill_manual(values = c("#4682B4", "#E64B35")) +
  theme_minimal() +
  theme(
    text = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14, face = "bold"),
    legend.position = "top",
    legend.title = element_blank(),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    plot.caption = element_text(size = 14, face = "bold", hjust = 0)
  ) +
  labs(
    x = "Tissue Type",
    y = "Expression Level",
    title = "MIF mRNA Expression Across Human Malignancies",
    caption = "(*) Mann-Whitney p<0.05 and expression >10 in tumor or normal"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2)))

# Add sample size labels
p <- p + geom_text(
  data = sample_labels,
  aes(x = Organ, y = y_min, label = label),
  vjust = 1.5,
  size = 2.8,
  fontface = "bold",
  inherit.aes = FALSE
)

# Add significance markers
p <- p + geom_text(
  data = stat_results,
  aes(x = Organ, y = y.position, label = "*"),
  inherit.aes = FALSE,
  size = 5,
  color = "red"
)

# Save plot
ggsave(
  "nature_style_boxplot.pdf",
  p,
  width = 14,
  height = 4,
  units = "in",
  dpi = 300,
  device = cairo_pdf
)






# Required libraries
library(ggplot2)
library(tidyr)
library(dplyr)
library(readr)
library(ggpubr)
library(writexl)

rm(list=ls())

# Extract group name from filename and column name
extract_group_name <- function(filename) {
  data <- read.delim(filename, sep="\t", header=TRUE, nrows=1)
  col_name <- names(data)[3]
  gseid <- sub(".*GSE([0-9]+)_.*", "\\1", filename)
  
  gpl_match <- regexpr("GPL[0-9]+", filename)
  gpl_suffix <- ""
  if (gpl_match > 0) {
    gpl_number <- regmatches(filename, gpl_match)
    gpl_suffix <- paste0("_", gpl_number)
  }
  
  return(paste0(col_name, "_", gseid, gpl_suffix))
}

# Process expression data with z-score normalization
read_expression_data <- function(file_path) {
  data <- read.delim(file_path, sep="\t", header=TRUE)
  group_name <- extract_group_name(file_path)
  
  result <- data.frame(
    Sample = 1:nrow(data),
    Group = group_name,
    Tumor = as.numeric(data[,3]),
    Normal = as.numeric(data[,4])
  )
  
  normal_mean <- mean(result$Normal)
  normal_sd <- sd(result$Normal)
  
  result$Tumor <- (result$Tumor - normal_mean) / normal_sd
  result$Normal <- (result$Normal - normal_mean) / normal_sd
  
  return(result)
}

# Process all clinical data files
txt_files <- list.files(pattern = ".*clinical_.*\\.txt$")
all_data <- do.call(rbind, lapply(txt_files, read_expression_data))

# Prepare data for plotting
long_data <- all_data %>%
  pivot_longer(cols = c(Tumor, Normal),
               names_to = "Type",
               values_to = "Expression") %>%
  mutate(Type = factor(Type, levels = c("Normal", "Tumor")),
         Group = factor(Group)) %>%
  group_by(Group) %>%
  mutate(x_position = as.numeric(Group) + if_else(Type == "Normal", -0.2, 0.2))

# Create paired boxplot
p <- ggplot(long_data, aes(x = Group, y = Expression)) +
  geom_boxplot(aes(fill = Type), 
               position = position_dodge(width = 0.8), 
               alpha = 0.5, 
               outlier.shape = NA, 
               width = 0.7) +
  geom_line(aes(x = x_position, group = interaction(Sample, Group)), 
            color = "grey", 
            alpha = 0.3) +
  geom_point(aes(x = x_position, color = Type),
             size = 1, 
             alpha = 0.5) +
  scale_fill_manual(values = c("Normal" = "#0B3D91", "Tumor" = "#B22222")) +
  scale_color_manual(values = c("Normal" = "#0B3D91", "Tumor" = "#B22222")) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 10),
    axis.text.y = element_text(face = "bold", size = 10),
    axis.title = element_text(face = "bold", size = 12),
    legend.text = element_text(face = "bold", size = 10),
    legend.title = element_text(face = "bold", size = 12),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.position = "top",
    plot.title = element_text(face = "bold", size = 14)
  ) +
  labs(
    y = "Expression Level (Z-score)",
    x = NULL,
    title = "MIF mRNA Expression in GEO Paired Samples"
  ) +
  stat_compare_means(aes(group = Type),
                     paired = TRUE,
                     method = "wilcox.test",
                     label = "p.signif",
                     vjust = 0.65)

# Save plot with dynamic width
plot_width <- max(10, length(txt_files) * 2.5) * 1/6
ggsave("paired_expression_comparison.pdf", p, width = plot_width, height = 3, limitsize = FALSE)

# Calculate and export expression ratios
median_ratios <- long_data %>%
  group_by(Group, Type) %>%
  summarise(median_expr = median(Expression), .groups = 'drop') %>%
  pivot_wider(names_from = Type, values_from = median_expr) %>%
  mutate(
    median_ratio = Tumor / Normal,
    formatted_ratio = sprintf("%.3f", median_ratio)
  ) %>%
  arrange(desc(abs(median_ratio)))

# Export results
excel_data <- median_ratios %>%
  select(
    Dataset = Group,
    `Normal Median` = Normal,
    `Tumor Median` = Tumor,
    `Tumor/Normal Ratio` = median_ratio
  )

write_xlsx(excel_data, "median_expression_ratios.xlsx")






# Required packages
library(ggplot2)
library(tidyr)
library(dplyr)
library(gridExtra)
library(ggsci)
library(scales)
library(rstatix)
library(stats)

rm(list=ls())

# Load and process data
data <- read.delim("MIF_PlotData.txt", header=TRUE, sep="\t")

# Statistical analysis with Bonferroni correction
stat_results <- data %>%
  group_by(CancerType) %>%
  wilcox_test(MIF ~ Type) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance() %>%
  mutate(
    effect_size = wilcox_effsize(data = data, MIF ~ Type, grouping.var = CancerType)$effsize,
    direction = case_when(
      p.adj < 0.05 & statistic > 0 ~ "Up",
      p.adj < 0.05 & statistic < 0 ~ "Down",
      TRUE ~ "not_sig"
    )
  )

# Calculate median comparisons
median_comparison <- data %>%
  group_by(CancerType, Type) %>%
  summarise(median_MIF = median(MIF), .groups = "drop") %>%
  pivot_wider(names_from = Type, values_from = median_MIF) %>%
  mutate(
    tumor_gt_normal = Tumor > Normal,
    y.position = max(Tumor) * 2.5
  )

# Merge statistical results
stat_results <- stat_results %>%
  left_join(median_comparison, by = "CancerType") %>%
  mutate(
    significance = case_when(
      p.adj < 0.001 ~ "***",
      p.adj < 0.01 ~ "**",
      p.adj < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  ) %>%
  filter(p.adj < 0.05 & tumor_gt_normal)

# Filter significant data
significant_cancertypes <- stat_results$CancerType
data_filtered <- data %>% 
  filter(CancerType %in% significant_cancertypes)

# Calculate sample sizes
sample_size <- data_filtered %>%
  group_by(CancerType, Type) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = Type, values_from = n) %>%
  mutate(label = paste0("n=", Normal, "/", Tumor))

# Set y-axis limits
y_min <- min(data_filtered$MIF)
y_max <- max(data_filtered$MIF) * 1.2

# Create visualization
p <- ggplot(data_filtered, aes(x=CancerType, y=MIF, fill=Type)) +
  geom_boxplot(
    position = position_dodge(0.8),
    width = 0.7,
    outlier.shape = NA
  ) +
  scale_fill_manual(values = c("#AF7AC5", "#B37B49")) +
  theme_minimal() +
  theme(
    text = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14, face = "bold"),
    legend.position = "top",
    legend.title = element_blank(),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    plot.caption = element_text(size = 14, face = "bold", hjust = 0)
  ) +
  labs(
    x = "Tissue Type",
    y = "Expression Level",
    title = "MIF Protein Levels Across CPTAC Cancer Cohorts",
    caption = "* p < 0.05, ** p < 0.01, *** p < 0.001\n(p-values adjusted using Bonferroni correction)"
  ) +
  scale_y_continuous(
    limits = c(y_min, y_max),
    expand = expansion(mult = c(0.1, 0.15))
  )

# Add sample size labels
p <- p + geom_text(
  data = sample_size,
  aes(x = CancerType, y = y_min, label = label),
  vjust = 1.5,
  size = 3.2,
  fontface = "bold",
  inherit.aes = FALSE
)

# Add significance markers
p <- p + geom_text(
  data = stat_results,
  aes(x = CancerType, y = y.position, label = significance),
  inherit.aes = FALSE,
  size = 5,
  vjust = -1,
  color = "red"
)

# Save plot
ggsave(
  "protein_expression_analysis.pdf",
  p,
  width = 6,
  height = 4.5,
  units = "in",
  dpi = 300,
  device = cairo_pdf
)






