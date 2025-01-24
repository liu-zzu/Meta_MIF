# Required libraries
library(ggplot2)
library(scales)
library(readr)

# Data processing
data <- read_lines("fig1c.txt")[-1]
years <- as.numeric(data)
pub_counts <- as.data.frame(table(years))
names(pub_counts) <- c("Year", "Count")
pub_counts$Year <- as.numeric(as.character(pub_counts$Year))

# Create decade labels
pub_counts$label <- ifelse(pub_counts$Year %% 10 == 0, pub_counts$Count, "")

# Visualization
ggplot(pub_counts, aes(x = Year, y = Count)) +
  geom_bar(stat = "identity", fill = "#4292c6", alpha = 0.8) +
  geom_text(aes(label = label), vjust = -0.5, size = 6) +
  geom_smooth(method = "loess", color = "#de2d26", size = 2.1, se = FALSE) +
  theme_minimal() +
  theme(
    text = element_text(size = 16, face = "bold"),
    axis.text = element_text(color = "black", size = 16, face = "bold"),
    axis.title = element_text(size = 18, face = "bold"),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 19, face = "bold"),
    plot.subtitle = element_text(size = 16, face = "bold")
  ) +
  labs(
    title = "Cancer Liquid Biopsy Literature Trends",
    subtitle = "PubMed Publications (1950-2024)", 
    x = "Year",
    y = "Number of Papers"
  ) +
  scale_y_continuous(labels = comma, expand = c(0.1, 0)) +
  scale_x_continuous(breaks = seq(1950, 2025, by = 10))

# Save plot
ggsave("fig1c.pdf", width = 5.5, height = 4, dpi = 300)










# Required libraries
library(readxl)
library(dplyr)
library(tidyr)
library(networkD3)
library(htmlwidgets)
library(webshot)

rm(list=ls())

# Read and process Excel data
data <- read_excel("04_all_year.processed_data_with_stats2_all_updated_number_format.14.15.16.17.18.19_1_2.xlsx", 
                   col_names = FALSE) %>%
  select(14, 15) %>%
  rename(primary_method = `...14`, specific_method = `...15`)

# Split multiple techniques
split_techniques <- function(x) {
  unlist(strsplit(as.character(x), ", "))
}

# Process data and calculate frequencies
processed_data <- data %>%
  mutate(specific_method = sapply(specific_method, split_techniques)) %>%
  unnest(specific_method) %>%
  group_by(primary_method, specific_method) %>%
  summarise(value = n()) %>%
  ungroup()

# Get top 15 categories for each column
top_main <- processed_data %>%
  group_by(primary_method) %>%
  summarise(total = sum(value)) %>%
  top_n(15, total) %>%
  pull(primary_method)

top_specific <- processed_data %>%
  group_by(specific_method) %>%
  summarise(total = sum(value)) %>%
  top_n(15, total) %>%
  pull(specific_method)

# Filter data for top categories
filtered_data <- processed_data %>%
  filter(primary_method %in% top_main & specific_method %in% top_specific)

# Calculate node totals and font sizes
source_totals <- filtered_data %>%
  group_by(primary_method) %>%
  summarise(total = sum(value))

target_totals <- filtered_data %>%
  group_by(specific_method) %>%
  summarise(total = sum(value))

base_font_size <- 17
max_increment <- 22

calculate_font_size <- function(total, min_total, max_total, base_size = base_font_size, max_inc = max_increment) {
  base_size + (total - min_total) / (max_total - min_total) * max_inc
}

# Calculate font sizes based on totals
all_totals <- c(source_totals$total, target_totals$total)
min_total <- min(all_totals)
max_total <- max(all_totals)

source_totals$fontSize <- calculate_font_size(source_totals$total, min_total, max_total)
target_totals$fontSize <- calculate_font_size(target_totals$total, min_total, max_total)

# Create nodes and links data
nodes <- data.frame(
  name = c(source_totals$primary_method, target_totals$specific_method),
  fontSize = c(source_totals$fontSize, target_totals$fontSize)
)

links <- filtered_data %>%
  mutate(
    source = match(primary_method, nodes$name) - 1,
    target = match(specific_method, nodes$name) - 1
  )

# Custom JavaScript for styling
custom_js <- sprintf('
  function(el, x) {
    var fontSizes = %s;
    d3.select(el)
      .selectAll(".node text")
      .style("font-size", function(d, i) { return fontSizes[i] + "px"; })
      .style("font-weight", "bold");
    
    d3.select(el)
      .selectAll(".link")
      .style("stroke", "#444444")
      .style("stroke-opacity", function(d) {
        var maxValue = d3.max(x.links, function(d) { return d.value; });
        var minValue = d3.min(x.links, function(d) { return d.value; });
        return 0.2 + (d.value - minValue) / (maxValue - minValue) * 0.6;
      });
  }
', jsonlite::toJSON(nodes$fontSize))

# Create and customize Sankey diagram
sankey_network <- sankeyNetwork(Links = links, 
                                Nodes = nodes,
                                Source = "source", 
                                Target = "target",
                                Value = "value", 
                                NodeID = "name",
                                sinksRight = FALSE, 
                                nodeWidth = 30)

sankey_network <- onRender(sankey_network, custom_js)

# Save outputs
saveWidget(sankey_network, "sankey_diagram.html", selfcontained = TRUE)
webshot("sankey_diagram.html", 
        "sankey_diagram.pdf",
        delay = 2,
        zoom = 2, 
        vwidth = 1060,
        vheight = 500)






# Required libraries
library(readxl)
library(dplyr)
library(stringr)
library(wordcloud2)
library(RColorBrewer)
library(htmlwidgets)

rm(list=ls())

# Data processing
data <- read_excel("04_all_year.processed_data_with_stats2_all_updated_number_format.14.15.16.17.18.19_1_2.xlsx", 
                   col_names = T)

text <- data[[16]]
words <- unlist(strsplit(paste(text, collapse = ","), ","))
words <- str_trim(words)

# Calculate frequencies
word_freq <- table(words)
df <- data.frame(word = names(word_freq), 
                 freq = as.numeric(word_freq)) %>% 
  arrange(desc(freq))

# Generate PDF output
pdf("wordcloud.pdf", width = 8.5, height = 6.5, pointsize = 12)
par(mar = c(0,0,0,0), bg = "white")

# Color scheme
colors <- c("#2E5A87", "#4682B4", "#6495ED", "#87CEEB", "#B0E0E6")

# Create word cloud
wordcloud2(data = df,
           size = 0.6,
           minSize = 0.3,
           color = rep(colors, length.out = nrow(df)),
           backgroundColor = "white",
           fontFamily = "Arial",
           fontWeight = "bold",
           shuffle = FALSE,
           rotateRatio = 0)

# Save outputs
saveWidget(wordcloud2(df), "wordcloud.html", selfcontained = TRUE)
write.csv(df, "word_frequencies.csv", row.names = FALSE)
dev.off()




# Required libraries
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
rm(list=ls())

# Data import and preprocessing
data <- read_excel("04_all_year.processed_data_with_stats2_all_updated_number_format.14.15.16.17.18.19_1_2.xlsx")

df <- data %>%
  select(17, 18) %>%
  rename(Cancer_Type = 1, Biomarker_Type = 2)

# Process and expand data
df <- df %>%
  mutate(Cancer_Type = ifelse(Cancer_Type == "Unspecified", "Unspecified", Cancer_Type)) %>%
  mutate(Cancer_Type = strsplit(as.character(Cancer_Type), ",\\s*")) %>%
  mutate(Biomarker_Type = strsplit(as.character(Biomarker_Type), ",\\s*")) %>%
  unnest(Cancer_Type) %>%
  unnest(Biomarker_Type) %>%
  filter(Biomarker_Type != "none" & Biomarker_Type != "") %>%
  mutate(across(everything(), trimws))

# Select top categories
top_12_cancer_types <- df %>%
  count(Cancer_Type) %>%
  top_n(12, n) %>%
  pull(Cancer_Type)

top_12_biomarker_types <- df %>%
  count(Biomarker_Type) %>%
  top_n(12, n) %>%
  pull(Biomarker_Type)

# Filter and calculate frequencies
filtered_df <- df %>%
  filter(Cancer_Type %in% top_12_cancer_types, 
         Biomarker_Type %in% top_12_biomarker_types)

freq_data <- filtered_df %>%
  group_by(Cancer_Type, Biomarker_Type) %>%
  summarise(Count = n(), .groups = "drop")

# Color scheme
colors <- c(
  "#6B8CC2", "#72B484", "#D86F73", "#9A89C9", 
  "#D9C364", "#52A3C9", "#6B8CC2", "#72B484", 
  "#D86F73", "#9A89C9", "#D9C364", "#52A3C9"
)

# Create visualization
pdf("cancer_biomarker_distribution.pdf", width = 6, height = 5)
ggplot(freq_data, aes(x = reorder(Cancer_Type, -Count), y = Count, fill = Biomarker_Type)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size = 13, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 10, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 16, face = "bold")
  ) +
  labs(title = "Distribution of Biomarker Types Across Cancer Types",
       x = "Cancer Type", 
       y = "Paper Count", 
       fill = "Biomarker Type") +
  scale_fill_manual(values = colors)
dev.off()

# Export processed data
write.csv(freq_data, "biomarker_cancer_frequency.csv", row.names = FALSE)





# Required libraries
library(readxl)
library(ggplot2)
library(dplyr)
library(scales)
library(ggrepel)
library(showtext)
library(extrafont)
library(viridis)
library(gridExtra)
library(ggthemes)

rm(list=ls())

# Font setup
font_add_google("Roboto", "Roboto")
showtext_auto()

# Data processing
data <- read_excel("04_all_year.processed_data_with_stats2_all_updated_number_format.14.15.16.17.18.19_1_2.xlsx")
df <- data.frame(sample_size = data[[20]]) %>% 
  mutate(sample_size = as.numeric(sample_size)) %>%
  filter(!is.na(sample_size) & sample_size > 0)

# Calculate statistics
stats <- list(
  n_samples = nrow(df),
  median = median(df$sample_size),
  mean = mean(df$sample_size),
  max = max(df$sample_size),
  min = min(df$sample_size),
  q1 = quantile(df$sample_size, 0.25),
  q3 = quantile(df$sample_size, 0.75)
)

# Custom theme definition
custom_theme <- function() {
  theme_minimal() %+replace%
    theme(
      text = element_text(family = "Roboto", color = "#000000", size = 14),
      plot.title = element_text(size = 19, face = "bold", hjust = 0.5, 
                                margin = margin(b = 20), color = "#000000"),
      plot.subtitle = element_text(size = 18, hjust = 0.5, 
                                   margin = margin(b = 10), color = "#000000"),
      axis.title.y = element_text(size = 18, margin = margin(r = 10), 
                                  angle = 90, face = "bold", color = "#000000"),
      axis.text.y = element_text(size = 16, angle = 90, 
                                 face = "bold", color = "#000000"),
      axis.text.x = element_blank(),
      panel.grid.major.y = element_line(color = "#ecf0f1", 
                                        linewidth = 0.2, linetype = "dashed"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      legend.position = "none",
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(30, 80, 30, 30)
    )
}

# Create visualization
p <- ggplot(df, aes(x = "", y = sample_size)) +
  stat_boxplot(geom = 'errorbar', width = 0.2, color = "#2980b9") +
  geom_boxplot(
    fill = "#3498db",
    alpha = 0.6,
    color = "#2980b9",
    outlier.shape = 21,
    outlier.fill = "white",
    outlier.color = "#2980b9",
    outlier.alpha = 0.6,
    outlier.size = 4,
    linewidth = 0.8
  ) +
  scale_y_log10(
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x))
  ) +
  geom_hline(
    yintercept = stats$mean,
    linetype = "dashed",
    color = "#e74c3c",
    linewidth = 0.8,
    alpha = 0.8
  ) +
  labs(
    title = "Distribution of Sample Sizes",
    subtitle = paste("n =", format(stats$n_samples, big.mark = ",")),
    y = "Sample Size (log scale)",
    x = NULL
  ) +
  custom_theme()

# Add annotations
annotations <- tibble(
  x = 4,
  y = c(stats$max, stats$q3, stats$median, stats$q1, stats$min, stats$mean),
  label = c(
    sprintf("Maximum: %s", format(stats$max, big.mark = ",")),
    sprintf("Q3: %s", format(round(stats$q3, 1), big.mark = ",")),
    sprintf("Median: %s", format(round(stats$median, 1), big.mark = ",")),
    sprintf("Q1: %s", format(round(stats$q1, 1), big.mark = ",")),
    sprintf("Minimum: %s", format(stats$min, big.mark = ",")),
    sprintf("Mean: %s", format(round(stats$mean, 1), big.mark = ","))
  ),
  color = c(rep("#000000", 5), "#e74c3c")
)

# Add annotation layer
p <- p + geom_text_repel(
  data = annotations,
  aes(x = x, y = y, label = label, color = I(color)),
  size = 5,
  family = "Roboto",
  fontface = "bold",
  direction = "y",
  hjust = 0,
  segment.size = 0.4,
  segment.color = "#95a5a6",
  segment.alpha = 0.6,
  box.padding = unit(0.5, "lines"),
  point.padding = unit(0.3, "lines"),
  force = 2
)

# Save output
ggsave(
  "sample_size_distribution.pdf",
  p,
  width = 4.5,
  height = 5,
  units = "in",
  dpi = 300,
  device = cairo_pdf
)






# Required libraries
library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(writexl)

rm(list=ls())

# Analysis of tumor markers and cancer types
analyze_marker_types <- function(file_path) {
  data <- read_excel(file_path, col_names = FALSE)
  
  column_data <- data[[19]]
  
  processed_data <- column_data %>%
    strsplit(",") %>%
    unlist() %>%
    .[str_detect(., "\\|")] %>%
    strsplit("\\|") %>%
    do.call(rbind, .) %>%
    as.data.frame() %>%
    setNames(c("Marker", "Cancer")) %>%
    mutate(across(everything(), trimws))
  
  result <- processed_data %>%
    group_by(Marker) %>%
    summarise(
      Unique_Cancer_Types = n_distinct(Cancer),
      Cancer_Types = paste(unique(Cancer), collapse = ", ")
    )
  
  write_xlsx(result, "marker_analysis_result.xlsx")
  return(result)
}

# Analysis of TCGA tumor types
analyze_tcga_types <- function(file_path) {
  data <- read_excel(file_path)
  
  process_cell <- function(cell) {
    groups <- unlist(strsplit(cell, ",\\s*"))
    extracted <- sapply(groups, function(group) {
      parts <- unlist(strsplit(group, ":"))
      if(length(parts) > 1) {
        return(trimws(gsub('"', '', parts[2])))
      } else {
        return(NA)
      }
    })
    return(length(unique(extracted)))
  }
  
  data$Type_Count <- sapply(data$API_result1, process_cell)
  
  write_xlsx(data, "tcga_analysis_result.xlsx")
  return(data)
}

# Execute analyses
marker_results <- analyze_marker_types("tumor_marker_data.xlsx")
tcga_results <- analyze_tcga_types("tcga_data.xlsx")

# Print summary statistics
print("Marker Analysis Results:")
print(head(marker_results))

print("TCGA Analysis Results:")
print(summary(tcga_results$Type_Count))
print(table(tcga_results$Type_Count))



# Required libraries
library(readxl)
library(ggplot2)
library(dplyr)
library(scales)
library(ggrepel)

rm(list=ls())

# Data processing
data <- read_excel("biomarker_analysis_result.xlsx")
colnames(data) <- c("Biomarker", "Associated_Cancer_Types", "Associated_Cancer_Type2s")

data <- data %>% 
  arrange(desc(Associated_Cancer_Types)) %>%
  mutate(
    percentile_rank = percent_rank(Associated_Cancer_Types),
    label_text = if_else(Associated_Cancer_Types >= quantile(Associated_Cancer_Types, 0.95),
                         paste0(Biomarker, "\n(", Associated_Cancer_Types, " cancer types)"),
                         NA_character_)
  )

top_140 <- head(data, 140)

# Calculate statistics
stats_summary <- summarise(top_140,
                           Total_Markers = n(),
                           Mean_Cancer_Types = mean(Associated_Cancer_Types),
                           Median_Cancer_Types = median(Associated_Cancer_Types),
                           SD_Cancer_Types = sd(Associated_Cancer_Types),
                           Max_Cancer_Types = max(Associated_Cancer_Types),
                           Min_Cancer_Types = min(Associated_Cancer_Types))

# Create visualization
p <- ggplot(top_140, aes(x = reorder(Biomarker, -Associated_Cancer_Types), 
                         y = Associated_Cancer_Types,
                         color = percentile_rank)) +
  geom_segment(aes(xend = Biomarker, 
                   y = 0, 
                   yend = Associated_Cancer_Types),
               color = "gray90", size = 0.3) +
  geom_point(size = 4, alpha = 0.9) +
  geom_label_repel(aes(label = label_text),
                   na.rm = TRUE,
                   size = 3,
                   box.padding = 0.5,
                   point.padding = 0.5,
                   force = 10,
                   segment.color = "gray50") +
  scale_color_gradient(low = "#BFD3E6", high = "#08519C") +
  theme_minimal() +
  theme(
    text = element_text(color = "black"),
    axis.text.x = element_text(size = 8, angle = 45, hjust = 1, vjust = 1, face = "bold"),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 20, face = "bold"),
    plot.subtitle = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 0.3),
    axis.ticks = element_line(color = "black", size = 0.3),
    plot.margin = margin(30, 30, 60, 30),
    legend.position = "top",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  ) +
  labs(
    x = "Biomarker",
    y = "Number of Associated Cancer Types",
    title = "Distribution of Biomarkers by Associated Cancer Types",
    subtitle = paste(nrow(top_140), "unique biomarkers filtered by presence in more than 15 cancer types"),
    color = "Percentile Rank"
  ) +
  scale_y_continuous(
    limits = c(15, max(top_140$Associated_Cancer_Types) + 1),
    expand = expansion(mult = c(0, 0.15)),
    breaks = pretty_breaks(n = 10)
  )

# Save output
ggsave("biomarker_distribution.pdf", 
       p, 
       width = 14,
       height = 6,
       units = "in",
       dpi = 300,
       limitsize = FALSE)

