library(pROC)

# Function to generate ROC curves for different comparison types
generate_roc_curves <- function(comparison_type, output_filename) {
  # Read and process data
  data <- read.table("Patent_all_data_English_select.txt", header = TRUE, sep = "\t", 
                     stringsAsFactors = FALSE)
  
  # Filter data based on comparison type
  if (comparison_type == "healthy_vs_malignant") {
    filtered_data <- subset(data, Benign_malignant %in% c("malignant", "healthy"))
    label_mapping <- c("malignant" = "Malignant", "healthy" = "Healthy")
    main_title <- "ROC Curves for Healthy vs. Malignant in Different Parts"
  } else if (comparison_type == "healthy_vs_benign") {
    filtered_data <- subset(data, Benign_malignant %in% c("benign", "healthy"))
    label_mapping <- c("benign" = "Benign", "healthy" = "Healthy")
    main_title <- "ROC Curves for Healthy vs. Benign in Different Parts"
  } else {
    filtered_data <- subset(data, Benign_malignant %in% c("benign", "malignant"))
    label_mapping <- c("benign" = "Benign", "malignant" = "Malignant")
    main_title <- "ROC Curves for Benign vs. Malignant in Different Parts"
  }
  
  # Define parts and colors
  parts <- c("uterus", "lung", "liver", "colon", "breast", "kidney", 
             "esophagus", "stomach", "rectum")
  colors <- c("blue", "yellow", "darkgreen", "green", "red", "pink", 
              "brown", "orange", "purple")
  
  # Create plot
  pdf(output_filename, width = 10, height = 10)
  par(mfrow = c(3, 3), mar = c(4, 4, 3, 2) + 0.1, oma = c(1, 1, 2, 1))
  
  for (i in seq_along(parts)) {
    part <- parts[i]
    color <- colors[i]
    
    if (comparison_type %in% c("healthy_vs_malignant", "healthy_vs_benign")) {
      data_part <- filtered_data[filtered_data$Part %in% c("healthy", part), ]
      data_part$Label <- factor(label_mapping[data_part$Benign_malignant], 
                                levels = c("Healthy", tail(label_mapping, 1)))
      roc_obj <- roc(data_part$Label, data_part$MIF_concentrate2)
    } else {
      data_part <- filtered_data[filtered_data$Part == part, ]
      roc_obj <- roc(data_part$Benign_malignant, data_part$MIF_concentrate2, 
                     levels = c("benign", "malignant"))
    }
    
    plot(1, type = "n", xlim = c(1, 0), ylim = c(0, 1),
         xlab = "1 - Specificity", ylab = "Sensitivity",
         main = paste("ROC Curve for", sub("_vs_", " vs. ", comparison_type), "in", part),
         cex.lab = 1.2, cex.main = 1.2, cex.axis = 1)
    
    abline(a = 0, b = 1, lty = 2, col = "gray")
    lines(roc_obj, col = color, lwd = 2)
    legend("bottomright", legend = paste("AUC =", round(auc(roc_obj), 3)), cex = 1)
    grid(lwd = 0.5, col = "lightgray")
  }
  
  title(main_title, outer = TRUE, line = -1, cex.main = 1.5)
  dev.off()
}

# Generate ROC curves for each comparison
generate_roc_curves("healthy_vs_malignant", "ROC_Curve_healthy_malignant_3x3.pdf")
generate_roc_curves("healthy_vs_benign", "ROC_Curve_healthy_benign_3x3.pdf")
generate_roc_curves("benign_vs_malignant", "ROC_Curve_benign_malignant_3x3.pdf")






# Required packages
library(ggplot2)
library(tidyr)
library(dplyr)
library(scales)

# Parse mean and standard deviation from string format
parse_mean_sd <- function(str) {
  parts <- strsplit(str, "±")[[1]]
  mean <- as.numeric(gsub("[^0-9.]", "", parts[1]))
  sd <- as.numeric(gsub("[^0-9.]", "", parts[2]))
  return(c(mean = mean, sd = sd))
}

# Generate bar plot comparing healthy vs malignant samples
generate_comparison_plot <- function(data_file) {
  data_raw <- read.table(data_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  
  # Process data into long format
  data_long <- data.frame(
    row_id = rep(1:nrow(data_raw), each = 2),
    study_id = rep(data_raw$Study_ID, each = 2),
    disease = rep(data_raw$group1_name_e, each = 2),
    group = factor(rep(c("Control", "Cancer"), times = nrow(data_raw)), 
                   levels = c("Control", "Cancer"))
  )
  
  # Add mean and sd columns
  for(i in 1:nrow(data_raw)) {
    g1 <- parse_mean_sd(data_raw$group1_mean_sd[i])
    g2 <- parse_mean_sd(data_raw$group2_mean_sd[i])
    data_long$mean[c(2*i-1, 2*i)] <- c(g2["mean"], g1["mean"])
    data_long$sd[c(2*i-1, 2*i)] <- c(g2["sd"], g1["sd"])
  }
  
  # Calculate z-scores and create labels
  data_long <- data_long %>%
    group_by(disease) %>%
    mutate(
      mean_z = (mean - min(mean)) / sd(mean),
      sd_z = sd / mean(mean)
    ) %>%
    ungroup()
  
  data_long$disease_label <- paste0(data_long$disease, " (", data_long$study_id, ")")
  
  # Set disease order and add sample sizes
  disease_order <- data_long %>%
    filter(group == "Cancer") %>%
    arrange(mean_z) %>%
    pull(disease_label) %>%
    unique()
  
  data_long$disease_label <- factor(data_long$disease_label, levels = disease_order)
  
  data_long <- data_long %>%
    group_by(disease_label) %>%
    mutate(
      sample_size = ifelse(
        group == "Cancer",
        paste0("n=", data_raw$group1_n[row_id[1]]),
        paste0("n=", data_raw$group2_n[row_id[1]])
      )
    )
  
  # Create visualization
  p <- ggplot(data_long, aes(y = mean_z, x = disease_label, fill = group)) +
    geom_bar(stat = "identity", position = position_dodge(0.9), alpha = 0.7) +
    geom_errorbar(aes(ymin = pmax(mean_z - sd_z, 0), ymax = mean_z + sd_z),
                  position = position_dodge(0.9), width = 0.25) +
    geom_text(aes(label = sample_size, y = -0.1), 
              position = position_dodge(0.9),
              angle = 45,
              vjust = 1,
              hjust = 1,
              size = 2.5) +
    labs(title = "Healthy-Malignant Comparison from External Studies",
         x = "Disease Type",
         y = "Standardized Concentration",
         fill = "Group") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9, face = "bold"),
      plot.margin = margin(t = 10, r = 10, b = 50, l = 10),
      axis.text = element_text(size = 12, face = "bold"),
      plot.title = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 12, face = "bold"),
      legend.title = element_text(size = 12, face = "bold"),
      legend.position = "top"
    ) +
    coord_cartesian(clip = "off") +
    scale_y_continuous(expand = expansion(mult = c(0.2, 0.1))) +
    scale_fill_manual(values = c("Control" = "#0D57A1", "Cancer" = "#BA3A2F"), 
                      breaks = c("Control", "Cancer"),
                      labels = c("Healthy", "Malignant"))
  
  ggsave("disease_comparisons.pdf", width = 6, height = 5)
  return(data_raw)
}

# Generate ROC curves
generate_roc_curves <- function(data_raw) {
  create_roc_curve <- function(auc, num_points = 100) {
    fpr <- seq(0, 1, length.out = num_points)
    tpr <- sapply(fpr, function(x) x^((1/auc) - 1))
    return(data.frame(FPR = fpr, TPR = tpr))
  }
  
  nature_theme <- theme_minimal() + theme(
    text = element_text(size = 9, face = "bold"),
    axis.title = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 6.5, face = "bold"),
    plot.title = element_text(size = 9, face = "bold"),
    legend.position = c(0.55, 0.2),
    legend.title = element_blank(),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key.size = unit(0.8, "lines"),
    legend.spacing.y = unit(0.1, "cm"),
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
    panel.grid.major = element_line(color = "grey90", size = 0.2),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
  )
  
  n_groups <- ceiling(nrow(data_raw) / 4)
  for(group in 1:n_groups) {
    start_idx <- (group-1)*4 + 1
    end_idx <- min(group*4, nrow(data_raw))
    current_data <- data_raw[start_idx:end_idx, ]
    roc_data <- data.frame()
    
    for(i in 1:nrow(current_data)) {
      auc <- as.numeric(current_data$MEDIAN5AUC_original[i])
      curve_data <- create_roc_curve(auc)
      curve_data$Group <- paste0(current_data$group1_name_e[i], 
                                 "_", current_data$Study_ID[i],
                                 " (AUC: ", round(auc, 3), ")")
      roc_data <- rbind(roc_data, curve_data)
    }
    
    p2 <- ggplot(roc_data, aes(x = FPR, y = TPR, color = Group)) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", 
                  color = "grey50", size = 0.3) +
      geom_line(size = 0.75) +
      labs(x = "1 - Specificity",
           y = "Sensitivity",
           title = "ROC Curves for Healthy vs. Malignant Parts") +
      scale_x_continuous(breaks = seq(0, 1, 0.2),
                         labels = number_format(accuracy = 0.1),
                         expand = c(0.01, 0)) +
      scale_y_continuous(breaks = seq(0, 1, 0.2),
                         labels = number_format(accuracy = 0.1),
                         expand = c(0.01, 0)) +
      scale_color_manual(values = c("#E64B35", "#4DBBD5", "#4DAF4A", "#984EA3")) +
      nature_theme +
      guides(color = guide_legend(override.aes = list(size = 1)))
    
    ggsave(paste0("roc_curves_group_", group, ".pdf"), 
           plot = p2, width = 3, height = 3, dpi = 300)
  }
}

# Run analysis
data_raw <- generate_comparison_plot("plot_data.txt")
generate_roc_curves(data_raw)


