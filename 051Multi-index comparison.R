# Required packages
library(data.table)
library(dplyr)
library(pROC)
library(showtext)

# Process and combine datasets
process_data <- function(file1, file2) {
  data1 <- as.data.frame(fread(file1, header = TRUE, sep = "\t", check.names = FALSE))
  data1 <- data1[!duplicated(data1[,1]),]
  write.table(data1, file = "DATA1_unique.txt", sep = "\t", row.names = FALSE, quote = FALSE)
  
  data2 <- as.data.frame(fread(file2, header = TRUE, sep = "\t", check.names = FALSE))
  data2$number_ID <- as.character(data2$number_ID)
  data1$number_ID <- as.character(data1$number_ID)
  
  combined_data <- data2 %>% inner_join(data1, by = "number_ID")
  write.table(combined_data, file = "combined_data.txt", sep = "\t", row.names = FALSE, quote = FALSE)
  
  return(combined_data)
}

# Generate ROC curves for multiple features
generate_roc_curves <- function(data_file) {
  data <- as.data.frame(fread(data_file, header = TRUE, sep = "\t", check.names = FALSE))
  
  # Select features and outcome
  X <- data[, 2:18]
  y <- data[, 24]
  
  features <- colnames(X)[1:16]
  main_feature <- colnames(X)[17]
  colors <- rainbow(length(features))
  
  # Setup font
  font_add("SimHei", "SimHei.ttf")
  showtext_auto()
  
  # Generate ROC curves for each feature
  for (i in 1:length(features)) {
    feature <- features[i]
    color <- colors[i]
    
    # Handle missing values
    valid_indices <- which(!is.na(y) & !is.na(X[, i]) & !is.na(X[, main_feature]))
    y_valid <- y[valid_indices]
    X_valid <- X[valid_indices, i]
    main_X_valid <- X[valid_indices, main_feature]
    
    valid_indices <- which(!is.na(X_valid) & !is.na(main_X_valid))
    y_valid <- y_valid[valid_indices]
    X_valid <- as.numeric(as.character(X_valid[valid_indices]))
    main_X_valid <- as.numeric(as.character(main_X_valid[valid_indices]))
    
    # Calculate sample sizes
    benign_count <- sum(y_valid == "benign")
    malignant_count <- sum(y_valid == "malignant")
    
    # Generate ROC curves
    roc_obj <- roc(y_valid, X_valid)
    main_roc_obj <- roc(y_valid, main_X_valid)
    
    # Create plot
    pdf(paste0(feature, "_ROC.pdf"), width = 4.5, height = 4.5)
    plot(roc_obj, col = color, print.auc = FALSE, 
         main = paste("ROC Curve for", feature))
    plot(main_roc_obj, col = "red", print.auc = FALSE, add = TRUE)
    
    legend("bottomright", 
           legend = c(
             paste("Benign, n=", benign_count, ";", "Malignant, n=", malignant_count),
             paste(feature, "AUC:", round(roc_obj$auc, 3)),
             paste(main_feature, "AUC:", round(main_roc_obj$auc, 3))
           ),
           col = c("black", color, "red"), 
           lty = 1, 
           cex = 0.8)
    dev.off()
  }
}

# Execute analysis
combined_data <- process_data("DATA1.txt", "DATA2.txt")
generate_roc_curves("DATA5_new2.txt")




# Required packages
library(ggplot2)
library(pROC)

# Generate theoretical ROC curve data for a given AUC
create_roc_curve <- function(auc, num_points = 100) {
  fpr <- seq(0, 1, length.out = num_points)
  tpr <- sapply(fpr, function(x) x^((1/auc) - 1))
  return(data.frame(FPR = fpr, TPR = tpr))
}

# Define study parameters
auc_values <- c(0.737, 0.696, 0.82, 0.682)
curve_labels <- c(
  "MIF (AUC = 0.737)",
  "Shimu Luo et al. (AUC = 0.696)",
  "Emanuela Flamini et al.(AUC = 0.82)",
  "Hai Luo et al.(AUC = 0.682)"
)

# Generate ROC curves data
roc_data <- do.call(rbind, lapply(seq_along(auc_values), function(i) {
  curve_data <- create_roc_curve(auc_values[i])
  curve_data$Group <- curve_labels[i]
  return(curve_data)
}))

# Create visualization
p <- ggplot(roc_data, aes(x = FPR, y = TPR, color = Group)) +
  geom_line(size = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  labs(
    title = "CEA_CRC ROC Curves",
    x = "False Positive Rate",
    y = "True Positive Rate"
  ) +
  theme_minimal() +
  scale_color_manual(values = c("blue", "red", "green", "purple")) +
  theme(
    legend.position = c(0.75, 0.25),
    legend.title = element_blank(),
    legend.background = element_rect(fill = alpha('white', 0.5)),
    legend.box.background = element_rect(colour = "black"),
    legend.margin = margin(t = 10, r = 10, b = 10, l = 10),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank()
  ) +
  guides(color = guide_legend(override.aes = list(size = 1)))

# Save plot
ggsave("comparative_roc_curves.pdf", plot = p, width = 4, height = 4)







