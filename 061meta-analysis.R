


# Required packages
library(mada)
library(meta)
library(metafor)
library(ggplot2)
library(dplyr)

# Perform publication bias analysis
analyze_publication_bias <- function(data_file, correction = 0.5) {
  # Read and process data
  data <- read.delim(data_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # Apply continuity correction
  data$TP <- data$TP + correction
  data$FP <- data$FP + correction
  data$TN <- data$TN + correction
  data$FN <- data$FN + correction
  
  # Calculate metrics
  data$ESS <- sqrt(data$TP + data$FP + data$TN + data$FN)
  data$DOR <- (data$TP * data$TN) / (data$FP * data$FN)
  data$lnDOR <- log(data$DOR)
  
  return(data)
}

# Generate Deeks' funnel plot
create_funnel_plot <- function(data, output_file) {
  plot <- ggplot(data, aes(x = lnDOR, y = ESS)) +
    theme_bw() +
    geom_point(size = 3, shape = 21, fill = "blue", alpha = 0.6) +
    geom_smooth(method = "lm", col = "red", se = TRUE, linetype = "dashed") +
    labs(
      title = "Deeks' Funnel Plot",
      x = "Log Diagnostic Odds Ratio",
      y = "Root of Effective Sample Size",
      subtitle = paste("Studies:", nrow(data))
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 12, face = "bold"),
      legend.position = "none",
      panel.grid.major = element_line(colour = "grey90"),
      panel.grid.minor = element_line(colour = "grey95")
    )
  
  # Perform Deeks' test
  model <- lm(ESS ~ lnDOR, data = data)
  p_value <- format.pval(summary(model)$coefficients[2,4], digits = 3)
  
  # Add p-value annotation
  plot <- plot +
    annotate(
      "text",
      x = min(data$lnDOR),
      y = max(data$ESS),
      label = paste("Deeks' test p-value =", p_value),
      hjust = 0,
      vjust = 1,
      size = 6,
      fontface = "bold"
    )
  
  pdf(output_file, width = 4, height = 4)
  print(plot)
  dev.off()
  
  return(model)
}

# Generate analysis report
generate_report <- function(data, model, output_file) {
  sink(output_file)
  cat("Publication Bias Analysis Report\n")
  cat("================================\n\n")
  
  # Study information
  cat("1. Study Information:\n")
  cat("   Number of studies:", nrow(data), "\n\n")
  
  # Summary statistics
  cat("2. Summary Statistics:\n")
  cat("   Mean lnDOR:", round(mean(data$lnDOR), 4), "\n")
  cat("   Median lnDOR:", round(median(data$lnDOR), 4), "\n")
  cat("   SD lnDOR:", round(sd(data$lnDOR), 4), "\n")
  cat("   Range lnDOR:", round(min(data$lnDOR), 4), "to", round(max(data$lnDOR), 4), "\n\n")
  
  # Test results
  deeks_test <- summary(model)
  p_value <- format.pval(deeks_test$coefficients[2,4], digits = 3)
  cat("3. Deeks' Test Results:\n")
  cat("   P-value:", p_value, "\n")
  cat("   Slope coefficient:", round(deeks_test$coefficients[2,1], 4), "\n")
  cat("   R-squared:", round(deeks_test$r.squared, 4), "\n\n")
  
  # Interpretation
  cat("4. Interpretation:\n")
  if(as.numeric(sub("<", "", p_value)) < 0.05) {
    cat("   Evidence of significant publication bias (p < 0.05)\n")
  } else {
    cat("   No significant evidence of publication bias (p ≥ 0.05)\n")
  }
  sink()
  
  # Export summary data
  summary_data <- data.frame(
    Study_ID = 1:nrow(data),
    ESS = round(data$ESS, 2),
    lnDOR = round(data$lnDOR, 2),
    DOR = round(data$DOR, 2)
  )
  write.csv(summary_data, "analysis_summary.csv", row.names = FALSE)
}

# Execute analysis
data <- analyze_publication_bias("deeks_data.txt")
model <- create_funnel_plot(data, "deeks_funnel_plot.pdf")
generate_report(data, model, "publication_bias_report.txt")






# Required packages
library(meta)
library(mada)
library(dplyr)
library(ggplot2)

# Function to create forest plots for sensitivity and specificity
generate_forest_plots <- function(data_file, output_file) {
  # Read and process data
  data <- read.table(data_file, header = TRUE, sep = "\t", 
                     quote = "", stringsAsFactors = FALSE)
  
  data <- data %>%
    mutate(studyid = paste0(Author, " ", Year, "_", row_number()))
  
  # Prepare sensitivity data
  sens_data <- metaprop(
    event = TP,
    n = TP + FN,
    studlab = studyid,
    data = data,
    sm = "PLOGIT",
    method.tau = "ML",
    common = FALSE,
    random = TRUE,
    prediction = TRUE
  )
  
  # Prepare specificity data
  spec_data <- metaprop(
    event = TN,
    n = TN + FP,
    studlab = studyid,
    data = data,
    sm = "PLOGIT",
    method.tau = "ML",
    common = FALSE,
    random = TRUE,
    prediction = TRUE
  )
  
  # Set common forest plot parameters
  forest_params <- list(
    prediction = TRUE,
    hetstat = TRUE,
    test.overall = TRUE,
    digits = 2,
    print.I2 = TRUE,
    fs.study = 12,
    fs.heading = 14,
    fs.random = 12,
    fs.hetstat = 12,
    fs.axis = 12,
    squaresize = 0.8,
    smlab = ""
  )
  
  # Generate plots
  pdf(output_file, width = 8.5, height = 8, useDingbats = FALSE)
  par(mfrow = c(1,2),
      mar = c(4, 8, 4, 2),
      family = "Arial",
      mgp = c(2, 0.6, 0),
      las = 1,
      cex = 1.5)
  
  # Sensitivity forest plot
  do.call(forest, c(
    list(sens_data,
         leftlabs = c("Study", "TP", "N(+)"),
         main = "Sensitivity Analysis",
         col.diamond = "#2171b5",
         col.predict = "#9ecae1",
         col.random = "#2171b5",
         col.study = "#252525",
         col.square = "#2171b5",
         col.study.lines = "gray40"),
    forest_params
  ))
  
  # Specificity forest plot
  do.call(forest, c(
    list(spec_data,
         leftlabs = c("Study", "TN", "N(-)"),
         main = "Specificity Analysis",
         col.diamond = "#ef3b2c",
         col.predict = "#fcae91",
         col.random = "#ef3b2c",
         col.study = "#252525",
         col.square = "#ef3b2c",
         col.study.lines = "gray40"),
    forest_params
  ))
  
  dev.off()
  
  return(list(
    sensitivity = sens_data,
    specificity = spec_data
  ))
}

# Execute analysis
results <- generate_forest_plots(
  "meta_analysis_data.txt",
  "forest_plots.pdf"
)






# Required packages
library(ggplot2)
library(dplyr)
library(scales)
library(ggrepel)
library(showtext)
library(extrafont)

# Custom scientific notation formatter
scientific_format_custom <- function(x) {
  text <- gsub("1e\\+0?", "10^", format(x, scientific = TRUE))
  text <- gsub("1e-0?", "10^-", text)
  parse(text = text)
}

# Custom theme for consistent styling
custom_theme <- function() {
  theme_minimal() %+replace%
    theme(
      text = element_text(color = "#000000", size = 14, family = "sans"),
      plot.title = element_text(
        size = 19,
        face = "bold",
        hjust = 0.5,
        margin = margin(b = 20),
        color = "#000000"
      ),
      plot.subtitle = element_text(
        size = 18,
        hjust = 0.5,
        margin = margin(b = 10),
        color = "#000000"
      ),
      axis.title.y = element_text(
        size = 18,
        margin = margin(r = 10),
        angle = 90,
        face = "bold"
      ),
      axis.text.y = element_text(
        size = 16,
        face = "bold"
      ),
      axis.text.x = element_blank(),
      panel.grid.major.y = element_line(
        color = "#ecf0f1",
        linewidth = 0.2,
        linetype = "dashed"
      ),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(30, 80, 30, 30)
    )
}

# Generate threshold analysis visualizations
generate_threshold_plots <- function(data_file) {
  # Read and process data
  raw_data <- read.table(data_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  data <- raw_data %>%
    mutate(row_num = row_number()) %>%
    filter(Threshold_ngml != "N_A" & 
             Threshold_ngml != "3.546 log(Raman intensity)") %>%
    mutate(
      Study = paste0(Author, " ", Year, "_", row_num),
      Threshold = as.numeric(Threshold_ngml)
    )
  
  # Generate ordered point plot
  p1 <- ggplot(data, aes(x = reorder(Study, Threshold), y = Threshold)) +
    geom_point(size = 3, color = "blue", alpha = 0.6) +
    theme_minimal() +
    coord_flip() +
    labs(
      title = "Diagnostic Thresholds Across Studies",
      subtitle = paste0("(Based on ", nrow(data), " of ", nrow(raw_data), " studies)"),
      x = "Study",
      y = "Diagnostic Threshold (ng/mL)"
    ) +
    scale_y_log10(labels = scientific_format_custom) +
    theme(
      text = element_text(size = 14, face = "bold"),
      plot.title = element_text(size = 16, face = "bold"),
      axis.title = element_text(size = 14, face = "bold"),
      axis.text.y = element_text(size = 12, face = "bold"),
      axis.text.x = element_text(size = 14, face = "bold")
    )
  
  # Generate distribution plot
  stats <- list(
    median = median(data$Threshold),
    mean = mean(data$Threshold),
    max = max(data$Threshold),
    min = min(data$Threshold),
    q1 = quantile(data$Threshold, 0.25),
    q3 = quantile(data$Threshold, 0.75)
  )
  
  annotations <- data.frame(
    x = rep(1, 6),
    y = c(stats$max, stats$q3, stats$median, stats$q1, stats$min, stats$mean),
    label = c(
      sprintf("Maximum: %.2f", stats$max),
      sprintf("Q3: %.2f", stats$q3),
      sprintf("Median: %.2f", stats$median),
      sprintf("Q1: %.2f", stats$q1),
      sprintf("Minimum: %.2f", stats$min),
      sprintf("Mean: %.2f", stats$mean)
    ),
    color = c(rep("#000000", 5), "#e74c3c")
  )
  
  p2 <- ggplot(data, aes(x = "", y = Threshold)) +
    geom_jitter(width = 0.2, shape = 21, fill = "#3498db", 
                color = "#2980b9", alpha = 0.5, size = 3) +
    stat_boxplot(geom = 'errorbar', width = 0.2, color = "#2980b9") +
    geom_boxplot(fill = "#3498db", alpha = 0.6, color = "#2980b9",
                 outlier.shape = NA, linewidth = 0.8) +
    geom_hline(yintercept = stats$mean, linetype = "dashed",
               color = "#e74c3c", linewidth = 0.8, alpha = 0.8) +
    labs(
      title = "Distribution of Diagnostic Thresholds",
      subtitle = paste0("n = ", nrow(data), " of ", nrow(raw_data), " studies"),
      y = "Diagnostic Threshold (ng/mL)",
      x = NULL
    ) +
    scale_y_log10(labels = scientific_format_custom) +
    custom_theme() +
    geom_text_repel(
      data = annotations,
      aes(x = x, y = y, label = label, color = I(color)),
      size = 5,
      fontface = "bold",
      direction = "y",
      hjust = 0,
      segment.size = 0.4,
      segment.color = "#95a5a6",
      segment.alpha = 0.6,
      box.padding = 0.5,
      point.padding = 0.3,
      force = 2
    )
  
  # Save plots
  ggsave("threshold_ordered_plot.pdf", p1, width = 3.9, height = 4.2)
  ggsave("threshold_distribution.pdf", p2, width = 4.5, height = 5,
         device = cairo_pdf)
  
  return(list(ordered_plot = p1, distribution_plot = p2))
}

# Execute analysis
plots <- generate_threshold_plots("threshold_data.txt")







# Required packages
library(meta)
library(mada)
library(ggplot2)

analyze_threshold_effect <- function(data_file, output_file) {
  # Read and process data
  data <- read.table(data_file, header = TRUE, sep = "\t")
  
  # Apply continuity correction
  TP <- data$TP + 0.5
  FP <- data$FP + 0.5
  TN <- data$TN + 0.5
  FN <- data$FN + 0.5
  
  # Calculate sensitivity and specificity
  new_sens <- TP/(TP + FN)
  new_spec <- TN/(TN + FP)
  
  # Calculate logit transformations
  logit_se <- log(new_sens/(1-new_sens))
  logit_sp <- log(new_spec/(1-new_spec))
  
  # Prepare plot data
  plot_data <- data.frame(
    Author = data$Author,
    Year = data$Year,
    logit_se = logit_se, 
    logit_sp = logit_sp
  )
  
  # Calculate Spearman correlation
  cor_test <- cor.test(logit_se, logit_sp, method = "spearman")
  
  # Create visualization
  plot <- ggplot(plot_data, aes(x = logit_se, y = logit_sp)) +
    theme_bw() +
    geom_point(size = 3, shape = 21, fill = "#1B5E20", alpha = 0.6) +
    geom_smooth(method = "lm", col = "#EC8D63", se = TRUE, linetype = "dashed") +
    labs(
      title = "Diagnostic Threshold Effect Analysis",
      x = "Logit Sensitivity",
      y = "Logit Specificity",
      subtitle = paste("Studies:", nrow(plot_data))
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title = element_text(size = 13, face = "bold"),
      axis.text = element_text(size = 11, face = "bold"),
      legend.position = "none",
      panel.grid.major = element_line(colour = "grey90"),
      panel.grid.minor = element_line(colour = "grey95")
    )
  
  # Add correlation statistics
  plot <- plot +
    annotate(
      "text",
      x = min(logit_se),
      y = max(logit_sp),
      label = paste("Spearman correlation =", round(cor_test$estimate, 3),
                    "\np-value =", round(cor_test$p.value, 3)),
      hjust = 0,
      vjust = 1,
      size = 5,
      fontface = "bold"
    )
  
  # Save plot
  pdf(output_file, width = 4, height = 4)
  print(plot)
  dev.off()
  
  return(list(
    correlation = cor_test$estimate,
    p_value = cor_test$p.value,
    plot = plot
  ))
}

# Execute analysis
results <- analyze_threshold_effect(
  "threshold_effect_data.txt",
  "threshold_effect_analysis.pdf"
)









# Required packages
library(meta)
library(metafor)
library(dplyr)
library(ggplot2)
library(gridExtra)

# Function to perform meta-analysis for a specific subgroup
perform_subgroup_meta <- function(data, subgroup_var, subgroup_label) {
  metagen(
    TE = lnDOR,
    seTE = se_lnDOR,
    studlab = Author,
    data = data,
    sm = "ROM",
    fixed = FALSE,
    random = TRUE,
    method.tau = "DL",
    byvar = data[[subgroup_var]]
  )
}

# Function to create forest plot
create_forest_plot <- function(meta_obj, subgroup_var, subgroup_label) {
  forest(meta_obj,
         leftcols = c("studlab", subgroup_var),
         leftlabs = c("Study", subgroup_label),
         smlab = "log Diagnostic Odds Ratio",
         prediction = TRUE,
         print.tau2 = TRUE,
         print.I2 = TRUE)
}

# Function to analyze heterogeneity
analyze_heterogeneity <- function(meta_obj, group_name) {
  overall_stats <- data.frame(
    Group = "Overall",
    Subgroup = "All studies",
    N_Studies = meta_obj$k,
    Q = round(meta_obj$Q, 2),
    Q_df = meta_obj$df.Q,
    Q_pval = round(meta_obj$pval.Q, 4),
    I2 = round(meta_obj$I2, 1),
    Tau2 = round(meta_obj$tau2, 4),
    H = round(sqrt(meta_obj$Q/meta_obj$df.Q), 2),
    stringsAsFactors = FALSE
  )
  
  results_list <- list(overall_stats)
  
  if (!is.null(meta_obj$bylevs) && length(meta_obj$bylevs) > 0) {
    subgroup_stats <- data.frame(
      Group = group_name,
      Subgroup = meta_obj$bylevs,
      N_Studies = meta_obj$k.w,
      Q = round(meta_obj$Q.w, 2),
      Q_df = meta_obj$df.Q.w,
      Q_pval = round(meta_obj$pval.Q.w, 4),
      I2 = round(meta_obj$I2.w, 1),
      Tau2 = round(meta_obj$tau2.w, 4),
      H = round(sqrt(meta_obj$Q.w/meta_obj$df.Q.w), 2),
      stringsAsFactors = FALSE
    )
    results_list[[2]] <- subgroup_stats
    
    if (!is.null(meta_obj$Q.b) && !is.na(meta_obj$Q.b)) {
      between_group <- data.frame(
        Group = paste0(group_name, " (Between)"),
        Subgroup = "Between subgroups",
        N_Studies = NA,
        Q = round(meta_obj$Q.b, 2),
        Q_df = meta_obj$df.Q.b,
        Q_pval = round(meta_obj$pval.Q.b, 4),
        I2 = NA,
        Tau2 = NA,
        H = NA,
        stringsAsFactors = FALSE
      )
      results_list[[3]] <- between_group
    }
  }
  
  do.call(rbind, results_list)
}

# Main analysis function
perform_meta_analysis <- function(data_file) {
  # Read and process data
  data <- read.delim(data_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # Apply continuity correction
  data$TP_adj <- data$TP + 0.5
  data$FP_adj <- data$FP + 0.5
  data$TN_adj <- data$TN + 0.5
  data$FN_adj <- data$FN + 0.5
  
  # Calculate DOR and its logarithm
  data$DOR_adj <- (data$TP_adj * data$TN_adj)/(data$FP_adj * data$FN_adj)
  data$lnDOR <- log(data$DOR_adj)
  data$se_lnDOR <- sqrt(1/data$TP_adj + 1/data$FP_adj + 1/data$TN_adj + 1/data$FN_adj)
  
  # Define subgroup analyses
  subgroups <- list(
    list(var = "Cancer_System", label = "Cancer System"),
    list(var = "Study_Type_Group", label = "Study Type"),
    list(var = "Region_Group", label = "Region"),
    list(var = "Control_Group", label = "Control Group"),
    list(var = "Sample_Type", label = "Sample Type"),
    list(var = "Sample_Size_Group", label = "Sample Size Group")
  )
  
  # Perform analyses and create plots
  pdf("meta_analysis_results.pdf", width = 12, height = 20)
  
  meta_results <- lapply(subgroups, function(sg) {
    meta_obj <- perform_subgroup_meta(data, sg$var, sg$label)
    create_forest_plot(meta_obj, sg$var, sg$label)
    return(list(
      meta = meta_obj,
      heterogeneity = analyze_heterogeneity(meta_obj, sg$label)
    ))
  })
  
  dev.off()
  
  # Perform meta-regression
  meta_reg_year <- try({
    metareg(meta_results[[1]]$meta, formula = ~as.numeric(Year))
  })
  
  # Create meta-regression plot
  pdf("meta_regression_plot.pdf", width = 5.5, height = 4.5)
  
  ggplot(data, aes(x = as.numeric(Year), y = lnDOR)) +
    theme_bw() +
    geom_point(aes(size = 1/se_lnDOR), fill = "black", shape = 21, alpha = 0.6) +
    geom_smooth(method = "lm", color = "blue", linetype = "dashed", se = TRUE) +
    labs(
      title = "Publication Year vs ln(DOR)",
      subtitle = paste("Studies:", nrow(data)),
      x = "Publication Year",
      y = "ln(DOR)",
      size = "Precision"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 19, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.title = element_text(size = 16, face = "bold"),
      axis.text = element_text(size = 13, face = "bold"),
      legend.position = "right",
      legend.title = element_text(size = 15, face = "bold"),
      legend.text = element_text(size = 14),
      panel.grid.major = element_line(colour = "grey90"),
      panel.grid.minor = element_line(colour = "grey95")
    )
  
  dev.off()
  
  return(list(
    meta_results = meta_results,
    meta_regression = meta_reg_year
  ))
}

# Execute analysis
results <- perform_meta_analysis("meta_analysis_data.txt")










# Required packages
library(mada)
library(ggplot2)
library(meta)

generate_sroc_analysis <- function(data_file, output_file, correction = 0.01) {
  # Read and process data
  data <- read.delim(data_file, sep = "\t", header = TRUE)
  
  # Apply continuity correction
  data$TP <- data$TP + correction
  data$FP <- data$FP + correction
  data$TN <- data$TN + correction
  data$FN <- data$FN + correction
  
  # Prepare study data
  study_data <- data.frame(
    TP = data$TP,
    FP = data$FP,
    FN = data$FN,
    TN = data$TN
  )
  
  # Calculate sensitivity and specificity
  sens <- data$TP/(data$TP + data$FN)
  spec <- data$TN/(data$TN + data$FP)
  n <- data$TP + data$FP + data$FN + data$TN
  
  # Perform SROC analysis
  fit <- reitsma(study_data)
  
  # Process AUC meta-analysis
  ci_values <- strsplit(data$AUC95CI, "-")
  ci_lower <- as.numeric(sapply(ci_values, function(x) x[1]))
  ci_upper <- as.numeric(sapply(ci_values, function(x) x[2]))
  se <- (ci_upper - ci_lower) / (2 * 1.96)
  
  meta_auc <- metagen(
    TE = data$AUC_value,
    seTE = se,
    studlab = paste(data$Author, data$Year),
    sm = "AUC",
    fixed = FALSE,
    random = TRUE,
    method.tau = "DL"
  )
  
  # Extract combined AUC results
  auc <- meta_auc$TE.random
  ci_lower <- meta_auc$lower.random
  ci_upper <- meta_auc$upper.random
  
  # Generate visualization
  pdf(output_file, width = 4.5, height = 4.5)
  
  par(mar = c(5, 5, 4, 2), 
      mgp = c(2.5, 1, 0),
      cex.lab = 1.2,
      cex.axis = 1.1,
      cex.main = 1.5)
  
  # Create base plot
  plot(0, 0, type = "n",
       xlim = c(0, 1), 
       ylim = c(0, 1),
       main = "",
       xlab = "1-Specificity",
       ylab = "Sensitivity",
       font.lab = 2)
  
  mtext("SROC Curve Analysis", side = 3, line = 2, cex = 1.5, font = 2)
  mtext(sprintf("Studies: %d", nrow(data)), side = 3, line = 0.5, cex = 1.2, font = 2)
  
  # Add plot elements
  grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
  abline(0, 1, lty = 2, col = "gray50")
  
  ROCellipse(fit, 
             add = TRUE,
             pch = 19,
             cex = 1.5,
             col = "pink",
             lwd = 2)
  
  # Generate and add SROC curve
  fpr_points <- seq(0, 1, length.out = 100)
  sroc_points <- sroc(fit, fpr = fpr_points)
  lines(fpr_points, sroc_points[,2], col = "red", lwd = 2)
  
  # Add study points and labels
  points(1-spec, sens, pch = 19, col = "darkblue", cex = 1.5)
  
  auc_text <- sprintf("AUC = %.3f (95%% CI: %.3f-%.3f)", auc, ci_lower, ci_upper)
  text(0.0051, 0.0051, auc_text, pos = 4, cex = 1.1, font = 2, col = "darkred")
  
  text(1-spec, sens, 
       labels = n-correction*4,
       pos = 1, 
       cex = 0.8, 
       col = "darkgray")
  
  # Add legend
  legend("bottomright",
         legend = c("SROC Curve", "Study Estimate", "95% CI", "Random Line"),
         col = c("red", "darkblue", "pink", "gray50"),
         pch = c(NA, 19, NA, NA),
         lty = c(1, NA, 1, 2),
         lwd = c(2, NA, 2, 1),
         pt.cex = 1.5,
         bty = "n",
         cex = 1.1,
         inset = 0.08)
  
  dev.off()
  
  return(list(
    fit = fit,
    meta_auc = meta_auc,
    auc_summary = auc_text
  ))
}

# Execute analysis
results <- generate_sroc_analysis(
  "analysis_data.txt",
  "sroc_analysis.pdf"
)











