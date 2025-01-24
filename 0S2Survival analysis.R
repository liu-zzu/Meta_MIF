# Load required packages
library(pacman)
p_load(
  data.table,
  magrittr,
  tidyverse,
  survival,
  survminer,
  limma
)

# Clear workspace
rm(list=ls())

# Load and preprocess survival data
rt <- as.data.frame(fread("all_survival3.txt", header=TRUE, sep="\t", check.names=FALSE))
rt$futime <- rt$futime / 30  # Convert time to months

# Determine optimal cutpoint
sur.cut <- surv_cutpoint(
  rt,
  time="futime",
  event="fustat",
  minprop=0.3,
  variables="concentrate2"
)

# Categorize based on cutpoint
sur.cat <- surv_categorize(sur.cut)

# Fit survival model
fit <- survfit(Surv(futime, fustat) ~ sur.cat[,3], data=sur.cat)

# Create survival plot
pdf(file="survival_analysis.pdf", width=4, height=5)
ggsurvplot(
  fit,
  data=sur.cat,
  conf.int=FALSE,
  pval=TRUE,
  palette="jco",
  legend.title="Group",
  legend.labs=c("MIF_high", "MIF_low"),
  risk.table=TRUE,
  xlab="Time (months)",
  font.main=c(15, "bold"),
  font.x=c(13, "bold"),
  font.y=c(13, "bold"),
  font.tickslab=c(12, "bold"),
  font.legend=c(13, "bold"),
  risk.table.fontsize=5,
  break.time.by=10,
  xlim=c(0, 40),
  tables.y.text=FALSE
)
dev.off()





# Load required packages
library(pacman)
p_load(
  forestplot,
  data.table,
  magrittr,
  tidyverse
)

# Clear workspace
rm(list=ls())

# Create data frame for meta-analysis
data <- data.frame(
  Author = c('Y. Mohri et al', 'Theresa H Wirtz et al', 'Yi-Ming Zhao et al', 
             'Mira Lanki et al', 'Azaz Ahmed et al', 'Faruk Tas et al', 'Yuhuan Zheng et al'),
  Cancer_Type = c('Stomach Cancer', 'Liver Cancer', 'Liver Cancer', 
                  'Pancreatic Cancer', 'Pancreatic Cancer', 'Ovarian Cancer', 'Multiple Myeloma'),
  HR = c(2.84, 1.957, 1.754, 1.743, 1.908, NA, NA),
  CI_Lower = c(1.27, 1.268, 1.103, 1.001, 1.084, NA, NA),
  CI_Upper = c(6.68, 3.022, 5.432, 3.036, 3.36, NA, NA),
  P_value = c(0.03, 0.002, 0.012, 0.05, 0.04, 0.01, 0.02),
  Sample_Size = c(105, 50, 442, 173, 48, 50, 173)
)

# Format table text for forest plot
tabletext <- cbind(
  c("Author", data$Author),
  c("Cancer Type", data$Cancer_Type),
  c("HR (95% CI)", ifelse(is.na(data$HR), 
                          "Data Missing", 
                          paste0(data$HR, " (", data$CI_Lower, "-", data$CI_Upper, ")"))),
  c("P-value", data$P_value),
  c("Sample Size", data$Sample_Size)
)

# Generate forest plot
pdf("meta_analysis_forest_plot.pdf", width=10, height=5)
forestplot(
  tabletext,
  mean = c(NA, data$HR),
  lower = c(NA, data$CI_Lower),
  upper = c(NA, data$CI_Upper),
  is.summary = c(TRUE, rep(FALSE, nrow(data))),
  xlab = "Hazard Ratio (HR)",
  col = fpColors(box="royalblue", lines="darkblue", summary="royalblue"),
  txt_gp = fpTxtGp(
    label = gpar(fontfamily="", cex=1.0, fontface="bold"),
    ticks = gpar(cex=0.8, fontface="bold"),
    xlab = gpar(cex=0.9, fontface="bold")
  ),
  boxsize = 0.3,
  ci.vertices = TRUE,
  ci.vertices.height = 0.1,
  zero = 1,
  graph.pos = 3,
  align = c("l", "l", "l", "l", "l")
)
dev.off()


















