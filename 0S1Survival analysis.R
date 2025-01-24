rm(list=ls())

# Load and preprocess data
data2 <- read.table('Coxb_log2_.txt', header=TRUE, sep='\t')
data2$Coxb_HR = log2(data2$Coxb_HR+1)

# Process confidence intervals
nameCI <- as.character(data2[,3])
a <- strsplit(nameCI, '-')
data2$'Lower' <- as.numeric(sapply(a, function(x) x[1]))
data2$'Upper' <- as.numeric(sapply(a, function(x) x[2]))

# Set annotation positions
data2$'annotation' <- rep(6.5, 36)
data2$'ltext' <- rep(-0.5, 36)
as.factor(data2$CA_name)

# Load required libraries
library(ggbreak)
library(ggplot2)
library(RColorBrewer)

# Set up color palette
colourCount = 36
getPalette = colorRampPalette(brewer.pal(8, "Set1"))

# Create forest plot
p <- ggplot(data2, aes(data2$Coxb_HR, data2$ID))
p + geom_point(size=0.8, color=getPalette(36)) +
  scale_fill_brewer() +
  geom_errorbarh(aes(xmax=Upper, xmin=Lower), height=0.5, color=getPalette(36), size=0.4) +
  scale_x_continuous(limits=c(-1.5, 7.5), breaks=seq(-1, 7.5, 2)) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour="black", size=0.1),
    axis.text.x = element_text(size=6),
    axis.text.y = element_text(size=6),
    axis.title.x = element_text(size=8, face="bold")
  ) +
  geom_vline(xintercept=1, lwd=0.3, col='black') +
  xlab('Hazard Ratio') + 
  ylab(' ') +
  geom_text(x=data2$annotation, y=data2$ID, label=data2[,4], size=2, col='black') +
  geom_text(x=6.5, y=37, label='P value', size=2.5, col='black') +
  geom_text(x=data2$ltext, y=data2$ID, label=data2$Coxb_CI, size=2, col='black') +
  geom_text(x=-0.1, y=37, label='95% CI', size=2.5, col='black') +
  scale_y_discrete(expand=expansion(mult=0, 1.5))

# Save plot
ggsave("forest_plot.pdf", width=4.5, height=3.5)





# Load required libraries
library(ggbreak)
library(ggplot2)
library(RColorBrewer)

# Load and preprocess data
data2 <- read.table('Coxb_log2_.txt', header=TRUE, sep='\t')
data2$Coxb_HR = log2(data2$Coxb_HR+1)

# Process confidence intervals
nameCI <- as.character(data2[,3])
a <- strsplit(nameCI, '-')
data2$'Lower' <- as.numeric(sapply(a, function(x) x[1]))
data2$'Upper' <- as.numeric(sapply(a, function(x) x[2]))

# Set annotation positions
data2$'annotation' <- rep(6.5, 33)
data2$'ltext' <- rep(-0.5, 33)
as.factor(data2$CA_name)

# Set up color palette
colourCount = 33
getPalette = colorRampPalette(brewer.pal(8, "Set1"))

# Create forest plot
p <- ggplot(data2, aes(data2$Coxb_HR, data2$ID))
p + geom_point(size=0.8, color=getPalette(33)) +
  scale_fill_brewer() +
  geom_errorbarh(aes(xmax=Upper, xmin=Lower), height=0.5, color=getPalette(33), size=0.4) +
  scale_x_continuous(limits=c(-1.5, 7.5), breaks=seq(-1, 7.5, 2)) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour="black", size=0.1),
    axis.text.x = element_text(size=6),
    axis.text.y = element_text(size=6),
    axis.title.x = element_text(size=8, face="bold")
  ) +
  geom_vline(xintercept=1, lwd=0.3, col='black') +
  xlab('Hazard Ratio') + 
  ylab(' ') +
  # Add annotations for P-values and CI
  geom_text(x=data2$annotation, y=data2$ID, label=data2[,4], size=2, col='black') +
  geom_text(x=6.5, y=34, label='P value', size=2.5, col='black') +
  geom_text(x=data2$ltext, y=data2$ID, label=data2$Coxb_CI, size=2, col='black') +
  geom_text(x=-0.5, y=34, label='95% CI', size=2.5, col='black') +
  scale_y_discrete(expand=expansion(mult=0, 1.5))

# Save plot
ggsave("forest_plot.pdf", width=4.5, height=3.5)



