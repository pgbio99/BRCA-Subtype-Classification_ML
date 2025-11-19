setwd("E:/Bioteck_Project/volcano_plot_degs")
getwd()

install.packages("ggplot2")
library(ggplot2)

data <- read.csv("sign.genes.csv", header = TRUE)
colnames(data)

logFC.threshold <- 1
padj.threshold <- 0.05

data$Regulation <- "Non-regulated"
data$Regulation[data$logFC > logFC.threshold & data$adj.P.Val < padj.threshold] <- "Upregulated"
data$Regulation[data$logFC < -logFC.threshold & data$adj.P.Val < padj.threshold] <- "Downregulated"

ggplot(data, aes(x = logFC, y = -log10(adj.P.Val), color=Regulation))+
  geom_point(alpha=0.8) +
  scale_color_manual(values = c("Upregulated"="red", "Downregulated"="blue", "Non-regulated"="yellow")) + 
  theme_minimal() +
  labs(titile="Volvano Plot",
       x="log2(Fold Change)",
       y="-log10(adj.P.Val)",
       color="Regulation") +
  geom_vline(xintercept = c(-logFC.threshold, logFC.threshold), linetype="dashed", color="black") +
  geom_hline(yintercept = -log10(padj.threshold), linetype="dashed", color="black")