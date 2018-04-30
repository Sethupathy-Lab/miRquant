library(DESeq2)
library(tidyverse)

# Get command line argument
args <- commandArgs(TRUE)

# Set working directory
setwd(dirname(args[1]))

# Create directory for output results
dir.create(file.path(getwd(), 'DESeq_output'), showWarnings = FALSE)

# Load Data
counts <- read.csv(args[1], header = T, sep = ",", check.names = F) %>%
  select(miR = 1, everything()) %>%
  mutate_if(is.numeric, funs(round))

row.names(counts) <- counts$miR
counts$miR <- NULL

conditions <- read.csv(args[2], header = T, sep = ",")
rownames(conditions) <- conditions$Sample

print('check to make sure names matching up')
print(all(rownames(si) == colnames(rc)))

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = conditions,
                              design = ~ Condition)


dds <- DESeq(dds)

normalized.counts <- as.data.frame(counts(dds, normalized=TRUE ))
write.csv(normalized.counts, 'DESeq_output/DESeq_normalized_miR_counts.csv', row.names = T)


###################################################################################################
#
#  Getting fold changes from direct comparisons with control
#
###################################################################################################
##Highlight genes that have an absolute fold change > 2 and a p-value < Bonferroni cut-off
volcanoPlot <- function(res_vol, fc_cut = 0, pv_cut = .05) {
  res_vol$threshold = as.factor(ifelse(res_vol$padj <= pv_cut  & res_vol$log2FoldChange < -fc_cut, -1,
                                       ifelse(res_vol$padj <= pv_cut  & res_vol$log2FoldChange > fc_cut, 1, 0)))
  
  ##Construct the plot object
  g = ggplot(data=as.data.frame(res_vol), aes(x=log2FoldChange, y=-log10(pvalue), colour=threshold, label = name)) +
    geom_point(alpha=0.4, size=1.75) +
    scale_colour_manual(values = c("blue", "gray", "red")) +
    geom_vline(xintercept = c(-fc_cut, fc_cut), linetype = 'dashed') +
    geom_hline(yintercept = -log10(pv_cut), linetype = 'dashed') +
    annotate('text', label = 'FDR', x = -8, y = -log10(pv_cut) + .6, vjust = 0, hjust = 0) +
    scale_x_continuous("log2 fold change",
                       limit = c(-8, 8),
                       expand = c(0.025,0)) +
    scale_y_continuous("-log10 p-value",
                       limit = c(0, 50),
                       expand = c(0, 0)) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_rect(colour = "black", fill = "white"), 
          legend.position = "none")
  return(g)
}

###################################################################################################
#
#  Getting fold changes from direct comparisons with control
#
###################################################################################################

comparisons <- read.csv(args[2], header = T, sep = ",")
print(comparisons)

res <- results( dds, contrast = c("Condition", "Active", "Quiescent-2") )
head(res)

## Add gene name to res file
res.5<-res[res$baseMean>5, ]
head(res.5)

## Adjust p-value according to Benjamini & Hochberg method (need to do this since we filtered out by base mean 5 above)
res.5$padj <- p.adjust(res.5$pvalue, method="BH")

res.5 <- res.5[!is.na(res.5$pvalue),]

## Write res.cont DESeq data to output file
write.csv(res.5, file="DESeq_output/DESeq_activeVSquiescent-2.csv", quote=F)

res.5$name <- rownames(res.5)
# Volcano plot
png('DESeq_output/FLCHCCvsAHEP_PVAL.05_VolcanoPlot.png', units = 'in', width = 8, height = 8, res = 250)
g <- volcanoPlot(res.5, 1.5)
g
dev.off()

res.5 %>% as.data.frame() %>% filter(padj < .05) %>% select(name, baseMean, padj)
