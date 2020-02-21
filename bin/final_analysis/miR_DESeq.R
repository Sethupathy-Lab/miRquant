### Load functions
# Function to check if packages exist
pkgLoad <- function(x)
{
  if (!require(x,character.only = TRUE, quietly = TRUE))
  {
    stop(paste0('Package ', x, ' not found, please install R package and run again'))
  }
  #now load library and suppress warnings
  suppressPackageStartupMessages(x)
}

# Load packages
suppressMessages(pkgLoad('DESeq2'))
suppressMessages(pkgLoad('tidyverse'))

# Get command line argument
args <- commandArgs(TRUE)

# Set working directory
setwd(dirname(args[1]))

# Create directory for output results
dir.create(file.path(getwd(), 'DESeq_output'), showWarnings = FALSE)

# Load raw miRNA counts
counts <- read.csv(args[1], header = T, sep = ",", check.names = F) %>%
  select(miR = 1, everything()) %>%
  mutate_if(is.numeric, funs(round))

row.names(counts) <- counts$miR
counts$miR <- NULL

# Load conditions file
conditions <- read.csv(args[2], header = T, sep = ",")
rownames(conditions) <- conditions$Sample

# Remove samples to not be included in final analysis
counts <- counts[,rownames(conditions)]

cat('\nConfirm sample names in raw counts and conditions file are the same.\n')
if (all(rownames(conditions) == colnames(counts)) == F) {
  stop('Samples in conditions not the same as in raw counts, exiting...')
} else{
  cat("Confirmed!\n\n")
}


# Make DESeq object and calculate adjusted counts based on size factors
cat('Executing DESeq...\n')
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

volcanoPlot <- function(res_vol, pv_cut = .05) {
  
  res_vol$threshold = as.factor(ifelse(res_vol$padj <= pv_cut & res_vol$log2FoldChange < 0, -1,
                                       ifelse(res_vol$padj <= pv_cut & res_vol$log2FoldChange > 0, 1, 0)))
  
  res_vol$padj <- -log10(res_vol$padj)
  
  ##Construct the plot object
  g = ggplot(data=as.data.frame(res_vol), aes(x=log2FoldChange, y=padj, colour=threshold, label = name)) +
    geom_point(alpha=0.4, size=2.75) +
    scale_colour_manual(values = c("blue", "gray", "red")) +
    geom_hline(yintercept = -log10(pv_cut), linetype = 'dashed') +
    annotate('text', label = 'FDR', x = -8, y = -log10(pv_cut) + .6, vjust = 0, hjust = 0) +
    scale_x_continuous("log2 fold change",
                       limit = c(-8, 8),
                       expand = c(0.025,0)) +
    scale_y_continuous("-log10 p-value",
                       limit = c(0, max(res_vol$padj) * 1.1),
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

comparisons <- read.csv(args[3], header = F, sep = ",")

for (row in c(1:nrow(comparisons))) {

  cond1 <- as.character(comparisons[row,1])
  cond2 <- as.character(comparisons[row,2])
  
  cat(paste0('\nRunning differential expression for ', cond1, ' VS ', cond2, '\n'))
  
  if (!(cond1 %in% conditions$Condition) | !(cond2 %in% conditions$Condition)) {
    cat(paste0('Either ', cond1, ' or ', cond2, ' not a condition in conditions file, skipping comparison.\n\n'))
    next()
  }
  
  res <- results( dds, contrast = c("Condition", cond1, cond2))
  
  ## Filter for a baseMean greater than 5
  res.5<-res[res$baseMean>5, ]
  
  ## Adjust p-value according to Benjamini & Hochberg method (need to do this since we filtered out by base mean 5 above)
  res.5$padj <- p.adjust(res.5$pvalue, method="BH")
  
  res.5 <- res.5[!is.na(res.5$pvalue),]

  ## Get average normalized counts for each condition and add to file
  cond1_s <- conditions %>%
  filter(Condition == cond1) %>%
  pull(Sample) %>%
  as.vector()

  cond2_s <- conditions %>%
  filter(Condition == cond2) %>%
  pull(Sample) %>%
  as.vector()

  avg_df <- as.data.frame(normalized.counts) %>%
  mutate(miR = rownames(.),
         cond1_avg = rowMeans(.[,cond1_s]),
         cond2_avg = rowMeans(.[,cond2_s]))

  avg_df <- avg_df[,c('miR', 'cond1_avg', 'cond2_avg')]

  names(avg_df) <- c('miR', paste0('avg_', cond1), paste0('avg_', cond2))

  ## Write res.cont DESeq data to output file
  res.5 <- res.5[order(res.5$pvalue),]
  res.5 <- res.5 %>%
  as.data.frame() %>%
  mutate(miR = rownames(.)) %>%
  select(miR, everything()) %>%
  left_join(avg_df)

  write_csv(res.5, paste0("DESeq_output/DESeq_", cond1, "vs", cond2, ".csv"))
  
  res.5$name <- rownames(res.5)
  
  ## Volcano plot
#  png(paste0("DESeq_output/DESeq_", cond1, "vs", cond2, "_PVAL.05_VolcanoPlot.png"), units = 'in', width = 6, height = 6, res = 250)
#  g <- volcanoPlot(res.5)
#  print(g)
#  dev.off()
#  print(g)
}
