# Get command line argument
args <- commandArgs(TRUE)

# Set working directory
setwd(dirname(args[1]))

# Load your normalized RRM file for your samples
data <- read.csv(args[1], header = T, row.names = 1)

print("Dimensions (# of miRs, # of samples):")
dim(data)

# Get correlation values
data <- data.matrix(data)

# get the correlation values for the data
data_corr <- cor(as.matrix(data))

# write correlation values to file
write.table(data_corr, file = "sample_correlation_values.csv", sep = ",", quote = F)

# make heat map using the correlation values
library(RColorBrewer)
library(pheatmap)

colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

if (length(args) == 2) {
  treat_ann <- read.table(args[2], row.names = 1, header = T)
  names(treat_ann) <- 'Condition'

  png("sample_correlation_values.png")
  pheatmap(data_corr,
           col=colors,
           annotation = treat_ann,
           main="Sample to Sample Correlations")
  invisible(dev.off())
} else {
  png("sample_correlation_values.png")
  pheatmap(data_corr,
           col=colors,
           main="Sample to Sample Correlations")
  invisible(dev.off())
}
