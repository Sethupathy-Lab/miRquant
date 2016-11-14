# Get command line argument
args <- commandArgs(TRUE)

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
#library(gplots)
library(pheatmap)

colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
colnames(data_corr) <- NULL

png("sample_correlation_heatmap.png")
pheatmap(data_corr,
         col=colors,
         main="Sample to Sample Distance")
dev.off()
