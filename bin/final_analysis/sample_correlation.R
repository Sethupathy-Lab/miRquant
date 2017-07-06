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

# Function to calculate row variance
rvar <- function (x,na.rm = TRUE) 
{
  sqr = function(x) x * x
  n = rowSums(!is.na(x))
  n[n <= 1] = NA
  return(rowSums(sqr(x - rowMeans(x,na.rm = na.rm)), na.rm = na.rm)/(n - 1))
}

# Load libraries
invisible(pkgLoad('RColorBrewer'))
invisible(pkgLoad('pheatmap'))

# Get command line argument
args <- commandArgs(TRUE)

# Set working directory
setwd(dirname(args[1]))

# Load your normalized RPMMM file for your samples
data <- read.csv(args[1], header = T, row.names = 1, check.names = F)

print("Dimensions (# of miRs over 50 RPMMM, # of samples):")
dim(data)

# Subset to top 50 most variable miRs (if 50 exist)
n_var_mirs <- nrow(data)

if (nrow(data) > 50) {
  data$var <- rvar(data)
  data <- data[order(-data$var)[1:50],]
  data$var <- NULL
  n_var_mirs <- 50
}

# Get correlation values
data <- data.matrix(data)

# get the correlation values for the data
data_corr <- cor(as.matrix(data))

# write correlation values to file
write.table(data_corr, file = "sample_correlation_values.csv", sep = ",", quote = F)

# make heat map using the correlation values
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

if (length(args) == 2) {
  treat_ann <- read.csv(args[2], row.names = 1, header = T)
  names(treat_ann) <- 'Condition'

  png("sample_correlation_values.png", units = 'in', width = 8, height = 8, res = 150)
  pheatmap(data_corr,
           col=colors,
           annotation = treat_ann,
           main=paste("Sample to Sample Correlations\n(Based on top ", n_var_mirs, "most variable miRs)"))
  invisible(dev.off())
} else {
  png("sample_correlation_values.png", units = 'in', width = 8, height = 8, res = 150)
  pheatmap(data_corr,
           col=colors,
           main=paste("Sample to Sample Correlations\n(Based on top ", n_var_mirs, "most variable miRs)"))
  invisible(dev.off())
}
