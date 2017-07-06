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
invisible(pkgLoad('ggplot2'))

# Get command line argument
args <- commandArgs(TRUE)

# Set working directory
setwd(dirname(args[1]))

# Create blank pca.csv file for assembly of excel output
invisible(file.create('PCA.csv'))

# Load your normalized RPMMM file for your samples
data <- read.csv(args[1], header = T, row.names = 1, check.names = F)

# Add small amount for log step
data <- data + .01

# Remove rows with no variance
data <- data[rvar(data) != 0,]

# Get pca of transposed and logged dataframe
pca <- prcomp(t(log(data)),
              center = TRUE,
              scale. = TRUE)

# Get the variance explained by each PC
li <- pca$sdev^2 / sum(pca$sdev^2)
pc1v <- round(li[1] * 100, 1)
pc2v <- round(li[2] * 100, 1)

# Reformate data for plotting
pca_s <- as.data.frame(pca$x)
pca_s <- pca_s[,c(1,2)]

# Plot data
if (length(args) == 2) {
  cond <- read.csv(args[2], row.names = 1, header = T)
  pca_s <- merge(pca_s, cond, by = 0)
  names(pca_s)[1] <- 'name'
  
  png('PCA.png', units = 'in', width = 12, height = 8, res = 250)
  print(
  ggplot(pca_s, aes(PC1, PC2)) +
    geom_point(aes(color = Condition), size=4.5) +
    geom_text(aes(label=name)) +
    theme_bw() +
    xlab(paste0("PC1: ",pc1v,"% variance")) +
    ylab(paste0("PC2: ",pc2v,"% variance")) +
    ggtitle("PCA Analysis") + 
    theme(text = element_text(size = 26),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0),
          legend.position = 'left',
          legend.text.align = 0,
          legend.key.size = unit(2, 'lines'))
  )
  invisible(dev.off())
} else {
  pca_s$name = row.names(pca_s)

  png('PCA.png', units = 'in', width = 12, height = 8, res = 250)
  print(
  ggplot(pca_s, aes(PC1, PC2)) +
    geom_point(size=4.5, color = 'grey50') +
    geom_text(aes(label=name)) +
    theme_bw() +
    xlab(paste0("PC1: ",pc1v,"% variance")) +
    ylab(paste0("PC2: ",pc2v,"% variance")) +
    ggtitle("PCA Analysis") + 
    theme(text = element_text(size = 26),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0),
          legend.position = 'left',
          legend.text.align = 0,
          legend.key.size = unit(2, 'lines'))
  )
  invisible(dev.off())
}
