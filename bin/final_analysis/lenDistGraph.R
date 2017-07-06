# Load libraries
library(ggplot2)
library(reshape2)

# Get command line argument
args <- commandArgs(TRUE)

# Set working directory
setwd(dirname(args[1]))

# Load Data
lenData <- read.csv(args[1], header = T, sep = ",")

# Change names to first 10 characters
new_names = c()
for (name in names(lenData)) {
  if (nchar(name) > 20) {
    new_names <- append(new_names, substr(name, 1, 20))
  } else {
    new_names <- append(new_names, name)
  }
}

# Confirm shortened names are unique before using
if (length(new_names) == length(unique(new_names))) {
  names(lenData) <- new_names
}

# Convert data for ggplot input
bit <- melt(lenData, id.vars = "Length")

# Determine number of samples per role in output (based on total # of samples)
wrap = round(length(lenData) * .33)
height_multi = ceiling(length(lenData) / 5)

# Write output to lenDistHistogram.png
png("length_distribution.png", width = 10, height = 2 * height_multi, units = "in", res = 150)
ggplot(bit, aes(x = Length, y = value)) +
  geom_bar(stat="identity", width = 1) +
  facet_wrap(~variable, ncol = 5) +
  ggtitle('Read Length Distribution') +
  scale_y_continuous('Proportion of reads',
                     breaks = seq(0, max(bit$value), 5),
                     limits = c(0, max(bit$value) * 1.1),
                     expand = c(0,0)) +
  scale_x_continuous('Read length') +
  theme_bw()
invisible(dev.off())
