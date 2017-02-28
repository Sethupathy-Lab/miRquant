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
names(lenData) <- new_names

# Convert data for ggplot input
bit <- melt(lenData, id.vars = "X")

# Determine number of samples per role in output (based on total # of samples)
wrap = round(length(lenData) * .33)

# Write output to lenDistHistogram.png
png("length_distribution.png", width = 10, height = 10, units = "in", res = 150)
ggplot(bit, aes(x = X, y = value, fill = variable)) +
  geom_bar(stat="identity") +
  facet_wrap(~variable, ncol = wrap) +
  scale_fill_discrete(guide=FALSE) +
  labs(list(title="Length Distribution", x="Read length", y="Ratio of total reads")) 
invisible(dev.off())
