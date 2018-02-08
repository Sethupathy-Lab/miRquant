### Load functions
# Function to check if packages exist
pkgLoad <- function(x)
{
  if (!require(x,character.only = TRUE))
  {
    stop(paste0('Package ', x, ' not found, please install R package and run again'))
  }
  #now load library and suppress warnings
  suppressPackageStartupMessages(x)
}

# Load libraries
print('Loading R libraries...')
pkgLoad('tidyverse')

# Get command line argument
args <- commandArgs(TRUE)

# Set working directory
setwd(dirname(args[1]))

# Load your mapping statistics file for your samples
map_stats <- read.csv(args[1]) %>%
  select(Sample_name, contains('Percent'))

# Filter for only the percent mapped columns
map_stats <- read.csv('Mapping_Statistics.csv') %>%
  select(Sample_name, contains('Percent'))

# Create mapping statistics image
png('Mapping_Statistics.png', units = 'in', width = 8, height = 6, res = 250)
map_stats %>%
  gather(category, percent, -Sample_name) %>%
  mutate(category = factor(category, levels = c("Percent.Trimmed.Reads",
                                                "Percent.Too.Short",
                                                "Percent.Exact.Matches",
                                                "Percent.Mismatched",
                                                "Percent.Mapped",
                                                "Percent.miR.Mapped",
                                                "Percent.tRNA.Mapped",
                                                "Percent.yRNA.Mapped"))) %>%
  ggplot(., aes(category, percent)) +
  geom_boxplot() +
  scale_y_continuous('Percent',
                     limits = c(0,100),
                     breaks = c(0,25,50,75,100)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(color = 'black'),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5, color = 'black'),
        axis.title.x = element_blank())
invisible(dev.off())
