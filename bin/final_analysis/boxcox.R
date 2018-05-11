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

# Box-Cox function
bc <- function(x) {
  lmbd <- BoxCox.lambda(x, method='loglik')
  return(as.vector(BoxCox(x, lmbd)))
}


# Load libraries
print('Loading R libraries...')
pkgLoad('tidyverse')
pkgLoad('forecast')

# Get command line argument
args <- commandArgs(TRUE)

# Set working directory
setwd(dirname(args[1]))


# Load data
df <- read.csv(args[1], row.names = 1, check.names = F)
df <- df + 1
       
df_out <- as.data.frame(t(apply(df, 1, bc)))
names(df_out) <- names(df)
write.csv(df_out, 'RPMMM_mirs_over_50_boxcox.csv', quote = F)
