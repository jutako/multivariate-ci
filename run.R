
# If packages are missing run this to install them:
# install.packages(c('grid','ggplot2','zoo','signal','RColorBrewer'))

require(ggplot2)
require(grid)
require(zoo)
require(signal)

start_time <- Sys.time()

# Set data locations
dirlst <- list(pubfig.dir = './figs',
               expres.dir = './rds')
filelst <- list(stockdata.file = './rds/stockdata/stockdata.rds')

# Recompute confidence bands
update.res = T

# Do not recompute confidence bands
#update.res = F

OPTS <- list(dir = dirlst,
             file = filelst,
             update.res = update.res)


if (!dir.exists(OPTS$dir$pubfig.dir)){ dir.create(OPTS$dir$pubfig.dir) }


# Load some functions to namespace
source("mvci.R")
source("mvci_min_interval.R")
source("mvci_tools.R")
source("experiments/mvci_data_loaders.R")
source("experiments/experiment_tools.R")

# Plot motivating example
source("experiments/motivating_example.R")

# Plot stock data example
source("experiments/stock_example.R")

as.numeric(difftime(Sys.time(), start_time, units = "min"))
