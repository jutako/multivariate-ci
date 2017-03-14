# Load always needed functions
source("mvci.R")
source("mvci_min_interval.R")
source("mvci_tools.R")
source("experiments/experiment_tools.R")

# Load toy data generation function
source("experiments/mvci_data_loaders.R")


# Create data to play with
dmat <- make.toy.data2(Nbase = 100, Nbump = 99, M = 80)

# Compute confidence band
alpha <- 0.1
k <- floor(alpha * nrow(dmat))
cb_l0 <- findcb_topdown(dmat, K = k, L = 0)

# Extract band in various formats
cbm0 <- cbr_extract_cb(dmat, cb_l0)
dimnames(cbm0) <- list(c('up','down'), sprintf("%d",1:ncol(cbm0)) )
pd0 <- make.ggplot.df(cbm0)
head(pd0)

# Visualize
plot.data.ci(dmat, cb_l0$Z, cb_l0$upmask, cb_l0$downmask)

plot.data.cb(dmat, cb_l0)

p <- ggplot(pd0, aes(x = variable, y = value, group = row))
p <- p + geom_line()
p <- p + geom_point()
p
