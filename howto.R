## Load always needed functions
source("mvci.R")
source("mvci_min_interval.R")
source("mvci_tools.R")
source("experiments/experiment_tools.R")

## Load toy data generation function
source("experiments/mvci_data_loaders.R")


## Create data to play with
dmat <- make.toy.data2(Nbase = 100, Nbump = 99, M = 80)

## Compute confidence band
alpha <- 0.1
k <- floor(alpha * nrow(dmat))
cb_l0 <- findcb_topdown(dmat, K = k, L = 0)
cb_naive <-  cb_quantile(dmat, k)

## Compute a calibrated condfidence band

# confidence band computation function
mycb <- function(dmat, k){
  cb <- findcb_topdown(dmat, K = k, L = 0)
  matrix(c(dmat[cb$downmask], dmat[cb$upmask]), nrow = 2, byrow = T)
}

#mycb <- function(dmat, k){
#  cb <- cb_quantile(dmat, k)
#  cb$cb
#}

# function to test if a row is within band
mycontain <- cb.violations

# compute profile
k.max <- floor(alpha * nrow(dmat))

# generic version - very slow for MWE
fp1 <- fold.profile(mycb, mycontain, dmat, B = 2, k.max = k.max)
pd1 <- data.frame(k = fp1$profile.k, coverage = fp1$profile)
pd1$id <- 'generic'

# MWE specific version - much faster
fp2 <- mwe.fold.profile(dmat, B = 2, k.max = k.max)
pd2 <- data.frame(k = fp2$profile.k, coverage = fp2$profile)
pd2$id <- 'MWE'

pd <- rbind(pd1, pd2)
p <- qplot(x = k, y = coverage, data = pd, geom = "point", xlim = c(-1, k.max),
           group = id, color = id)
p


## Extract band in various formats
cbm0 <- cbr_extract_cb(dmat, cb_l0)
dimnames(cbm0) <- list(c('up','down'), sprintf("%d",1:ncol(cbm0)) )
pd0 <- make.ggplot.df(cbm0)
head(pd0)

## Visualize
plot.data.ci(dmat, cb_l0$Z, cb_l0$upmask, cb_l0$downmask)

plot.data.cb(dmat, cb_l0)

p <- ggplot(pd0, aes(x = variable, y = value, group = row))
p <- p + geom_line()
p <- p + geom_point()
p
