
## Set some options
k <- 40
l <- 125
blw <- 5
update.cb <- OPTS$update.res #T
res.save.path <- file.path(OPTS$dir$expres.dir, 'stockdata')
if(!dir.exists(res.save.path)){ dir.create(res.save.path) }

fig.save.path <- OPTS$dir$pubfig.dir 


stockdl <- load_data_stock(rdsfile = OPTS$file$stockdata.file,
                      type = 'raw')
myfun <- function(x){x - mean(x[1:blw])}
stockd <- t(apply(stockdl$data, 1, myfun))

cbfile <- file.path(res.save.path,
              sprintf('stockdata_blw%d_K%d_L%d.rds', blw, k, l))
if (update.cb){
  cb <-  findcb_topdown(stockd, K = k, L = l)
  saveRDS(cb, cbfile)
} else {
  cb <- readRDS(cbfile)
}

plot.stock.data <- function(stockdl, stockd, cb){
  N <- nrow(stockd)
  
  xlim <- c(stockdl$time[1], stockdl$time[ncol(stockd)])
  ylim <- c(-25, 100)
  stockdl$time = stockdl$time
  
  xlab = ''
  ylab = 'Price (USD)'
  title = ''
  yintercept = 0
  xaxt = 'n'
  samp.frac = 0.5
  colors <- RColorBrewer::brewer.pal(3, 'Dark2')
  
  stockd.out <- stockd[setdiff(1:N, cb$row.inc.idx),]
  stockd.in <- stockd[cb$row.inc.idx, ]
  
  # Crate plot
  plot(stockdl$time, colMeans(stockd), lwd=0.75, col='white',
       ylab = ylab, xlab = xlab, main = title,
       xlim = xlim,
       ylim = ylim,
       xaxt = xaxt)
  
  # rows inside bands
  # subsetting rows and columns to make smaller pdf
  set.seed(42)
  idx <- sample(1:nrow(stockd.in), floor(samp.frac * nrow(stockd.in)), replace = F)
  tidx <- seq(1, ncol(stockd.in), by = 2)
  for ( i in idx ) {
    lines(stockdl$time[tidx], stockd.in[i,tidx], lwd = 0.75, col = '#DBD6FF' )
  }
  
  # points outside bands
  # subsetting rows and columns to make smaller pdf
  if (nrow(stockd.out) > 0){
    set.seed(42)
    idx <- sample(1:nrow(stockd.out), floor(samp.frac * nrow(stockd.out)), replace = F) 
    tidx <- seq(1, ncol(stockd.out), by = 2)
    for ( i in idx ) {
      lines(stockdl$time[tidx], stockd.out[i,tidx], lwd = 0.75, col = colors[2])
    }

  }
  
  # confidence bands
  lines(stockdl$time, stockd[cb$upmask], lwd=2.5, col = colors[3])
  lines(stockdl$time, stockd[cb$downmask], lwd=2.5, col = colors[3])

  # "All inlier"
  in.idx <- which(rownames(stockd) == 'MORN')
  lines(stockdl$time, stockd[in.idx,], col = 'black', lwd = 2.5)
  
  # Code to browser "all inlier"s
  # idx1 <- which(rowSums(cb$Z)==1258)
  # lines(stockdl$time, stockd[idx1[1],], col = 'red', lwd = 2.5)
  # # browse observations that are completely within bands
  # stockid.arr <- rownames(stockd)[idx1]
  # tmp <- subset(stockdl$info, Symbol %in% stockid.arr)
  # View(tmp)
  
  # global outlier
  # idx2 <- which(rowSums(cb$Z)==0)
  # lines(stockdl$time, stockd[idx2[30],], col = 'red', lwd = 2.5)
  
  # local outlier
  tind <- which(stockdl$time == as.Date('2012-07-02'))
  rind <- which.max(stockd[cb$row.inc.idx, tind])
  orind <- cb$row.inc.idx[rind]
  # Note: stockdl$info is not in the same order as stockdl$data !
  stockid <- rownames(stockd)[orind]
  stockdl$info[stockdl$info$Symbol == stockid,]
  lines(stockdl$time, stockd[orind,], col = '#D04545', lwd = 2.5)

  # Set x-axis
  axis.Date(1, at = c(seq(xlim[1], xlim[2], "quarter"), xlim[2]),
            labels = c('11/Q1','','11/Q3','',
                       '12/Q1','','12/Q3','',
                       '13/Q1','','13/Q3','',
                       '14/Q1','','14/Q3','',
                       '15/Q1','','15/Q3','',
                       '16/Q1') )
  
}

savefile <- file.path(fig.save.path,
                      sprintf('stockdata_blw%d_K%d_L%d.png', blw, k, l))
png(savefile, width = 15, height = 5, units = 'in', res = 150)
plot.stock.data(stockdl, stockd, cb)
dev.off()

savefile <- file.path(fig.save.path,
                      sprintf('stockdata_blw%d_K%d_L%d.pdf', blw, k, l))
pdf(savefile, width = 15, height = 5)
plot.stock.data(stockdl, stockd, cb)
dev.off()
