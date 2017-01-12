# Implied packages:
# require(ElemStatLearn) #makeX()

#' Generate data from the ozone dataset
makeX <- function( n,m ) {
  require(ElemStatLearn)
  ozone <- ElemStatLearn::ozone
  xpts <- seq( min(ozone$radiation), max(ozone$radiation), length.out=m )
  t(sapply( 1:n, function(i) {
    idx <- sample( nrow(ozone), replace=TRUE )
    stats::ksmooth( x=ozone$radiation[idx], y=ozone$ozone[idx],
             kernel='normal', bandwidth=100, x.points=xpts )$y
  } ))
}


#' Load stock data from an rds file
load_data_stock <- function(rdsfile = OPTS$file$stockdata.file, 
                            type = 'prc'){
  
  stockd <- readRDS(rdsfile)
  
  if (type == 'prc'){
    list(data = stockd$closingprice_pc,
         time = stockd$time_pc,
         info = stockd$meta)
  } else {
    list(data = stockd$closingprice,
         time = stockd$time,
         info = stockd$meta)
  }
}
