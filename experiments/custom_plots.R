

motivating.plot <- function(pdn, pd, pd0, pd1, xvar, yvar){
  
  
  lw <- 0.8
  lw2 <- 0.8
  colors <- RColorBrewer::brewer.pal(3, 'Dark2')
  p <- ggplot(data = pd, mapping = aes(x = variable, y = value, group = row))
  
  #p <- p + geom_line(alpha = 0.1)
  p <- p + geom_line(color = "#DBDBDB")
  
  # p <- p + geom_line(data = pdcb, aes(x = variable, y = value, group = row, color = type),
  #                    linetype = 'solid' , size = lw2)
  
  p <- p + geom_line(data = pdn, color = colors[[3]], linetype = 'solid', size = lw2)
  p <- p + geom_line(data = pd0, color = colors[[1]], linetype = 'solid', size = lw2)
  p <- p + geom_line(data = pdl, color = colors[[2]], linetype = 'solid', size = lw2)
  
  # p <- p + geom_line(data = pdn, color = 'red', linetype = 'solid' , size = lw2)
  # p <- p + geom_line(data = pd0, color = 'blue', linetype = 'solid', size = lw2)
  # p <- p + geom_line(data = pdl, color = 'green',, linetype = 'solid', size = lw2)
  
  
  p <- p + geom_line(data = subset(pd, row == 1), size = lw, linetype = 'dotdash') #normal
  p <- p + geom_line(data = subset(pd, row == 2), size = lw, linetype = 'solid') #local outlier
  p <- p + geom_line(data = subset(pd, row == 3), size = lw, linetype = 'dashed') #global outier
  
  p <- p + theme_light()
  p <- p + theme( axis.ticks = element_blank(),
                  axis.text.x = element_blank(),
                  axis.text.y = element_blank(), 
                  panel.grid.minor = element_blank(),
                  panel.grid.major = element_blank(),
                  legend.position="bottom", 
                  panel.border = element_blank(),
                  plot.margin = unit(c(1,1,1,1.5), "lines"))
  p <- p + labs(x = NULL, y = NULL)
  
  # annotate confidence bands
  datatext <- data.frame(x = c('V1','V1','V1','V1','V1','V1'),
                         y = c(0.08, 0.05, 0.035, -0.03, -0.05, -0.07),
                         label = c('L=0','Naive','L=25','L=25','Naive','L=0'),
                         stringsAsFactors = F)
  p <- p + geom_text(mapping = aes(x=x, y=y, label=label, group=NULL, hjust = "right"),
                     nudge_x = -1,
                     data =  datatext,
                     size = 3)
  
  # annotate confidence bands
  datatext <- data.frame(x = c( xvar, yvar, xvar, yvar),
                         y = c( -0.11, -0.11, -0.09, -0.09),
                         label = c('x','y', '|', '|'),
                         stringsAsFactors = F)
  p <- p + geom_text(mapping = aes(x=x, y=y, label=label,
                      group=NULL, hjust = "center"),
                      nudge_x = -1,
                      data =  datatext,
                      size = 3)
  
  # some code to override clipping
  # From: http://stackoverflow.com/questions/12409960/ggplot2-annotate-outside-of-plot
  gt <- ggplot_gtable(ggplot_build(p))
  gt$layout$clip[gt$layout$name == "panel"] <- "off"
  #grid.draw(gt)
  gt
}



cross.plot <- function(x, y, K, lwd = 2, lwd.bracket = 3){
  require(RColorBrewer)
  colors <- RColorBrewer::brewer.pal(3, 'Dark2')
  
  qx <- smallest_ci(x, K)
  qy <- smallest_ci(y, K)
  
  ## Greedy
  dm <- matrix(c(x, y), ncol=2, byrow = F)
  cb0 <- findcb_topdown(dm, K = K, L = 0)
  qx0 <- c(dm[cb0$downmask][1], dm[cb0$upmask][1])
  qy0 <- c(dm[cb0$downmask][2], dm[cb0$upmask][2])
  
  cb1 <- findcb_topdown(matrix(c(x, y), ncol=2, byrow = F),
                        K = K, L = 1)
  qx1 <- c(dm[cb1$downmask][1], dm[cb1$upmask][1])
  qy1 <- c(dm[cb1$downmask][2], dm[cb1$upmask][2])
  
  par(mar=rep(1.5,4))
  
  lll <- 0.25
  
  plot(c(min(x)-3*lll,max(x)),c(min(y)-3*lll,max(y)),
       type="n",main="",xlab="x",ylab="y",
       xaxt="n",yaxt="n",bty="n",mgp=c(0.3,0.3,0.3))
  
  
  # naive
  lines(c(qx[1],qx[1],qx[2],qx[2]),
        c(min(y)+0.6*lll-3*lll,min(y)-3*lll,min(y)-3*lll,min(y)+0.6*lll-3*lll), col = colors[3], lwd = lwd.bracket)
  lines(c(qx[1],qx[1]),
        c(min(y)-3*lll+0.4*lll,max(y)),lty="dashed", col = colors[3], lwd = lwd)
  lines(c(qx[2],qx[2]),
        c(min(y)-3*lll+0.4*lll,max(y)),lty="dashed", col = colors[3], lwd = lwd)
  rug(x,side=1)
  
  lines(c(min(x)-3*lll+0.4*lll,max(x)),c(qy[1],qy[1]),lty="dashed", col = colors[3], lwd = lwd)
  lines(c(min(x)-3*lll+0.4*lll,max(x)),c(qy[2],qy[2]),lty="dashed", col = colors[3], lwd = lwd)
  lines(c(min(x)+0.6*lll-3*lll,min(x)-3*lll,min(x)-3*lll,min(x)+0.6*lll-3*lll),
        c(qy[1],qy[1],qy[2],qy[2]), col = colors[3], lwd = lwd.bracket)
  rug(y,side=2)
  
  text(qx[2],min(y)-3*lll+0.2*lll,pos=4,cex=1,labels="naive", col = colors[3])
  text(qx0[2],min(y)-2*lll+0.2*lll,pos=4,cex=1,labels="L = 0", col = colors[1])
  text(qx1[2]+0.1,min(y)-lll+0.2*lll,pos=4,cex=1,labels="L = 1", col = colors[2])
  
  
  polygon(c(min(x),qx1[1],qx1[1],qx1[2],qx1[2],max(x),max(x),qx1[2],qx1[2],qx1[1],qx1[1],min(x)),
          c(qy1[1],qy1[1],min(y),min(y),qy1[1],qy1[1],qy1[2],qy1[2],max(y),max(y),qy1[2],qy1[2]),
          border=NA,density=10,angle=-45,col="gray90")
  
  
  # L = 0
  polygon(c(qx0[1],qx0[2],qx0[2],qx0[1]),
          c(qy0[1],qy0[1],qy0[2],qy0[2]),
          border=colors[1],density=10,angle=45,col="gray90")
  # brackets
  lines(c(qx0[1],qx0[1],qx0[2],qx0[2]),
        c(min(y)+0.6*lll-2*lll,min(y)-2*lll,min(y)-2*lll,min(y)+0.6*lll-2*lll),col=colors[1],lwd = lwd.bracket)
  lines(c(min(x)+0.6*lll-2*lll,min(x)-2*lll,min(x)-2*lll,min(x)+0.6*lll-2*lll),
        c(qy0[1],qy0[1],qy0[2],qy0[2]),col=colors[1], lwd = lwd.bracket)
  
  # L = 1
  lines(c(min(x),qx1[1],qx1[1]),
        c(qy1[1],qy1[1],min(y)),col=colors[2], lwd = lwd)
  lines(c(qx1[2],qx1[2],max(x)),
        c(min(y),qy1[1],qy1[1]),col=colors[2], lwd = lwd)
  lines(c(max(x),qx1[2],qx1[2]),
        c(qy1[2],qy1[2],max(y)),col=colors[2], lwd = lwd)
  lines(c(qx1[1],qx1[1],min(x)),
        c(max(y),qy1[2],qy1[2]),col=colors[2], lwd = lwd)
  #brackets
  lines(c(qx1[1],qx1[1],qx1[2],qx1[2]),
        c(min(y)+0.6*lll-lll,min(y)-lll,min(y)-lll,min(y)+0.6*lll-lll),
        col=colors[2],lwd = lwd.bracket)
  lines(c(min(x)+0.6*lll-lll,min(x)-lll,min(x)-lll,min(x)+0.6*lll-lll),
        c(qy1[1],qy1[1],qy1[2],qy1[2]),
        col=colors[2],lwd = lwd.bracket)
  
  
  points(x,y,pch=20)
}
