 figure.timeseries.conditional = function( B, main="", xlab="", ylab="", span=0.4 ) {
 
    names(B) = c( "y", "x", "f" ) # value, value, factor
    nf = length( unique(B$f) )  # a factor
    #browser()

    #ylim=range( B$y, na.rm=T); dy=(ylim[2]-ylim[1])/10; ylim=ylim+c(-dy, +dy)
    ylim='NULL'
    rg = range(B$y, finite=TRUE)
    ylim[1] = rg[1] - (rg[2]*0.1)
    ylim[2] = rg[2]*1.1
    #ylim =range(B$y, na.rm=T)
    xlim=range( B$x, na.rm=T); dx=1; xlim=xlim+c(-dx, +dx)

    xlabels = seq(min(xlim), max(xlim), 1)
   # xlabels = c(xlim[1], xlim[1]+4, xlim[1]+8, xlim[1]+12, xlim[1]+16, xlim[1]+18) 
    ylabels = round(seq(as.numeric(rg[1]), round(as.numeric(rg[2]), 2), length.out=6),2)

    setup.lattice.options()
    pl = xyplot( y ~ x | f, data=B, 
      layout=c(1,nf), 
      par.strip.text = list(cex=1),
      par.settings=list(  
        axis.text=list(cex=0.75), 
        par.main.text = list(cex=1),
        layout.heights=list(strip=1, panel=1, main=0.3 ),
        strip.background=list(col="lightgrey"), 
        box.umbrella = list(col = "black"), 
        box.rectangle = list(col="black"), 
        box.dot = list(col = "black", pch = 3, cex=2),
        plot.symbol=list(col='black', fill='darkgrey', cex=0.75, pch=21)
      ),
      #xlim=xlim, 
      ylim=(c(as.numeric(ylim[1]), as.numeric(ylim[2]))),
      scales = list(y =list(at=ylabels, labels=ylabels), x=list(at=xlabels, labels=xlabels, rot=45)),
      #scales = list(y = "free", x=list(at=xlabels, labels=xlabels, rot=45)),
      main=main,
      xlab=xlab, ylab=ylab,
      panel = function( x, y, ...) {
        panel.bwplot( x, y,  box.ratio=0.5 , horizontal=F, notch=F, pch="|", col="black", lwd=1, ...)
        panel.xyplot( x, y, type=c("a"), col="black", lwd=2, ...)
        #panel.xyplot( x, y, type=c("smooth", "a"), span=span, col="red", lty=2, lwd=1, ...)
      }
    )
     
      return(pl)
  }

