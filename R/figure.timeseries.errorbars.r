 

  figure.timeseries.errorbars = function( dat, outdir, fname="tmp" ) {
    
    x = dat[,1]
    y = dat[,2]
    
    if ( ncol(dat) == 4 ) {
      lb = dat[,3]
      ub = dat[,4]
    } else if( ncol(dat) ==3 ) {
      e = dat[,3]
      ub = y + e
      lb = y - e
    }
   
    sdy = sd( y )

    xlim = range( x, na.rm=T )
    xlim[1] = xlim[1]-0.5
    xlim[2] = xlim[2]+0.5
    
    ylim = range( c( ub, lb, 0) , na.rm=T)
    
    par( cex = 2, cex.lab=4 )

    plot( x, y, type="p", pch=19, xlim=xlim, ylim=ylim, ylab="", xlab="")
    lines (  x,y,  col="gray50" , lwd = 2 )
    lines ( lowess( x,y, f=0.25 ), col="blue" , lwd = 3 )
    arrows( x, ub, x, lb,  col="gray50", length=0.1, angle=90, code=3, lwd=1.5)
    
    Pr(dev="png", dname=outdir, fname=fname, width=16, height=12)

    return("Done")  
  } 



