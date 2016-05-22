
  figure.timeseries.larvae = function( outdir ) {
    larvae = read.table(file=file.path(  project.datadirectory("snowcrab"), "data", "larvae", "ts.brachyura.csv"),  sep = "\t")
    colnames(larvae) = c("yr", "mean.n.m3", "sdev", "n")
    larvae$se = larvae$sdev / sqrt(larvae$n-1)

    ts.plotandsave( x=larvae$yr, y=larvae$mean.n.m3, lb=larvae$mean.n.m3-larvae$se, ub=larvae$mean.n.m3+larvae$se,
      w=larvae$n, outdir=outdir, action="save", title="Brachyura", smooth=0.75 ) 

    l.monthly = read.table(file=file.path(  project.datadirectory("snowcrab"), "data", "larvae", "ts.brachyura.monthly.csv"),  sep = "\t")
    colnames(l.monthly) = c("month", "mean.n.m3", "n")

    x = l.monthly$month
    y = l.monthly$mean.n.m3
    w = l.monthy$n

    xrange = range(-0.5, 12.5)
   
  fn = file.path(outdir, "BrachyuraMonthly")
   Cairo( file=fn, type="pdf", bg="white", units="in", width=8, height=4 )
  
    plot( x, y, type="n", axes=F, xlab="Month", xlim=xrange, ylab=toupper("Brachyura"), cex=2 )
    lines( loess( y ~ x, weights=w, control=loess.control(surface = "direct")), col="orange", lty="solid", lwd=4 )

    points(x, y, pch=10)
    axis( 1 )
    axis( 2 )
   dev.off()
    cmd( "convert   -trim -quality 9  -geometry 200% -frame 2% -mattecolor white -antialias ", paste(fn, "pdf", sep="."),  paste(fn, "png", sep=".") )
  
  }



