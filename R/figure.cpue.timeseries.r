  figure.cpue.timeseries = function( yearmax, outdir=NULL, outfile=NULL, outfile2=NULL ) {
   
    regions = c("cfanorth", "cfasouth", "cfa4x")
    k = NULL
    for (r in regions) {
      res = get.fishery.stats.by.region(Reg=r)
      if(r=='cfanorth') {ii=which(res$yr==2014); res[ii,'cpue'] <- 104.5; ii=which(res$yr==2013); res[ii,'cpue'] <- 106}
      k = cbind( k, res$cpue )
    }
    
    k = as.data.frame( k )
    colnames(k) = regions
    rownames(k) = res$yr
    
    k = k[ which( as.numeric(rownames(k)) <= yearmax ), ] 
    uyrs = as.numeric(rownames(k) ) 

    pts = c(19, 22, 24)
    lns = c(1, 1, 1)
    cols = c("grey10", "grey10",  "grey20")

    yrange = range (k, na.rm=T)
    yrange[1] = 0
    xrange = range(uyrs)
    xrange[1] = xrange[1]
    xrange[2] = xrange[2]
#    xlabels = c(xrange[1], xrange[1]+8, xrange[1]+18, xrange[1]+28, xrange[2])
    xlabels = seq(xrange[1]+1, xrange[2], 2)
    
    dir.create( outdir, recursive=T, showWarnings=F  )
    fn = file.path( outdir, paste( outfile, "pdf", sep="." ) )

    #Cairo( file=fn, type="png", bg="white", , pointsize=30, units="in", width=6, height=4, dpi=300 )
    pdf(file=fn, width=7, height=7, bg='white')
    #png( file=fn,units='in', width=7,height=7,pointsize=10, res=350,type='cairo')
      m=1; plot( uyrs, k[,m],  type="b", ylab="Catch rate (kg/trap)", xlab="Year", col=cols[m], lwd=3, lty=lns[m], pch=pts[m], xaxt="n", xlim=xrange, ylim=yrange)
      m=2; points(uyrs, k[,m], type="b", col=cols[m], lwd=3, lty=lns[m], pch=pts[m])
      m=3; points(uyrs, k[,m], type="b", col=cols[m], lwd=3, lty=lns[m], pch=pts[m])
      axis( 1, at=xlabels, labels=FALSE )
      text(x=xlabels+1, y=par('usr')[3], labels=xlabels, srt=45, adj=c(1.5,1), xpd=TRUE)

      axis( 2 )
      legend(x=1980, y=100, c("N-ENS", "S-ENS", "4X"), bty="n", lty=lns, lwd=2, pch=pts, col=cols, cex=1.2 )
    
    dev.off()

   # cmd( "convert -trim -frame 10x10 -mattecolor white ", fn, fn )
    #table.view(k)
    return( fn )
  }


