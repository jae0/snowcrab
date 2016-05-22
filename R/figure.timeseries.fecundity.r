
  figure.timeseries.fecundity = function( p, outdir, all.areas=T ) {
  
    set = snowcrab.db("set.merge.det")

    if (all.areas) {
      areas = c("cfa4x", "cfasouth", "cfanorth" )
      regions = c("4X", "S-ENS", "N-ENS")
    } else {
      areas = c("cfasouth", "cfanorth" )
      regions = c("S-ENS", "N-ENS")
    }

    n.regions = length(regions)
    n.areas = length(areas)

    v = "fecundity"

    eps = 1e-6

    p$vars.to.model = v
    p$years.to.model = 1998:p$current.assessment.year
    p = make.list( list(y=p$years.to.model, v=p$vars.to.model ), Y=p )
    K = interpolation.db( DS="interpolation.simulation", p=p )

    td = K
 
    xlim = range( td$yr, na.rm=T )
    xlim[1] = xlim[1]-0.5
    xlim[2] = xlim[2]+0.5
 
    yy =  "ln( Total potential fecundity ) "
    convert = 10^6  # convert from  

    varstocheck = c("total", "ubound", "lbound")
    for (vs in varstocheck) {
        td[,vs] = td[,vs] / convert
        kk = which(td[,vs] <= eps)
        if (length(kk)>0) td[kk,vs] = 0
    }

    td = td[ which( is.finite(td$total)) ,]
    td = td[order(td$yr),]
    td$region = factor(td$region, levels=areas, labels=regions)

    fn = file.path( outdir, v )
    Cairo( file=fn, type="pdf", bg="white", units="in", width=6, height=8 )
 
    setup.lattice.options()
    pl = xyplot( total~yr|region, data=td,  lbound=td$lbound, ubound=td$ubound,
        layout=c(1,3), xlim=xlim, scales = list(y = "free"),
           main="Potential Egg Production", xlab="Year", 
           ylab=expression(paste("Potential Egg Production ( x", 10^6, ")")),
           panel = function(x, y, subscripts, lbound, ubound, ...) {
             larrows(x, lbound[subscripts], x, ubound[subscripts],  angle = 90, code = 3, length=0.02, lwd=3)
             panel.abline(h=mean(y, na.rm=T), col="gray40", lwd=1.5,...)
             panel.xyplot(x, y, type="b", pch=19, lwd=1.5, lty="11", col="black", ...)
#            panel.loess(x,y, span=0.15, lwd=2.5, col="darkblue", ... )
       }
    )
    print(pl)
    dev.off()
 
    cmd( "convert   -trim -quality 9  -geometry 200% -frame 2% -mattecolor white -antialias ", paste(fn, "pdf", sep="."),  paste(fn, "png", sep=".") )
   
     table.view ( td )  
     return(td)

  }  

 

