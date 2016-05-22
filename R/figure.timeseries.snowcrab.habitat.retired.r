
  figure.timeseries.snowcrab.habitat.retired = function(p ) {
    
    td = prediction.surface(p, vclass="R0.mass" )
    areas = c("cfa4x", "cfasouth", "cfanorth" )
    regions = c("4X", "S-ENS", "N-ENS")
    td$region = factor(td$cfa, levels=areas, labels=regions)

    td = td[is.finite(td$sa) ,]
    td$sa = td$sa / 1000
    td = td[order(td$region,td$yr),]
    xlim=range(td$yr); xlim[1]=xlim[1]-0.5; xlim[2]=xlim[2]+0.5
    
    fn = file.path( p$annual.results, "timeseries",  "kriged", "snowcrab.habitat.sa" )
    outdir = dirname( fn ) 
    dir.create( outdir, recursive=T, showWarnings=F )
    Cairo( file=fn, type="svg", bg="white", units="in", width=6, height=8, dpi=75)

    setup.lattice.options()
    pl = xyplot( sa~yr|region, data=td, 
        layout=c(1,3), xlim=xlim, scales = list(y = "free"),
            main="Potential snow crab habitat", xlab="Year", ylab=expression(paste("Surface area ( X", 10^3, km^2, ")")),
            panel = function(x, y, subscripts, ...) {
            panel.abline(h=mean(y, na.rm=T), col="gray40", lwd=1.5,...)
            panel.xyplot(x, y, type="b", pch=19, lwd=1.5, lty="11", col="black", ...)
#            panel.loess(x,y, span=0.15, lwd=2.5, col="darkblue", ... )
       }
    )
    print(pl)
    dev.off()
    cmd( "convert   -trim -quality 9  -geometry 100% -frame 2% -mattecolor white -antialias ", paste(fn, "svg", sep="."),  paste(fn, "png", sep=".") )
 
    means = tapply(td$sa, td$region, mean, na.rm=T)

    print("mean annual SA X 1000 km^2")
    print(means)
    print("SD:")
    print( tapply(td$sa, td$region, sd, na.rm=T)) 
    print( "latest year:" )
    print( td$sa[ td$yr==p$current.assessment.year ])
    table.view( td )

    return( fn )  
  } 


