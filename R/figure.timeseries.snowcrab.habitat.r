 
  figure.timeseries.snowcrab.habitat = function(p ) {

    td = bio.snowcrab::interpolation.db( p=p, DS="timeseries" )
    
    areas = c("cfa4x", "cfasouth", "cfanorth" )
    regions = c("4X", "S-ENS", "N-ENS")
    td$region = factor(td$region, levels=areas, labels=regions)
    td$sa = td$sa.region

    td = td[is.finite(td$sa) ,]
    td$sa = td$sa / 1000
    td = td[order(td$region, td$yr),]
    xlim=range(td$yr); xlim[1]=xlim[1]-0.5; xlim[2]=xlim[2]+0.5
    
    fn = file.path( p$annual.results, "timeseries",  "interpolated", "snowcrab.habitat.sa" )
    outdir = dirname( fn ) 
    dir.create( outdir, recursive=T, showWarnings=F )
    png( filename=paste(fn, "png", sep="."), width=3072, height=2304, pointsize=40, res=300 )

    setup.lattice.options()
    pl = xyplot( sa~yr|region, data=td, 
        layout=c(1,3), xlim=xlim, scales = list(y = "free"),
            main="Potential snow crab habitat", xlab="Year", ylab=expression("Surface area; X 10^3 km^2"),
            panel = function(x, y, subscripts, ...) {
            panel.abline(h=mean(y, na.rm=T), col="gray40", lwd=1.5,...)
            panel.xyplot(x, y, type="b", pch=19, lwd=1.5, lty="11", col="black", ...)
#            panel.loess(x,y, span=0.15, lwd=2.5, col="darkblue", ... )
       }
    )
    print(pl)
    dev.off()
 
    means = tapply(td$sa, td$region, mean, na.rm=T)

    print("mean annual SA X 1000 km^2")
    print(means)
    print("SD:")
    print( tapply(td$sa, td$region, sd, na.rm=T)) 
    print( "latest year:" )
    print( td$sa[ td$yr==p$year.assessment ])
   # table.view( td )

    return( fn )  
  } 


