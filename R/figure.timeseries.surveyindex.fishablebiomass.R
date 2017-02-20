
  figure.timeseries.surveyindex.fishablebiomass = function(p ) {

    td = bio.snowcrab::interpolation.db( p=p, DS="timeseries" )
    
    areas = c("cfa4x", "cfasouth", "cfanorth" )
    regions = c("4X", "S-ENS", "N-ENS")
    td$region = factor(td$region, levels=areas, labels=regions)

    td = td[is.finite(td$total) ,]
    td = td[order(td$region, td$yr),]
    xlim=range(td$yr); xlim[1]=xlim[1]-0.5; xlim[2]=xlim[2]+0.5
    
    fn = file.path( p$annual.results, "timeseries",  "interpolated", "snowcrab.biomass.index" )
    outdir = dirname( fn ) 
    dir.create( outdir, recursive=T, showWarnings=F )
    Cairo( file=paste(fn, "svg", sep="."), type="svg", bg="white", units="in", width=6, height=8, dpi=75)

    setup.lattice.options()
    pl = xyplot( total~yr|region, data=td, 
        layout=c(1,3), xlim=xlim, scales = list(y = "free"),
            main="Potential snow crab habitat", xlab="Year", ylab=expression(paste("Surface area ( X", 10^3, km^2, ")")),
            panel = function(x, y, subscripts, ...) {
            panel.abline(h=mean(y, na.rm=T), col="gray40", lwd=1.5,...)
            panel.xyplot(x, y, type="b", pch=19, lwd=1.5, lty="11", col="black", ...)
            panel.loess(x,y, span=0.15, lwd=2.5, col="darkblue", ... )
       }
    )
    print(pl)
    dev.off()
    cmd( "convert   -trim -quality 9  -geometry 100% -frame 2% -mattecolor white -antialias ", paste(fn, "svg", sep="."),  paste(fn, "png", sep=".") )
 
    means = tapply(td$total, td$region, mean, na.rm=T)

    print("t/km^2")
    print(means)
    print("SD:")
    print( tapply(td$total, td$region, sd, na.rm=T)) 
    print( "latest year:" )
    print( td$total[ td$yr==p$year.assessment ])
    table.view( td )

    return( fn )  
  } 


