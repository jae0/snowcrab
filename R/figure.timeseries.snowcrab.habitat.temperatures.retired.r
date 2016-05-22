 figure.timeseries.snowcrab.habitat.temperatures.retired = function(p) {
        
    td = prediction.surface(p, vclass="R0.mass" )
    areas = c("cfa4x", "cfasouth", "cfanorth" )
    regions = c("4X", "S-ENS", "N-ENS")
    td$region = factor(td$cfa, levels=areas, labels=regions)
    td$ubound = td$t + td$t.sd 
    td$lbound = td$t - td$t.sd 
    td = td[is.finite(td$t) ,]
    td = td[order(td$region,td$yr),]
    xlim=range(td$yr); xlim[1]=xlim[1]-0.5; xlim[2]=xlim[2]+0.5
    ylim= range( c(td$t, td$ubound, td$lbound)  ); ylim[1]=ylim[1]-0.5; ylim[2]=ylim[2]+0.5
    
    fn = file.path( p$annual.results, "timeseries",  "kriged", "mean.bottom.temp.snowcrab.habitat" )
    
    dir.create( dirname(fn), recursive=T, showWarnings=F  )
    
    Cairo( file=fn, type="svg", bg="white", units="in", width=6, height=8, dpi=75 )

    setup.lattice.options()
    pl = xyplot( t~yr|region, data=td, ub=td$ubound, lb=td$lbound,
        layout=c(1,3), xlim=xlim, ylim=ylim, cex=3, # scales = list(y = "free"),
            main="Temperature in potential habitats", xlab="Year", ylab="Celcius",
            panel = function(x, y, subscripts, ub, lb, ...) {
            panel.abline(h=mean(y, na.rm=T), col="gray40", lwd=2, ...)
            larrows(x, lb[subscripts],
                    x, ub[subscripts],
                   angle = 90, code = 3, length=0.005, lwd=1, col="darkgray")
            panel.xyplot(x, y, type="b", lwd=4, pch=".", col="black", ...)
#            panel.loess(x,y, span=0.15, lwd=2)
       }
    )
    print( pl )
    dev.off()
    cmd( "convert   -trim -quality 9  -geometry 100% -frame 2% -mattecolor white -antialias ", paste(fn, "svg", sep="."),  paste(fn, "png", sep=".") )
    
    means = tapply(td$t, td$region, mean, na.rm=T)
    print("mean temperatures:")
    print(means)
    print("SD:")
    print( tapply(td$t, td$region, sd, na.rm=T) ) 
    print( "latest year:" )
    print( td$t[ td$yr==p$current.assessment.year ])
    table.view( td )
    return(fn)  
  }


