
  figure.timeseries.R0.reduced.stations = function( outdir, all.areas=T ) {
 
    set = snowcrab.db( DS="set.merge.det")
  
    if (all.areas) {
      areas = c("cfa4x", "cfasouth", "cfanorth" )
      regions = c("4X", "S-ENS", "N-ENS")
    } else {
      areas = c("cfasouth", "cfanorth" )
      regions = c("S-ENS", "N-ENS")
    }

    n.regions = length(regions)
    n.areas = length(areas)

    v = "R0.mass"

    td =  get.time.series ( from.file=T,reduced.stations=T )
    td = td[ which( td$variable == v) ,]
    td = td[ order(td$region, td$year) , ]
    td$region = factor(td$region, levels=areas, labels =regions)
    #   td[ which(td$region=="4X" & td$year < 2004), c("mean", "se", "ub", "lb", "n")] = NA

    ylim=range(c(td$mean), na.rm=T); ylim[1]=ylim[1]-0.1*ylim[2]; ylim[2] = ylim[2]+ylim[2]*0.2
    xlim=range(td$year); xlim[1]=xlim[1]-0.5; xlim[2]=xlim[2]+0.5

    dir.create( outdir, recursive=T, showWarnings=F )
    fn = file.path( outdir, paste( v, "reduced.stations.combined",  sep="." ) )
        cex.main = 1.4
    cex.lab = 1.3
    cex.axis = 1.3

 Cairo( file=fn, type="pdf", bg="white",  units="in", width=6, height=8 )
    setup.lattice.options()
    pl = xyplot( mean~year|region, data=td, ub=td$ub, lb=td$lb,
          layout=c(1,n.regions), xlim=xlim, ylim=ylim,
              main="Fishable Biomass Reduced Stations", xlab="Year", ylab=list("Geometric mean t / km^2"),
              cex.lab=cex.lab, cex.axis=cex.axis, cex.main = cex.main,
              panel = function(x, y, subscripts, ub, lb, ...) {
             larrows(x, lb[subscripts],
                     x, ub[subscripts],
                     angle = 90, code = 3, length=0.05)
             panel.xyplot(x, y, type="b", lty=1, lwd=2, pch=20, col="black", ...)
             panel.abline(h=median(y), col="gray75", ...)
         }
      )
    
    print(pl)
    dev.off()
   cmd( "convert   -trim -quality 9  -geometry 200% -frame 2% -mattecolor white -antialias ", paste(fn, "pdf", sep="."),  paste(fn, "png", sep=".") )
     return("Done")
    }



