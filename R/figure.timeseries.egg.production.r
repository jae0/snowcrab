
  figure.timeseries.egg.production= function( outdir=file.path(project.datadirectory('snowcrab'), "assessments","2013"), vars='fecundity', all.areas=T ) {
 
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
    varnames = names(set)
#    vars = varnames[ grep ( vars, varnames) ]

    for (v in vars ) {

      td =  get.time.series ( from.file=T )
      ii =  which( td$variable == v)
      if (length(ii) < 5) next()
      td = td[ii ,]
      td = td[ order(td$region, td$year) , ]
      td$region = factor(td$region, levels=areas, labels =regions)

      ylim=range(c(td$mean), na.rm=T); ylim[1]=ylim[1]-0.1*ylim[2]; ylim[2] = ylim[2]+ylim[2]*0.2
      xlim=range(td$year); xlim[1]=xlim[1]-0.5; xlim[2]=xlim[2]+0.5

      dir.create( outdir, recursive=T, showWarnings=F )
      fn = file.path( outdir, paste( v, "combined", "png", sep="." ) )
      Cairo( file=fn, type="png", bg="white",  units="in", width=5, height=6, dpi=300)
      setup.lattice.options()
      pl = xyplot( mean~year|region, data=td, ub=td$ub, lb=td$lb,
        layout=c(1,n.areas), 
        par.strip.text = list(cex=3.0),
        par.settings=list(  
          axis.text=list(cex=2.5), 
          par.main.text = list(cex=4),
          layout.heights=list(strip=0.2, panel=1, main=0.5 ) 
        ),
        xlim=xlim, scales = list(y = "free"),
        main=v, xlab=list("Year", cex=3), ylab=list("Geometric mean '000 / km^2", cex=3),
        panel = function(x, y, subscripts, ub, lb, ...) {
          larrows(x, lb[subscripts], x, ub[subscripts], angle = 90, code = 3, length=0.05,lwd=3)
          panel.xyplot(x, y, type="b", lty=1, lwd=6, pch=20, col="black", ...)
          panel.abline(h=0, col="gray75", lwd=5, ...)
        }
      )
       
      print(pl)
      dev.off()
      cmd( "convert -trim -frame 10x10 -mattecolor white ", fn, fn )
    
    }

    return("Done")
    }



