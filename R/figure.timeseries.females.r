
  figure.timeseries.females = function( outdir, vars="totno.female",  all.areas=T ) {
 
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
    vars = varnames[ grep ( vars, varnames) ]

    for (v in vars ) {

      td =  get.time.series ( from.file=T )
      ii =  which( td$variable == v)
      if (length(ii) < 5) next()
      td = td[ii ,]
      td = td[ order(td$region, td$year) , ]
      td$region = factor(td$region, levels=areas, labels =regions)

      #ylim=range(c(td$mean), na.rm=T); ylim[1]=ylim[1]-0.1*ylim[2]; ylim[2] = ylim[2]+ylim[2]*0.2
      xlim=range(td$year); xlim[1]=xlim[1]-1; xlim[2]=xlim[2]+1
      ylim=NULL
      ylim[2]=max(td$ub) *0.8
      ylim[1]=round(min(td$lb), 0) - (max(td$ub)*0.05)
      
      #xlabels = c(xlim[1]+1, xlim[1]+5, xlim[1]+9, xlim[1]+13, xlim[1]+17, xlim[1]+20) 
      xlabels = seq(min(xlim), max(xlim), 2)

      ylabels = round(seq(0, round(ylim[2], -2), length.out=6),-1)

      dir.create( outdir, recursive=T, showWarnings=F )
      fn = file.path( outdir, paste( v, "combined", "pdf", sep=".") )
      cex.main = 1
      cex.lab = 1
      cex.axis = 1

 pdf(file=fn, width=4, height=6, bg='white')
 #Cairo( file=fn, type="png", bg="white",  units="in", width=4, height=6, dpi=350 )
    setup.lattice.options()
    pl = xyplot( mean~year|region, data=td, ub=td$ub, lb=td$lb,
          layout=c(1,n.regions), xlim=xlim, ylim=ylim,
              main=v, xlab="Year", ylab=list("Geometric mean No. / km^2"),
              cex.lab=cex.lab, cex.axis=cex.axis, cex.main = cex.main,
              par.settings=list(plot.symbol=list(col='black', fill='darkgrey', cex=0.75, pch=21),
                strip.background=list(col="lightgrey")),
              scales = list(y =list(at=ylabels, labels=ylabels), x=list(at=xlabels, labels=xlabels, rot=45)),
              panel = function(x, y, subscripts, ub, lb, ...) {
             larrows(x, lb[subscripts],
                     x, ub[subscripts],
                     angle = 90, code = 3, length=0.01)
             panel.xyplot(x, y, type="b", lty=1, lwd=1, col="black", ...)
             #panel.abline(h=0.5, col="gray75", ...)
         }
      )
    
    print(pl)
    dev.off()
   #cmd( "convert   -trim -quality 9  -geometry 200% -frame 2% -mattecolor white -antialias ", paste(fn, "pdf", sep="."),  paste(fn, "png", sep=".") )
    }

    return("Done")
  }



