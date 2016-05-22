
  figure.timeseries.males = function( outdir, vars="totno.male",  all.areas=T ) {
 
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
      
      ylim=NULL

      #ylim=range(c(td$mean), na.rm=T); ylim[1]=ylim[1]-0.1*ylim[2]; ylim[2] = ylim[2]+ylim[2]*0.4 +200
      ylim[2]=max(td$ub)
      ylim[1]=round(min(td$lb), 0) - (max(td$ub)*0.05)
      xlim=range(td$year); xlim[1]=xlim[1]; xlim[2]=xlim[2]+1

      #xlabels = c(xlim[1], xlim[1]+2, xlim[1]+4, xlim[1]+6, xlim[1]+8, xlim[1]+10, xlim[1]+12, xlim[1]+14, xlim[1]+16, xlim[1]+18) 
      #xlabels = c(xlim[1], xlim[1]+4, xlim[1]+8, xlim[1]+12, xlim[1]+16, xlim[1]+18) 
      xlabels = seq(min(xlim), max(xlim), 2)
      ylabels = round(seq(0, round(ylim[2], -2), length.out=6),-1)


      dir.create( outdir, recursive=T, showWarnings=F )
      fn = file.path( outdir, paste( v, "combined", "pdf", sep="." ) )
      pdf(file=fn, width=4, height=6, bg='white')
      #Cairo( file=fn, type="png", bg="white",  units="in", width=4, height=6, dpi=450)
      setup.lattice.options()
      pl = xyplot( mean~year|region, data=td, ub=td$ub, lb=td$lb,
        layout=c(1,n.areas), 
        par.strip.text = list(cex=1.0),
        par.settings=list(  
          plot.symbol=list(col='black', fill='darkgrey', cex=0.75, pch=21),
          axis.text=list(cex=0.7), 
          par.main.text = list(cex=1),
          layout.heights=list(strip=1, panel=1, main=0.5 ), 
          strip.background=list(col="lightgrey")
        ),
        #ylim=ylim,
        #xlim = xlabels, 
        ylim=(c(ylim[1], ylim[2])),
        scales = list(y =list(at=ylabels, labels=ylabels), x=list(at=xlabels, labels=xlabels, rot=45)),
        main=v, xlab=list("Year", cex=1), ylab=list("Geometric mean No. / km^2", cex=1),
        panel = function(x, y, subscripts, ub, lb, ...) {
          larrows(x, lb[subscripts], x, ub[subscripts], angle = 90, code = 3, length=0.01,lwd=1)
          panel.xyplot(x, y, type = "b", lty=1, lwd=1, col="black", ...)
          #panel.abline(h=0, col="gray75", lwd=1, ...)
        }
      )
      print(pl)
      
      dev.off()
      #cmd( "convert -trim -frame 10x10 -mattecolor white ", fn, fn )
    }
    return("Done")
  }

  
