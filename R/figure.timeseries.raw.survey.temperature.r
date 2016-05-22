
  figure.timeseries.raw.survey.temperature = function( outdir, vars="t",  all.areas=T ) {

    
    set = snowcrab.db( DS="set.merge.det")
  
    if (all.areas) {
      areas = c("cfa4x", "cfasouth", "cfanorth" )
      regions = c("4X", "S-ENS", "N-ENS")
    } else {
      areas = c("cfasouth", "cfanorth" )
      regions = c("S-ENS", "N-ENS")
    }
    cex.main = 1.4
    cex.lab =1.5
    cex.axis = 1.5

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
      
      ylim='NULL'
      ylim[2]=max(td$ub)
      ylim[1]=round(min(td$lb), 0) - (max(td$ub)*0.05)
      xlim=range(td$year); xlim[1]=xlim[1]; xlim[2]=xlim[2]
      
      xlabels = seq(min(xlim), max(xlim), 2)
      ylabels = round(seq(0, round(as.numeric(ylim[2]), 1), length.out=8), 1)

      #ylim=range(c(td$mean), na.rm=T); ylim[1]=ylim[1]-0.1*ylim[2]; ylim[2] = ylim[2]+ylim[2]*0.2
      #xlim=range(td$year); xlim[1]=xlim[1]-0.5; xlim[2]=xlim[2]+0.5

      dir.create( outdir, recursive=T, showWarnings=F )
      fn = file.path( outdir, paste( 'temperature', "combined.pdf",  sep="." ) )
      #Cairo( file=fn, type="png", bg="white",  units="in", width=5, height=8, dpi=450 )
      pdf(file=fn, width=5, height=8, bg='white')
      setup.lattice.options()
      pl = xyplot( mean~year|region, data=td, ub=td$ub, lb=td$lb,
          layout=c(1,n.regions), 
          par.strip.text=list(
            plot.symbol=list(col='black', fill='darkgrey', cex=0.75, pch=21),
            axis.text=list(cex=0.7),
            par.main.text=list(cex=1),
            layout.heights=list(strip=1, panel=1, main=0.5),
            strip.background=list(col='lightgrey')),
          #xlim=xlim, 
          ylim=(c(as.numeric(ylim[1]), as.numeric(ylim[2]))),
          scales=list(y=list(at=ylabels, labels=ylabels, cex=0.65), x=list(at=xlabels, labels=xlabels, rot=50, cex=0.8)),
          main="Snow Crab Survey Temperature", xlab="Year", ylab=list("Geometric mean Temperature"),
          #cex.lab=cex.lab, 
          cex.axis=cex.axis, cex.main = cex.main,
          panel = function(x, y, subscripts, ub, lb, ...) {
          larrows(x, lb[subscripts], x, ub[subscripts], angle = 90, code = 3, length=0.05)
             panel.xyplot(x, y, type="b", lty=1, lwd=1.5, pch=21, fill= 'darkgrey', col="black", ...)
             panel.abline(median(y), col="gray75", ...)
         }
      )
    
    print(pl)
    dev.off()
  # cmd( "convert   -trim -quality 9  -geometry 200% -frame 2% -mattecolor white -antialias ", paste(fn, "pdf", sep="."),  paste(fn, "png", sep=".") )
   
    }

    return("Done")
  }
 