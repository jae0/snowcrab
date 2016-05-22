
  figure.timeseries.CW = function( outdir, all.areas=T, type="trawl" ) {
    
    if (all.areas) {
      areas = c("cfa4x", "cfasouth", "cfanorth" )
      regions = c("4X", "S-ENS", "N-ENS")
    } else {
      areas = c("cfasouth", "cfanorth" )
      regions = c("S-ENS", "N-ENS")
    }

    n.regions = length(regions)
    n.areas = length(areas)

    cex.main = 1.4
    cex.lab = 1.3
    cex.axis = 1.3

    if (type == "trawl") {
      set = snowcrab.db("set.merge.det")
      v = "cw.comm.mean"
      fn = file.path( outdir, paste( v, "combined.pdf", sep="." ) )
      td =  get.time.series ( from.file=T )
      td = td[ which( td$variable == v) ,]
      td$mean =  10^td$mean
      td$ub = 10^td$ub
      td$lb = 10^td$lb

    }
    if (type == "observer") {
      odb = observer.db( DS="odb" )
      odb = odb[ which(odb$cw>=95),]
      v = "cw"
      fn = file.path( outdir, paste( v, "combined.observer.pdf", sep="." ) )
      td =  get.time.series (odb, areas, v, outfile=file.path(tempdir(), "ts.tmp.csv") )
    }
 
    td = td[ order(td$region, td$year) , ]
    td$region = factor(td$region, levels=areas, labels =regions)
    #   td[ which(td$region=="4X" & td$year < 2004), c("mean", "se", "ub", "lb", "n")] = NA

    #ylim=range(c(td$mean), na.rm=T); ylim[1]=ylim[1]-0.1*ylim[2]; ylim[2] = ylim[2]+ylim[2]*0.2
    xlim=range(td$year); xlim[1]=xlim[1]; xlim[2]=xlim[2]+1
    ylim='NULL'
    rg= range(td$mean, finite=TRUE)
    ylim[2]= rg[2]*1.02
    ylim[1]=rg[1] - rg[2]*0.02

    xlabels = seq(min(xlim), max(xlim), 1)
   # xlabels = c(xlim[1], xlim[1]+4, xlim[1]+8, xlim[1]+12, xlim[1]+16, xlim[1]+18) 
    ylabels = round(seq(as.numeric(rg[1]), round(as.numeric(rg[2]), 0), length.out=6),0)


    dir.create( outdir, recursive=T, showWarnings=F  )
    pdf(file=fn, width=5, height=8, bg='white')
   # Cairo( file=fn, type="png", bg="white",  units="in", width=5, height=8, dpi=350)
    setup.lattice.options()
    pl = xyplot( mean~year|region, data=td, ub=td$ub, lb=td$lb,
    #    layout=c(1,n.areas), xlim=xlim, scales = list(y = "free"),
        layout=c(1,n.areas), 
        par.strip.text = list(cex=1),
        par.settings=list(
        axis.text=list(cex=0.75), 
        par.main.text = list(cex=1),
        layout.heights=list(strip=1, panel=1, main=0.3 ),
        strip.background=list(col="lightgrey"), 
        box.umbrella = list(col = "black"), 
        box.rectangle = list(col="black"), 
        box.dot = list(col = "black", pch = 3, cex=2),
        plot.symbol=list(col='black', fill='darkgrey', cex=0.75, pch=21)
        ),
        #xlim=xlim,
        ylim=(c(as.numeric(ylim[1]), as.numeric(ylim[2]))),
        scales = list(y =list(at=ylabels, labels=ylabels), x=list(at=xlabels, labels=xlabels, rot=45)),
        main="Carapace width of the fishable biomass", xlab="Year", ylab="Carapace width (mm)",
        cex.lab=cex.lab, cex.axis=cex.axis, cex.main = cex.main,
        panel = function(x, y, subscripts, ub, lb, ...) {
           larrows(x, lb[subscripts], x, ub[subscripts], angle = 90, code = 3, length=0.05)
           panel.xyplot(x, y, type="b", lty=1, lwd=1.5, pch=21, fill='darkgrey', col="black", ...)
           #panel.abline(h=0, col="gray75", ...)
       }
    )
   
    print(pl) 
    dev.off()
   # cmd( "convert   -trim -quality 9  -geometry 200% -frame 2% -mattecolor white -antialias ", paste(fn, "pdf", sep="."),  paste(fn, "png", sep=".") )
     return("Done")
  }



