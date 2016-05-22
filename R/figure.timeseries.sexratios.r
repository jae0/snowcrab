
  figure.timeseries.sexratios = function( outdir, all.areas=T, type="mature" ) {
 #browser()
    set = snowcrab.db("set.merge.det")
    #set2015 = set[which(set$yr == 2015),]
    #print(head(set2015))
  
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

    if (type=="mature") {
      v = "sexratio.mat"
      outfile =  paste(v, ".pdf", sep="")
      figure.title = "Sex ratios -- mature"
      td =  get.time.series ( from.file=T )
      td = td[ which( td$variable == v) ,]
      # td = td[ - which(td$region=="cfa4x" ),]
      td = td[ order(td$region, td$year) , ]
      td$region = factor(td$region, levels=areas,labels=regions)
      #   td[ which(td$region=="4X" & td$year < 2004), c("mean", "se", "ub", "lb", "n")] = NA
      #ylim=c(0,1)
      ylim= 'NULL'
      ylim[2]=max(is.finite(td$ub))
      ylim[1]=round(min(is.finite(td$lb)), 0) - (max(is.finite(td$ub))*0.05)
      xlim=range(td$year); xlim[1]=xlim[1]; xlim[2]=xlim[2]+1
      # ylim=range(c(td$mean), na.rm=T); ylim[1]=ylim[1]-0.2*ylim[2]; ylim[2] = ylim[2]+ylim[2]*0.2
            #xlim=range(td$year); xlim[1]=xlim[1]-0.5; xlim[2]=xlim[2]+0.5
      xlabels = seq(min(xlim), max(xlim), 2)
      ylabels = seq(0, as.numeric(ylim[2]), length.out=6)
    }

    if (type=="immature") {
      v = "sexratio.imm"
      outfile = paste (v, ".pdf", sep="")
      figure.title = "Sex ratios -- immature"
      td =  get.time.series ( from.file=T )
      td = td[ which( td$variable == v) ,]
      td = td[ order(td$region, td$year) , ]
      td$region = factor(td$region, levels=areas, labels =regions )
      #   td[ which(td$region=="4X" & td$year < 2004), c("mean", "se", "ub", "lb", "n")] = NA
      # ylim=c(0,1)
      # ylim=range(c(td$mean), na.rm=T); ylim[1]=ylim[1]-0.2*ylim[2]; ylim[2] = ylim[2]+ylim[2]*0.2
      #xlim=range(td$year); xlim[1]=xlim[1]-0.5; xlim[2]=xlim[2]+0.5
      ylim= 'NULL'
      ylim[2]=max(is.finite(td$ub))
      ylim[1]=round(min(is.finite(td$lb)), 0) - (max(is.finite(td$ub))*0.05)
      xlim=range(td$year); xlim[1]=xlim[1]; xlim[2]=xlim[2]+1
      # ylim=range(c(td$mean), na.rm=T); ylim[1]=ylim[1]-0.2*ylim[2]; ylim[2] = ylim[2]+ylim[2]*0.2
            #xlim=range(td$year); xlim[1]=xlim[1]-0.5; xlim[2]=xlim[2]+0.5
      xlabels = seq(min(xlim), max(xlim), 2)
      ylabels = seq(0, as.numeric(ylim[2]), length.out=6)

    }

    dir.create( outdir, recursive=T, showWarnings=F  )
    fn = file.path( outdir, outfile )
    pdf(file=fn, width=4, height=6, bg='white')
    #Cairo( file=fn, type="png", bg="white",  units="in", width=4, height=6, dpi=350)
    setup.lattice.options()
    pl = xyplot( mean~year|region, data=td, ub=td$ub, lb=td$lb,
          layout=c(1,n.regions), 
          par.strip.text = list(cex=1.0),
          par.settings=list(
            plot.symbol=list(col='black', fill='darkgrey', cex=0.75, pch=21),
            axis.text=list(cex=0.7),
            par.main.text = list(cex=1),
            layout.heights=list(strip=1, panel=1, main=0.5),
            strip.background=list(col="lightgrey")),
              ylim=(c(as.numeric(ylim[1]), as.numeric(ylim[2]))),
              scales=list(y=list(at=ylabels, labels=ylabels), x=list(at=xlabels, labels=xlabels, rot=45)),
              main=figure.title, xlab=list("Year", cex=1), ylab=list("Proportion female", cex=1),
              cex.axis=cex.axis, cex.main = cex.main,
              panel = function(x, y, subscripts, ub, lb, ...) {
                larrows(x, lb[subscripts], x, ub[subscripts], angle = 90, code = 3, length=0.05)
                panel.xyplot(x, y, type="b", lty=1, lwd=1.5, pch=21, fill='darkgrey', col="black", ...)
                #panel.abline(h=0.5, col="gray75", ...)
         }
      )
    
    print(pl)
    dev.off()
   #cmd( "convert   -trim -quality 9  -geometry 200% -frame 2% -mattecolor white -antialias ", paste(fn, "pdf", sep="."),  paste(fn, "png", sep=".") )
     return("Done")
  }