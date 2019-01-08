
  figure.timeseries.bycatch = function(p, outdir, species,  plotyears, type='mass', all.areas=T, graphic='pdf' ) {
 #browser()
    outdir=file.path(p$annual.results, "timeseries", "survey")

    if (all.areas) {
      areas = c("cfa4x", "cfasouth", "cfanorth" )
      regions = c("4X", "S-ENS", "N-ENS")
    } else {
      areas = c("cfasouth", "cfanorth" )
      regions = c("S-ENS", "N-ENS")
    }

    n.regions = length(regions)
    n.areas = length(areas)

    spcd = groundfish.db( DS="spcodes")
    tdb = snowcrab.timeseries.db( DS="biologicals", p=p )
    if(missing(species)){
      cat = snowcrab.db( DS="cat.initial" )
      species = unique(cat$spec)
    }
    v = paste("ms",type,species,sep='.')
    spcd = subset(spcd,code%in%species)
    if(missing(plotyears))plotyears = unique(tdb$year)



    for(i in 1:length(v)){

      td = tdb[ which( tdb$variable == v[i] & tdb$year%in%plotyears) ,]
      td = td[ order(td$region, td$year) , ]
      td$region = factor(td$region, levels=areas, labels =regions)
      #   td[ which(td$region=="4X" & td$year < 2004), c("mean", "se", "ub", "lb", "n")] = NA


      common = spcd$comm[which(spcd$code==species[i])]

      ylim='NULL'
      ylim[2]=max(td$ub,na.rm=T)*1.05
      ylim[1]=0
      xlim=range(td$year); xlim[1]=xlim[1]; xlim[2]=xlim[2]

      xlabels = seq(min(xlim), max(xlim), 1)
      ylabels = pretty(ylim,8)

      dir.create( outdir, recursive=T, showWarnings=F )
      fn = file.path( outdir, paste( v[i], graphic,  sep="." ) )
      cex.main = 1.4
      cex.lab = 1
      cex.axis = 0.2

   if(graphic=='png') Cairo( file=fn, type="png", bg="white",  units="in", width=5, height=6, dpi=350 )
   if(graphic=='pdf') pdf(file=fn, width=5, height=6, bg='white')
   if(graphic=='R') plot.new()
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
              scales=list(y=list(at=ylabels, labels=ylabels, cex=0.65), x=list(at=xlabels, labels=xlabels, rot=50, cex=0.65)),
                main=paste(capwords(common,T,F),"Biomass"), xlab=list("Year", cex=1), ylab=list("Geometric mean kg / km^2", cex=1),
                #cex.lab=cex.lab,
                cex.axis=cex.axis,
                cex.main = cex.main,
                panel = function(x, y, subscripts, ub, lb, ...) {
               larrows(x, lb[subscripts], x, ub[subscripts], angle = 90, code = 3, length=0.05)
               panel.xyplot(x, y, type="b", lty=1, lwd=1.5, pch=21, fill='darkgrey', col="black", ...)
               panel.abline(h=median(y), col="gray", ...)
          }
        )

      print(pl)
      if(graphic!='R') dev.off()
      print(fn)
    }
   #cmd( "convert   -trim -quality 9  -geometry 200% -frame 2% -mattecolor white -antialias ", paste(fn, "pdf", sep="."),  paste(fn, "png", sep=".") )
     
    return("Done")
    }



