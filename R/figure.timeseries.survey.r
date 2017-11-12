
  figure.timeseries.survey = function(outdir, variables, plotyears,type="biologicals", all.areas=T, minN=10, u=NULL, graphic='pdf',...) {


    if (all.areas) {
      areas = c("cfa4x", "cfasouth", "cfanorth" )
      regions = c("4X", "S-ENS", "N-ENS")
    } else {
      areas = c("cfasouth", "cfanorth" )
      regions = c("S-ENS", "N-ENS")
    }

    n.regions = length(regions)
    n.areas = length(areas)

    # base data
    tdb = snowcrab.timeseries.db( DS=type )

    if(missing(variables)){
      variables =  c( 
         emaf::variable.list.expand("all.to.model"), 
         emaf::variable.list.expand("snowcrab.cw"), 
         emaf::variable.list.expand("physical"),
         'sexratio.mat','sexratio.imm','sexratio.all' 
      )
      variables = intersect( variables, unique(tdb$variable))
    }

    # base data
    tdb = snowcrab.timeseries.db( DS=type )

    if(missing(plotyears))plotyears = unique(tdb$year)
    
    tdb = subset(tdb,variable%in%variables&year%in%plotyears)
    tdb$region = factor(tdb$region, levels=areas, labels =regions)

    #  load transformation tables associated with a given variable
    tl = emaf::lookup.datatransformation('snowcrab')

    if (file.exists( tl$repository) ) {
      load (tl$repository)
    } else {
      REPOS = emaf::recode.variable.initiate.db ( db )
    }
    tvars = REPOS$varname[which(REPOS$transform=='log10')]

    for ( v in variables ) {
      td = tdb[ which( tdb$variable == v) ,] 
      
      oo = which( td$mean > 0 | is.finite( td$mean ) )
      if (length(oo) < minN ) next()

      if(is.null(u))u=NULL # for variable specific units, needs a lookup table 

      ylim=c(0,max(c(td$ub,td$mean),na.rm=T))
      #browser()
      if(length(grep('ratio',v))==1)ylim=c(0,1)

      xlim=range(td$year)
      if(v %in% tvars){
        ylab = list( paste("Geometric mean", u) , cex=1)
      } else {
        ylab = list( paste("Mean", u) , cex=1)
        ylim[1] = min(c(td$lb,td$mean),na.rm=T)
      }
      xlab = list("Year", cex=1)
      ylim[1] = ylim[1]-diff(ylim)*0.04
      ylim[2] = ylim[2]+diff(ylim)*0.04

      main = capwords(gsub("."," ",v,fixed=T),F,F)

      xlabels = seq(xlim[1], xlim[2], 1)
      ylabels = pretty(ylim,7)

      dir.create( outdir, recursive=T, showWarnings=F )

      fn = file.path( outdir, paste( v, graphic,  sep="." ) )
      if(type=='groundfish.t'){
        fn = file.path( outdir, paste( type, graphic,  sep="." ) )
        main = "Groundfish Survey Temperature"
        xlabels = seq(xlim[1], xlim[2], 2)
      }
      dline = ifelse(length(grep('ratio',v))==1,0.5,NA)
      if(graphic=='png')Cairo( file=fn, type="png", bg="white",  units="in",dpi=350,... )
      if(graphic=='pdf')pdf(file=fn, bg='white', ...)
      if(graphic=='R')plot.new()
      setup.lattice.options()
      pl = xyplot( mean~year|region, data=td, ub=td$ub, lb=td$lb, dline=dline,
            layout=c(1,n.regions),
            par.strip.text=list(
              plot.symbol=list(col='black', fill='darkgrey', cex=0.75, pch=21),
              axis.text=list(cex=0.7),
              par.main.text=list(cex=1),
              layout.heights=list(strip=1, panel=1, main=0.5),
              strip.background=list(col='lightgrey')),
              #xlim=xlim,
              ylim=ylim,
              scales=list(y=list(at=ylabels, labels=ylabels, cex=0.65), x=list(at=xlabels, labels=xlabels, rot=50, cex=0.65)),
                main=main, xlab=xlab, ylab=ylab, 
                cex.axis=0.2,
                cex.main = 1.4,
                panel = function(x, y, subscripts, ub, lb, ...) {
               larrows(x, lb[subscripts], x, ub[subscripts], angle = 90, code = 3, length=0.05)
               panel.abline(h=median(y,na.rm=T), col="gray", ...)
               panel.abline(h=dline, col="gray", lty=2,...)
               panel.xyplot(x, y, type="b", lty=1, lwd=1.5, pch=21, fill='darkgrey', col="black", ...)
          }
        )

      print(pl)
      if(graphic!='R')dev.off()
    }
    return("Done")
  }


