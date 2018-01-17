
  map.logbook.locations = function(p, basedir, newyear=T, map.method="lattice"  ) {

    x = logbook.db( DS="logbook" )
    x = x[aegis::polygon_inside(x, region="isobath1000m"),]
    years = sort( unique( x$yr ) )
    if (newyear) years = p$year.assessment
    x = x[, c("yr", "lon", "lat")]
    x = x[ is.finite( rowSums(x) ) ,]

    if (map.method=="lattice" ) {
      for (y in years) {
        ii =  which(x$yr==y)
        if ( length(ii)  < 10 ) next()
        toplot = x[ ii, c("lon", "lat")]
        annot = paste ("Logbook locations", y)
        outfn = paste("logbook.locations", y, sep=".")
        print(outfn)        
        dir.create (basedir, showWarnings=FALSE, recursive =TRUE)
        fn = file.path( basedir, paste(outfn, "png", sep="." ) )
        png( filename=fn, width=3072, height=2304, pointsize=40, res=300 )
        lp = aegis::aegis_map( toplot, cfa.regions=TRUE, depthcontours=TRUE, annot=annot, corners=p$corners )
        print(lp)
        dev.off()
      }
    }
  }



