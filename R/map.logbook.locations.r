
  map.logbook.locations = function(p, basedir, newyear=T, map.method="lattice"  ) {

    x = logbook.db( DS="logbook" )
    x = x[emgis::polygon_inside(x, region="isobath1000m"),]
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
        fn = paste("logbook.locations", y, sep=".")
        print(fn)
        emei::emei_map( toplot, cfa.regions=T, depthcontours=T, annot=annot, fn=fn, loc=basedir, corners=p$corners )
      }

    }


  }



