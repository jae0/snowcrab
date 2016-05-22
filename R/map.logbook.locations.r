
  map.logbook.locations = function(p, basedir, newyear=T, map.method="GMT"  ) {
    
    x = logbook.db( DS="logbook" )
    x = x[filter.region.polygon(x, region="isobath1000m"),]
    years = sort( unique( x$yr ) )
    if (newyear) years = p$current.assessment.year
    x = x[, c("yr", "lon", "lat")]
    x = x[ is.finite( rowSums(x) ) ,]
    
    if ( map.method=="GMT" ) {
      p$psxyoptions = "-Sc0.1c -G20"  # Sc = circle with size x cm, G is color/grayscale
      p$basedir =  basedir
      
      for (y in years) {
        ii =  which(x$yr==y)
        if ( length(ii)  < 10 ) next()
        toplot = x[ ii, c("lon", "lat")]
        p$outfile.basename = file.path(p$basedir, paste("logbook.locations", y, sep=".") )
        gmt.xyplot ( p, toplot, y, conversions=p$conversions )
      }

      pause(30)
      files.to.delete = list.files( p$basedir, "^logbook.locations.*.ps$", all.files=T, full.names=T, recursive=T )
      remove.files (files.to.delete) 
    }
    if (map.method=="lattice" ) {
      
      for (y in years) {
        ii =  which(x$yr==y)
        if ( length(ii)  < 10 ) next()
        toplot = x[ ii, c("lon", "lat")]
        annot = paste ("Logbook locations", y)
        fn = paste("logbook.locations", y, sep=".") 
        print(fn)
        map( toplot, cfa.regions=T, depthcontours=T, annot=annot, fn=fn, loc=basedir, corners=planar.corners )
      }
     
    }


  }



