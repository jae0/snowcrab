
  map.observer.locations = function(p, basedir, newyear=T , map.method="GMT" ) {
    
    odb = observer.db( DS="odb")
    odb$yr = odb$fishyr  # use fishyr and not the real year ###################
    years = sort( unique( odb$yr ) )
    if (newyear) years = p$current.assessment.year
     
    odb = odb[, c("yr", "lon", "lat")]
    odb = odb[ is.finite( rowSums(odb) ) ,]

    if ( map.method=="GMT" ) {
      
      p$psxyoptions = "-Sc0.1c -G20"  # Sc = circle with size 0.1cm, G is color/grayscale
      p$basedir =  basedir
       
      for (y in years) { 
        ii =  which(odb$yr==y)
        if ( length(ii)  < 10)  next()
        toplot = odb[ii, c("lon", "lat")]
        p$outfile.basename = file.path(p$basedir, paste("observer.locations", y, sep=".") )
        gmt.xyplot ( p, toplot, y, conversions=p$conversions )
      }
      pause(30)
      files.to.delete =  list.files( p$basedir, "^observer.locations.*.ps$", all.files=T, full.names=T, recursive=T )
      remove.files ( files.to.delete ) 
    }

    if (map.method=="lattice" ) {
 
      for (y in years) {
        ii =  which(odb$yr==y)
        if ( length(ii)  < 10 ) next()
        toplot = odb[ ii, c("lon", "lat")]
        annot = paste ("Observer locations", y)
        fn = paste("observer.locations", y, sep=".") 
        print(fn)
        map( xyz=toplot,  cfa.regions=T, depthcontours=T, annot=annot, fn=fn, loc=basedir, corners=planar.corners )
       
      }
     
    }

    return ("Done" )


  }


