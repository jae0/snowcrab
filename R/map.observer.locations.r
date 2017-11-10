
  map.observer.locations = function(p, basedir, newyear=T , map.method="lattice" ) {

    odb = observer.db( DS="odb")
    odb$yr = odb$fishyr  # use fishyr and not the real year ###################
    years = sort( unique( odb$yr ) )
    if (newyear) years = p$year.assessment

    odb = odb[, c("yr", "lon", "lat")]
    odb = odb[ is.finite( rowSums(odb) ) ,]

    if (map.method=="lattice" ) {

      corners = data.frame(rbind( cbind( plon=c(220, 990), plat=c(4750, 5270) )))
      for (y in years) {
        ii =  which(odb$yr==y)
        if ( length(ii)  < 10 ) next()
        toplot = odb[ ii, c("lon", "lat")]
        annot = paste ("Observer locations", y)
        fn = paste("observer.locations", y, sep=".")
        print(fn)
        emei::emei_map( xyz=toplot,  cfa.regions=T, depthcontours=T, annot=annot, fn=fn, loc=basedir, corners=corners )

      }

    }

    return ("Done" )


  }


