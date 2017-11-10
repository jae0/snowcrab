 map.survey.locations = function(p, basedir, newyear=T, map.method="lattice" ) {


    set = snowcrab.db( DS="set.clean")
    years = sort( unique( set$yr ) )
    if (newyear) years = p$year.assessment

    if (map.method=="lattice" ) {

      set = set[, c("yr", "plon", "plat")]
      set = set[ is.finite( rowSums(set) ) ,]
      corners = data.frame(rbind( cbind( plon=c(220, 990), plat=c(4750, 5270) )))

      for (y in years) {
        toplot = set[ which(set$yr==y), c("plon", "plat")]
        annot = paste ("Survey locations", y)
        fn = paste("survey.locations", y, sep=".")
        print(fn)
        emei::emei_map( toplot, cfa.regions=T, depthcontours=T, annot=annot, fn=fn, loc=basedir, corners=corners )

      }
    }

 
    return ("Done" )
  }


