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
        annot = paste (y)
        fn = paste("survey.locations", y, sep=".")
        print(fn)
        dir.create (basedir, showWarnings=FALSE, recursive =TRUE)
        png( filename=file.path(basedir, paste(fn, "png", sep=".")), width=3072, height=2304, pointsize=40, res=300 )
        print(
          aegis::aegis_map( xyz=toplot, cfa.regions=TRUE, depthcontours=TRUE, annot=annot, annot.cex=2.8, corners=corners)
        )
        dev.off()
      }
    }

 
    return ("Done" )
  }


