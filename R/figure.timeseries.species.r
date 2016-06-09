
  figure.timeseries.survey = function(p, areas="", from.file=F) {
     
    set = snowcrab.db( DS ="set.logbook" ) 
     
    if ( areas=="small.areas" ) { # default
      areas = c( "cfa20", "cfa21", "cfa22", "cfa23", "cfa24", "cfa4x", 
      "cfa23slope", "cfa24slope")  # i.e., must be mutually exclusive categories
      for (r in areas) {  
        qr = filter.region.polygon(set, r)
        set$region[qr] = r
      }
    } else if (areas=="cfall" ) {
      set$region = factor("cfall", labels="Scotian.Shelf")
    } else {  # default
        set$region = factor( set$cfa, 
        levels = c("cfa4x", "cfasouth", "cfanorth" ), 
        labels = c("4X", "S-ENS", "N-ENS") 
      )
    }

    #na.omit()
    #variables =  c( variable.list.expand("all.to.model"), variable.list.expand("snowcrab.cw"), variable.list.expand("physical"),"landings", "cpue", "notraps", "effort" )
    #variables = c('RO.mass', 'R1.no', 'totno.female')
    #species = c(10, 30, 201, 50, 2521, 2511, 202, 204, 2211) 
    #cod, halibut, thornyskate, wolfish, lessertoadcrab, jonahcrab, smoothskate, winterskate, northernshrimp, 
    variables = c('ms.mass.10', 'ms.mass.30')

    names.set = names(set)
    for ( v in variables ) {
      cex=1
      main = list(v, cex=cex) 
  
      if (!( v %in% names.set))  next() 
      B = set[, c( v, "yr", "region" )]
      B[,1] = variable.recode( B[,1], v, direction="forward", db="snowcrab" )
      oo = which( B[,1] > 0 | is.finite( B[,1] ) )
      if (length(oo) < 10 ) next()

      ylab = list( paste("Geometric mean", v) , cex=cex)
      xlab = list("Year", cex=cex)

      fn = file.path( project.datadirectory("snowcrab"), "R", "timeseries", "survey", paste( v, "png", sep="." ) )
      dir.create( dirname(fn), recursive=T, showWarnings=F )
      ftype="png"
      if (grepl( "png$", fn)) ftype="png"
      if (grepl( "pdf$", fn)) ftype="pdf"
     Cairo( file=fn, type=ftype, bg="white", units="in", width=5, height=6.5, dpi=350)
        print (figure.timeseries.conditional(B, main=main, xlab=xlab, ylab=ylab ) )
      dev.off()
      #cmd( "convert -trim -frame 10x10 -mattecolor white ", fn, fn )
    }
    return("Done")
  }


