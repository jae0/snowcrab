
  habitat.usage = function( usevar, covariate, outdir ) {

    b = snowcrab.db("set.complete")
    logbook = logbook.db(DS="fisheries.complete")
    
    logbook$landings = logbook$landings / 1000  # convert kg to t
    L = logbook[ which(logbook$depth<500 & logbook$depth > 20),]

    if (usevar == "totno.all")     ylab=expression(paste("Number of all crab (1000 / ",km^2, ")"))
    if (usevar == "totno.male.com")   ylab=expression(paste("Number of commercial males (1000/ ",km^2, ")"))
    if (usevar == "totmass.male.com")   ylab=expression(paste("Mass of commercial males (kg / ",km^2, ")"))
    if (usevar == "totno.male.mat")     ylab=expression(paste("Number of mature males (1000 / ",km^2, ")"))
    if (usevar == "totno.male.imm")   ylab=expression(paste("Number of immature males (1000 / ",km^2, ")"))
    if (usevar == "totmass.male.mat")     ylab=expression(paste("Number of mature males (1000 / ",km^2, ")"))
    

    if (covariate=="depth") {
      # survey densities
      freq1 = errorbars( b=b, v=c("z", usevar), br=30, nfilter=10, lowess=0.25, trim=0.1,
          xlab="Depth (m)", ylab=ylab, db="snowcrab")
      Pr(dev="png", dname=outdir, fname=paste("survey.depth", usevar, sep=".") )
     
      # dockside landings
      usevar = "landings"
      freq1 = errorbars( b=L, v= c("depth", usevar), br=50, nfilter=10, lowess=0.25,  trim=0.1,
          xlab="Depth (m)", ylab="Landings (t)", db="snowcrab")
      Pr(dev="png", dname=outdir, fname=paste("landings.depth", usevar, sep="."))

      # logbook cpue
      usevar = "cpue"
      freq1 = errorbars( b=L, v= c("depth", usevar), br=50, nfilter=10, lowess=0.25,  trim=0.1,
          xlab="Depth (m)", ylab="CPUE (kg/trap haul)", db="snowcrab")
      Pr(dev="png", dname=outdir, fname=paste("cpue.depth", usevar, sep="."))
    }

    if (covariate=="temperature") {
      # survey densities
      freq1 = errorbars( b=b, v=c("t", usevar), br=50, nfilter=10, lowess=0.25, trim=0.1,
          xlab="Temperature (C)", ylab=ylab, db="snowcrab" )
      Pr(dev="png", dname=outdir, fname=paste("survey.temperatures", usevar, sep=".") )
    }

    if (covariate=="bottom.slope") {     
      # survey densities
      freq1 = errorbars( b=b, v=c("dZ", usevar), br=50, nfilter=10, lowess=0.25, trim=0.1,
          xlab="ln (abs( Slope; m/m) )", ylab=ylab, db="snowcrab")
      Pr(dev="png", dname=outdir, fname=paste("survey.slope", usevar, sep=".") )

      # dockside landings
      usevar = "landings"
      freq1 = errorbars( b=L, v= c("dZ", usevar), br=50, nfilter=10, lowess=0.25,  trim=0.1,
          xlab="ln (abs( Slope; m/m) )", ylab="Landings (t)", db="snowcrab" )
      Pr(dev="png", dname=outdir, fname=paste("landings.slope", usevar, sep="."))

      # logbook cpue
      usevar = "cpue"
      freq1 = errorbars( b=L, v= c("dZ", usevar), br=50, nfilter=10, lowess=0.25,  trim=0.1,
          xlab="ln (abs( Slope; m/m) )", ylab="CPUE (kg/trap haul)", db="snowcrab")
      Pr(dev="png", dname=outdir, fname=paste("cpue.slope", usevar, sep="."))
    }


    if (covariate=="bottom.curvature") {
      # survey densities
      freq1 = errorbars( b=b, v=c("ddZ", usevar), br=50, nfilter=10, lowess=0.25, trim=0.1,
          xlab="ln (abs ( Curvature; m/m2) ) )", ylab=ylab, db="snowcrab")
      Pr(dev="png", dname=outdir, fname=paste("survey.curvature", usevar, sep=".") )
    
      # dockside landings
      usevar = "landings"
      freq1 = errorbars( b=L, v= c("ddZ", usevar), br=50, nfilter=10, lowess=0.25,  trim=0.1,
          xlab="ln (abs ( Curvature; m/m2) ) )", ylab="Landings (t)", db="snowcrab")
      Pr(dev="png", dname=outdir, fname=paste("landings.curvature", usevar, sep="."))

      # logbook cpue
      usevar = "cpue"
      freq1 = errorbars( b=L, v= c("ddZ", usevar), br=30, nfilter=3, lowess=0.3,  trim=0.1,
          xlab="ln (abs ( Curvature; m/m2) ) )", ylab="CPUE (kg/trap haul)", db="snowcrab")
      Pr(dev="png", dname=outdir, fname=paste("cpue.curvature", usevar, sep="."))

    }

    if (covariate=="substrate") {
      # survey densities
      freq1 = errorbars( b=b, v=c("substrate.mean", usevar), br=50, nfilter=10, lowess=0.25, trim=0.1,
          xlab="ln ( substrate grainsize; mm )", ylab=ylab, db="snowcrab")
      Pr(dev="png", dname=outdir, fname=paste("survey.curvature", usevar, sep=".") )
    
      # dockside landings
      usevar = "landings"
      freq1 = errorbars( b=L, v= c("substrate.mean", usevar), br=50, nfilter=10, lowess=0.25,  trim=0.1,
          xlab="ln ( substrate grainsize; mm )", ylab="Landings (t)", db="snowcrab")
      Pr(dev="png", dname=outdir, fname=paste("landings.curvature", usevar, sep="."))

      # logbook cpue
      usevar = "cpue"
      freq1 = errorbars( b=L, v= c("substrate.mean", usevar), br=30, nfilter=3, lowess=0.3,  trim=0.1,
          xlab="ln ( substrate grainsize; mm )", ylab="CPUE (kg/trap haul)", db="snowcrab")
      Pr(dev="png", dname=outdir, fname=paste("cpue.curvature", usevar, sep="."))

 
    }

  }
   

