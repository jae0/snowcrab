snowcrab_interpolations_db = function (p, DS, voi=NULL) {

  # over-ride default dependent variable name if it exists
  if (is.null(voi)) if (exists("variables",p)) if(exists("Y", p$variables)) voi=p$variables$Y


  if (DS=="lbm_data") {
    oo = snowcrab.db( ...)


  }


  if (DS %in% c("prediction.surface", "prediction.surface.redo") ) {

    outdir = file.path( project.datadirectory("bio.indicators"), "PS", p$spatial.domain )
    dir.create(outdir, recursive=T, showWarnings=F)

    dyear_index = 1
    if (exists("dyears", p) & exists("prediction.dyear", p))  dyear_index = which.min( abs( p$prediction.dyear - p$dyears))

    outfile =  file.path( outdir, paste("PS", dyear_index, "rdata", sep=".") )

    if ( DS=="prediction.surface" ) {
      PS = NULL
      if (file.exists(outfile)) load( outfile )
      return (PS)
    }

    # this is the same as the p`rediction surface as used for indicators.db but the domain is smaller
    # .. more will be added to it belwo
    PS = indicators.db(p=p, DS="spatial")
    names(PS)[which(names(PS)=="tmean")] = "tmean.climatology"
    names(PS)[which(names(PS)=="tsd")] = "tsd.climatology"
    names(PS)[which(names(PS)=="tmin")] = "tmin.climatology"
    names(PS)[which(names(PS)=="tmax")] = "tmax.climatology"
    names(PS)[which(names(PS)=="amplitude")] = "tamplitude.climatology"

    nPS = nrow( PS )
    PS = as.list(PS)

    p0 = bio.temperature::temperature.parameters(p=p, current.year=p$current.year )
    p0 = bio.temperature::temperature.parameters( DS="lbm", p=p0 )
    p0 = bio.spacetime::spatial_parameters( p=p0, type=p$spatial.domain ) # return to correct domain

    yr_index = match( p$yrs, p0$tyears )
    u = indicators.db(p=p, DS="spatial.annual")
    for ( vn in names(u) ){
      u[[vn]] = u[[vn]][,yr_index]
    }
    PS = c( PS, u)

    # now we add the other covariate fields for modelling and prediction
    for ( vn in p$indicators.variables ) {
      u = bio.indicators::indicators.db( p=p, DS="baseline", voi=vn )

    }

    save (PS, file=outfile, compress=T )
    return( outfile )

  }



}
