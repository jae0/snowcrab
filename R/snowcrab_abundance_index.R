snowcrab_abundance_index = function( p=NULL, RES=NULL, filename=NULL, operation="load_RES", ... ) {

  # require areal_units_fn,

  # deal with additional passed parameters
  p_add = list(...)
  if ( is.null(p) ) p=list()
  if (length(p_add) > 0 ) p = c(p, p_add)
  i = which(duplicated(names(p), fromLast = TRUE ) )
  if ( length(i) > 0 ) p = p[-i] # give any passed parameters a higher priority, overwriting pre-existing variable

  # same file naming as in carstm ..

  if ( is.null(filename) ) {
    fn = filename
  } else {
    areal_units_fns = p$areal_units_fn
    if (exists( "inputdata_spatial_discretization_planar_km", p )) areal_units_fns = paste( areal_units_fns, round(p$inputdata_spatial_discretization_planar_km, 6),   sep="_" )
    if (exists( "inputdata_temporal_discretization_yr", p )) areal_units_fns = paste( areal_units_fns, round(p$inputdata_temporal_discretization_yr, 6),   sep="_" )
    areal_units_fns_suffix = paste( areal_units_fns, p$variabletomodel, p$carstm_modelengine, "aggregated_timeseries",  "rdata", sep="." )
    outputdir = file.path(p$modeldir, p$carstm_model_label)
    if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )
    fn = file.path( outputdir, paste("carstm_modelled_results", areal_units_fns_suffix, sep="." ) )
  }

  if (is.null(RES)) {
    save( RES, file=fn, compress=TRUE )
    return(fn)
  }

  if ( operation="load_timeseries" ) {
    RES = NA
    if (file.exists(fn)) load( fn)
    return( RES )
  }

  if ( operation="load_spatial" ) {
    out = NA
    if (file.exists(fn_out)) load( fn_out )
    return( out )
  }

  if (operation=="compute") {
    # construct meanweights matrix used to convert number to weight
    M = snowcrab_carstm( p=p, DS="carstm_inputs" )
    M$yr = M$year  # req for meanweights
    sppoly = areal_units( p=p )
    weight_year = meanweights_by_arealunit(
      set=M[M$tag=="observations",],
      AUID=as.character( sppoly$AUID ),
      yrs=p$yrs,
      fillall=TRUE,
      annual_breakdown=TRUE
    )
    # convert numerical density to total number and convert to biomass:  / 10^6  # 10^6 kg -> kt .. kg/km * km
    out = res[[ paste( p$variabletomodel, "predicted", sep=".")]]
    out[!is.finite(out)] = NA
    out[out > 1e10] = NA
    save( out, file=fn_out, compress=TRUE )

    RES = data.frame( yrs = p$yrs )
    RES$cfaall    = colSums( out * weight_year * sppoly$au_sa_km2/ 10^6, na.rm=TRUE )
    RES$cfanorth  = colSums( out * weight_year * sppoly$cfanorth_surfacearea/ 10^6, na.rm=TRUE )
    RES$cfasouth  = colSums( out * weight_year * sppoly$cfasouth_surfacearea/ 10^6, na.rm=TRUE )
    RES$cfa23     = colSums( out * weight_year * sppoly$cfa23_surfacearea/ 10^6, na.rm=TRUE )
    RES$cfa24     = colSums( out * weight_year * sppoly$cfa24_surfacearea/ 10^6, na.rm=TRUE )
    RES$cfa4x     = colSums( out * weight_year * sppoly$cfa4x_surfacearea/ 10^6, na.rm=TRUE
    save( RES, file=fn, compress=TRUE )
    return( fn )
  }

}

