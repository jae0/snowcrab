snowcrab_abundance_index = function( p=NULL, operation="load_RES", ... ) {

  # require areal_units_fn,

  # deal with additional passed parameters
  p_add = list(...)
  if ( is.null(p) ) p=list()
  if (length(p_add) > 0 ) p = c(p, p_add)
  i = which(duplicated(names(p), fromLast = TRUE ) )
  if ( length(i) > 0 ) p = p[-i] # give any passed parameters a higher priority, overwriting pre-existing variable


  required.vars = c("areal_units_fn", "inputdata_spatial_discretization_planar_km", "inputdata_temporal_discretization_yr", "variabletomodel",
    "carstm_modelengine", "modeldir", "carstm_model_label", "yrs", "variabletomodel" )

  for (i in required.vars) {
    if (!exists(i, p)) {
      message( "Missing parameter" )
      message( i )
      stop()
    }
  }


  # same file naming as in carstm ..
  outputdir = file.path(p$modeldir, p$carstm_model_label)
  areal_units_fns = p$areal_units_fn
  if (exists( "inputdata_spatial_discretization_planar_km", p )) areal_units_fns = paste( areal_units_fns, round(p$inputdata_spatial_discretization_planar_km, 6),   sep="_" )
  if (exists( "inputdata_temporal_discretization_yr", p )) areal_units_fns = paste( areal_units_fns, round(p$inputdata_temporal_discretization_yr, 6),   sep="_" )
  if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )
  fn     = file.path( outputdir, paste("carstm_modelled_results", paste( areal_units_fns, p$variabletomodel, p$carstm_modelengine, "aggregated_timeseries",  "rdata", sep="." ), sep="." ) )
  fn_no = file.path( outputdir, paste("carstm_modelled_results", paste( areal_units_fns, p$variabletomodel, p$carstm_modelengine, "space_timeseries_number",  "rdata", sep="." ), sep="." ) )
  fn_bio = file.path( outputdir, paste("carstm_modelled_results", paste( areal_units_fns, p$variabletomodel, p$carstm_modelengine, "space_timeseries_biomass",  "rdata", sep="." ), sep="." ) )
  fn_wgts = file.path( outputdir, paste("carstm_modelled_results", paste( areal_units_fns, p$variabletomodel, p$carstm_modelengine, "space_timeseries_weights",  "rdata", sep="." ), sep="." ) )

  if ( operation=="load_timeseries" ) {
    RES = NA
    if (file.exists(fn)) load( fn)
    return( RES )
  }

  if ( operation=="load_spacetime_number" ) {
    nums = NA
    if (file.exists(fn_no)) load( fn_no )
    return( nums )
  }

  if ( operation=="load_spacetime_biomass" ) {
    biom = NA
    if (file.exists(fn_bio)) load( fn_bio )
    return( biom )
  }


  if ( operation=="load_spacetime_weights" ) {
    weight_year = NA
    if (file.exists(fn_wgts)) load( fn_wgts )
    return( weight_year )
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
      annual_breakdown=TRUE,
      robustify_quantiles=c(0, 0.95)  # high upper bounds are more dangerous
    )

    save (weight_year, file=fn_wgts, compress=TRUE)

    res = carstm_model( p=p, DS="carstm_modelled", carstm_model_label=p$carstm_model_label ) # to load currently saved res

    # convert numerical density to total number and convert to biomass:  / 10^6  # 10^6 kg -> kt .. kg/km * km
    nums = res[[ paste( p$variabletomodel, "predicted", sep=".")]]
    nums[!is.finite(nums)] = NA
    qnt = quantiles( nums, probs=0.95, na.rm=TRUE)

    nums[nums > qnt] = qnt
    save( nums, file=fn_no, compress=TRUE )

    biom = nums * weight_year
    save( biom, file=fn_bio, compress=TRUE )



    RES = data.frame( yrs = p$yrs )
    RES$cfaall    = colSums( biom * sppoly$au_sa_km2/ 10^6, na.rm=TRUE )
    RES$cfanorth  = colSums( biom * sppoly$cfanorth_surfacearea/ 10^6, na.rm=TRUE )
    RES$cfasouth  = colSums( biom * sppoly$cfasouth_surfacearea/ 10^6, na.rm=TRUE )
    RES$cfa23     = colSums( biom * sppoly$cfa23_surfacearea/ 10^6, na.rm=TRUE )
    RES$cfa24     = colSums( biom * sppoly$cfa24_surfacearea/ 10^6, na.rm=TRUE )
    RES$cfa4x     = colSums( biom * sppoly$cfa4x_surfacearea/ 10^6, na.rm=TRUE )
    save( RES, file=fn, compress=TRUE )
    return( fn )
  }

}

