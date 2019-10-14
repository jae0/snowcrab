
snowcrab_carstm = function( p=NULL, DS=NULL, sppoly=NULL, redo=FALSE, ... ) {

  #\\ Note inverted convention: depths are positive valued
  #\\ i.e., negative valued for above sea level and positive valued for below sea level

  # deal with additional passed parameters
  if ( is.null(p) ) p=list()
  p_add = list(...)
  if (length(p_add) > 0 ) p = c(p, p_add)
  i = which(duplicated(names(p), fromLast = TRUE ))
  if ( length(i) > 0 ) p = p[-i] # give any passed parameters a higher priority, overwriting pre-existing variable


  # -----------------------


  if ( DS=="carstm_inputs") {

    fn = file.path( p$modeldir, paste( "snowcrab", "carstm_inputs", p$auid,
      p$inputdata_spatial_discretization_planar_km,
      round(p$inputdata_temporal_discretization_yr, 6),
      "rdata", sep=".") )

    if (!redo)  {
      if (file.exists(fn)) {
        load( fn)
        return( M )
      }
    }
    message( "Generating carstm_inputs ... ")

    # prediction surface
    sppoly = areal_units( p=p )  # will redo if not found
    crs_lonlat = sp::CRS(projection_proj4string("lonlat_wgs84"))


    # do this immediately to reduce storage for sppoly (before adding other variables)
    M = snowcrab.db( p=p, DS="biological_data", redo=TRUE )  # will redo if not found .. not used here but used for data matching/lookup in other aegis projects that use bathymetry

    # reduce size
    M = M[ which( M$lon > p$corners$lon[1] & M$lon < p$corners$lon[2]  & M$lat > p$corners$lat[1] & M$lat < p$corners$lat[2] ), ]
    # levelplot(z.mean~plon+plat, data=M, aspect="iso")

    M$StrataID = over( SpatialPoints( M[, c("lon", "lat")], crs_lonlat ), spTransform(sppoly, crs_lonlat ) )$StrataID # match each datum to an area
    M$lon = NULL
    M$lat = NULL
    M$plon = NULL
    M$plat = NULL
    M = M[ which(is.finite(M$StrataID)),]
    M$StrataID = as.character( M$StrataID )  # match each datum to an area

    M$snowcrab = M$snowcrab.mean
    M$tag = "observations"

    APS = as.data.frame(sppoly)
    APS$StrataID = as.character( APS$StrataID )
    APS$tag ="predictions"
    APS$t = NA
    APS$z = NA

    pb = aegis.bathymetry::bathymetry_parameters( p=p, project_class =="carstm_auid" ) # transcribes relevant parts of p to load bathymetry
    BI = bathymetry_carstm ( p=pb, DS="carstm_modelled" )  # unmodeled!

    jj = match( as.character( APS$StrataID), as.character( BI$StrataID) )
    APS$z = BI$z.predicted[jj]
    jj =NULL
    BI = NULL

    vn = c("t", "tag", "StrataID", "z")
    APS = APS[, vn]

    # expand APS to all time slices
    n_aps = nrow(APS)
    APS = cbind( APS[ rep.int(1:n_aps, p$nt), ], rep.int( p$prediction_ts, rep(n_aps, p$nt )) )
    names(APS) = c(vn, "tiyr")

    M$snowcrab = M$snowcrab.mean
    M$tiyr = M$yr + M$dyear
    M = rbind( M[, names(APS)], APS )
    APS = NULL

    M$StrataID  = factor( as.character(M$StrataID), levels=levels( sppoly$StrataID ) ) # revert to factors
    sppoly = NULL


  # Get/create predictions surface

    APS = bathymetry_carstm( p=p, sppoly=sppoly, DS="carstm_modelled", redo=FALSE )

    APS = snowcrab_carstm( p=p, sppoly=sppoly, DS="carstm_modelled", redo=FALSE )



  ## ----------------------------------
  # covariate estimates for prediction in strata and year
  # collapse PS vars with time into APS (and regrid via raster)
    APS = aegis_db_extract(
      vars=p$lookupvars,
      yrs=p$yrs,
      spatial_domain=p$spatial_domain,
      dyear=p$prediction_dyear,
      areal_units_resolution_km=p$areal_units_resolution_km,
      aegis_proj4string_planar_km=sp::CRS(p$aegis_proj4string_planar_km),
      returntype="data.frame",
      redo = FALSE
    )

    APS$yr = as.numeric( APS$year)
    APS[,p$variabletomodel] = NA
    APS$data_offset = 1  # force to be density n/km^2
    APS$tag = "predictions"


    # StrataID reset to be consistent in both data and prediction areal units
    o = over( SpatialPoints( set[,c("lon", "lat")], sp::CRS(projection_proj4string("lonlat_wgs84")) ), spTransform(sppoly, sp::CRS(projection_proj4string("lonlat_wgs84")) ) ) # match each datum to an area
    set$StrataID = o$StrataID


    o = over( SpatialPoints( APS[,c("lon", "lat")], sp::CRS(projection_proj4string("lonlat_wgs84")) ), spTransform(sppoly, sp::CRS(projection_proj4string("lonlat_wgs84")) ) ) # match each datum to an area
    APS$StrataID = o$StrataID

    o = NULL

    #  good data
    ok = which(
      is.finite(set[,p$variabletomodel]) &   # INLA can impute Y-data
      is.finite(set$data_offset) &
      is.finite(set$StrataID)
    )


  # construct meanweights matrix
  weight_year = meanweights_by_strata( set=set, StrataID=as.character( sppoly$StrataID ), yrs=p$yrs, fillall=TRUE, annual_breakdown=TRUE )
  # weight_year = weight_year[, match(as.character(p$yrs), colnames(weight_year) )]
  # weight_year = weight_year[ match(as.character(sppoly$StrataID), rownames(weight_year) )]


  varstokeep = unique( c( p$variabletomodel, "StrataID", "yr", "data_offset", "tag", p$lookupvars) )

  M = rbind( set[ok, varstokeep], APS[,varstokeep] )

  M = M[ which(
        is.finite(M$data_offset) &
        is.finite(M$StrataID)
      ) , ]

  M$yr_factor = factor( as.character(M$yr) )
  M$StrataID  = factor( as.character(M$StrataID), levels=levels( sppoly$StrataID ) )
  M$strata  = as.numeric( M$StrataID)
  M$year  = as.numeric( M$yr_factor)
  M$iid_error = 1:nrow(M) # for inla indexing for set level variation


  M$t[!is.finite(M$t)] = median(M$t, na.rm=TRUE )  # missing data .. quick fix .. do something better
  M$tsd[!is.finite(M$tsd)] = median(M$tsd, na.rm=TRUE )  # missing data .. quick fix .. do something better
  M$tmin[!is.finite(M$tmin)] = median(M$tmin, na.rm=TRUE )  # missing data .. quick fix .. do something better
  M$tmax[!is.finite(M$tmax)] = median(M$tmax, na.rm=TRUE )  # missing data .. quick fix .. do something better
  M$degreedays[!is.finite(M$degreedays)] = median(M$degreedays, na.rm=TRUE )  # missing data .. quick fix .. do something better
  M$z[!is.finite(M$z)] = median(M$z, na.rm=TRUE )  # missing data .. quick fix .. do something better
  M$dZ[!is.finite(M$dZ)] = median(M$dZ, na.rm=TRUE )  # missing data .. quick fix .. do something better
  M$ddZ[!is.finite(M$ddZ)] = median(M$ddZ, na.rm=TRUE )  # missing data .. quick fix .. do something better
  M$substrate.grainsize[!is.finite(M$substrate.grainsize)] = median(M$substrate.grainsize, na.rm=TRUE )  # missing data .. quick fix .. do something better


  M$ti = discretize_data( M$t, p$discretization$t )
  M$tisd = discretize_data( M$tsd, p$discretization$tsd )
  M$timin = discretize_data( M$tmin, p$discretization$tmin )
  M$timax = discretize_data( M$tmax, p$discretization$tmax )
  M$di = discretize_data( M$t, p$discretization$degreedays )
  M$zi = discretize_data( M$t, p$discretization$z )
  M$zid = discretize_data( M$t, p$discretization$dZ )
  M$zidd = discretize_data( M$t, p$discretization$ddZ )
  M$si = discretize_data( M$t, p$discretization$substrate.grainsize )


  totest = setdiff(1:ncol(M), which(names(M) %in% c( p$variabletomodel , "StrataID", "tag", "yr_factor") ))
  ii = which(is.finite(rowSums(M[,totest])))

  M = M[ii,]

  # ---------------------
  # generic PC priors
  m = log( {set[,p$variabletomodel] / set$data_offset}[ok] )
  m[!is.finite(m)] = min(m[is.finite(m)])





    save( M, file=fn, compress=TRUE )
    return( M )
  }


  # ------------



  if ( DS=="carstm_modelled") {

    auids = paste(  p$auid, p$inputdata_spatial_discretization_planar_km,
      round(p$inputdata_temporal_discretization_yr, 6),   sep="_" )

    fn = file.path( p$modeldir, paste("snowcrab", "carstm_modelled", p$carstm_modelengine, auids, "rdata", sep="." ) )
    fn_fit = file.path( p$modeldir, paste( "snowcrab", "carstm_modelled_fit", p$carstm_modelengine, auids,  "rdata", sep=".") )

    if (!redo)  {
         if (file.exists(fn)) {
        load( fn)
        return( res )
      }
      if (DS=="carstm_modelled_fit") {
        if (file.exists(fn_fit)) {
          load( fn_fit )
          return( fit )
        }
      }
    }

    # prediction surface
    sppoly = areal_units( p=p )  # will redo if not found
    res = list(StrataID = sppoly[["StrataID"]])  # init results list
    res$strata = as.numeric(res$StrataID)

    M = snowcrab_carstm ( p=p, DS="carstm_inputs" )  # 16 GB in RAM just to store!
    M$strata  = as.numeric( M$StrataID)

    fit  = NULL

    if ( grepl("glm", p$carstm_modelengine) ) {
      assign("fit", eval(parse(text=paste( "try(", p$carstm_modelcall, ")" ) ) ))
      if (is.null(fit)) warning("model fit error")
      if ("try-error" %in% class(fit) ) warning("model fit error")
      save( fit, file=fn_fit, compress=TRUE )
      ii = which( M$tag=="predictions" & M$strata %in% M[ which(M$tag=="observations"), "strata"] )
      jj = match( M$strata[ii], res$strata )
      preds = predict( fit, newdata=M[ii,], type="link", na.action=na.omit, se.fit=TRUE )  # no/km2
      res[,"snowcrab.predicted"] = exp( preds$fit[jj])
      res[,"snowcrab.predicted_se"] = exp( preds$se.fit[jj])
      res[,"snowcrab.predicted_lb"] = exp( preds$fit[jj] - preds$se.fit[jj] )
      res[,"snowcrab.predicted_ub"] = exp( preds$fit[jj] + preds$se.fit[jj] )
      save( res, file=fn, compress=TRUE )
    }

    if ( grepl("gam", p$carstm_modelengine) ) {
      assign("fit", eval(parse(text=paste( "try(", p$carstm_modelcall, ")" ) ) ))
      if (is.null(fit)) warning("model fit error")
      if ("try-error" %in% class(fit) ) warning("model fit error")
      save( fit, file=fn_fit, compress=TRUE )
      ii = which( M$tag=="predictions" & M$strata %in% M[ which(M$tag=="observations"), "strata"] )
      jj = match( M$strata[ii], res$strata )
      preds = predict( fit, newdata=M[ii,], type="link", na.action=na.omit, se.fit=TRUE )  # no/km2
      res[,"snowcrab.predicted"] = exp( preds$fit[jj])
      res[,"snowcrab.predicted_se"] = exp( preds$se.fit[jj])
      res[,"snowcrab.predicted_lb"] = exp( preds$fit[jj] - preds$se.fit[jj] )
      res[,"snowcrab.predicted_ub"] = exp( preds$fit[jj] + preds$se.fit[jj] )
      save( res, file=fn, compress=TRUE )
    }


    if ( grepl("inla", p$carstm_modelengine) ) {
      H = carstm_hyperparameters( sd(M$snowcrab, na.rm=TRUE), alpha=0.5, median( M$snowcrab, na.rm=TRUE) )

      M$zi = discretize_data( M$z, p$discretization$z )
      M$tiyr2 = M$tiyr  # use a copy for "seasonal" models
      M$year = floor(M$tiyr)
      M$dyear  =  M$tiyr - M$year
      M$iid_error = 1:nrow(M) # for inla indexing for set level variation

      assign("fit", eval(parse(text=paste( "try(", p$carstm_modelcall, ")" ) ) ))
      if (is.null(fit)) warning("model fit error")
      if ("try-error" %in% class(fit) ) warning("model fit error")
      save( fit, file=fn_fit, compress=TRUE )

      # reformat predictions into matrix form
      ii = which(
        M$tag=="predictions" &
        M$strata %in% res$strata &
        M$year %in% p$yrs
      )  # filter by strata and years in case additional data in other areas and times are used in the input data

      # res[ res>1e10] = NA

      matchfrom = list( strata=M$strata[ii], year=as.character(M$year[ii]), dyear=M$dyear[ii] )
      matchto   = list( strata=res$strata, year=as.character(p$yrs), dyear=factor(p$dyears) )

      res$snowcrab.predicted = reformat_to_array(
        input = fit$summary.fitted.values[ ii, "mean" ],
        matchfrom=matchfrom, matchto=matchto
      )

      res$snowcrab.predicted_lb = reformat_to_array(
        input = fit$summary.fitted.values[ ii, "0.025quant" ],
        matchfrom=matchfrom, matchto=matchto
      )

      res$snowcrab.predicted_ub = reformat_to_array(
        input =  fit$summary.fitted.values[ ii, "0.975quant" ],
        matchfrom=matchfrom, matchto=matchto
      )

      # random effects results ..
      if (exists("summary.random", fit)) {

        nstrata = length(res$StrataID)

        # reformat predictions into matrix form
        ii = which(
          M$tag=="predictions" &
          M$strata %in% res$strata &
          M$year %in% p$yrs
        )  # filter by strata and years in case additional data in other areas and times are used in the input data


        if (exists("iid_error", fit$summary.random)) {
          # IID random effects
          matchfrom = list( strata=M$strata[ii], year=M$year[ii], dyear=M$dyear[ii] )
          matchto   = list( strata=res$strata, year=p$yrs, dyear=factor(p$dyears) )
          res$snowcrab.random_sample_iid = reformat_to_array(
            input =  fit$summary.random$iid_error[ ii, "mean" ],
            matchfrom=matchfrom, matchto=matchto
          )
          # carstm_plot( p=p, res=res, vn="snowcrab.random_sample_iid", time_match=list(year="1950", dyear="0") )
        }

        if (exists("strata", fit$summary.random)) {

          if (nrow(fit$summary.random$strata) == nstrata*2) {
            # CAR random effects can be of variable length depending upon model construct:

            # a single spatial and nonspatial effect (no grouping across time)
            jj = 1:nstrata
            matchfrom = list( strata=fit$summary.random$strata$ID[jj]  )
            matchto   = list( strata=res$strata  )
            res$snowcrab.random_strata_nonspatial = reformat_to_array(
              fit$summary.random$strata[ jj, "mean" ],
              matchfrom=matchfrom, matchto=matchto
            )
            res$snowcrab.random_strata_spatial =reformat_to_array(
              fit$summary.random$strata[ jj+nstrata, "mean" ],
              matchfrom=matchfrom, matchto=matchto
            )
            # carstm_plot( p=p, res=res, vn="snowcrab.random_strata_nonspatial"  )
            # carstm_plot( p=p, res=res, vn="snowcrab.random_strata_spatial" )

          } else if (nrow(fit$summary.random$strata) == nstrata*2 * p$ny ) {
            # spatial and nonspatial effects grouped by year
            matchfrom = list( strata=M$strata[ii], year=M$year[ii] )
            matchto   = list( strata=res$strata, year=p$yrs )

            res$snowcrab.random_strata_nonspatial = reformat_to_array(
              input =  fit$summary.random$strata[ ii, "mean" ],
              matchfrom = matchfrom, matchto = matchto
            )
            res$snowcrab.random_strata_spatial = reformat_to_array(
              input = fit$summary.random$strata[ ii+max(ii), "mean" ],
              matchfrom = matchfrom, matchto = matchto
            )
            # carstm_plot( p=p, res=res, vn="snowcrab.random_strata_nonspatial", time_match=list(year="2000" ) )
            # carstm_plot( p=p, res=res, vn="snowcrab.random_strata_spatial", time_match=list(year="2000" ) )

          } else if (nrow(fit$summary.random$strata) == nstrata*2 * p$nt ) {

            # need to test/fix ...
            matchfrom = list( StrataID=M$StrataID[ii], year=M$year[ii], dyear=M$dyear[ii] ),
            matchto   = list( StrataID=res$StrataID, year=p$yrs, dyear=factor(p$dyears) )

            res$snowcrab.random_strata_nonspatial = reformat_to_array(
              input = fit$summary.random$strata[ jj, "mean" ],
              matchfrom = matchfrom,  matchto   = matchto
            )
            res$snowcrab.random_strata_spatial = reformat_to_array(
              input = fit$summary.random$strata[ ii+max(ii), "mean" ],
              matchfrom = matchfrom, matchto   = matchto
            )

            # carstm_plot( p=p, res=res, vn="snowcrab.random_strata_nonspatial", time_match=list(year="2000", dyear="0.8" ) )
            # carstm_plot( p=p, res=res, vn="snowcrab.random_strata_spatial", time_match=list(year="2000", dyear="0.8" ) )

          }
        }
      save( res, file=fn, compress=TRUE )
    }


    return( res )
  }

}
