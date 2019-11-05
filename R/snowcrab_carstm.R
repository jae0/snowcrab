
snowcrab_carstm = function( p=NULL, DS="parameters", redo=FALSE, ...) {

	# handles all basic data tables, etc. ...

  # over-ride default dependent variable name if it exists

    if ( is.null(p)) {
      p = bio.snowcrab::snowcrab_parameters(...)
    } else {
      p = bio.snowcrab::snowcrab_parameters(p=p, ...)
    }


  # ---------------------
  # create/update library list
  p$libs = c( p$libs, RLibrary ( "colorspace",  "fields", "geosphere", "lubridate",  "lattice",
    "maps", "mapdata", "maptools", "parallel",  "rgdal", "rgeos",  "sp", "splancs", "GADMTools" ) )
  p$libs = c( p$libs, project.library ( "aegis", "aegis.bathymetry", "aegis.polygons", "aegis.coastline", "aegis.survey", "netmensuration", "bio.snowcrab", "carstm" ) )

  # p$taxa.of.interest = aegis.survey::groundfish.variablelist("catch.summary")
  # p$season = "summer"
  # p$taxa =  "maxresolved"


  # ---------------------

  if (DS =="parameters_override") {
    # translate param values from one project to a unified representation
    # must be first to catch p
    pc = snowcrab_carstm(
      DS = "parameters",
      project_class = "carstm", # defines which parameter class / set to load
      project_name = "snowcrab",
      variabletomodel = "totno",
      spatial_domain = p$spatial_domain,  # defines spatial area, currenty: "snowcrab" or "SSE"
      areal_units_overlay = p$areal_units_overlay, # currently: "snowcrab_managementareas",  "groundfish_strata" .. additional polygon layers for subsequent analysis for now ..
      areal_units_resolution_km = p$areal_units_resolution_km, # km dim of lattice ~ 1 hr
      areal_units_proj4string_planar_km = p$areal_units_proj4string_planar_km,  # coord system to use for areal estimation and gridding for carstm
      inputdata_spatial_discretization_planar_km = p$inputdata_spatial_discretization_planar_km,  # 1 km .. some thinning .. requires 32 GB RAM and limit of speed -- controls resolution of data prior to modelling to reduce data set and speed up modelling
      inputdata_temporal_discretization_yr = p$inputdata_temporal_discretization_yr,  # ie., weekly .. controls resolution of data prior to modelling to reduce data set and speed up modelling
      modeldir = p$modeldir,  # outputs all go the the main project's model output directory
      auid = p$auid
    )
    return(pc)
  }


  # ---------------------



  if (DS=="parameters") {


    if ( !exists("project_name", p)) p$project_name = "snowcrab"

    if ( !exists("areal_units_strata_type", p)) p$areal_units_strata_type = "lattice" # "stmv_lattice" to use ageis fields instead of carstm fields ... note variables are not the same

    if ( !exists("areal_units_overlay", p)) p$areal_units_overlay = "snowcrab_managementareas" # currently: "snowcrab_managementareas",  "groundfish_strata" .. additional polygon layers for subsequent analysis for now ..
    if ( !exists("areal_units_resolution_km", p)) p$areal_units_resolution_km = 25 # km dim of lattice ~ 1 hr
    if ( !exists("areal_units_proj4string_planar_km", p)) p$areal_units_proj4string_planar_km = projection_proj4string("utm20")  # coord system to use for areal estimation and gridding for carstm
    if ( !exists("areal_units_constrain_to_data", p)) p$areal_units_constrain_to_data = TRUE

    # if ( !exists("areal_units_proj4string_planar_km", p)) p$areal_units_proj4string_planar_km = projection_proj4string("omerc_nova_scotia")  # coord system to use for areal estimation and gridding for carstm
    p$inputdata_spatial_discretization_planar_km = p$pres  # 1 km .. requires 32 GB RAM and limit of speed -- controls resolution of data prior to modelling to reduce data set and speed up modelling
    p$inputdata_temporal_discretization_yr = 1/12  # ie., monthly .. controls resolution of data prior to modelling to reduce data set and speed up modelling }

    if ( !exists("carstm_modelengine", p)) p$carstm_modelengine = "inla.default"  # {model engine}.{label to use to store}

    if ( !exists("carstm_modelcall", p)) {
      if ( grepl("inla", p$carstm_modelengine) ) {
        p$libs = unique( c( p$libs, project.library ( "INLA" ) ) )

        # use strata which is a numeric representation of the StrataID factor
        p$carstm_model_label = "default_inla"
        p$carstm_modelcall = paste(
          'inla( formula =', p$variabletomodel,
          ' ~ 1
            + f(tiyr2, model="seasonal", season.length=10 )
            + f(ti, model="rw2", scale.model=TRUE, diagonal=1e-6, hyper=H$rw2)
            + f(zi, model="rw2", scale.model=TRUE, diagonal=1e-6, hyper=H$rw2)
            + f(gsi, model="rw2", scale.model=TRUE, diagonal=1e-6, hyper=H$rw2)
            + f(strata, model="bym2", graph=sppoly@nb, scale.model=TRUE, constr=TRUE, hyper=H$bym2)
            + f(iid_error, model="iid", hyper=H$iid),
            family = "normal",
            data= M,
            control.compute=list(dic=TRUE, config=TRUE),
            control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
            control.predictor=list(compute=FALSE, link=1 ),
            control.fixed=H$fixed,  # priors for fixed effects, generic is ok
            # control.inla=list(int.strategy="eb") ,# to get empirical Bayes results much faster.
            # control.inla=list( strategy="laplace", cutoff=1e-6, correct=TRUE, correct.verbose=FALSE ),
            num.threads=4,
            #blas.num.threads=4,
            verbose=TRUE
          )'
        )
      }
        #    + f(tiyr, model="ar1", hyper=H$ar1 )
        # + f(year,  model="ar1", hyper=H$ar1 )

      if ( grepl("glm", p$carstm_modelengine) ) {
        p$carstm_model_label = "default_glm"
        p$carstm_modelcall = paste(
          'glm( formula =',  p$variabletomodel,
          ' ~ 1 + StrataID + t + z + substrate.grainsize +tiyr,  # StrataID is a factor
            data= M[ which(M$tag=="observations"), ],
            family=gaussian(link="identity")
          )'
        )
      }

      if ( grepl("gam", p$carstm_modelengine) ) {
        p$libs = unique( c( p$libs, project.library ( "mgcv"  ) ) )

        p$carstm_model_label = "default_gam"
        p$carstm_modelcall = paste(
          'gam( formula =',  p$variabletomodel,
          ' ~ 1 + StrataID + s(t) + s(z) + s(substrate.grainsize) + s(yr) + s(dyear),
            data= M[ which(M$tag=="observations"), ],
            family=gaussian(link="identity")
          )'
        )
      }
    }

    p = carstm_parameters( p=p, DS="basic" )  # fill in anything missing and some checks

  #  boundingbox = list( xlim = c(-70.5, -56.5), ylim=c(39.5, 47.5)), # bounding box for plots using spplot

    return(p)
  }

  # -------------------------

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

    pB = bathymetry_carstm( p=p, DS="parameters_override" ) # transcribes relevant parts of p to load bathymetry
    BI = carstm_model ( p=pB, DS="carstm_modelled" )  # unmodeled!

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

    APS = carstm_model( p=p$p_bathymetry, DS="carstm_modelled" )
    APS2 = carstm_model( p=p$p_snowcrab, DS="carstm_modelled" )



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
   vn = c( p$variabletomodel, pB$variabletomodel,  ps$variabletomodel,  pt$variabletomodel, "tag", "StrataID" )


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


      M$strata  = as.numeric( M$StrataID)
      M$iid_error = 1:nrow(M) # for inla indexing for set level variation

      M$zi = discretize_data( M[, pB$variabletomodel], p$discretization$z )

      M$tiyr  = trunc( M$tiyr / p$tres )*p$tres    # discretize for inla .. midpoints
      M$year = floor(M$tiyr)
      M$dyear  =  factor( as.character( trunc(  (M$tiyr - M$year )/ p$tres )*p$tres), levels=p$dyears)


    save( M, file=fn, compress=TRUE )
    return( M )
  }


}  ## end snowcrab.db
