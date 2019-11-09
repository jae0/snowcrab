
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
      project_name = "snowcrab",
      yrs = p$yrs,
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
        p$carstm_model_label = "production"
        p$carstm_modelcall = paste(
          'inla( formula =', p$variabletomodel,
          ' ~ 1
            + offset( log(data_offset))
            + f(year_factor, model="ar1", hyper=H$ar1 )
            + f(dyri, model="rw2", scale.model=TRUE, diagonal=1e-6, hyper=H$rw2 )
            + f(ti, model="rw2", scale.model=TRUE, diagonal=1e-6, hyper=H$rw2)
            + f(zi, model="rw2", scale.model=TRUE, diagonal=1e-6, hyper=H$rw2)
            + f(gsi, model="rw2", scale.model=TRUE, diagonal=1e-6, hyper=H$rw2)
            + f(pca1i, model="rw2", scale.model=TRUE, diagonal=1e-6, hyper=H$rw2)
            + f(pca2i, model="rw2", scale.model=TRUE, diagonal=1e-6, hyper=H$rw2)
            + f(strata, model="bym2", graph=sppoly@nb, group=year_factor, scale.model=TRUE, constr=TRUE, hyper=H$bym2),
            family = "poisson",
            data= M,
            control.compute=list(dic=TRUE, config=TRUE),
            control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
            control.predictor=list(compute=FALSE, link=1 ),
            control.fixed=H$fixed,  # priors for fixed effects, generic is ok
            # control.inla=list(int.strategy="eb") ,# to get empirical Bayes results much faster.
            # control.inla=list( strategy="laplace", cutoff=1e-6, correct=TRUE, correct.verbose=FALSE ),
            # num.threads=4,
            # blas.num.threads=4,
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
          ' ~ 1 + StrataID + s(t) + s(z) + s(substrate.grainsize) + s(year) + s(dyear),
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
      p$variabletomodel,
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
    M = snowcrab.db( p=p, DS="biological_data" )  # will redo if not found .. not used here but used for data matching/lookup in other aegis projects that use bathymetry

    # globally remove all unrealistic data
    # p$quantile_bounds_data = c(0.0005, 0.9995)
    if (exists("quantile_bounds_data", p)) {
      TR = quantile(M[,p$variabletomodel], probs=p$quantile_bounds_data, na.rm=TRUE ) # this was -1.7, 21.8 in 2015
      keep = which( M[,p$variabletomodel] >=  TR[1] & M[,p$variabletomodel] <=  TR[2] )
      if (length(keep) > 0 ) M = M[ keep, ]
        # this was -1.7, 21.8 in 2015
    }

    # reduce size
    M = M[ which( M$lon > p$corners$lon[1] & M$lon < p$corners$lon[2]  & M$lat > p$corners$lat[1] & M$lat < p$corners$lat[2] ), ]
    # levelplot(z.mean~plon+plat, data=M, aspect="iso")

    M$StrataID = over( SpatialPoints( M[, c("lon", "lat")], crs_lonlat ), spTransform(sppoly, crs_lonlat ) )$StrataID # match each datum to an area

    names(M)[which(names(M)=="yr") ] = "year"
    # M = M[ which(M$year %in% p$yrs), ]
    # M$tiyr = lubridate::decimal_date ( M$timestamp )
    # M$dyear = M$tiyr - M$year


    # tailored specifically for snow crab resolution -- override defaults
    p$discretization = list()
    p$discretization$z = c(0, 5, 10, 20, 40, 60, 80, 100, 150, 200, 250, 300, 350, 400, 500, 1000, 2000, 5000 )  # depth cut points
    p$discretization$substrate.grainsize = c( 0, 1, 2, 4, 6, 8, 10, 12, 16, 24, 32 )
    p$discretization$pca1 = c( -1, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 1 )
    p$discretization$pca2 = c( -1, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 1 )
    p$discretization$t = c( -4, -2, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 28 )
    p$discretization$dyear = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
    p$n.season = length(p$discretization[["dyear"]]) - 1


    pB = bathymetry_carstm( p=p, DS="parameters_override" ) # transcribes relevant parts of p to load bathymetry
    pS = substrate_carstm( p=p, DS="parameters_override" ) # transcribes relevant parts of p to load bathymetry
    pT = temperature_carstm( p=p, DS="parameters_override" ) # transcribes relevant parts of p to load
    pPC1 = speciescomposition_carstm( p=p, DS="parameters_override", varnametomodel="pca1" ) # transcribes relevant parts of p to load

    pPC2 = speciescomposition_carstm( p=p, DS="parameters_override", varnametomodel="pca2" ) # transcribes relevant parts of p to load

    if (!(exists(pB$variabletomodel, M ))) M[,pB$variabletomodel] = NA
    if (!(exists(pS$variabletomodel, M ))) M[,pS$variabletomodel] = NA
    if (!(exists(pT$variabletomodel, M ))) M[,pT$variabletomodel] = NA
    if (!(exists(pPC1$variabletomodel, M ))) M[,pPC1$variabletomodel] = NA
    if (!(exists(pPC2$variabletomodel, M ))) M[,pPC2$variabletomodel] = NA


    kk = which(!is.finite(M[, pB$variabletomodel]))
    if (length(kk) > 0 ) M[kk, pB$variabletomodel] = lookup_bathymetry_from_surveys( p=pB, locs=M[kk, c("lon", "lat")] )

    kk = which(!is.finite(M[, pS$variabletomodel]))
    if (length(kk) > 0 ) M[kk, pS$variabletomodel] = lookup_substrate_from_surveys(  p=pS, locs=M[kk, c("lon", "lat")] )

    kk = which(!is.finite(M[, pT$variabletomodel]))
    if (length(kk) > 0 ) M[kk, pT$variabletomodel] = lookup_temperature_from_surveys(  p=pT, locs=M[kk, c("lon", "lat")], timestamp=M$timestamp )

    kk = which(!is.finite(M[, pPC1$variabletomodel]))
    if (length(kk) > 0 ) {
      pc1 = speciescomposition.db( p=pPC1, DS="speciescomposition"  )
      ii = match( M$id, pc1$id)
      M[kk, pPC1$variabletomodel] = pc1[ii, pPC1$variabletomodel]
    }

    kk = which(!is.finite(M[, pPC2$variabletomodel]))
    if (length(kk) > 0 ) {
      pc2 = speciescomposition.db( p=pPC2, DS="speciescomposition"  )
      ii = match( M$id, pc2$id)
      M[kk, pPC2$variabletomodel] = pc2[ii, pPC2$variabletomodel]
    }


    # if any still missing then use a mean depth by StrataID
    kk =  which( !is.finite(M[, pB$variabletomodel]))
    if (length(kk) > 0) {
      AD = bathymetry.db ( p=pB, DS="aggregated_data"  )  # 16 GB in RAM just to store!
      AD = AD[ which( AD$lon > p$corners$lon[1] & AD$lon < p$corners$lon[2]  & AD$lat > p$corners$lat[1] & AD$lat < p$corners$lat[2] ), ]
      # levelplot( eval(paste(p$variabletomodel, "mean", sep="."))~plon+plat, data=M, aspect="iso")
      AD$StrataID = over( SpatialPoints( AD[, c("lon", "lat")], crs_lonlat ), spTransform(sppoly, crs_lonlat ) )$StrataID # match each datum to an area
      oo = tapply( AD[, paste(pB$variabletomodel, "mean", sep="." )], AD$StrataID, FUN=median, na.rm=TRUE )
      jj = match( as.character( M$StrataID[kk]), as.character( names(oo )) )
      M[kk, pB$variabletomodel] = oo[jj ]
    }

    # if any still missing then use a mean substrate by StrataID
    kk =  which( !is.finite(M[, pS$variabletomodel]))
    if (length(kk) > 0) {
      AD = substrate.db ( p=pS, DS="aggregated_data"  )  # 16 GB in RAM just to store!
      AD = AD[ which( AD$lon > p$corners$lon[1] & AD$lon < p$corners$lon[2]  & AD$lat > p$corners$lat[1] & AD$lat < p$corners$lat[2] ), ]
      # levelplot( eval(paste(p$variabletomodel, "mean", sep="."))~plon+plat, data=M, aspect="iso")
      AD$StrataID = over( SpatialPoints( AD[, c("lon", "lat")], crs_lonlat ), spTransform(sppoly, crs_lonlat ) )$StrataID # match each datum to an area
      oo = tapply( AD[, paste(pS$variabletomodel, "mean", sep="." )], AD$StrataID, FUN=median, na.rm=TRUE )
      jj = match( as.character( M$StrataID[kk]), as.character( names(oo )) )
      M[kk, pS$variabletomodel] = oo[jj ]
    }

    # substrate coverage poor .. add from modelled results
    kk =  which( !is.finite(M[, pS$variabletomodel]))
    if (length(kk) > 0) {
      SI = carstm_model ( p=pS, DS="carstm_modelled" )
      jj = match( as.character( M$StrataID[kk]), as.character( SI$StrataID) )
      M[kk, pS$variabletomodel] = SI[[ paste(pS$variabletomodel,"predicted",sep="." )]] [jj]
    }


    # if any still missing then use a mean temp by StrataID
    kk =  which( !is.finite(M[, pT$variabletomodel]))
    if (length(kk) > 0) {
      AD = temperature.db ( p=pT, DS="aggregated_data"  )  # 16 GB in RAM just to store!
      AD = AD[ which( AD$lon > p$corners$lon[1] & AD$lon < p$corners$lon[2]  & AD$lat > p$corners$lat[1] & AD$lat < p$corners$lat[2] ), ]
      # levelplot( eval(paste(p$variabletomodel, "mean", sep="."))~plon+plat, data=M, aspect="iso")

      AD$StrataID = over( SpatialPoints( AD[, c("lon", "lat")], crs_lonlat ), spTransform(sppoly, crs_lonlat ) )$StrataID # match each datum to an area
      AD$uid = paste(AD$StrataID, AD$year, AD$dyear, sep=".")

      M_dyear_discret = discretize_data( M$dyear, p$discretization$dyear )  # AD$dyear is discretized. . match discretization
      M$uid =  paste(M$StrataID, M$year, M_dyear_discret, sep=".")

      oo = tapply( AD[, paste(pT$variabletomodel, "mean", sep="." )], AD$uid, FUN=median, na.rm=TRUE )

      jj = match( as.character( M$uid[kk]), as.character( names(oo )) )
      M[kk, pT$variabletomodel] = oo[jj ]
    }




    kk =  which( !is.finite(M[, pPC1$variabletomodel]))
    if (length(kk) > 0) {
      pc1 = speciescomposition.db( p=pPC1, DS="speciescomposition"  )
      pc1 = planar2lonlat( pc1, proj.type=p$aegis_proj4string_planar_km )
      pc1$StrataID = over( SpatialPoints( pc1[, c("lon", "lat")], crs_lonlat ), spTransform(sppoly, crs_lonlat ) )$StrataID # match each datum to an area
      oo = tapply( pc1[, pPC1$variabletomodel ], pc1$StrataID, FUN=median, na.rm=TRUE )
      jj = match( as.character( M$StrataID[kk]), as.character( names(oo )) )
      if (length(jj) > 0) M[kk, pPC1$variabletomodel] = oo[jj ]
    }
    kk =  which( !is.finite(M[, pPC1$variabletomodel]))
    if (length(kk) > 0) {
      PI = carstm_model ( p=pPC1, DS="carstm_modelled" )
      strata_map = match( as.numeric(M$StrataID[kk]), levels(sppoly$StrataID[as.numeric(dimnames(PI)$strata)]  ) )
      year_map = match( as.character(M$year[kk]), dimnames(PI)$year )
      dindex = cbind(strata_map, year_map )
      M[kk, pPC1$variabletomodel] = PI [dindex]
    }



    kk =  which( !is.finite(M[, pPC2$variabletomodel]))
    if (length(kk) > 0) {
      pc2 = speciescomposition.db( p=pPC2, DS="speciescomposition"  )
      pc2 = planar2lonlat( pc2, proj.type=p$aegis_proj4string_planar_km )
      pc2$StrataID = over( SpatialPoints( pc2[, c("lon", "lat")], crs_lonlat ), spTransform(sppoly, crs_lonlat ) )$StrataID # match each datum to an area
      oo = tapply( pc2[, pPC2$variabletomodel ], pc2$StrataID, FUN=median, na.rm=TRUE )
      jj = match( as.character( M$StrataID[kk]), as.character( names(oo )) )
      if (length(jj) > 0) M[kk, pPC2$variabletomodel] = oo[jj ]
    }
    kk =  which( !is.finite(M[, pPC2$variabletomodel]))
    if (length(kk) > 0) {
      PI = carstm_model ( p=pPC2, DS="carstm_modelled" )
      strata_map = match( as.numeric(M$StrataID[kk]), levels(sppoly$StrataID[as.numeric(dimnames(PI)$strata)]  ) )
      year_map = match( as.character(M$year[kk]), dimnames(PI)$year )
      dindex = cbind(strata_map, year_map )
      M[kk, pPC2$variabletomodel] = PI [dindex]
    }


    PI = NULL
    M$plon = NULL
    M$plat = NULL
    M$lon = NULL
    M$lat = NULL

    M = M[ which(is.finite(M[, pB$variabletomodel] )),]
    M = M[ which(is.finite(M[, pS$variabletomodel] )),]
    M = M[ which(is.finite(M[, pT$variabletomodel] )),]
    M = M[ which(is.finite(M$StrataID)),]
    M$StrataID = as.character( M$StrataID )  # match each datum to an area

    M$tag = "observations"


    APS = as.data.frame(sppoly)
    APS$StrataID = as.character( APS$StrataID )
    APS$tag ="predictions"
    APS[,p$variabletomodel] = NA


    BI = carstm_model ( p=pB, DS="carstm_modelled" )
    jj = match( as.character( APS$StrataID), as.character( BI$StrataID) )
    APS[, pB$variabletomodel] = BI[[ paste(pB$variabletomodel,"predicted",sep="." ) ]] [jj]
    jj =NULL
    BI = NULL

    SI = carstm_model ( p=pS, DS="carstm_modelled" )
    jj = match( as.character( APS$StrataID), as.character( SI$StrataID) )
    APS[, pS$variabletomodel] = SI[[ paste(pS$variabletomodel,"predicted",sep="." )]] [jj]
    jj =NULL
    SI = NULL

    # to this point APS is static, now add time dynamics (teperature)
    # ---------------------

    vn = c( p$variabletomodel, pB$variabletomodel,  pS$variabletomodel, "tag", "StrataID" )
    APS = APS[, vn]

    # expand APS to all time slices
    n_aps = nrow(APS)
    APS = cbind( APS[ rep.int(1:n_aps, p$nt), ], rep.int( p$prediction_ts, rep(n_aps, p$nt )) )
    names(APS) = c(vn, "tiyr")
    APS$year = floor( APS$tiyr)
    APS$dyear = APS$tiyr - APS$year


    TI = carstm_model ( p=pT, DS="carstm_modelled" )
    TI = TI[[ paste(pT$variabletomodel,"predicted",sep="." )]]
    strata_map = match( as.numeric(APS$StrataID),levels(sppoly$StrataID[as.numeric(dimnames(TI)$strata)]  ) )
    year_map = match( as.character(APS$year), dimnames(TI)$year )
    dyear_breaks = c(p$dyears, p$dyears[length(p$dyears)]+ diff(p$dyears)[1] )
    dyear_map = as.numeric( cut( APS$dyear, breaks=dyear_breaks, include.lowest=TRUE, ordered_result=TRUE, right=FALSE ) )
    dindex = cbind(strata_map, year_map, dyear_map )
    APS[, pT$variabletomodel] = TI[ dindex]
    TI = NULL


    PI = carstm_model ( p=pPC1, DS="carstm_modelled" )
    PI = PI[[ paste(pPC1$variabletomodel,"predicted",sep="." )]]
    strata_map = match( as.numeric(APS$StrataID), levels(sppoly$StrataID[as.numeric(dimnames(PI)$strata)]  ) )
    year_map = match( as.character(APS$year), dimnames(PI)$year )
    dindex = cbind(strata_map, year_map )
    APS[, pPC1$variabletomodel] = PI [dindex]
    PI = NULL

    PI = carstm_model ( p=pPC2, DS="carstm_modelled" )
    PI = PI[[ paste(pPC2$variabletomodel,"predicted",sep="." )]]
    strata_map = match( as.numeric(APS$StrataID), levels(sppoly$StrataID[as.numeric(dimnames(PI)$strata)]  ) )
    year_map = match( as.character(APS$year), dimnames(PI)$year )
    dindex = cbind(strata_map, year_map  )
    APS[, pPC2$variabletomodel] = PI [dindex]
    PI = NULL

    # useful vars to have for analyses outside of carstm_model
    varstoadd = c( "totwgt", "totno", "sa", "data_offset",  "zn", "qn" )

    for (vn in varstoadd) if (!exists( vn, APS)) APS[,vn] = NA
    APS$data_offset = 1  # force to solve for unit area

    M = rbind( M[, names(APS)], APS )
    APS = NULL


    M$StrataID  = factor( as.character(M$StrataID), levels=levels( sppoly$StrataID ) ) # revert to factors
    M$strata  = as.numeric( M$StrataID)

    M$zi  = discretize_data( M[, pB$variabletomodel], p$discretization[[pB$variabletomodel]] )
    M$ti  = discretize_data( M[, pT$variabletomodel], p$discretization[[pT$variabletomodel]] )
    M$gsi = discretize_data( M[, pS$variabletomodel], p$discretization[[pS$variabletomodel]] )

    M$pca1i = discretize_data( M[, pPC1$variabletomodel], p$discretization[[pPC1$variabletomodel]] )
    M$pca2i = discretize_data( M[, pPC2$variabletomodel], p$discretization[[pPC2$variabletomodel]] )

    M$tiyr  = trunc( M$tiyr / p$tres )*p$tres    # discretize for inla .. midpoints

    M$year = trunc( M$tiyr)
    M$year_factor = as.numeric( factor( M$year, levels=p$yrs))
    M$dyear =  M$tiyr - M$year   # revert dyear to non-discretized form

    M$dyri = discretize_data( M[, "dyear"], p$discretization[["dyear"]] )

    # M$seasonal = (as.numeric(M$year_factor) - 1) * length(p$dyears)  + as.numeric(M$dyear)


    save( M, file=fn, compress=TRUE )
    return( M )
  }


}  ## end snowcrab.db
