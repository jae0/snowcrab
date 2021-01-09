
snowcrab_carstm = function( p=NULL, DS="parameters", redo=FALSE, extrapolation_limit=NA, extrapolation_replacement=NA, ...) {

	# handles all basic data tables, etc. ...



  if (DS=="parameters") {

    # over-ride default dependent variable name if it exists

      if ( is.null(p)) {
        p = bio.snowcrab::snowcrab_parameters(...)
      } else {

        # if p is passed, assume it is a secondary call ... overwrite nonrelevent params to force use of project defaults
        p$data_root = NULL
        p$datadir  = NULL
        p$data_transformation= NULL
        p$carstm_model_call = NULL  # defaults to generic

        p = bio.snowcrab::snowcrab_parameters(p=p, ...)
      }

    # ---------------------
    # create/update library list2
    p$libs = unique( c( p$libs, RLibrary (
      "colorspace",  "fields", "geosphere", "lubridate",  "lattice",
      "maps", "mapdata", "maptools", "parallel",  "rgdal", "rgeos",
      "sp", "sf", "splancs", "GADMTools",
      "aegis",
      "aegis.bathymetry",
      "aegis.polygons",
      "aegis.coastline",
      "aegis.survey",
      "aegis.temperature",
      "aegis.substrate",
      "aegis.speciescomposition",
      "netmensuration",
      "bio.snowcrab",
      "carstm"
    ) ) )


    # p$taxa.of.interest = aegis.survey::groundfish_variablelist("catch.summary")
    # p$season = "summer"
    # p$taxa =  "maxresolved"

    if ( !exists("assessment.years", p)) stop( "must define assessment.years")

    p = parameters_add_without_overwriting( p,
      groundfish_species_code=2526,
      speciesname = "Snow crab", 
      spatial_domain = "snowcrab", 
      yrs = p$assessment.years,
      inputdata_spatial_discretization_planar_km = 1 ,  # 1 km .. some thinning .. requires 32 GB RAM and limit of speed -- controls resolution of data prior to modelling to reduce data set and speed up modelling
      inputdata_temporal_discretization_yr = 1/12,
      trawlable_units = "sweptarea",  # <<<<<<<<<<<<<<<<<< also:  "standardtow", "sweptarea" (for groundfish surveys)
      areal_units_type = "lattice", # "stmv_lattice" to use ageis fields instead of carstm fields ... note variables are not the same
      areal_units_overlay = "snowcrab_managementareas", # currently: "snowcrab_managementareas",  "groundfish_strata" .. additional polygon layers for subsequent analysis for now ..
      areal_units_resolution_km = 25, # km dim of lattice ~ 1 hr
      areal_units_constraint = "snowcrab",  # locations of data as constraint .. "snowcrab" loads these automatically, otherwise a xy matrix of positions

      areal_units_constraint_nmin = 3, 

      areal_units_proj4string_planar_km = aegis::projection_proj4string("utm20"),  # coord system to use for areal estimation and gridding for carstm
      quantile_bounds =c(0, 0.99) # trim upper bounds
    )

    
    if ( !exists("selection", p)) p$selection=list()
    if ( !exists("type", p$selection) ) p$selection$type = "number"
    if ( !exists("biologicals", p$selection) ) p$selection$biologicals = list(
      spec_bio=bio.taxonomy::taxonomy.recode( from="spec", to="parsimonious", tolookup=p$groundfish_species_code ),
      sex=0, # male
      mat=1, # do not use maturity status in groundfish data as it is suspect ..
      len= c( 95, 200 )/10, #  mm -> cm ; aegis_db in cm
      ranged_data="len"
    )

    if ( !exists("survey", p$selection) ) p$selection$survey = list(
      data.source = ifelse (p$selection$type=="number", c("snowcrab"), c("snowcrab", "groundfish")),
      yr = p$assessment.years,      # time frame for comparison specified above
      settype = 1, # same as geartype in groundfish_survey_db
      polygon_enforce=TRUE,  # make sure mis-classified stations or incorrectly entered positions get filtered out
      strata_toremove = NULL #,  # emphasize that all data enters analysis initially ..
      # ranged_data = c("dyear")  # not used .. just to show how to use range_data
    )


    if ( !exists("variabletomodel", p)) p$variabletomodel = "totno"

    if ( !exists("carstm_modelengine", p)) p$carstm_modelengine = "inla"  # {model engine}.{label to use to store}

    if ( !exists("carstm_model_call", p)) {
      if ( grepl("inla", p$carstm_modelengine) ) {
        p$libs = unique( c( p$libs, project.library ( "INLA" ) ) )
        if ( !exists("carstm_model_label", p))  p$carstm_model_label = "production"
        # 281 configs .. 8hrs

        #number of nodes can be adjusted below to smooth (n=9 instead of n=13)
        p$carstm_model_call = paste(
          'inla( formula =', p$variabletomodel,
          ' ~ 1
            + offset( log(data_offset))
            + f( dyri, model="ar1", hyper=H$ar1 )
            + f( inla.group( t, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2)
            + f( inla.group( z, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2)
            + f( inla.group( substrate.grainsize, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2)
            + f( inla.group( pca1, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2)
            + f( inla.group( pca2, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2)
            + f( auid, model="bym2", graph=slot(sppoly, "nb"), group=year_factor, scale.model=TRUE, constr=TRUE, hyper=H$bym2, control.group=list(model="ar1", hyper=H$ar1_group)),
            family = "poisson",
            data= M,
            control.compute=list(cpo=TRUE, waic=TRUE, dic=TRUE, config=TRUE),
            control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
            control.predictor=list(compute=FALSE, link=1 ),
            #control.fixed = list(prec.intercept = 0.1),
            control.fixed = H$fixed,  # priors for fixed effects, generic is ok
            control.inla = list(cmin = 0, h=1e-4, tolerance=1e-9, strategy="adaptive", optimise.strategy="smart"),
            # control.inla = list( h=1e-6, tolerance=1e-12), # increase in case values are too close to zero
            #control.inla=list( strategy="laplace", cutoff=1e-6, correct=TRUE, correct.verbose=FALSE ),
            # control.inla = list( h=1e-6, tolerance=1e-12), # increase in case values are too close to zero
            # control.inla = list(h=1e-3, tolerance=1e-9, cmin=0), # restart=3), # restart a few times in case posteriors are poorly defined
            # control.mode = list( restart=TRUE, result=RES ), # restart from previous estimates
          verbose=TRUE
          )'
        )
      }
        #    + f(tiyr, model="ar1", hyper=H$ar1 )
        # + f(year,  model="ar1", hyper=H$ar1 )

      if ( grepl("glm", p$carstm_modelengine) ) {
        if ( !exists("carstm_model_label", p))  p$carstm_model_label = "default_glm"
        p$carstm_model_call = paste(
          'glm( formula =',  p$variabletomodel,
          ' ~ 1 + factor(AUID) + t + z + substrate.grainsize +tiyr,
            data= M[ which(M$tag=="observations"), ],
            family=gaussian(link="identity")
          )'
        )
      }

      if ( grepl("gam", p$carstm_modelengine) ) {
        p$libs = unique( c( p$libs, project.library ( "mgcv"  ) ) )

        if ( !exists("carstm_model_label", p))  p$carstm_model_label = "default_gam"
        p$carstm_model_call = paste(
          'gam( formula =',  p$variabletomodel,
          ' ~ 1 + factor(AUID) + s(t) + s(z) + s(substrate.grainsize) + s(yr) + s(dyear),
            data= M[ which(M$tag=="observations"), ],
            family=gaussian(link="identity")
          )'
        )
      }
    }

    if (!exists("nsims", p) ) p$nsims = 10000

    p = carstm_parameters( p=p)  # fill in anything missing and some checks
    p$carstm_inputs_aggregated = FALSE

    return(p)
  }

  # -------------------------


  if ( DS=="carstm_inputs") {

    # prediction surface
    crs_lonlat = st_crs(projection_proj4string("lonlat_wgs84"))
    sppoly = areal_units( p=p )  # will redo if not found
    sppoly = st_transform(sppoly, crs=crs_lonlat )
    areal_units_fn = attributes(sppoly)[["areal_units_fn"]]

    # shared accross various secneario using the same polys
    #.. store at the modeldir level as default
    # outputdir = file.path(p$modeldir, p$carstm_model_label)
    outputdir = file.path(p$modeldir )
    if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )

    fn = file.path( outputdir,
      paste( "snowcrab", "carstm_inputs", areal_units_fn,
        p$variabletomodel, paste0(p$selection$survey$data.source, collapse=""),
        p$inputdata_spatial_discretization_planar_km,
        round(p$inputdata_temporal_discretization_yr, 6),
        "rdata",
        sep="."
      )
    )

    if (!redo)  {
      if (file.exists(fn)) {
        load( fn )
        return( M )
      }
    }
    message( "Generating carstm_inputs ... ")


    # do this immediately to reduce storage for sppoly (before adding other variables)
    M = snowcrab.db( p=p, DS="biological_data" )  # will redo if not found .. not used here but used for data matching/lookup in other aegis projects that use bathymetry
    # some survey timestamps extend into January (e.g., 2020) force them to be part of the correct "survey year", i.e., "yr"
    i = which(lubridate::month(M$timestamp)==1)
    if (length(i) > 0) M$timestamp[i] = M$timestamp[i] - lubridate::duration(month=1)

    M$tiyr = lubridate::decimal_date(M$timestamp)

    # M$totno = M$totno_adjusted / M$cf_set_no   # convert density to counts
    # M$totwgt = M$totwgt_adjusted / M$cf_set_mass # convert density to total wgt

    # M$data_offset = 1 / M$cf_set_no  ## offset only used in poisson model


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

    M$AUID = st_points_in_polygons(
      pts = st_as_sf( M, coords=c("lon","lat"), crs=crs_lonlat ),
      polys = sppoly[, "AUID"],
      varname = "AUID"
    )
    M = M[!is.na(M$AUID),]

    names(M)[which(names(M)=="yr") ] = "year"
    # M = M[ which(M$year %in% p$yrs), ]
    # M$tiyr = lubridate::decimal_date ( M$timestamp )
    # M$dyear = M$tiyr - M$year

    pB = bathymetry_parameters( p=parameters_reset(p), project_class="carstm"  )
    vnmod = pB$variabletomodel
    if (!(exists(vnmod, M ))) M[,vnmod] = NA
  
      iM = which(!is.finite( M[, vnmod] ))
      if (length(iM > 0)) {
        LU = bathymetry_db ( p=bathymetry_parameters( spatial_domain=p$spatial_domain, project_class="core"  ), DS="aggregated_data" )  # raw data
        LU = lonlat2planar(LU, proj.type=p$aegis_proj4string_planar_km)
        LU_map = array_map( "xy->1", LU[,c("plon","plat")], gridparams=p$gridparams )
        M_map = array_map( "xy->1", M[iM,c("plon","plat")], gridparams=p$gridparams )
        M[iM, vnmod] = LU[match( M_map, LU_map ), paste(vnmod, "mean", sep=".") ]

        # if any still missing then use a mean depth by AUID
        iM = NULL
        iM =  which( !is.finite(M[, vnmod]))
        if (length(iM) > 0) {
          LU$AUID = st_points_in_polygons(
            pts = st_as_sf( LU, coords=c("lon","lat"), crs=crs_lonlat ),
            polys = sppoly[, "AUID"],
            varname="AUID"
          )
          oo = tapply( LU[, paste(vnmod, "mean", sep="." )], LU$AUID, FUN=median, na.rm=TRUE )
          jj = match( as.character( M$AUID[iM]), as.character( names(oo )) )
          M[iM, vnmod] = oo[jj ]
        }
      }
      M = M[ is.finite(M[ , vnmod]  ) , ]

      if (p$carstm_inputs_aggregated) {
        if ( exists("spatial_domain", p)) {
          M = M[ geo_subset( spatial_domain=p$spatial_domain, Z=M ) , ] # need to be careful with extrapolation ...  filter depths
        }
      }



      pS = substrate_parameters( p=parameters_reset(p), project_class="carstm"  )
      vnmod = pS$variabletomodel
      if (!(exists(vnmod, M ))) M[,vnmod] = NA
  
      iM = which(!is.finite(M[, vnmod]))
      if (length(iM) > 0 ) {
        LU = substrate_db ( p=substrate_parameters( spatial_domain=p$spatial_domain, project_class="core"  ), DS="aggregated_data" )  # raw data
        LU = LU[ which( LU$lon > p$corners$lon[1] & LU$lon < p$corners$lon[2]  & LU$lat > p$corners$lat[1] & LU$lat < p$corners$lat[2] ), ]
        LU = lonlat2planar(LU, proj.type=p$aegis_proj4string_planar_km)
        # levelplot( eval(paste(vnmod, "mean", sep="."))~plon+plat, data=M, aspect="iso")
        LU_map = array_map( "xy->1", LU[,c("plon","plat")], gridparams=p$gridparams )
        M_map  = array_map( "xy->1", M[iM, c("plon","plat")], gridparams=p$gridparams )
        M[iM, vnmod] = LU[ match( M_map, LU_map ), paste(vnmod, "mean", sep=".") ]
        LU_map = NULL
        M_map = NULL
        iM = NULL

        # if any still missing then use a mean substrate by AUID
        iM =  which( !is.finite(M[, vnmod]))
        if (length(iM) > 0) {
          LU$AUID = st_points_in_polygons(
            pts = st_as_sf( LU, coords=c("lon","lat"), crs=crs_lonlat ),
            polys = sppoly[, "AUID"],
            varname="AUID"
          )
          LU = tapply( LU[, paste(vnmod, "mean", sep="." )], LU$AUID, FUN=median, na.rm=TRUE )
          iML = match( as.character( M$AUID[iM]), as.character( names(LU )) )
          M[iM, vnmod] = LU[iML ]
        }
        iM = NULL
        LU = NULL
        iML = NULL
      }


      pT = temperature_parameters( p=parameters_reset(p), project_class="carstm"  )
      vnmod = pT$variabletomodel
      if (!(exists(vnmod, M ))) M[,vnmod] = NA

      iM = which(!is.finite(M[, vnmod]))
      if (length(iM) > 0 ) {

        tz = "America/Halifax"
        T = data.frame( timestamp = M$timestamp[iM] )
        if (! "POSIXct" %in% class(T$timestamp)  ) T$timestamp = as.POSIXct( T$timestamp, tz=tz, origin=lubridate::origin  )
        T$yr = lubridate::year(T$timestamp)        
        T$dyear = lubridate::decimal_date( T$timestamp ) - T$yr
        LU = temperature_db ( p=temperature_parameters( spatial_domain=p$spatial_domain, project_class="core" ), year.assessment=max(p$yrs), DS="aggregated_data" )  # raw data
        names(LU)[ which(names(LU) =="temperature.mean") ] = vnmod
        LU = LU[ which( LU$lon > p$corners$lon[1] & LU$lon < p$corners$lon[2]  & LU$lat > p$corners$lat[1] & LU$lat < p$corners$lat[2] ), ]
        LU = lonlat2planar(LU, proj.type=p$aegis_proj4string_planar_km)
        LUT_map = array_map( "ts->1", LU[,c("yr", "dyear")], dims=c(p$ny, p$nw), res=c( 1, 1/p$nw ), origin=c( min(p$yrs), 0) )
        LUS_map = array_map( "xy->1", LU[,c("plon","plat")], gridparams=p$gridparams )
        T_map = array_map( "ts->1", T[, c("yr", "dyear")], dims=c(p$ny, p$nw), res=c( 1, 1/p$nw ), origin=c( min(p$yrs), 0) )
        M_map = array_map( "xy->1", M[iM, c("plon","plat")], gridparams=p$gridparams )
        iLM = match( paste(M_map, T_map, sep="_"), paste(LUS_map, LUT_map, sep="_") )
        M[ iM, vnmod ] = LU[ iLM, paste(vnmod, "mean", sep="." ) ]
        LU = LUT_map = LUS_map = T = T_map = M_map = iLM = NULL

        iM =  which( !is.finite(M[, vnmod]))
        if (length(iM) > 0) {
          LU$AUID = st_points_in_polygons(
            pts = st_as_sf( LU, coords=c("lon","lat"), crs=crs_lonlat ),
            polys = sppoly[, "AUID"],
            varname="AUID"
          )
          LU_dyear_discret = discretize_data( LU$dyear, p$discretization$dyear ) 
          M_dyear_discret = discretize_data( M$dyear, p$discretization$dyear )  # LU$dyear is discretized. . match discretization
          LU$uid = paste(LU$AUID, LU$yr, LU_dyear_discret, sep=".")
          M$uid =  paste(M$AUID, M$yr, M_dyear_discret, sep=".")
          LU = tapply( LU[, paste(vnmod, "mean", sep="." )], LU$uid, FUN=median, na.rm=TRUE )
          iLM = match( as.character( M$uid[iM]), as.character( names(LU )) )
          M[iM, vnmod] = LU[iLM ]
        }

      }


    pPC1 = speciescomposition_parameters( p=parameters_reset(p), project_class="carstm", variabletomodel="pca1" )
    vnmod = pPC1$variabletomodel
    if (!(exists(vnmod, M ))) M[,vnmod] = NA

    iM = which(!is.finite(M[, vnmod]))
    if (length(iM) > 0 ) {
      LU = speciescomposition_db( p=pPC1, DS="speciescomposition"  )
      LU = planar2lonlat(LU, proj.type=p$aegis_proj4string_planar_km )
      iM = match( M$id, LU$id)
      M[iM, vnmod] = LU[iM, vnmod]
      iM =  which( !is.finite(M[, vnmod]))
      if (length(iM) > 0) {
        LU$AUID = st_points_in_polygons(
          pts = st_as_sf( LU, coords=c("lon","lat"), crs=crs_lonlat ),
          polys = sppoly[, "AUID"],
          varname="AUID"
        )
        LU = tapply( LU[, vnmod], LU$AUID, FUN=median, na.rm=TRUE )
        iLM = match( as.character( M$AUID[iM]), as.character( names(LU )) )
        if (length(jj) > 0) M[iM, vnmod] = LU[ iLM ]
        iM =  which( !is.finite(M[, vnmod]))
        if (length(iM) > 0) {
          LU = carstm_summary ( p=pPC1 )
          au_map = match( M$AUID[iM], dimnames(LU)$AUID )
          year_map = match( as.character(M$yr[iM]), dimnames(LU)$yr )
          dindex = cbind(au_map, year_map )
          M[iM, vnmod] = LU [dindex]
        }
      }
    }



    pPC2 = speciescomposition_parameters( p=parameters_reset(p), project_class="carstm", variabletomodel="pca2" )
    vnmod = pPC2$variabletomodel
    if (!(exists(vnmod, M ))) M[,vnmod] = NA

    iM = which(!is.finite(M[, vnmod]))
    if (length(iM) > 0 ) {
      LU = speciescomposition_db( p=pPC2, DS="speciescomposition"  )
      LU = planar2lonlat(LU, proj.type=p$aegis_proj4string_planar_km )
      iM = match( M$id, LU$id)
      M[iM, vnmod] = LU[iM, vnmod]
      iM =  which( !is.finite(M[, vnmod]))
      if (length(iM) > 0) {
        LU$AUID = st_points_in_polygons(
          pts = st_as_sf( LU, coords=c("lon","lat"), crs=crs_lonlat ),
          polys = sppoly[, "AUID"],
          varname="AUID"
        )
        LU = tapply( LU[, vnmod], LU$AUID, FUN=median, na.rm=TRUE )
        iLM = match( as.character( M$AUID[iM]), as.character( names(LU )) )
        if (length(jj) > 0) M[iM, vnmod] = LU[ iLM ]
        iM =  which( !is.finite(M[, vnmod]))
        if (length(iM) > 0) {
          LU = carstm_summary ( p=pPC2 )
          au_map = match( M$AUID[iM], dimnames(LU)$AUID )
          year_map = match( as.character(M$yr[iM]), dimnames(LU)$yr )
          dindex = cbind(au_map, year_map )
          M[iM, vnmod] = LU [dindex]
        }
      }
    }


    PI = NULL
    M$plon = NULL
    M$plat = NULL
    M$lon = NULL
    M$lat = NULL

    M = M[ which(is.finite(M[, pB$variabletomodel] )),]
    M = M[ which(is.finite(M[, pS$variabletomodel] )),]
    M = M[ which(is.finite(M[, pT$variabletomodel] )),]
    M = M[ which(!is.na(M$AUID)),]
    M$AUID = as.character( M$AUID )  # match each datum to an area

    M$tag = "observations"



      # end observations
      # ----------

      # predicted locations (APS)


    region.id = slot( slot(sppoly, "nb"), "region.id" )
    APS = st_drop_geometry(sppoly)
  
    APS$AUID = as.character( APS$AUID )
    APS$tag ="predictions"
    APS[,p$variabletomodel] = NA


    BI = carstm_summary ( p=pB )
    jj = match( as.character( APS$AUID), as.character( BI$AUID) )
    APS[, pB$variabletomodel] = BI[[ paste(pB$variabletomodel,"predicted",sep="." ) ]] [jj]
    jj =NULL
    BI = NULL

    SI = carstm_summary ( p=pS )
    jj = match( as.character( APS$AUID), as.character( SI$AUID) )
    APS[, pS$variabletomodel] = SI[[ paste(pS$variabletomodel,"predicted",sep="." )]] [jj]
    jj =NULL
    SI = NULL

    # to this point APS is static, now add time dynamics (teperature)
    # ---------------------

    vn = c( p$variabletomodel, pB$variabletomodel,  pS$variabletomodel, "tag", "AUID" )
    APS = APS[, vn]

    # expand APS to all time slices
    n_aps = nrow(APS)
    APS = cbind( APS[ rep.int(1:n_aps, p$nt), ], rep.int( p$prediction_ts, rep(n_aps, p$nt )) )
    names(APS) = c(vn, "tiyr")
    APS$yr = aegis_floor( APS$tiyr)
    APS$dyear = APS$tiyr - APS$yr


    TI = carstm_summary ( p=pT )
    TI = TI[[ paste(pT$variabletomodel,"predicted",sep="." )]]
    au_map = match( APS$AUID, dimnames(TI)$AUID )
    year_map = match( as.character(APS$yr), dimnames(TI)$yr )
    dyear_breaks = c(p$dyears, p$dyears[length(p$dyears)]+ diff(p$dyears)[1] )
    dyear_map = as.numeric( cut( APS$dyear, breaks=dyear_breaks, include.lowest=TRUE, ordered_result=TRUE, right=FALSE ) )
    dindex = cbind(au_map, year_map, dyear_map )
    APS[, pT$variabletomodel] = TI[ dindex]
    TI = NULL


    PI = carstm_summary ( p=pPC1 )
    PI = PI[[ paste(pPC1$variabletomodel,"predicted",sep="." )]]
    au_map = match( APS$AUID, dimnames(PI)$AUID )
    year_map = match( as.character(APS$yr), dimnames(PI)$yr )
    dindex = cbind(au_map, year_map )
    APS[, pPC1$variabletomodel] = PI [dindex]
    PI = NULL

    PI = carstm_summary ( p=pPC2 )
    PI = PI[[ paste(pPC2$variabletomodel,"predicted",sep="." )]]
    au_map = match( APS$AUID, dimnames(PI)$AUID )
    year_map = match( as.character(APS$yr), dimnames(PI)$yr )
    dindex = cbind(au_map, year_map  )
    APS[, pPC2$variabletomodel] = PI [dindex]
    PI = NULL

    # useful vars to have for analyses outside of carstm_summary
    varstoadd = c( "totwgt", "totno", "sa", "data_offset",  "zn", "qn" )

    for (vn in varstoadd) if (!exists( vn, APS)) APS[,vn] = NA
    APS$data_offset = 1  # force to solve for unit area

    M = rbind( M[, names(APS)], APS )
    APS = NULL


    M$AUID  = as.character(M$AUID)  # revert to factors
    M$auid  = as.numeric( factor(M$AUID) )

    M$zi  = discretize_data( M[, pB$variabletomodel], p$discretization[[pB$variabletomodel]] )
    M$ti  = discretize_data( M[, pT$variabletomodel], p$discretization[[pT$variabletomodel]] )
    M$gsi = discretize_data( M[, pS$variabletomodel], p$discretization[[pS$variabletomodel]] )

    M$pca1i = discretize_data( M[, pPC1$variabletomodel], p$discretization[[pPC1$variabletomodel]] )
    M$pca2i = discretize_data( M[, pPC2$variabletomodel], p$discretization[[pPC2$variabletomodel]] )

    M$tiyr  = aegis_floor( M$tiyr / p$tres )*p$tres    # discretize for inla .. midpoints

    M$yr = aegis_floor( M$tiyr)
    M$year_factor = as.numeric( factor( M$yr, levels=p$yrs))

    M$dyear =  M$tiyr - M$yr   # revert dyear to non-discretized form

    M$dyri = discretize_data( M[, "dyear"], p$discretization[["dyear"]] )

    # M$seasonal = (as.numeric(M$year_factor) - 1) * length(p$dyears)  + as.numeric(M$dyear)

    save( M, file=fn, compress=TRUE )

    if (redo) M=fn  # when redo'ing .. return file name aftert the save

    return( M )
  }



  # -------------------------


  if ( DS=="carstm_inputs_hybrid") {

    # various global data sources

    # prediction surface
    crs_lonlat = st_crs(projection_proj4string("lonlat_wgs84"))
    sppoly = areal_units( p=p )  # will redo if not found
    sppoly = st_transform(sppoly, crs=crs_lonlat )
    areal_units_fn = attributes(sppoly)[["areal_units_fn"]]

    # shared accross various secneario using the same polys
    #.. store at the modeldir level as default
    # outputdir = file.path(p$modeldir, p$carstm_model_label)
    outputdir = file.path(p$modeldir )
    if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )

    fn = file.path( outputdir,
      paste( "snowcrab", "carstm_inputs", areal_units_fn,
        p$variabletomodel, paste0(p$selection$survey$data.source, collapse=""),
        p$inputdata_spatial_discretization_planar_km,
        round(p$inputdata_temporal_discretization_yr, 6),
        "rdata",
        sep="."
      )
    )

    if (!redo)  {
      if (file.exists(fn)) {
        load( fn )
        return( M )
      }
    }
    message( "Generating carstm_inputs ... ")


    # do this immediately to reduce storage for sppoly (before adding other variables)
    M = snowcrab.db( p=p, DS="biological_data" )  # will redo if not found .. not used here but used for data matching/lookup in other aegis projects that use bathymetry

    # some survey timestamps extend into January (e.g., 2020) force them to be part of the correct "survey year", i.e., "yr"
    i = which(lubridate::month(M$timestamp)==1)
    if (length(i) > 0) M$timestamp[i] = M$timestamp[i] - lubridate::duration(month=1)


    M$tiyr=lubridate::decimal_date(M$timestamp)

    # M$totno = M$totno_adjusted / M$cf_set_no   # convert density to counts
    # M$totwgt = M$totwgt_adjusted / M$cf_set_mass # convert density to total wgt

    # M$data_offset = 1 / M$cf_set_no  ## offset only used in poisson model


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

    M$AUID = st_points_in_polygons(
      pts = st_as_sf( M, coords=c("lon","lat"), crs=crs_lonlat ),
      polys = sppoly[, "AUID"],
      varname = "AUID"
    )
    M = M[!is.na(M$AUID),]

    # names(M)[which(names(M)=="yr") ] = "year"
    # M = M[ which(M$year %in% p$yrs), ]
    # M$tiyr = lubridate::decimal_date ( M$timestamp )
    # M$dyear = M$tiyr - M$year



    pB = bathymetry_parameters( p=parameters_reset(p), project_class="carstm"  )
    vnmod = pB$variabletomodel
    if (!(exists(vnmod, M ))) M[,vnmod] = NA
        # already has depth .. but in case some are missing data
 
    iM = which(!is.finite( M[, vnmod] ))
      if (length(iM > 0)) {
        LU = bathymetry_db ( p=bathymetry_parameters( spatial_domain=p$spatial_domain, project_class="core"  ), DS="aggregated_data" )  # raw data
        LU = lonlat2planar(LU, proj.type=p$aegis_proj4string_planar_km)
        LU_map = array_map( "xy->1", LU[,c("plon","plat")], gridparams=p$gridparams )
        M_map = array_map( "xy->1", M[,c("plon","plat")], gridparams=p$gridparams )
        M[iM, vnmod] = LU[match( M_map, LU_map ), paste(vnmod, "mean", sep=".")  ]

        # if any still missing then use a mean depth by AUID
        iM = NULL
        iM =  which( !is.finite(M[, vnmod]))
        if (length(iM) > 0) {
          LU$AUID = st_points_in_polygons(
            pts = st_as_sf( LU, coords=c("lon","lat"), crs=crs_lonlat ),
            polys = sppoly[, "AUID"],
            varname="AUID"
          )
          oo = tapply( LU[, paste(vnmod, "mean", sep="." )], LU$AUID, FUN=median, na.rm=TRUE )
          jj = match( as.character( M$AUID[iM]), as.character( names(oo )) )
          M[iM, vnmod] = oo[jj ]
        }
      }
      M = M[ is.finite(M[ , vnmod]  ) , ]

      if (p$carstm_inputs_aggregated) {
        if ( exists("spatial_domain", p)) {
          M = M[ geo_subset( spatial_domain=p$spatial_domain, Z=M ) , ] # need to be careful with extrapolation ...  filter depths
        }
      }


      # substrate observations
      
      pS = substrate_parameters( p=parameters_reset(p), project_class="carstm"  )
      vnmod = pS$variabletomodel
      if (!(exists(vnmod, M ))) M[,vnmod] = NA
    
      iM = which(!is.finite(M[, vnmod]))
      if (length(iM) > 0 ) {
        LU = substrate_db ( p=substrate_parameters( spatial_domain=p$spatial_domain, project_class="core"  ), DS="aggregated_data" )  # raw data
        LU = LU[ which( LU$lon > p$corners$lon[1] & LU$lon < p$corners$lon[2]  & LU$lat > p$corners$lat[1] & LU$lat < p$corners$lat[2] ), ]
        LU = lonlat2planar(LU, proj.type=p$aegis_proj4string_planar_km)
        # levelplot( eval(paste(vnmod, "mean", sep="."))~plon+plat, data=M, aspect="iso")
        LU_map = array_map( "xy->1", LU[,c("plon","plat")], gridparams=p$gridparams )
        M_map  = array_map( "xy->1", M[iM, c("plon","plat")], gridparams=p$gridparams )
        M[iM, vnmod] = LU[ match( M_map, LU_map ), paste(vnmod, "mean", sep=".") ]
        LU_map = NULL
        M_map = NULL
        iM = NULL

        # if any still missing then use a mean substrate by AUID
        iM =  which( !is.finite(M[, vnmod]))
        if (length(iM) > 0) {
          LU$AUID = st_points_in_polygons(
            pts = st_as_sf( LU, coords=c("lon","lat"), crs=crs_lonlat ),
            polys = sppoly[, "AUID"],
            varname="AUID"
          )
          LU = tapply( LU[, paste(vnmod, "mean", sep="." )], LU$AUID, FUN=median, na.rm=TRUE )
          iML = match( as.character( M$AUID[iM]), as.character( names(LU )) )
          M[iM, vnmod] = LU[iML ]
        }
        iM = NULL
        LU = NULL
        iML = NULL
      }


      # temperature observations
      pT = temperature_parameters( p=parameters_reset(p), project_class="carstm"  )
      if (!(exists(pT$variabletomodel, M ))) M[,pT$variabletomodel] = NA

      iM = which(!is.finite(M[, pT$variabletomodel]))
      if (length(iM) > 0 ) {

        tz = "America/Halifax"
        T = data.frame( timestamp = M$timestamp[iM] )
        if (! "POSIXct" %in% class(T$timestamp)  ) T$timestamp = as.POSIXct( T$timestamp, tz=tz, origin=lubridate::origin  )
        T$yr = lubridate::year(T$timestamp)        
        T$dyear = lubridate::decimal_date( T$timestamp ) - T$yr
        LU = temperature_db ( p=temperature_parameters( spatial_domain=p$spatial_domain, project_class="core" ), year.assessment=max(p$yrs), DS="aggregated_data" )  # raw data
        names(LU)[ which(names(LU) =="temperature.mean") ] = vnmod
        LU = LU[ which( LU$lon > p$corners$lon[1] & LU$lon < p$corners$lon[2]  & LU$lat > p$corners$lat[1] & LU$lat < p$corners$lat[2] ), ]
        LU = lonlat2planar(LU, proj.type=p$aegis_proj4string_planar_km)
        LUT_map = array_map( "ts->1", LU[,c("yr", "dyear")], dims=c(p$ny, p$nw), res=c( 1, 1/p$nw ), origin=c( min(p$yrs), 0) )
        LUS_map = array_map( "xy->1", LU[,c("plon","plat")], gridparams=p$gridparams )
        T_map = array_map( "ts->1", T[, c("yr", "dyear")], dims=c(p$ny, p$nw), res=c( 1, 1/p$nw ), origin=c( min(p$yrs), 0) )
        M_map = array_map( "xy->1", M[iM, c("plon","plat")], gridparams=p$gridparams )
        iLM = match( paste(M_map, T_map, sep="_"), paste(LUS_map, LUT_map, sep="_") )
        M[ iM, vnmod ] = LU[ iLM, paste(vnmod, "mean", sep="." ) ]
        LU = LUT_map = LUS_map = T = T_map = M_map = iLM = NULL

        iM =  which( !is.finite(M[, vnmod]))
        if (length(iM) > 0) {
          LU$AUID = st_points_in_polygons(
            pts = st_as_sf( LU, coords=c("lon","lat"), crs=crs_lonlat ),
            polys = sppoly[, "AUID"],
            varname="AUID"
          )
          LU_dyear_discret = discretize_data( LU$dyear, p$discretization$dyear ) 
          M_dyear_discret = discretize_data( M$dyear, p$discretization$dyear )  # LU$dyear is discretized. . match discretization
          LU$uid = paste(LU$AUID, LU$yr, LU_dyear_discret, sep=".")
          M$uid =  paste(M$AUID, M$yr, M_dyear_discret, sep=".")
          LU = tapply( LU[, paste(vnmod, "mean", sep="." )], LU$uid, FUN=median, na.rm=TRUE )
          iLM = match( as.character( M$uid[iM]), as.character( names(LU )) )
          M[iM, vnmod] = LU[iLM ]
        }


      }




    pPC1 = speciescomposition_parameters( p=parameters_reset(p), project_class="carstm", variabletomodel="pca1" )
    if (!(exists(pPC1$variabletomodel, M ))) M[,pPC1$variabletomodel] = NA
  
    iM = which(!is.finite(M[, pPC1$variabletomodel]))
    if (length(iM) > 0 ) {
      pc1 = speciescomposition_db( p=pPC1, DS="speciescomposition"  )
      ii = match( M$id, pc1$id)
      M[iM, pPC1$variabletomodel] = pc1[ii, pPC1$variabletomodel]
    }
  
    iM =  which( !is.finite(M[, pPC1$variabletomodel]))
    if (length(iM) > 0) {
      pc1 = speciescomposition_db( p=pPC1, DS="speciescomposition"  )
      pc1 = planar2lonlat( pc1, proj.type=p$aegis_proj4string_planar_km )
      pc1$AUID = st_points_in_polygons(
        pts = st_as_sf( pc1, coords=c("lon","lat"), crs=crs_lonlat ),
        polys = sppoly[, "AUID"],
        varname = "AUID"
      )
      oo = tapply( pc1[, pPC1$variabletomodel ], pc1$AUID, FUN=median, na.rm=TRUE )
      jj = match( as.character( M$AUID[iM]), as.character( names(oo )) )
      if (length(jj) > 0) M[iM, pPC1$variabletomodel] = oo[jj ]
    }
    iM =  which( !is.finite(M[, pPC1$variabletomodel]))
    if (length(iM) > 0) {
      PI = carstm_summary ( p=pPC1 )
      au_map = match( M$AUID[iM], dimnames(PI)$AUID )
      year_map = match( as.character(M$yr[iM]), dimnames(PI)$yr )
      dindex = cbind(au_map, year_map )
      M[iM, pPC1$variabletomodel] = PI [dindex]
    }



    pPC2 = speciescomposition_parameters( p=parameters_reset(p), project_class="carstm", variabletomodel="pca2")
    if (!(exists(pPC2$variabletomodel, M ))) M[,pPC2$variabletomodel] = NA

    iM = which(!is.finite(M[, pPC2$variabletomodel]))
    if (length(iM) > 0 ) {
      pc2 = speciescomposition_db( p=pPC2, DS="speciescomposition"  )
      ii = match( M$id, pc2$id)
      M[iM, pPC2$variabletomodel] = pc2[ii, pPC2$variabletomodel]
    }



    iM =  which( !is.finite(M[, pPC2$variabletomodel]))
    if (length(iM) > 0) {
      pc2 = speciescomposition_db( p=pPC2, DS="speciescomposition"  )
      pc2 = planar2lonlat( pc2, proj.type=p$aegis_proj4string_planar_km )
      pc2$AUID = st_points_in_polygons(
        pts = st_as_sf( pc2, coords=c("lon","lat"), crs=crs_lonlat ),
        polys = sppoly[, "AUID"],
        varname = "AUID"
      )

      oo = tapply( pc2[, pPC2$variabletomodel ], pc2$AUID, FUN=median, na.rm=TRUE )
      jj = match( as.character( M$AUID[iM]), as.character( names(oo )) )
      if (length(jj) > 0) M[iM, pPC2$variabletomodel] = oo[jj ]
    }
    iM =  which( !is.finite(M[, pPC2$variabletomodel]))
    if (length(iM) > 0) {
      PI = carstm_summary ( p=pPC2 )
      au_map = match( M$AUID[iM], dimnames(PI)$AUID )
      year_map = match( as.character(M$yr[iM]), dimnames(PI)$yr )
      dindex = cbind(au_map, year_map )
      M[iM, pPC2$variabletomodel] = PI [dindex]
    }


    PI = NULL
    M$plon = NULL
    M$plat = NULL
    M$lon = NULL
    M$lat = NULL

    M = M[ which(is.finite(M[, pB$variabletomodel] )),]
    M = M[ which(is.finite(M[, pS$variabletomodel] )),]
    M = M[ which(is.finite(M[, pT$variabletomodel] )),]
    M = M[ which(!is.na(M$AUID)),]
    M$AUID = as.character( M$AUID )  # match each datum to an area

    M$tag = "observations"


      # end observations
      # ----------

      # predicted locations (APS)




    region.id = slot( slot(sppoly, "nb"), "region.id" )
    APS = st_drop_geometry(sppoly)

    APS$AUID = as.character( APS$AUID )
    APS$tag ="predictions"
    APS[,p$variabletomodel] = NA



    BI = carstm_summary ( p=pB )
    jj = match( as.character( APS$AUID), as.character( BI$AUID) )
    APS[, pB$variabletomodel] = BI[[ paste(pB$variabletomodel,"predicted",sep="." ) ]] [jj]
    jj =NULL
    BI = NULL

    SI = carstm_summary ( p=pS )
    jj = match( as.character( APS$AUID), as.character( SI$AUID) )
    APS[, pS$variabletomodel] = SI[[ paste(pS$variabletomodel,"predicted",sep="." )]] [jj]
    jj =NULL
    SI = NULL

    # to this point APS is static, now add time dynamics (teperature)
    # ---------------------

    vn = c( p$variabletomodel, pB$variabletomodel,  pS$variabletomodel, "tag", "AUID" )
    APS = APS[, vn]

    # expand APS to all time slices
    n_aps = nrow(APS)
    APS = cbind( APS[ rep.int(1:n_aps, p$nt), ], rep.int( p$prediction_ts, rep(n_aps, p$nt )) )
    names(APS) = c(vn, "tiyr")
    APS$yr = aegis_floor( APS$tiyr)
    APS$dyear = APS$tiyr - APS$yr


    TI = carstm_summary ( p=pT )
    TI = TI[[ paste(pT$variabletomodel,"predicted",sep="." )]]
    au_map = match( APS$AUID, dimnames(TI)$AUID )
    year_map = match( as.character(APS$yr), dimnames(TI)$yr )
    dyear_breaks = c(p$dyears, p$dyears[length(p$dyears)]+ diff(p$dyears)[1] )
    dyear_map = as.numeric( cut( APS$dyear, breaks=dyear_breaks, include.lowest=TRUE, ordered_result=TRUE, right=FALSE ) )
    dindex = cbind(au_map, year_map, dyear_map )
    APS[, pT$variabletomodel] = TI[ dindex]
    TI = NULL


    PI = carstm_summary ( p=pPC1 )
    PI = PI[[ paste(pPC1$variabletomodel,"predicted",sep="." )]]
    au_map = match( APS$AUID, dimnames(PI)$AUID )
    year_map = match( as.character(APS$yr), dimnames(PI)$yr )
    dindex = cbind(au_map, year_map )
    APS[, pPC1$variabletomodel] = PI [dindex]
    PI = NULL

    PI = carstm_summary ( p=pPC2 )
    PI = PI[[ paste(pPC2$variabletomodel,"predicted",sep="." )]]
    au_map = match( APS$AUID, dimnames(PI)$AUID )
    year_map = match( as.character(APS$yr), dimnames(PI)$yr )
    dindex = cbind(au_map, year_map  )
    APS[, pPC2$variabletomodel] = PI [dindex]
    PI = NULL

    # useful vars to have for analyses outside of carstm_summary
    varstoadd = c( "totwgt", "totno", "sa", "data_offset",  "zn", "qn" )

    for (vn in varstoadd) if (!exists( vn, APS)) APS[,vn] = NA
    APS$data_offset = 1  # force to solve for unit area

    M = rbind( M[, names(APS)], APS )
    APS = NULL


    M$AUID  = as.character(M$AUID)  # revert to factors
    M$auid = match( M$AUID, region.id )

    M$zi  = discretize_data( M[, pB$variabletomodel], p$discretization[[pB$variabletomodel]] )
    M$ti  = discretize_data( M[, pT$variabletomodel], p$discretization[[pT$variabletomodel]] )
    M$gsi = discretize_data( M[, pS$variabletomodel], p$discretization[[pS$variabletomodel]] )

    M$pca1i = discretize_data( M[, pPC1$variabletomodel], p$discretization[[pPC1$variabletomodel]] )
    M$pca2i = discretize_data( M[, pPC2$variabletomodel], p$discretization[[pPC2$variabletomodel]] )

    M$tiyr  = aegis_floor( M$tiyr / p$tres )*p$tres    # discretize for inla .. midpoints

    M$yr = aegis_floor( M$tiyr)
    M$year_factor = as.numeric( factor( M$yr, levels=p$yrs))

    M$dyear =  M$tiyr - M$yr   # revert dyear to non-discretized form

    M$dyri = discretize_data( M[, "dyear"], p$discretization[["dyear"]] )

    # M$seasonal = (as.numeric(M$year_factor) - 1) * length(p$dyears)  + as.numeric(M$dyear)

    save( M, file=fn, compress=TRUE )

    if (redo) M=fn  # when redo'ing .. return file name aftert the save

    return( M )
  }


  # --------------------------


  if ( any( grepl("carstm_output", DS) ) ) {

    # ie. usually run by "carstm_output_compute" which will bring you here

    sppoly = areal_units( p=p )
    areal_units_fn = attributes(sppoly)[["areal_units_fn"]]

    aufns = carstm_filenames( p=p, projecttype="carstm_outputs", areal_units_fn=areal_units_fn, dropextension=TRUE )

    # same file naming as in carstm ..
    outputdir = file.path(p$modeldir, p$carstm_model_label)

    if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )

    fn     = file.path( outputdir, paste("carstm_modelled_results", aufns, "aggregated_timeseries", "rdata", sep="." )  )
    fn_no  = file.path( outputdir, paste("carstm_modelled_results", aufns, "space_timeseries_number", "rdata", sep="." ) )
    fn_bio = file.path( outputdir, paste("carstm_modelled_results", aufns, "space_timeseries_biomass", "rdata", sep="." ) )
    fn_pa  = file.path( outputdir, paste("carstm_modelled_results", aufns, "space_timeseries_pa", "rdata", sep="." )  )

    if ( DS=="carstm_output_timeseries" ) {
      RES = NA
      if (file.exists(fn)) load( fn)
      return( RES )
    }

    if ( DS=="carstm_output_spacetime_number" ) {
      nums = NA
      if (file.exists(fn_no)) load( fn_no )
      return( nums )
    }

    if ( DS=="carstm_output_spacetime_biomass" ) {
      biom = NA
      if (file.exists(fn_bio)) load( fn_bio )
      return( biom )
    }

    if ( DS=="carstm_output_spacetime_pa")  {
      pa = NA
      if (file.exists(fn_pa)) load( fn_pa )
      return( pa )
    }

    # construct meanweights matrix used to convert number to weight

    fit = carstm_model( p=p, DS="carstm_modelled_fit" ) # to load currently saved res
    res = carstm_summary( p=p  )

    if (p$carstm_modelengine == "inla") {
      ps = inla.posterior.sample(n=p$nsims, fit, selection=list(Predictor=0))  # only predictions
    }

    if (p$carstm_modelengine %in% c( "glm", "gam") ) {
      # sample from marginal distributions as iid assumed
      mu = res[[ paste( p$variabletomodel, "predicted", sep=".")]]
      sigma = res[[ paste( p$variabletomodel, "predicted_se", sep=".")]]
      n = length( c(mu) )
      ncolres = ncol( mu )
      nrowres = nrow( mu )
      ps = tapply( 1:p$nsims, INDEX=1:p$nsims, FUN = function(x) { rnorm( n, mean=c(mu), sd=c(sigma) ) } )
   }


    if (p$selection$type %in% c("presence_absence") ) {
      pa =  res[[ paste( p$variabletomodel, "predicted", sep=".")]]
      pa[!is.finite(pa)] = NA
      pa = inverse.logit(pa)
      pa[!is.finite(pa)] = NA
      # if (is.na(extrapolation_limit)) extrapolation_limit = c(0,1)
      save( pa, file=fn_pa, compress=TRUE )

      sims = sapply( ps,
        function(x) {

          if (p$carstm_modelengine %in% c("glm", "gam") ) {
            pa = matrix(x, nrow=nrowres, ncol=ncolres)
          }

          if (p$carstm_modelengine == "inla") {
            input = x$latent[res$i_preds]
            input[!is.finite(input)] = NA
            input = inverse.logit( input )
            pa = reformat_to_array( input=input , matchfrom=res$matchfrom, matchto=res$matchto )
          }

          pa[!is.finite(pa)] = NA

          o = list()
          o$cfaall    = colSums( pa * sppoly$au_sa_km2/ sum(sppoly$au_sa_km2), na.rm=TRUE )
          o$cfanorth  = colSums( pa * sppoly$cfanorth_surfacearea/ sum(sppoly$cfanorth_surfacearea), na.rm=TRUE )
          o$cfasouth  = colSums( pa * sppoly$cfasouth_surfacearea/ sum(sppoly$cfasouth_surfacearea), na.rm=TRUE )
          o$cfa23     = colSums( pa * sppoly$cfa23_surfacearea/ sum(sppoly$cfa23_surfacearea), na.rm=TRUE )
          o$cfa24     = colSums( pa * sppoly$cfa24_surfacearea/ sum(sppoly$cfa24_surfacearea), na.rm=TRUE )
          o$cfa4x     = colSums( pa * sppoly$cfa4x_surfacearea/ sum(sppoly$cfa4x_surfacearea), na.rm=TRUE )
          return(o)
        }, simplify=TRUE
      )

    }


    if (p$selection$type %in% c("biomass", "number") ) {

      # M = snowcrab_carstm( p=p, DS="carstm_inputs" )

      M = res$M
      M$yr = M$year  # req for meanweights

      wgts = meanweights_by_arealunit(
        set=M[M$tag=="observations",],
        AUID=as.character( sppoly$AUID ),
        yrs=p$yrs,
        fillall=TRUE,
        annual_breakdown=TRUE,
        robustify_quantiles=c(0, 0.99)  # high upper bounds are more dangerous
      )

      if (p$selection$type == "biomass") {
        biom = res[[ paste( p$variabletomodel, "predicted", sep=".")]]
        biom[!is.finite(biom)] = NA

        NA_mask = NULL
        nnn = which( !is.finite(biom ))
        if (length(nnn)>0 ) NA_mask = nnn

        if (is.na(extrapolation_limit)) extrapolation_limit = max(M$totwgt/M$data_offset, na.rm=T) # 28921.8426
        uu = which( biom > extrapolation_limit )
        if (length(uu) > 0 ) {
          if (is.character(extrapolation_replacement)) if (extrapolation_replacement=="extrapolation_limit" ) extrapolation_replacement = extrapolation_limit
          biom[ uu] = extrapolation_replacement
          warning("\n Extreme-valued predictions were found, capping them to max observed rates .. \n you might want to have more informed priors, or otherwise set extrapolation_replacement=NA to replacement value \n")
        }
        biom[biom > extrapolation_limit[2]] = extrapolation_limit[2]
        nums = biom / wgts
        save( biom, file=fn_bio, compress=TRUE )
        save( nums, file=fn_no, compress=TRUE )


        sims = sapply( ps,
          function(x) {
            if (p$carstm_modelengine %in% c("glm", "gam") ) {
              biom = matrix(x, nrow=nrowres, ncol=ncolres)
            }

            if (p$carstm_modelengine == "inla") {
              input = exp( x$latent[res$i_preds])
              biom = reformat_to_array( input=input , matchfrom=res$matchfrom, matchto=res$matchto )
              if (!is.null(NA_mask)) biom[NA_mask] = NA
            }

            biom[!is.finite(biom)] = NA
            o = list()
            o$cfaall    = colSums( biom * sppoly$au_sa_km2/ 10^6, na.rm=TRUE )
            o$cfanorth  = colSums( biom * sppoly$cfanorth_surfacearea/ 10^6, na.rm=TRUE )
            o$cfasouth  = colSums( biom * sppoly$cfasouth_surfacearea/ 10^6, na.rm=TRUE )
            o$cfa23     = colSums( biom * sppoly$cfa23_surfacearea/ 10^6, na.rm=TRUE )
            o$cfa24     = colSums( biom * sppoly$cfa24_surfacearea/ 10^6, na.rm=TRUE )
            o$cfa4x     = colSums( biom * sppoly$cfa4x_surfacearea/ 10^6, na.rm=TRUE )
            return(o)
          }, simplify=TRUE
        )
      }


      if (p$selection$type == "number") {
        nums = res[[ paste( p$variabletomodel, "predicted", sep=".")]]
        nums[!is.finite(nums)] = NA
        NA_mask = NULL
        nnn = which( !is.finite(nums ))
        if (length(nnn)>0 ) NA_mask = nnn

        if (is.na(extrapolation_limit)) extrapolation_limit = max(M$totno/M$data_offset, na.rm=T) # 28921.8426
        uu = which( nums > extrapolation_limit )
        if (length(uu) > 0 ) {
          if (is.character(extrapolation_replacement)) if (extrapolation_replacement=="extrapolation_limit" ) extrapolation_replacement = extrapolation_limit
          nums[ uu] = extrapolation_replacement
          warning("\n Extreme-valued predictions were found, capping them to max observed rates .. \n you might want to have more informed priors, or otherwise set extrapolation=NA to replacement value \n")
        }
        nums[nums > extrapolation_limit[2]] = extrapolation_limit[2]
        biom = nums * wgts
        save( biom, file=fn_bio, compress=TRUE )
        save( nums, file=fn_no, compress=TRUE )

        sims = sapply( ps,
          function(x) {

            if (p$carstm_modelengine %in% c("glm", "gam") ) {
              nums = matrix(x, nrow=nrowres, ncol=ncolres  )
            }

            if (p$carstm_modelengine == "inla") {
              input = exp( x$latent[res$i_preds])
              nums = reformat_to_array( input=input , matchfrom=res$matchfrom, matchto=res$matchto )
              if (!is.null(NA_mask)) nums[NA_mask] = NA
            }

            nums[!is.finite(nums)] = NA
            biom = nums * wgts

            o = list()
            o$cfaall    = colSums( biom * sppoly$au_sa_km2/ 10^6, na.rm=TRUE )
            o$cfanorth  = colSums( biom * sppoly$cfanorth_surfacearea/ 10^6, na.rm=TRUE )
            o$cfasouth  = colSums( biom * sppoly$cfasouth_surfacearea/ 10^6, na.rm=TRUE )
            o$cfa23     = colSums( biom * sppoly$cfa23_surfacearea/ 10^6, na.rm=TRUE )
            o$cfa24     = colSums( biom * sppoly$cfa24_surfacearea/ 10^6, na.rm=TRUE )
            o$cfa4x     = colSums( biom * sppoly$cfa4x_surfacearea/ 10^6, na.rm=TRUE )
            return(o)
          }, simplify=TRUE
        )

      }
    }

    RES = data.frame( yrs = p$yrs )
    RES$cfaall = apply( simplify2array(sims["cfaall",]), 1, mean )
    RES$cfaall_sd = apply( simplify2array(sims["cfaall",]), 1, sd )
    RES$cfaall_median = apply( simplify2array(sims["cfaall",]), 1, median )
    RES$cfaall_lb = apply( simplify2array(sims["cfaall",]), 1, quantile, probs=0.025 )
    RES$cfaall_ub = apply( simplify2array(sims["cfaall",]), 1, quantile, probs=0.975 )

    RES$cfanorth = apply( simplify2array(sims["cfanorth",]), 1, mean )
    RES$cfanorth_sd = apply( simplify2array(sims["cfanorth",]), 1, sd )
    RES$cfanorth_median = apply( simplify2array(sims["cfanorth",]), 1, median )
    RES$cfanorth_lb = apply( simplify2array(sims["cfanorth",]), 1, quantile, probs=0.025 )
    RES$cfanorth_ub = apply( simplify2array(sims["cfanorth",]), 1, quantile, probs=0.975 )

    RES$cfasouth = apply( simplify2array(sims["cfasouth",]), 1, mean )
    RES$cfasouth_sd = apply( simplify2array(sims["cfasouth",]), 1, sd )
    RES$cfasouth_median = apply( simplify2array(sims["cfasouth",]), 1, median )
    RES$cfasouth_lb = apply( simplify2array(sims["cfasouth",]), 1, quantile, probs=0.025 )
    RES$cfasouth_ub = apply( simplify2array(sims["cfasouth",]), 1, quantile, probs=0.975 )

    RES$cfa23 = apply( simplify2array(sims["cfa23",]), 1, mean )
    RES$cfa23_sd = apply( simplify2array(sims["cfa23",]), 1, sd )
    RES$cfa23_median = apply( simplify2array(sims["cfa23",]), 1, median )
    RES$cfa23_lb = apply( simplify2array(sims["cfa23",]), 1, quantile, probs=0.025 )
    RES$cfa23_ub = apply( simplify2array(sims["cfa23",]), 1, quantile, probs=0.975 )

    RES$cfa24 = apply( simplify2array(sims["cfa24",]), 1, mean )
    RES$cfa24_sd = apply( simplify2array(sims["cfa24",]), 1, sd )
    RES$cfa24_median = apply( simplify2array(sims["cfa24",]), 1, median )
    RES$cfa24_lb = apply( simplify2array(sims["cfa24",]), 1, quantile, probs=0.025 )
    RES$cfa24_ub = apply( simplify2array(sims["cfa24",]), 1, quantile, probs=0.975 )

    RES$cfa4x = apply( simplify2array(sims["cfa4x",]), 1, mean )
    RES$cfa4x_sd = apply( simplify2array(sims["cfa4x",]), 1, sd )
    RES$cfa4x_median = apply( simplify2array(sims["cfa4x",]), 1, median )
    RES$cfa4x_lb = apply( simplify2array(sims["cfa4x",]), 1, quantile, probs=0.025 )
    RES$cfa4x_ub = apply( simplify2array(sims["cfa4x",]), 1, quantile, probs=0.975 )

    save( RES, file=fn, compress=TRUE )

    # repeat file extraction here  in case computation was required
    if ( DS=="carstm_output_timeseries" ) {
      RES = NA
      if (file.exists(fn)) load( fn)
      return( RES )
    }

    if ( DS=="carstm_output_spacetime_number" ) {
      nums = NA
      if (file.exists(fn_no)) load( fn_no )
      return( nums )
    }

    if ( DS=="carstm_output_spacetime_biomass" ) {
      biom = NA
      if (file.exists(fn_bio)) load( fn_bio )
      return( biom )
    }

    if ( DS=="carstm_output_spacetime_pa")  {
      pa = NA
      if (file.exists(fn_pa)) load( fn_pa )
      return( pa )
    }

    return(fn)
  }


}  ## end snowcrab.db
