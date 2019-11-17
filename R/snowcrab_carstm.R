
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
      modeldir = p$modeldir,  # outputs all go the the main project's model output directory
      yrs = p$yrs,
      variabletomodel = "totno",
      spatial_domain = p$spatial_domain,  # defines spatial area, currenty: "snowcrab" or "SSE"
      areal_units_overlay = p$areal_units_overlay, # currently: "snowcrab_managementareas",  "groundfish_strata" .. additional polygon layers for subsequent analysis for now ..
      areal_units_resolution_km = p$areal_units_resolution_km, # km dim of lattice ~ 1 hr
      areal_units_proj4string_planar_km = p$areal_units_proj4string_planar_km,  # coord system to use for areal estimation and gridding for carstm
      inputdata_spatial_discretization_planar_km = p$inputdata_spatial_discretization_planar_km,  # 1 km .. some thinning .. requires 32 GB RAM and limit of speed -- controls resolution of data prior to modelling to reduce data set and speed up modelling
      inputdata_temporal_discretization_yr = p$inputdata_temporal_discretization_yr,  # ie., weekly .. controls resolution of data prior to modelling to reduce data set and speed up modelling
      areal_units_fn = p$areal_units_fn,
      inla_num.threads= p$inla_num.threads,
      inla_blas.num.threads= p$inla_blas.num.threads
    )
    return(pc)
  }


  # ---------------------



  if (DS=="parameters") {


    if ( !exists("project_name", p)) p$project_name = "snowcrab"

    if ( !exists("groundfish_species_code", p)) p$groundfish_species_code = 2526
    if ( !exists("speciesname", p)) p$p$speciesname = "Snow crab"
    if ( !exists("runtype", p)) p$runtype = "number"  # "biomass", "presence_absence", "number"
    if ( !exists("spatial_domain", p)) p$spatial_domain = "snowcrab"  # defines spatial area, currenty: "snowcrab" or "SSE"

    if ( !exists("assessment.years", p)) stop( "must ddefine assessment.years")

    if ( !exists("yrs", p)) p$yrs = p$assessment.years
    if ( !exists("inputdata_spatial_discretization_planar_km", p)) p$inputdata_spatial_discretization_planar_km = 1  # 1 km .. some thinning .. requires 32 GB RAM and limit of speed -- controls resolution of data prior to modelling to reduce data set and speed up modelling
    if ( !exists("inputdata_temporal_discretization_yr", p)) p$inputdata_temporal_discretization_yr = 1/12
    if ( !exists("trawlable_units", p)) p$trawlable_units = "sweptarea"  # <<<<<<<<<<<<<<<<<< also:  "standardtow", "sweptarea" (for groundfish surveys)
    if ( !exists("areal_units_source", p)) p$areal_units_source = "lattice" # "stmv_lattice" to use ageis fields instead of carstm fields ... note variables are not the same
    if ( !exists("areal_units_overlay", p)) p$areal_units_overlay = "snowcrab_managementareas" # currently: "snowcrab_managementareas",  "groundfish_strata" .. additional polygon layers for subsequent analysis for now ..
    if ( !exists("areal_units_resolution_km", p)) p$areal_units_resolution_km = 25 # km dim of lattice ~ 1 hr
    if ( !exists("areal_units_constrain_to_data", p)) p$areal_units_constrain_to_data = TRUE
    if ( !exists("areal_units_fn", p)) p$areal_units_fn = "snowcrab_assessment_25"  # identifyer for areal units polygon filename
    if ( !exists("areal_units_proj4string_planar_km", p)) p$areal_units_proj4string_planar_km = aegis::projection_proj4string("utm20")  # coord system to use for areal estimation and gridding for carstm
    if ( !exists("quantile_bounds", p)) p$quantile_bounds =c(0, 0.99) # trim upper bounds
    if ( !exists("selection", p)) p$selection=list(
      type = p$runtype,
      biologicals=list(
        spec_bio=bio.taxonomy::taxonomy.recode( from="spec", to="parsimonious", tolookup=p$groundfish_species_code ),
        sex=0, # male
        mat=1, # do not use maturity status in groundfish data as it is suspect ..
        len= c( 95, 200 )/10, #  mm -> cm ; aegis_db in cm
        ranged_data="len"
      ),
      survey=list(
        data.source = ifelse (p$runtype=="number", c("snowcrab"), c("snowcrab", "groundfish")),
        yr = p$assessment.years,      # time frame for comparison specified above
        settype = 1, # same as geartype in groundfish db
        polygon_enforce=TRUE,  # make sure mis-classified stations or incorrectly entered positions get filtered out
        strata_toremove = NULL,  # emphasize that all data enters analysis initially ..
        ranged_data = c("dyear")  # not used .. just to show how to use range_data
      )
    )
    if ( !exists("variables", p)) p$variables = list(Y="totno")  # name to give (using stmv access methods)  .. redundant .. to remove (needed for now)
    if ( !exists("variabletomodel", p)) p$variabletomodel = "totno"

    if ( !exists("carstm_modelengine", p)) p$carstm_modelengine = "inla.default"  # {model engine}.{label to use to store}

    if ( !exists("carstm_modelcall", p)) {
      if ( grepl("inla", p$carstm_modelengine) ) {
        p$libs = unique( c( p$libs, project.library ( "INLA" ) ) )
        p$carstm_model_label = "production"
        p$carstm_modelcall = paste(
          'inla( formula =', p$variabletomodel,
          ' ~ 1
            + offset( log(data_offset))
            + f(year_factor, model="ar1", hyper=H$ar1 )
            + f(dyri, model="rw2", scale.model=TRUE, diagonal=1e-4, hyper=H$rw2 )
            + f(ti, model="rw2", scale.model=TRUE, diagonal=1e-4, hyper=H$rw2)
            + f(zi, model="rw2", scale.model=TRUE, diagonal=1e-4, hyper=H$rw2)
            + f(gsi, model="rw2", scale.model=TRUE, diagonal=1e-4, hyper=H$rw2)
            + f(pca1i, model="rw2", scale.model=TRUE, diagonal=1e-4, hyper=H$rw2)
            + f(pca2i, model="rw2", scale.model=TRUE, diagonal=1e-4, hyper=H$rw2)
            + f(auid, model="bym2", graph=sppoly@nb, group=year_factor, scale.model=TRUE, constr=TRUE, hyper=H$bym2),
            family = "poisson",
            data= M,
            control.compute=list(dic=TRUE, config=TRUE),
            control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
            control.predictor=list(compute=FALSE, link=1 ),
            control.fixed=H$fixed,  # priors for fixed effects, generic is ok
            # control.inla = list( h=1e-6, tolerance=1e-12), # increase in case values are too close to zero
            # control.mode = list( restart=TRUE, result=RES ), # restart from previous estimates
            # control.inla = list(cmin = 0 ),
            # control.inla=list( strategy="laplace", cutoff=1e-6, correct=TRUE, correct.verbose=FALSE ),
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
          ' ~ 1 + factor(AUID) + t + z + substrate.grainsize +tiyr,
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
          ' ~ 1 + factor(AUID) + s(t) + s(z) + s(substrate.grainsize) + s(year) + s(dyear),
            data= M[ which(M$tag=="observations"), ],
            family=gaussian(link="identity")
          )'
        )
      }
    }



if (0) {
  # formulafor binomials

  pa = presence.absence( X={set$totno / set$data_offset}, px=0.05 )  # determine presence absence and weighting
  set[, "pa"] = pa$pa
  set[, "wt"] = pa$probs
  pa = NULL

  # generic PC priors
  H = carstm_hyperparameters( sd(set$pa) )


  fit = inla(
    formula = pa ~ 1
    + f(ti, model="rw2", scale.model=TRUE, hyper=H$rw2)
    + f(zi, model="rw2", scale.model=TRUE, hyper=H$rw2)
    + f(di, model="rw2", scale.model=TRUE, hyper=H$rw2)
    + f(year, model="iid", hyper=H$iid)
    + f(auid, model="bym2", graph=sppoly@nb, scale.model=TRUE, constr=TRUE, hyper=H$bym2),
    family="binomial",  # alternates family="zeroinflatedbinomial0", family="zeroinflatedbinomial1",
    data=M,
    control.family=list(control.link=list(model="logit")),
    control.compute=list(dic=TRUE, config=TRUE),
    control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
    control.predictor=list(compute=TRUE, link=1 ), # compute=TRUE on each data location
    control.fixed=H$fixed,  # priors for fixed effects
    control.inla=list(  correct=TRUE, correct.verbose=FALSE ), # strategy="laplace", cutoff=1e-6,
    verbose=TRUE
  )


  APS = cbind( APS, fit$summary.fitted.values[ which(M$tag=="predictions"), ] )

  APS$iyr = match(APS$yr_factor, p$yrs)
  APS$iauid = match( APS$AUID, sppoly$AUID )

  # reformat predictions into matrix form
  out = matrix(NA, nrow=length(sppoly$AUID), ncol=length(p$yrs), dimnames=list( sppoly$AUID, p$yrs) )
  out[ cbind(APS$iauid, APS$iyr) ] = APS$mean
  RES$habitat_strata_CAR.yr_iid = colSums( {out * sppoly$au_sa_km2 }[sppoly$strata_to_keep,], na.rm=TRUE ) /sum(sppoly$au_sa_km2[sppoly$strata_to_keep]) # sa weighted average prob habitat


}


    p = carstm_parameters( p=p, DS="basic" )  # fill in anything missing and some checks

  #  boundingbox = list( xlim = c(-70.5, -56.5), ylim=c(39.5, 47.5)), # bounding box for plots using spplot

  MS = snowcrab.db( p=p, DS="biological_data"  )  # domain is  sse     # the underlying observations/data
  p$range_limit = list(
    z = range( MS$z, na.rm=TRUE ),
    t = range( MS$t, na.rm=TRUE)
  )


    return(p)
  }

  # -------------------------

  if ( DS=="carstm_inputs") {

    fn = file.path( p$modeldir, paste( "snowcrab", "carstm_inputs", p$areal_units_fn,
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

    M$AUID = over( SpatialPoints( M[, c("lon", "lat")], crs_lonlat ), spTransform(sppoly, crs_lonlat ) )$AUID # match each datum to an area

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


    # if any still missing then use a mean depth by AUID
    kk =  which( !is.finite(M[, pB$variabletomodel]))
    if (length(kk) > 0) {
      AD = bathymetry.db ( p=pB, DS="aggregated_data"  )  # 16 GB in RAM just to store!
      AD = AD[ which( AD$lon > p$corners$lon[1] & AD$lon < p$corners$lon[2]  & AD$lat > p$corners$lat[1] & AD$lat < p$corners$lat[2] ), ]
      # levelplot( eval(paste(p$variabletomodel, "mean", sep="."))~plon+plat, data=M, aspect="iso")
      AD$AUID = over( SpatialPoints( AD[, c("lon", "lat")], crs_lonlat ), spTransform(sppoly, crs_lonlat ) )$AUID # match each datum to an area
      oo = tapply( AD[, paste(pB$variabletomodel, "mean", sep="." )], AD$AUID, FUN=median, na.rm=TRUE )
      jj = match( as.character( M$AUID[kk]), as.character( names(oo )) )
      M[kk, pB$variabletomodel] = oo[jj ]
    }

    # if any still missing then use a mean substrate by AUID
    kk =  which( !is.finite(M[, pS$variabletomodel]))
    if (length(kk) > 0) {
      AD = substrate.db ( p=pS, DS="aggregated_data"  )  # 16 GB in RAM just to store!
      AD = AD[ which( AD$lon > p$corners$lon[1] & AD$lon < p$corners$lon[2]  & AD$lat > p$corners$lat[1] & AD$lat < p$corners$lat[2] ), ]
      # levelplot( eval(paste(p$variabletomodel, "mean", sep="."))~plon+plat, data=M, aspect="iso")
      AD$AUID = over( SpatialPoints( AD[, c("lon", "lat")], crs_lonlat ), spTransform(sppoly, crs_lonlat ) )$AUID # match each datum to an area
      oo = tapply( AD[, paste(pS$variabletomodel, "mean", sep="." )], AD$AUID, FUN=median, na.rm=TRUE )
      jj = match( as.character( M$AUID[kk]), as.character( names(oo )) )
      M[kk, pS$variabletomodel] = oo[jj ]
    }

    # substrate coverage poor .. add from modelled results
    kk =  which( !is.finite(M[, pS$variabletomodel]))
    if (length(kk) > 0) {
      SI = carstm_model ( p=pS, DS="carstm_modelled" )
      jj = match( as.character( M$AUID[kk]), as.character( SI$AUID) )
      M[kk, pS$variabletomodel] = SI[[ paste(pS$variabletomodel,"predicted",sep="." )]] [jj]
    }


    # if any still missing then use a mean temp by AUID
    kk =  which( !is.finite(M[, pT$variabletomodel]))
    if (length(kk) > 0) {
      AD = temperature.db ( p=pT, DS="aggregated_data"  )  # 16 GB in RAM just to store!
      AD = AD[ which( AD$lon > p$corners$lon[1] & AD$lon < p$corners$lon[2]  & AD$lat > p$corners$lat[1] & AD$lat < p$corners$lat[2] ), ]
      # levelplot( eval(paste(p$variabletomodel, "mean", sep="."))~plon+plat, data=M, aspect="iso")

      AD$AUID = over( SpatialPoints( AD[, c("lon", "lat")], crs_lonlat ), spTransform(sppoly, crs_lonlat ) )$AUID # match each datum to an area
      AD$uid = paste(AD$AUID, AD$year, AD$dyear, sep=".")

      M_dyear_discret = discretize_data( M$dyear, p$discretization$dyear )  # AD$dyear is discretized. . match discretization
      M$uid =  paste(M$AUID, M$year, M_dyear_discret, sep=".")

      oo = tapply( AD[, paste(pT$variabletomodel, "mean", sep="." )], AD$uid, FUN=median, na.rm=TRUE )

      jj = match( as.character( M$uid[kk]), as.character( names(oo )) )
      M[kk, pT$variabletomodel] = oo[jj ]
    }




    kk =  which( !is.finite(M[, pPC1$variabletomodel]))
    if (length(kk) > 0) {
      pc1 = speciescomposition.db( p=pPC1, DS="speciescomposition"  )
      pc1 = planar2lonlat( pc1, proj.type=p$aegis_proj4string_planar_km )
      pc1$AUID = over( SpatialPoints( pc1[, c("lon", "lat")], crs_lonlat ), spTransform(sppoly, crs_lonlat ) )$AUID # match each datum to an area
      oo = tapply( pc1[, pPC1$variabletomodel ], pc1$AUID, FUN=median, na.rm=TRUE )
      jj = match( as.character( M$AUID[kk]), as.character( names(oo )) )
      if (length(jj) > 0) M[kk, pPC1$variabletomodel] = oo[jj ]
    }
    kk =  which( !is.finite(M[, pPC1$variabletomodel]))
    if (length(kk) > 0) {
      PI = carstm_model ( p=pPC1, DS="carstm_modelled" )
      au_map = match( as.numeric(M$AUID[kk]), levels(sppoly$AUID[as.numeric(dimnames(PI)$auid)]  ) )
      year_map = match( as.character(M$year[kk]), dimnames(PI)$year )
      dindex = cbind(au_map, year_map )
      M[kk, pPC1$variabletomodel] = PI [dindex]
    }



    kk =  which( !is.finite(M[, pPC2$variabletomodel]))
    if (length(kk) > 0) {
      pc2 = speciescomposition.db( p=pPC2, DS="speciescomposition"  )
      pc2 = planar2lonlat( pc2, proj.type=p$aegis_proj4string_planar_km )
      pc2$AUID = over( SpatialPoints( pc2[, c("lon", "lat")], crs_lonlat ), spTransform(sppoly, crs_lonlat ) )$AUID # match each datum to an area
      oo = tapply( pc2[, pPC2$variabletomodel ], pc2$AUID, FUN=median, na.rm=TRUE )
      jj = match( as.character( M$AUID[kk]), as.character( names(oo )) )
      if (length(jj) > 0) M[kk, pPC2$variabletomodel] = oo[jj ]
    }
    kk =  which( !is.finite(M[, pPC2$variabletomodel]))
    if (length(kk) > 0) {
      PI = carstm_model ( p=pPC2, DS="carstm_modelled" )
      au_map = match( as.numeric(M$AUID[kk]), levels(sppoly$AUID[as.numeric(dimnames(PI)$auid)]  ) )
      year_map = match( as.character(M$year[kk]), dimnames(PI)$year )
      dindex = cbind(au_map, year_map )
      M[kk, pPC2$variabletomodel] = PI [dindex]
    }


    if( exists("spatial_domain", p)) M = geo_subset( spatial_domain=p$spatial_domain, Z=M ) # need to be careful with extrapolation ...  filter depths


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


    APS = as.data.frame(sppoly)
    APS$AUID = as.character( APS$AUID )
    APS$tag ="predictions"
    APS[,p$variabletomodel] = NA


    BI = carstm_model ( p=pB, DS="carstm_modelled" )
    jj = match( as.character( APS$AUID), as.character( BI$AUID) )
    APS[, pB$variabletomodel] = BI[[ paste(pB$variabletomodel,"predicted",sep="." ) ]] [jj]
    jj =NULL
    BI = NULL

    SI = carstm_model ( p=pS, DS="carstm_modelled" )
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
    APS$year = floor( APS$tiyr)
    APS$dyear = APS$tiyr - APS$year


    TI = carstm_model ( p=pT, DS="carstm_modelled" )
    TI = TI[[ paste(pT$variabletomodel,"predicted",sep="." )]]
    au_map = match( as.numeric(APS$AUID),levels(sppoly$AUID[as.numeric(dimnames(TI)$auid)]  ) )
    year_map = match( as.character(APS$year), dimnames(TI)$year )
    dyear_breaks = c(p$dyears, p$dyears[length(p$dyears)]+ diff(p$dyears)[1] )
    dyear_map = as.numeric( cut( APS$dyear, breaks=dyear_breaks, include.lowest=TRUE, ordered_result=TRUE, right=FALSE ) )
    dindex = cbind(au_map, year_map, dyear_map )
    APS[, pT$variabletomodel] = TI[ dindex]
    TI = NULL


    PI = carstm_model ( p=pPC1, DS="carstm_modelled" )
    PI = PI[[ paste(pPC1$variabletomodel,"predicted",sep="." )]]
    au_map = match( as.numeric(APS$AUID), levels(sppoly$AUID[as.numeric(dimnames(PI)$auid)]  ) )
    year_map = match( as.character(APS$year), dimnames(PI)$year )
    dindex = cbind(au_map, year_map )
    APS[, pPC1$variabletomodel] = PI [dindex]
    PI = NULL

    PI = carstm_model ( p=pPC2, DS="carstm_modelled" )
    PI = PI[[ paste(pPC2$variabletomodel,"predicted",sep="." )]]
    au_map = match( as.numeric(APS$AUID), levels(sppoly$AUID[as.numeric(dimnames(PI)$auid)]  ) )
    year_map = match( as.character(APS$year), dimnames(PI)$year )
    dindex = cbind(au_map, year_map  )
    APS[, pPC2$variabletomodel] = PI [dindex]
    PI = NULL

    # useful vars to have for analyses outside of carstm_model
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
