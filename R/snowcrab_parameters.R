
snowcrab_parameters = function( p=list(), year.assessment=NULL, project_name="bio.snowcrab", project_class="core", ... ) {

  # ---------------------
  # deal with additional passed parameters
  p = parameters_add(p, list(...)) # add passed args to parameter list, priority to args


  # ---------------------
  # create/update library list
  p$libs = unique( c( p$libs, RLibrary ( "colorspace",  "geosphere", "lattice", "GADMTools",
    "maps", "mapdata", "maptools", "parallel",  "rgdal", "rgeos",  "sp", "spdep", "sf" , "term", 
    "rgeos", "bigmemory", "numDeriv", "lubridate", "parallel", "fields", "mgcv" ) ) )
  p$libs = unique( c( p$libs, project.library (
    "aegis", "bio.taxonomy", "stmv",
    "aegis.bathymetry", "aegis.polygons", "aegis.coastline",
    "aegis.survey", "aegis.temperature",
    "aegis.substrate", "aegis.speciescomposition",
    "netmensuration", "bio.snowcrab" ) ) )


  p = parameters_add_without_overwriting( p, project_name = project_name )
  p = parameters_add_without_overwriting( p, data_root = project.datadirectory( p$project_name  ) )
  p = parameters_add_without_overwriting( p, datadir  = p$data_root )  # all unprocessed inputs (and simple manipulations) ..   #  usually the datadir is a subdirectory: "data" of data_root as in snowcrab.db, .. might cause problems
  p = parameters_add_without_overwriting( p, modeldir = file.path( p$data_root, "modelled" ) )  # all model outputs

  if ( !file.exists(p$datadir) ) dir.create( p$datadir, showWarnings=FALSE, recursive=TRUE )
  if ( !file.exists(p$modeldir) ) dir.create( p$modeldir, showWarnings=FALSE, recursive=TRUE )


  p$project.outputdir = project.datadirectory( p$project_name, "output" ) #required for interpolations and mapping
  p$transform_lookup = file.path( p$project.outputdir, "transform.lookup.rdata" ) # local storage of transforms for timeseries plots



  # ---------------------
  # define focal year. not required for pure spatial models but ignored by the spatial methods anyways
  if ( exists( "assessment.years", p)) {
    if ( !exists("year.assessment", p)) p$year.assessment = max( p$assessment.years )
    if ( !exists("yrs", p)) p$yrs = p$assessment.years
  }
  if (!is.null(year.assessment) ) {
    if ( !exists("year.assessment", p)) p$year.assessment = year.assessment  # over-ride
    if ( !exists("yrs", p)) p$yrs = c(1999:p$year.assessment)
  }
  if (!exists( "year.assessment", p)) {
    if ( exists("yrs", p)) p$year.assessment = max(p$yrs)
  }
  if (!exists( "year.assessment", p)) stop("year.assessment not defined" )
  if (!exists("yrs", p)) p$yrs = c(1999:p$year.assessment)


  # ---------------------
  # define years to map when mapping various stmv_variables, defaults to last four years
  if (!exists("mapyears", p)) p$mapyears = (as.numeric(p$year.assessment)-3):p$year.assessment


  p$seabird.yToload = intersect( p$yrs, 2012:p$year.assessment)
  p$minilog.yToload = intersect( p$yrs, 1999:p$year.assessment)
  p$netmind.yToload = intersect( p$yrs, 1999:p$year.assessment)
  p$esonar.yToload  = intersect( p$yrs, 2014:p$year.assessment)
  p$netmensuration.problems = c()

  p$aegis_dimensionality="space-year"  # a single solution each year at a given dyear (vs temperature)

  p$ny = length(p$yrs)
  p$nt = p$ny # must specify, else assumed = 1 (1= no time)  ## nt=ny annual time steps, nt = ny*nw is seassonal
  p$nw = 10 # default value of 10 time steps for all temp and indicators
  p$tres = 1/ p$nw # time resolution .. predictions are made with models that use seasonal components
  p$dyears = (c(1:p$nw)-1)  / p$nw # intervals of decimal years... fractional year breaks
  p$dyear_centre = p$dyears[ round(p$nw/2) ] + p$tres/2
  # used for creating timeslices and predictions  .. needs to match the values in aegis_parameters()
  p$prediction_dyear = lubridate::decimal_date( lubridate::ymd("0000/Sep/01"))
  # output timeslices for predictions in decimla years, yes all of them here
  p$prediction_ts = p$yrs + p$prediction_dyear

  p = temporal_parameters(p=p, aegis_dimensionality="space-year", timezone="America/Halifax")

  p$quantile_bounds =c(0, 0.99) # trim upper bounds

  p$spatial_domain="snowcrab"
  p$spatial_domain_subareas = NULL # add cfa's as subareas .. TODO
  p = spatial_parameters( p=p )  # data are from this domain .. so far

  p$data_sources = c("groundfish", "snowcrab")

  # output location for year-specific results
  p$annual.results = project.datadirectory("bio.snowcrab", "assessments", p$year.assessment )
  p$ofname = file.path(p$annual.results, paste("TSresults", p$year.assessment, "rdata", sep=".") )

  p$fisheries.grid.resolution = 2

  # required for lookups
  p = parameters_add_without_overwriting( p,
      inputdata_spatial_discretization_planar_km = p$pres ,  # 1 km .. some thinning .. requires 32 GB RAM and limit of speed -- controls resolution of data prior to modelling to reduce data set and speed up modelling
      inputdata_temporal_discretization_yr = 1/12
  )

  p$regions.to.model = c( "cfanorth", "cfasouth", "cfa4x", "cfaall" )
  p$plottimes=c("annual", "globalaverage")
  p$conversions=c("ps2png")
  p$recode.data = TRUE

  if (!exists("clusters", p)) p$clusters = rep("localhost", detectCores() )
  if (!exists("vars.to.model", p))  p$vars.to.model = bio.snowcrab::snowcrab.variablelist("all.to.model")

  p$habitat.threshold.quantile = 0.05 # quantile at which to consider zero-valued abundance
  p$threshold.distance = 5 # predict no farther than this distance km from survey stations


  # ---------------------

  if (project_class=="core") {
    p$project_class = "core"
    if (!exists( "variabletomodel", p)) p$variabletomodel = "totno"
    return(p)  # minimal specifications
  }

  # ---------------------

  if (project_class %in% c("stmv") ) {

    p$libs = unique( c( p$libs, project.library ( "stmv" ) ) )
    p$stmv_model_label="default"

    if (!exists( "variabletomodel", p)) p$variabletomodel = "totmass"
  
    if (!exists("stmv_variables", p)) p$stmv_variables = list()
    if (!exists("LOCS", p$stmv_variables)) p$stmv_variables$LOCS=c("plon", "plat")
    if (!exists("TIME", p$stmv_variables)) p$stmv_variables$TIME="tiyr"

    if (!exists("storage_backend", p)) p$storage_backend="bigmemory.ram"

    if (!exists("boundary", p)) p$boundary = FALSE
    if (!exists("stmv_filter_depth_m", p)) p$stmv_filter_depth_m = 0 # depth (m) stats locations with elevation > 0 m as being on land (and so ignore)

    if (!exists("stmv_rsquared_threshold", p)) p$stmv_rsquared_threshold = 0.25 # lower threshold
    if (!exists("stmv_distance_statsgrid", p)) p$stmv_distance_statsgrid = 4 # resolution (km) of data aggregation (i.e. generation of the ** statistics ** )

    if (!exists("stmv_distance_scale", p)) p$stmv_distance_scale = c(25, 35, 45) # km ... approx guess of 95% AC range

    if (!exists("stmv_nmin", p)) p$stmv_nmin = 250 # stmv_nmin/stmv_nmax changes with resolution must be more than the number of knots/edf
    # min number of data points req before attempting to model timeseries in a localized space
    if (!exists("stmv_nmax", p)) p$stmv_nmax = 6000 # actually can have a lot of data from logbooks ... this keeps things reasonable in terms of run-time


    # due to formulae being potentially created on the fly, these are required params

    if (!exists("Y", p$stmv_variables)) {
      if (exists("variabletomodel", p)) p$stmv_variables$Y = p$variabletomodel
    }

    if (!exists("Y", p$stmv_variables)) {
      if (exists("stmv_local_modelformula", p))  {
        if (!is.null(p$stmv_local_modelformula)) {
          if (p$stmv_local_modelformula != "none") {
            oo = all.vars( p$stmv_local_modelformula[[2]] )
            if (length(oo) > 0) p$stmv_variables$Y = oo
          }
        }
      }
    }

    if (!exists("Y", p$stmv_variables)) {
      if (exists("stmv_global_modelformula", p))  {
        if (!is.null(p$stmv_global_modelformula)) {
          if (p$stmv_global_modelformula != "none") {
            oo = all.vars( p$stmv_global_modelformula[[2]] )
            if (length(oo) > 0) p$stmv_variables$Y = oo
          }
        }
      }
    }

    if (!exists("Y", p$stmv_variables)) p$stmv_variables$Y = "not_defined" # this can be called to get covars.. do not stop


    # additional variable to extract from aegis_db for inputs
    p$aegis_variables = list()
    # p$aegis_project_datasources = c("speciescomposition", "speciesarea", "sizespectrum", "condition", "metabolism", "biochem")
    if (!exists("aegis_project_datasources", p)) p$aegis_project_datasources = "speciescomposition"
    for (id in p$aegis_project_datasources ) {

      pz = aegis_parameters( p=p, DS=id )
      pz_vars = intersect( pz$varstomodel, p$stmv_variables$COV )  # these are aegis vars to model
      if (length(pz_vars) > 0) p$aegis_variables[[id]] = pz_vars
    }

    if (!exists("stmv_local_modelengine", p)) p$stmv_local_modelengine ="gam"
    if (!exists("stmv_global_modelengine", p)) p$stmv_global_modelengine ="gam"
    if (!exists("stmv_global_family", p)) p$stmv_global_family = gaussian(link="log")

    # using covariates as a first pass essentially makes it ~ kriging with external drift .. no time or space here
    if (!exists("stmv_global_modelformula", p)) {
      p$stmv_global_modelformula = formula( paste(
        p$stmv_variables$Y,
        ' ~ s( t, k=3, bs="ts") + s( tsd, k=3, bs="ts") + s( tmax, k=3, bs="ts") + s( degreedays, k=3, bs="ts") ',
        ' + s( log(z), k=3, bs="ts") + s( log(dZ), k=3, bs="ts") + s( log(ddZ), k=3, bs="ts") ',
        ' + s( log(substrate.grainsize), k=3, bs="ts") + s(pca1, k=3, bs="ts") + s(pca2, k=3, bs="ts")  '
      ))  # no space
    }

    if (p$stmv_local_modelengine =="twostep") {
      # this is the time component (mostly) .. space enters as a rough constraint
      # if (!exists("stmv_local_modelformula", p))  p$stmv_local_modelformula = formula( paste(
      #   p$stmv_variables$Y, '~ s(yr, k=10, bs="ts") + s(cos.w, k=3, bs="ts") + s(sin.w, k=3, bs="ts") ',
      #     ' + s(cos.w, sin.w, yr, bs="ts", k=20) ',
      #     ' + s(plon, k=3, bs="ts") + s(plat, k=3, bs="ts") + s(plon, plat, k=20, bs="ts") ' ) )
      if (!exists("stmv_local_model_distanceweighted", p)) p$stmv_local_model_distanceweighted = TRUE
      if (!exists("stmv_gam_optimizer", p)) p$stmv_gam_optimizer=c("outer", "bfgs")

      if (!exists("stmv_twostep_time", p))  p$stmv_twostep_time = "gam"

      if (p$stmv_twostep_time == "gam") {
        if (!exists("stmv_local_modelformula_time", p)) {
          p$stmv_local_modelformula_time = formula( paste(
            p$stmv_variables$Y,
            ' ~ s(yr, k=12, bs="ts") + s(cos.w, k=3, bs="ts") + s(sin.w, k=3, bs="ts") ',
            ' + s(cos.w, sin.w, yr, bs="ts", k=30) ',
            ' + s( log(z), k=3, bs="ts" ) + s(plon, k=3, bs="ts") + s(plat, k=3, bs="ts") + s(log(z), plon, plat, k=30, bs="ts") '
          ) )
        }
      }

      # this is the spatial component
      if (!exists("stmv_twostep_space", p))  p$stmv_twostep_space = "fft"
      if (p$stmv_twostep_space == "gam") {
        if (!exists("stmv_local_modelformula_space", p))  p$stmv_local_modelformula_space = formula( paste(
        p$stmv_variables$Y, '~ s(log(z), k=3, bs="ts") + s(plon, k=3, bs="ts") + s(plat, k=3, bs="ts") + s( log(z), plon, plat, k=27, bs="ts")  ') )
      }
      if (!exists("stmv_fft_filter", p)) p$stmv_fft_filter="matern" #  matern, krige (very slow), lowpass, lowpass_matern

    }  else if (p$stmv_local_modelengine == "gam") {

      if (!exists("stmv_local_modelformula", p))  p$stmv_local_modelformula = formula( paste(
        p$stmv_variables$Y, '~ s(yr, bs="ts") + s(cos.w, k=3, bs="ts") + s(sin.w, k=3, bs="ts") ',
          ' + s(cos.w, sin.w, yr, bs="ts", k=25)  ',
          ' + s(plon, k=3, bs="ts") + s(plat, k=3, bs="ts") + s(plon, plat, k=25, bs="ts") ' ) )
      if (!exists("stmv_local_model_distanceweighted", p)) p$stmv_local_model_distanceweighted = TRUE
      # if (!exists("stmv_gam_optimizer", p)) p$stmv_gam_optimizer="perf"
      if (!exists("stmv_gam_optimizer", p)) p$stmv_gam_optimizer=c("outer", "bfgs")

    }  else if (p$stmv_local_modelengine == "bayesx") {
      # bayesx families are specified as characters, this forces it to pass as is and
      # then the next does the transformation internal to the "stmv__bayesx"
      # alternative models .. testing .. problem is that SE of fit is not accessible?
      p$stmv_local_modelformula = formula( paste(
        p$stmv_variables$Y, ' ~ sx(yr, bs="ps") + sx(cos.w, bs="ps") + s(sin.w, bs="ps") +s(z, bs="ps") + sx(plon, bs="ps") + sx(plat,  bs="ps")',
          ' + sx(plon, plat, cos.w, sin.w, yr, bs="te") ' )
          # te is tensor spline
      )
      p$stmv_local_model_bayesxmethod="MCMC"
      p$stmv_local_model_distanceweighted = FALSE
    } else {
      message( "The specified stmv_local_modelengine is not tested/supported ... you are on your own ;) ..." )
    }

    p = aegis_parameters(p=p, DS="stmv" ) # generics:

    return(p)
  }



  if (project_class %in% c("hybrid", "default") ) {

    p$libs = unique( c( p$libs, project.library ( "stmv" ) ) )
    p$stmv_model_label="default"

    if (!exists("varstomodel", p) ) p$varstomodel = c( "pca1", "pca2", "ca1", "ca2" )
    if (!exists("stmv_variables", p)) p$stmv_variables = list()
    if (!exists("LOCS", p$stmv_variables)) p$stmv_variables$LOCS=c("plon", "plat")
    if (!exists("TIME", p$stmv_variables)) p$stmv_variables$TIME="tiyr"


    if (!exists("storage_backend", p)) p$storage_backend="bigmemory.ram"

    if (!exists("boundary", p)) p$boundary = FALSE
    if (!exists("stmv_filter_depth_m", p)) p$stmv_filter_depth_m = 0 # depth (m) stats locations with elevation > 0 m as being on land (and so ignore)

    if (!exists("stmv_rsquared_threshold", p)) p$stmv_rsquared_threshold = 0.25 # lower threshold
    if (!exists("stmv_distance_statsgrid", p)) p$stmv_distance_statsgrid = 4 # resolution (km) of data aggregation (i.e. generation of the ** statistics ** )

    if (!exists("stmv_distance_scale", p)) p$stmv_distance_scale = c(25, 35, 45) # km ... approx guess of 95% AC range

    if (!exists("stmv_nmin", p)) p$stmv_nmin = 250 # stmv_nmin/stmv_nmax changes with resolution must be more than the number of knots/edf
    # min number of data points req before attempting to model timeseries in a localized space
    if (!exists("stmv_nmax", p)) p$stmv_nmax = 6000 # actually can have a lot of data from logbooks ... this keeps things reasonable in terms of run-time


    # due to formulae being potentially created on the fly, these are required params

    if (!exists("Y", p$stmv_variables)) {
      if (exists("variabletomodel", p)) p$stmv_variables$Y = p$variabletomodel
    }

    if (!exists("Y", p$stmv_variables)) {
      if (exists("stmv_local_modelformula", p))  {
        if (!is.null(p$stmv_local_modelformula)) {
          if (p$stmv_local_modelformula != "none") {
            oo = all.vars( p$stmv_local_modelformula[[2]] )
            if (length(oo) > 0) p$stmv_variables$Y = oo
          }
        }
      }
    }

    if (!exists("Y", p$stmv_variables)) {
      if (exists("stmv_global_modelformula", p))  {
        if (!is.null(p$stmv_global_modelformula)) {
          if (p$stmv_global_modelformula != "none") {
            oo = all.vars( p$stmv_global_modelformula[[2]] )
            if (length(oo) > 0) p$stmv_variables$Y = oo
          }
        }
      }
    }

    if (!exists("Y", p$stmv_variables)) p$stmv_variables$Y = "not_defined" # this can be called to get covars.. do not stop


    # additional variable to extract from aegis_db for inputs
    p$aegis_variables = list()
    # p$aegis_project_datasources = c("speciescomposition", "speciesarea", "sizespectrum", "condition", "metabolism", "biochem")
    if (!exists("aegis_project_datasources", p)) p$aegis_project_datasources = "speciescomposition"
    for (id in p$aegis_project_datasources ) {

      pz = aegis_parameters( p=p, DS=id )
      pz_vars = intersect( pz$varstomodel, p$stmv_variables$COV )  # these are aegis vars to model
      if (length(pz_vars) > 0) p$aegis_variables[[id]] = pz_vars
    }

    if (!exists("stmv_local_modelengine", p)) p$stmv_local_modelengine ="carstm"
    if (!exists("stmv_global_modelengine", p)) p$stmv_global_modelengine ="gam"
    if (!exists("stmv_global_family", p)) p$stmv_global_family = gaussian(link="log")

    # using covariates as a first pass essentially makes it ~ kriging with external drift .. no time or space here
    if (!exists("stmv_global_modelformula", p)) {
      p$stmv_global_modelformula = formula( paste(
        p$stmv_variables$Y,
        ' ~ s( t, k=3, bs="ts") + s( tsd, k=3, bs="ts") + s( tmax, k=3, bs="ts") + s( degreedays, k=3, bs="ts") ',
        ' + s( log(z), k=3, bs="ts") + s( log(dZ), k=3, bs="ts") + s( log(ddZ), k=3, bs="ts") ',
        ' + s( log(substrate.grainsize), k=3, bs="ts") + s(pca1, k=3, bs="ts") + s(pca2, k=3, bs="ts")  '
      ))  # no space
    }

    p = aegis_parameters(p=p, DS="stmv" ) # generics:

    return(p)
  }


  # ------------------------------------------------------


  if (project_class %in% c("carstm") ) {

    p$libs = c( p$libs, project.library ( "carstm", "INLA"  ) )

    p$project_class = "carstm"

    p = parameters_add_without_overwriting( p,
        groundfish_species_code=2526,
        speciesname = "Snow crab",
        spatial_domain = "snowcrab",
        yrs = p$yrs,
        tus = "yr",
        carstm_modelengine = "inla",  # {model engine}.{label to use to store}
        carstm_inputs_prefilter = FALSE,
        trawlable_units = "sweptarea"  # <<<<<<<<<<<<<<<<<< also:  "standardtow", "sweptarea" (for groundfish surveys)
    )

    # p$taxa.of.interest = aegis.survey::groundfish_variablelist("catch.summary")
    # p$season = "summer"
    # p$taxa =  "maxresolved"


    p = parameters_add_without_overwriting( p,
      areal_units_xydata = "snowcrab.db(p=p, DS='areal_units_input')",
      areal_units_overlay = "snowcrab_managementareas", # currently: "snowcrab_managementareas",  "groundfish_strata" .. additional polygon layers for subsequent analysis for now ..
      areal_units_constraint = "snowcrab",  # locations of data as constraint .. "snowcrab" loads these automatically, otherwise a xy matrix of positions
      areal_units_proj4string_planar_km = aegis::projection_proj4string("utm20"),  # coord system to use for areal estimation and gridding for carstm
      areal_units_timeperiod = "none",
      nAU_min = 30
    )
    
    if ( !p$areal_units_type %in% c("lattice", "tesselation")) stop("areal_units_type not defined")

    if ( p$areal_units_type == "lattice" ) {
      p = parameters_add_without_overwriting( p,
        areal_units_resolution_km = 25 ,
        areal_units_constraint_ntarget = 10,  
        areal_units_constraint_nmin = 3 
      )
      p = parameters_add_without_overwriting( p, sa_threshold_km2 = 0.1 * p$areal_units_resolution_km^2 )    # drop AU's smaller than 10% of this in km2
    }
 
  
    if ( p$areal_units_type =="tesselation" ) {
      p = parameters_add_without_overwriting( p,
        areal_units_resolution_km = 1, # km  
        areal_units_constraint_ntarget = 15,
        areal_units_constraint_nmin = 3,
        sa_threshold_km2 = 5,
        fraction_cv = 0.9,   # ie. stop if essentially a poisson distribution
        fraction_todrop = 1/7  # control tesselation
      )
    
    }

    if ( !exists("selection", p)) p$selection=list()
    if ( !exists("type", p$selection) ) p$selection$type = "number"
    if ( !exists("biologicals", p$selection) ) p$selection$biologicals = list(
      spec_bio=bio.taxonomy::taxonomy.recode( from="spec", to="parsimonious", tolookup=p$groundfish_species_code ),
      sex=0, # male
      mat=1, # do not use maturity status in groundfish data as it is suspect ..
      len= c( 95, 200 )/10, #  mm -> cm ; aegis_db in cm
      ranged_data="len"
    )

    ### survey data source can be variable depending upon what is being modelled
    if (p$selection$type %in% c("presence_absence") ) {
      p$data.source = c("snowcrab", "groundfish", "logbook" )
    } else if (p$selection$type %in% c("biomass", "number") ) {
      p$data.source =  "snowcrab" 
    }

    if ( !exists("survey", p$selection) ) p$selection$survey = list(
      data.source = "snowcrab",
      yr = p$yrs,      # time frame for comparison specified above
      settype = 1, # same as geartype in groundfish_survey_db
      polygon_enforce=TRUE,  # make sure mis-classified stations or incorrectly entered positions get filtered out
      strata_toremove = NULL #,  # emphasize that all data enters analysis initially ..
      # ranged_data = c("dyear")  # not used .. just to show how to use range_data
    )


    if ( !exists("carstm_inputdata_model_source", p))  p$carstm_inputdata_model_source = list()
    p$carstm_inputdata_model_source = parameters_add_without_overwriting( p$carstm_inputdata_model_source,
      bathymetry = "stmv",  # "stmv", "hybrid", "carstm"
      substrate = "stmv",  # "stmv", "hybrid", "carstm"
      temperature = "carstm",  # "stmv", "hybrid", "carstm"
      speciescomposition = "carstm" # "stmv", "hybrid", "carstm"
    )


    if ( !exists("carstm_modelengine", p)) p$carstm_modelengine = "inla"  # {model engine}.{label to use to store}



    if ( grepl("inla", p$carstm_modelengine) ) {


      if (p$selection$type =="number") {
        if ( !exists("variabletomodel", p)) p$variabletomodel = "totno"
        if ( !exists("carstm_model_label", p)) p$carstm_model_label = paste( p$variabletomodel, p$areal_units_type, p$selection$type, sep="_")

        if ( !exists("carstm_model_formula", p)  ) {
          
          p$carstm_model_formula = as.formula( paste(
          p$variabletomodel, ' ~ 1',
              ' + offset( log(data_offset)) ',
              ' + f( uid, model="iid" ) ',
              ' + f( season, model="rw2", hyper=H$rw2, cyclic=TRUE ) ',
              ' + f( time, model="ar1",  hyper=H$ar1 ) ',
              ' + f( inla.group( t, method="quantile", n=11 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
              ' + f( inla.group( z, method="quantile", n=11 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
              ' + f( inla.group( substrate.grainsize, method="quantile", n=11 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
              ' + f( inla.group( pca1, method="quantile", n=11 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
              ' + f( inla.group( pca2, method="quantile", n=11 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
              ' + f( space, model="bym2", graph=slot(sppoly, "nb"), scale.model=TRUE, constr=TRUE, hyper=H$bym2 ) ',
              ' + f( space_time, model="bym2", graph=slot(sppoly, "nb"), group=time_space, scale.model=TRUE, constr=TRUE, hyper=H$bym2, control.group=list(model="ar1", hyper=H$ar1_group)) '
          ) )
        }
        
        if ( !exists("carstm_model_family", p)  )  {
          p$carstm_model_family =  "poisson"  
          # p$carstm_model_label = "tesselation_overdispersed"   # default is the name of areal_units_type  
          # p$carstm_model_family  = "zeroinflatedpoisson0" #  "binomial",  # "nbinomial", "betabinomial", "zeroinflatedbinomial0" , "zeroinflatednbinomial0"
          # p$carstm_model_inla_control_familiy = NULL
        }

      } 

      if (p$selection$type =="biomass") {
        if ( !exists("variabletomodel", p)) p$variabletomodel = "totmass"
        if ( !exists("carstm_model_label", p)) p$carstm_model_label = paste( p$variabletomodel, p$areal_units_type, p$selection$type, sep="_")

        if ( !exists("carstm_model_formula", p)  ) {
          
          p$carstm_model_formula = as.formula( paste(
          p$variabletomodel, ' ~ 1',
              ' + offset( log(data_offset)) ',
              ' + f( uid, model="iid" ) ',
              ' + f( season, model="rw2", hyper=H$rw2, cyclic =TRUE ) ',
              ' + f( time, model="ar1",  hyper=H$ar1 ) ',
              ' + f( inla.group( t, method="quantile", n=11 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
              ' + f( inla.group( z, method="quantile", n=11 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
              ' + f( inla.group( substrate.grainsize, method="quantile", n=11 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
              ' + f( inla.group( pca1, method="quantile", n=11 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
              ' + f( inla.group( pca2, method="quantile", n=11 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
              ' + f( space, model="bym2", graph=slot(sppoly, "nb"), scale.model=TRUE, constr=TRUE, hyper=H$bym2 ) ',
              ' + f( space_time, model="bym2", graph=slot(sppoly, "nb"), group=time_space, scale.model=TRUE, constr=TRUE, hyper=H$bym2, control.group=list(model="ar1", hyper=H$ar1_group)) '
          ) )
        }
        if ( !exists("carstm_model_family", p)  )  p$carstm_model_family =  "gaussian"  
      } 

      if (p$selection$type =="presence_absence") {
        if ( !exists("variabletomodel", p)) p$variabletomodel = "pa"
        if ( !exists("carstm_model_label", p)) p$carstm_model_label = paste( p$variabletomodel, p$areal_units_type, p$selection$type, sep="_")

        if ( !exists("carstm_model_formula", p)  ) {

          p$carstm_model_formula = as.formula( paste(
            p$variabletomodel, ' ~ 1 ',
              ' + f( season, model="rw2", hyper=H$rw2, cyclic =TRUE ) ',
              ' + f( uid, model="iid" ) ',
              ' + f( time, model="ar1",  hyper=H$ar1 ) ',
              ' + f( inla.group( t, method="quantile", n=11 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
              ' + f( inla.group( z, method="quantile", n=11 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
              ' + f( inla.group( substrate.grainsize, method="quantile", n=11 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
              ' + f( inla.group( pca1, method="quantile", n=11 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
              ' + f( inla.group( pca2, method="quantile", n=11 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
              ' + f( space, model="bym2", graph=slot(sppoly, "nb"), scale.model=TRUE, constr=TRUE, hyper=H$bym2 ) ',
              ' + f( space_time, model="bym2", graph=slot(sppoly, "nb"), group=time_space, scale.model=TRUE, constr=TRUE, hyper=H$bym2, control.group=list(model="ar1", hyper=H$ar1_group)) '
          ) )
        }

        if ( !exists("carstm_model_family", p)  )  p$carstm_model_family = "binomial"  # "nbinomial", "betabinomial", "zeroinflatedbinomial0" , "zeroinflatednbinomial0"
        if ( !exists("carstm_model_inla_control_familiy", p)  )  p$carstm_model_inla_control_familiy = list(control.link=list(model='logit'))

      #  p$carstm_model_family  = "zeroinflatedbinomial1", #  "binomial",  # "nbinomial", "betabinomial", "zeroinflatedbinomial0" , "zeroinflatednbinomial0"
      #  p$carstm_model_inla_control_familiy = NULL


      } 
    
    }  # end carstm-based methods
 
    if ( grepl("glm", p$carstm_modelengine) ) {
      if ( !exists("carstm_model_label", p))  p$carstm_model_label = "default_glm"
      p$carstm_model_formula = as.formula( paste(
        p$variabletomodel, ' ~ 1 + factor(space) + t + z + substrate.grainsize + tiyr' ))
      p$carstm_model_formula = gaussian(link="identity")
      # data = M[ which(M$tag=="observations"), ]
    }

    if ( grepl("gam", p$carstm_modelengine) ) {
      if ( !exists("carstm_model_label", p))  p$carstm_model_label = "default_gam"
      p$carstm_model_formula = as.formula( paste(
        p$variabletomodel, ' ~ 1 + factor(space) + s(t) + s(z) + s(substrate.grainsize) + s(yr) + s(dyear) '))
      # data= M[ which(M$tag=="observations"), ],
      p$carstm_model_formula = gaussian(link="identity")
    }

    if (!exists("nsims", p) ) p$nsims = 10000

    p = carstm_parameters( p=p)  # fill in anything missing and some checks

    return(p)
  }

}
