
snowcrab_parameters = function( p=NULL, year.assessment=NULL, project_class="default", ... ) {

  # ---------------------
  # deal with additional passed parameters
  if ( is.null(p) ) p=list()
  p_add = list(...)
  if (length(p_add) > 0 ) p = c(p, p_add)
  i = which(duplicated(names(p), fromLast=TRUE))
  if ( length(i) > 0 ) p = p[-i] # give any passed parameters a higher priority, overwriting pre-existing variable


  # ---------------------
  # create/update library list
  p$libs = unique( c( p$libs, RLibrary ( "colorspace",  "geosphere", "lattice",
    "maps", "mapdata", "maptools", "parallel",  "rgdal", "rgeos",  "sp", "splancs", "rgeos", "bigmemory", "numDeriv", "lubridate", "parallel", "fields", "mgcv" ) ) )
  p$libs = unique( c( p$libs, project.library (
    "aegis", "bio.taxonomy", "stmv", "aegis.polygons", "aegis.coastline", "netmensuration", "bio.snowcrab" ) ) )


  # ---------------------
  if (!exists("project_name", p)) p$project_name = "bio.snowcrab"
  if (!exists("data_root", p)) p$data_root = project.datadirectory( p$project_name )

  p$project.outputdir = project.datadirectory( p$project_name, "output" ) #required for interpolations and mapping
  p$transform_lookup = file.path( p$project.outputdir, "transform.lookup.rdata" ) # local storage of transforms for timeseries plots

  if (!exists( "variabletomodel", p)) p$variabletomodel = "totno"


  # ---------------------
  # define focal year. not required for pure spatial models but ignored by the spatial methods anyways
  if (!is.null(year.assessment) ) p$year.assessment = year.assessment
  if (!exists( "year.assessment", p)) if ( exists("yrs", p)) p$year.assessment = max(p$yrs)
  if (!exists( "year.assessment", p)) p$year.assessment = lubridate::year(lubridate::now())

  # ---------------------
  # define years to map when mapping various variables, defaults to last four years
  if (!exists("mapyears", p)) p$mapyears = (as.numeric(p$year.assessment)-3):p$year.assessment

  if (!exists("yrs", p)) p$yrs = c(1999:p$year.assessment)

  p$seabird.yToload = intersect( p$yrs, 2012:p$year.assessment)
  p$minilog.yToload = intersect( p$yrs, 1999:p$year.assessment)
  p$netmind.yToload = intersect( p$yrs, 1999:p$year.assessment)
  p$esonar.yToload  = intersect( p$yrs, 2014:p$year.assessment)
  p$netmensuration.problems = c()

  p$aegis_dimensionality="space-year"

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

  p = temporal_parameters(p=p, aegis_dimensionality="space-year")

  p$spatial_domain="snowcrab"
  p$spatial_domain_subareas = NULL # add cfa's as subareas .. TODO
  p = spatial_parameters( p=p )  # data are from this domain .. so far

  p$data_sources = c("groundfish", "snowcrab")

  # output location for year-specific results
  p$annual.results = project.datadirectory("bio.snowcrab", "assessments", p$year.assessment )
  p$ofname = file.path(p$annual.results, paste("TSresults", p$year.assessment, "rdata", sep=".") )

  p$fisheries.grid.resolution = 2

  p$regions.to.model = c( "cfanorth", "cfasouth", "cfa4x", "cfaall" )
  p$plottimes=c("annual", "globalaverage")
  p$conversions=c("ps2png")
  p$recode.data = TRUE

  if (!exists("clusters", p)) p$clusters = rep("localhost", detectCores() )
  if (!exists("vars.to.model", p))  p$vars.to.model = bio.snowcrab::snowcrab.variablelist("all.to.model")

  p$habitat.threshold.quantile = 0.05 # quantile at which to consider zero-valued abundance
  p$threshold.distance = 5 # predict no farther than this distance km from survey stations


  if (project_class=="default") {
    return(p)
  }


  if (project_class=="stmv") {

    p$libs = unique( c( p$libs, project.library ( "stmv" ) ) )
    if (!exists("varstomodel", p) ) p$varstomodel = c( "pca1", "pca2", "ca1", "ca2" )
    if (!exists("variables", p)) p$variables = list()
    if (!exists("LOCS", p$variables)) p$variables$LOCS=c("plon", "plat")
    if (!exists("TIME", p$variables)) p$variables$TIME="tiyr"

    p$inputdata_spatial_discretization_planar_km = p$pres  # 1 km .. requires 32 GB RAM and limit of speed -- controls resolution of data prior to modelling to reduce data set and speed up modelling
    p$inputdata_temporal_discretization_yr = 1/12  # ie., monthly .. controls resolution of data prior to modelling to reduce data set and speed up modelling }

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
    if (!exists("Y", p$variables)) {
      if (exists("stmv_local_modelformula", p))  {
        if (!is.null(p$stmv_local_modelformula)) {
          if (p$stmv_local_modelformula != "none") {
            oo = all.vars( p$stmv_local_modelformula[[2]] )
            if (length(oo) > 0) p$variables$Y = oo
          }
        }
      }
      if (exists("stmv_global_modelformula", p))  {
        if (!is.null(p$stmv_global_modelformula)) {
          if (p$stmv_global_modelformula != "none") {
            oo = all.vars( p$stmv_global_modelformula[[2]] )
            if (length(oo) > 0) p$variables$Y = oo
          }
        }
      }
    }
    if (!exists("Y", p$variables)) p$variables$Y = "not_defined" # this can be called to get covars.. do not stop


    # additional variable to extract from aegis_db for inputs
    p$aegis_variables = list()
    # p$aegis_project_datasources = c("speciescomposition", "speciesarea", "sizespectrum", "condition", "metabolism", "biochem")
    if (!exists("aegis_project_datasources", p)) p$aegis_project_datasources = "speciescomposition"
    for (id in p$aegis_project_datasources ) {

      pz = aegis_parameters( p=p, DS=id )
      pz_vars = intersect( pz$varstomodel, p$variables$COV )  # these are aegis vars to model
      if (length(pz_vars) > 0) p$aegis_variables[[id]] = pz_vars
    }

    if (!exists("stmv_local_modelengine", p)) p$stmv_local_modelengine ="gam"
    if (!exists("stmv_global_modelengine", p)) p$stmv_global_modelengine ="gam"
    if (!exists("stmv_global_family", p)) p$stmv_global_family = gaussian(link="log")

    # using covariates as a first pass essentially makes it ~ kriging with external drift .. no time or space here
    if (!exists("stmv_global_modelformula", p)) {
      p$stmv_global_modelformula = formula( paste(
        p$variables$Y,
        ' ~ s( t, k=3, bs="ts") + s( tsd, k=3, bs="ts") + s( tmax, k=3, bs="ts") + s( degreedays, k=3, bs="ts") ',
        ' + s( log(z), k=3, bs="ts") + s( log(dZ), k=3, bs="ts") + s( log(ddZ), k=3, bs="ts") ',
        ' + s( log(substrate.grainsize), k=3, bs="ts") + s(pca1, k=3, bs="ts") + s(pca2, k=3, bs="ts")  '
      ))  # no space
    }

    if (p$stmv_local_modelengine =="twostep") {
      # this is the time component (mostly) .. space enters as a rough constraint
      # if (!exists("stmv_local_modelformula", p))  p$stmv_local_modelformula = formula( paste(
      #   p$variables$Y, '~ s(yr, k=10, bs="ts") + s(cos.w, k=3, bs="ts") + s(sin.w, k=3, bs="ts") ',
      #     ' + s(cos.w, sin.w, yr, bs="ts", k=20) ',
      #     ' + s(plon, k=3, bs="ts") + s(plat, k=3, bs="ts") + s(plon, plat, k=20, bs="ts") ' ) )
      if (!exists("stmv_local_model_distanceweighted", p)) p$stmv_local_model_distanceweighted = TRUE
      if (!exists("stmv_gam_optimizer", p)) p$stmv_gam_optimizer=c("outer", "bfgs")

      if (!exists("stmv_twostep_time", p))  p$stmv_twostep_time = "gam"

      if (p$stmv_twostep_time == "gam") {
        if (!exists("stmv_local_modelformula_time", p)) {
          p$stmv_local_modelformula_time = formula( paste(
            p$variables$Y,
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
        p$variables$Y, '~ s(log(z), k=3, bs="ts") + s(plon, k=3, bs="ts") + s(plat, k=3, bs="ts") + s( log(z), plon, plat, k=27, bs="ts")  ') )
      }
      if (!exists("stmv_fft_filter", p)) p$stmv_fft_filter="matern" #  matern, krige (very slow), lowpass, lowpass_matern

    }  else if (p$stmv_local_modelengine == "gam") {

      if (!exists("stmv_local_modelformula", p))  p$stmv_local_modelformula = formula( paste(
        p$variables$Y, '~ s(yr, bs="ts") + s(cos.w, k=3, bs="ts") + s(sin.w, k=3, bs="ts") ',
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
        p$variables$Y, ' ~ sx(yr, bs="ps") + sx(cos.w, bs="ps") + s(sin.w, bs="ps") +s(z, bs="ps") + sx(plon, bs="ps") + sx(plat,  bs="ps")',
          ' + sx(plon, plat, cos.w, sin.w, yr, bs="te") ' )
          # te is tensor spline
      )
      p$stmv_local_model_bayesxmethod="MCMC"
      p$stmv_local_model_distanceweighted = FALSE
    } else {
      message( "The specified stmv_local_modelengine is not tested/supported ... you are on your own ;) ..." )
    }

    # p = stmv_variablelist(p=p)  # decompose into covariates, etc

    p = aegis_parameters(p=p, DS="stmv" ) # generics:


    return(p)
  }



}
