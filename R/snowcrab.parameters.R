
snowcrab.parameters = function( p=NULL, DS="default", year.assessment=NULL, varname=NULL ) {

  if ( is.null(p) ) p=list()

  if ( DS=="default") {

    p$project.name = "bio.snowcrab"
    p$data_root = file.path( project.datadirectory( p$project.name ) )
    p$project.outputdir = project.datadirectory( p$project.name, "output" ) #required for interpolations and mapping
    
    rlibs = c( "lubridate", "rgdal", "parallel", "sp", "lattice", "fields", "mgcv" , 
               "geosphere",  "DBI", "Cairo", "Hmisc", "vegan", "akima",   "latticeExtra",  "maptools",  
               "boot", "grid", "RColorBrewer",  "spatstat", "rgeos", "bigmemory" ,"numDeriv")

    p$libs = c( p$libs, suppressMessages( RLibrary( rlibs ) ) )
    p$libs = c( p$libs, suppressMessages( project.library (
      "aegis.env", "bio.taxonomy", "stm", "aegis",  "netmensuration", 
      "bio.groundfish", 
      "bio.snowcrab" ) ) ) 
    p$libs = unique( p$libs )

    if (is.null(year.assessment) ) {
      if ( exists("year.assessment", p) ) {
        year.assessment=p$year.assessment
      } else {
        stop( "year.assessment was not defined" )
      }
    }

    p$year.assessment = year.assessment  

    p$seabird.yToload = 2012:p$year.assessment
    p$minilog.yToload = 1999:p$year.assessment
    p$netmind.yToload = 1999:p$year.assessment
    p$esonar.yToload  = 2014:p$year.assessment
    p$netmensuration.problems = c()

    p$yrs = c(1999:p$year.assessment)  
    p$ny = length(p$yrs)
    p$nt = p$ny # must specify, else assumed = 1 (1= no time)  ## nt=ny annual time steps, nt = ny*nw is seassonal
    p$nw = 10 # default value of 10 time steps for all temp and indicators
    p$tres = 1/ p$nw # time resolution .. predictions are made with models that use seasonal components
    p$dyears = (c(1:p$nw)-1)  / p$nw # intervals of decimal years... fractional year breaks
    p$dyear_centre = p$dyears[ round(p$nw/2) ] + p$tres/2
    # used for creating timeslices and predictions  .. needs to match the values in aegis_parameters()
    p$prediction.dyear = lubridate::decimal_date( lubridate::ymd("0000/Sep/01")) 
    # output timeslices for predictions in decimla years, yes all of them here
    p$prediction.ts = p$yrs + p$prediction.dyear 


    p$data.sources = c("groundfish", "snowcrab")
    p$spatial.domain = "snowcrab"
    p$spatial.domain.subareas = NULL # add cfa's as subareas .. TODO
    p = spatial_parameters( p=p )  # data are from this domain .. so far

    # output location for year-specific results
    p$annual.results = file.path( project.datadirectory("bio.snowcrab"), "assessments", p$year.assessment ) 
    p$ofname = file.path(p$annual.results, paste("TSresults", p$year.assessment, "rdata", sep=".") )
    

    p$fisheries.grid.resolution = 2
   
    p$regions.to.model = c( "cfanorth", "cfasouth", "cfa4x", "cfaall" )
    p$plottimes=c("annual", "globalaverage")
    p$conversions=c("ps2png")
    p$recode.data = TRUE

    if (!exists("clusters", p)) p$clusters = rep("localhost", detectCores() )
    if (!exists("vars.to.model", p))  p$vars.to.model = bio.snowcrab::variable.list.expand("all.to.model") # not sure why we have vars.to.model and vartomodel ... clean this up :: TODO

    p$habitat.threshold.quantile = 0.05 # quantile at which to consider zero-valued abundance
    p$threshold.distance = 5 # predict no farther than this distance km from survey stations
   
  }


  if (DS=="stm") {

    if (!("stm" %in% p$libs)) {
      RLibrary( "stm" )
      p$libs = c( p$libs, "stm" )  # required for parallel processing
    }
    
    if (!exists("storage.backend", p)) p$storage.backend="bigmemory.ram"
    if (!exists("clusters", p)) p$clusters = rep("localhost", detectCores() )

    if (!exists("boundary", p)) p$boundary = FALSE 
    if (!exists("depth.filter", p)) p$depth.filter = 0 # depth (m) stats locations with elevation > 0 m as being on land (and so ignore)
    if (!exists("stm_quantile_bounds", p)) p$stm_quantile_bounds = c(0.025, 0.975) # remove these extremes in interpolations
    
    if (!exists("stm_rsquared_threshold", p)) p$stm_rsquared_threshold = 0.2 # lower threshold
    if (!exists("stm_distance_statsgrid", p)) p$stm_distance_statsgrid = 3 # resolution (km) of data aggregation (i.e. generation of the ** statistics ** )
    if (!exists("stm_distance_prediction", p)) p$stm_distance_prediction = p$stm_distance_statsgrid*0.75  # this is a half window km
    if (!exists("stm_distance_scale", p)) p$stm_distance_scale = 50 # km ... approx guess of 95% AC range 
    if (!exists("stm_distance_min", p)) p$stm_distance_min = 2 
    if (!exists("stm_distance_max", p)) p$stm_distance_max = 65
  
    if (!exists("n.min", p)) p$n.min = 120 # n.min/n.max changes with resolution must be more than the number of knots/edf
    # min number of data points req before attempting to model timeseries in a localized space
    if (!exists("n.max", p)) p$n.max = 8000 # no real upper bound
    p$sampling = c( 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.1, 1.25 )  # 

    if (!exists("variables", p)) p$variables = list( 
      Y = varname, 
      LOCS = c("plon", "plat"), 
      TIME = "tiyr", 
      COV = c("z", "dZ", "ddZ", "log.substrate.grainsize", "t", "tmean.climatology", "tsd.climatology", "ca1", "ca2", "mr", "Npred", "smr" ) )
    p$varnames = c( p$variables$LOCS, p$variables$COV ) 
 
    # additional variable to extract from aegis_db for inputs
    p$aegis_variables = list()
    for (id in c("speciescomposition", "speciesarea", "sizespectrum", "condition", "metabolism", "biochem") ) {
      pz = aegis::aegis_parameters( p=p, DS=id )
      pz_vars = intersect( pz$varstomodel, p$variables$COV )  # these are aegis vars to model
      if (length(pz_vars) > 0) p$aegis_variables[[id]] = pz_vars 
    }

    if (!exists("stm_variogram_method", p)) p$stm_variogram_method = "fast"
    if (!exists("stm_local_modelengine", p)) p$stm_local_modelengine ="gam"
    if (!exists("stm_global_modelengine", p)) p$stm_global_modelengine ="gam"

    if (!exists("stm_global_family", p)) p$stm_global_family = gaussian(link=log)

    # using covariates as a first pass essentially makes it ~ kriging with external drift .. no time or space here
    if (!exists("stm_global_modelformula", p)) p$stm_global_modelformula = formula( paste( 
      varname, ' ~ s(t, k=3, bs="ts") + s(tmean.climatology, k=3, bs="ts") + s(tsd.climatology, k=3, bs="ts")  ', 
      ' + s( log(z), k=3, bs="ts") + s( log(dZ), k=3, bs="ts") + s( log(ddZ), k=3, bs="ts") ',
      ' + s( log(mr), k=3, bs="ts") + s( Npred, k=3, bs="ts") + s( smr, k=3, bs="ts")  ',
      ' + s(log.substrate.grainsize, k=3, bs="ts") + s(ca1, k=3, bs="ts") + s(ca2, k=3, bs="ts")   ' ))  # no space 

    if (p$stm_local_modelengine =="twostep") {

      # this is the time component (mostly) .. space enters as a rough constraint 
      if (!exists("stm_local_modelformula", p))  p$stm_local_modelformula = formula( paste(
        varname, '~ s(yr, k=10, bs="ts") + s(cos.w, k=3, bs="ts") + s(sin.w, k=3, bs="ts") ', 
          ' + s(cos.w, sin.w, yr, bs="ts", k=20) ',
          ' + s(plon, k=3, bs="ts") + s(plat, k=3, bs="ts") + s(plon, plat, k=20, bs="ts") ' ) )
      if (!exists("stm_local_model_distanceweighted", p)) p$stm_local_model_distanceweighted = TRUE

      # this is the spatial component
      # p$stm_twostep_space = "spatial.process"
      # p$stm_twostep_space = "fft"
      # p$stm_twostep_space = "tps"
      if (!exists("stm_twostep_space", p))  p$stm_twostep_space = "krige"
      # if (!exists("stm_twostep_space", p))  p$stm_twostep_space = "tps"
      if (!exists("stm_gam_optimizer", p)) p$stm_gam_optimizer=c("outer", "bfgs") 

    }  else if (p$stm_local_modelengine == "habitat") {

      p$stm_global_family = binomial( link=log )
      
      if (!exists("stm_local_modelformula", p))  p$stm_local_modelformula = formula( paste(
        varname, '~ s(yr, k=10, bs="ts") + s(cos.w, k=3, bs="ts") + s(sin.w, k=3, bs="ts") ', 
          ' + s(cos.w, sin.w, yr, bs="ts", k=10)  ',
          ' + s(plon, k=3, bs="ts") + s(plat, k=3, bs="ts") + s(plon, plat, k=10, bs="ts") ' ) )

      if (!exists("stm_local_model_distanceweighted", p)) p$stm_local_model_distanceweighted = TRUE
      # if (!exists("stm_gam_optimizer", p)) p$stm_gam_optimizer="perf"
      if (!exists("stm_gam_optimizer", p)) p$stm_gam_optimizer=c("outer", "bfgs") 
    
    }  else if (p$stm_local_modelengine == "gam") {

      if (!exists("stm_local_modelformula", p))  p$stm_local_modelformula = formula( paste(
        varname, '~ s(yr, bs="ts") + s(cos.w, k=3, bs="ts") + s(sin.w, k=3, bs="ts") ', 
          ' + s(cos.w, sin.w, yr, bs="ts", k=25)  ',
          ' + s(plon, k=3, bs="ts") + s(plat, k=3, bs="ts") + s(plon, plat, k=25, bs="ts") ' ) )

      if (!exists("stm_local_model_distanceweighted", p)) p$stm_local_model_distanceweighted = TRUE
      # if (!exists("stm_gam_optimizer", p)) p$stm_gam_optimizer="perf"
      if (!exists("stm_gam_optimizer", p)) p$stm_gam_optimizer=c("outer", "bfgs") 
    
    }  else if (p$stm_local_modelengine == "bayesx") {
 
      # bayesx families are specified as characters, this forces it to pass as is and 
      # then the next does the transformation internal to the "stm__bayesx"

      # alternative models .. testing .. problem is that SE of fit is not accessible?
      p$stm_local_modelformula = formula( paste( 
        varname, ' ~ sx(yr, bs="ps") + sx(cos.w, bs="ps") + s(sin.w, bs="ps") +s(z, bs="ps") + sx(plon, bs="ps") + sx(plat,  bs="ps")', 
          ' + sx(plon, plat, cos.w, sin.w, yr, bs="te") ' )
          # te is tensor spline
      )
      p$stm_local_model_bayesxmethod="MCMC"
      p$stm_local_model_distanceweighted = FALSE
    
    } else {
    
      message( "The specified stm_local_modelengine is not tested/supported ... you are on your own ;) ..." )

    }

  }

  return(p)




}
