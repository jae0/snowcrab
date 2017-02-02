
snowcrab.parameters = function( p=NULL, DS="default", current.year=NULL, varname=NULL ) {

  if ( is.null(p) ) p=list()

  if ( DS=="default") {

    p$project.name = "bio.snowcrab"
    p$project.root = file.path( project.datadirectory( "bio.indicators"), p$project.name )
    p$project.outdir.root = project.datadirectory( p$project.name, "R" ) #required for interpolations and mapping

    rlibs = c( "lubridate", "rgdal", "parallel", "sp", "lattice", "fields", "mgcv" ,
               "geosphere",  "DBI", "Cairo", "Hmisc", "vegan", "akima",   "latticeExtra",  "maptools",  
               "boot", "grid", "RColorBrewer",  "spatstat", "rgeos", "bigmemory" ,"numDeriv")

    p$libs = c( p$libs, RLibrary ( rlibs ) )
    p$libs = c( p$libs, bioLibrary (
      "bio.base", "bio.utilities", "bio.taxonomy", "bio.spacetime", "bio.polygons",  "netmensuration", 
      "bio.coastline",  "bio.bathymetry", "bio.temperature", "bio.substrate", "bio.groundfish", 
      "bio.snowcrab", "bio.indicators" ) )
    p$libs = unique( p$libs )

    if (is.null(current.year) ) stop( "must define current.year" )

    p$current.year = current.year  # this is a synonym for year.assessment ... will eventually remove one of these
    p$year.assessment = current.year

    p$yearstomodel = 1998:p$current.year
    p$seabird.yToload = 2012:p$current.year
    p$minilog.yToload = 1999:p$current.year
    p$netmind.yToload = 1999:p$current.year
    p$esonar.yToload  = 2014:p$current.year
    p$netmensuration.problems = c()

    p$yrs = c(1998:p$current.year)  # same as yearstomodel ... need to clean this up:: TODO
    p$ny = length(p$yrs)
    p$nt = p$ny # must specify, else assumed = 1 (1= no time)  ## nt=ny annual time steps, nt = ny*nw is seassonal
    p$nw = 10 # default value of 10 time steps for all temp and indicators
    p$tres = 1/ p$nw # time resolution .. predictions are made with models that use seasonal components
    p$dyears = (c(1:p$nw)-1)  / p$nw # intervals of decimal years... fractional year breaks
    p$dyear_centre = p$dyears[ round(p$nw/2) ] + p$tres/2
    # used for creating timeslices and predictions  .. needs to match the values in indicators.parameters()
    p$prediction.dyear = lubridate::decimal_date( lubridate::ymd("0000/Sep/01")) 
    # output timeslices for predictions in decimla years, yes all of them here
    p$prediction.ts = p$yrs + p$prediction.dyear 


    p$data.sources = c("groundfish", "snowcrab")
    p$spatial.domain = "snowcrab"
    p$spatial.domain.subareas = NULL
    p = spatial_parameters( p=p )  # data are from this domain .. so far

    # output location for year-specific results
    p$annual.results = file.path( project.datadirectory("bio.snowcrab"), "assessments", p$year.assessment ) 
    p$ofname = file.path(p$annual.results, paste("TSresults", p$year.assessment, "rdata", sep=".") )
    
    p$annot.cex=2
    p$do.parallel = TRUE

    p$fisheries.grid.resolution = 2
   
    p$regions.to.model = c( "cfanorth", "cfasouth", "cfa4x", "cfaall" )
    p$plottimes=c("annual", "globalaverage")
    p$conversions=c("ps2png")
    p$recode.data = TRUE

  
    if (1) {
      #     message( "This is marked for removal ..." ) # TODO

      p$ext2 = extent(matrix(c(-66.4, 42.2, -57.2, 47.4), nrow=2, ncol=2)) #MG extent of mapping frame
      p$extUTM = extent(matrix(c(219287.2, 4677581, 937584, 5265946), nrow=2, ncol=2)) #MG UTM extent of mapping frame
      p$geog.proj = "+proj=longlat +ellps=WGS84"

      ## these are kriging related parameters:: the method is deprecated
      p$years.to.model = c(1998:p$year.assessment)
      p$threshold.distance = 5  # in km for merging fisheries data into the trawl data for external drift kriging
      p$optimizers = c(  "bfgs", "nlm", "perf", "newton", "Nelder-Mead" )  # used by GAM
      p = gmt.parameters( p=p )

    }

  }


  if (DS=="lbm") {

    p$libs = RLibrary( c( p$libs, "lbm" ) ) # required for parallel processing
    p$storage.backend="bigmemory.ram"
    if (!exists("clusters", p)) p$clusters = rep("localhost", detectCores() )

    p$boundary = FALSE 
    p$depth.filter = 0 # depth (m) stats locations with elevation > 0 m as being on land (and so ignore)
    p$lbm_quantile_bounds = c(0.01, 0.99) # remove these extremes in interpolations
    
    p$lbm_rsquared_threshold = 0.1 # lower threshold
    p$lbm_distance_prediction = 7.5 # this is a half window km
    p$lbm_distance_statsgrid = 5 # resolution (km) of data aggregation (i.e. generation of the ** statistics ** )
    p$lbm_distance_scale = 50 # km ... approx guess of 95% AC range 
    p$lbm_distance_min = p$lbm_distance_statsgrid 
    p$lbm_distance_max = 75
  
    p$n.min = 100 # n.min/n.max changes with resolution must be more than the number of knots/edf
    # min number of data points req before attempting to model timeseries in a localized space
    p$n.max = 8000 # no real upper bound
    p$sampling = c( 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.1, 1.2, 1.5 )  # 

    p$variables = list( 
      Y = varname, 
      LOCS = c("plon", "plat"), 
      TIME = "tiyr", 
      COV = c("z", "dZ", "ddZ", "log.substrate.grainsize", "t", "tmean", "tamplitude", "ca1", "ca2" ) )
    p$varnames = c( p$variables$LOCS, p$variables$COV ) 

    p$indicators.variables = c("ca1", "ca2")  # additional variable to extract from indicators.db for inputs

    if (!exists("lbm_variogram_method", p)) p$lbm_variogram_method = "fast"
    
    p$lbm_local_modelengine ="twostep"
    # if (!exists("lbm_local_modelengine", p)) p$lbm_local_modelengine = "krige" # "twostep" might be interesting to follow up

    # using covariates as a first pass essentially makes it ~ kriging with external drift
    p$lbm_global_modelengine = "gam"
    p$lbm_global_modelformula = formula( paste( 
      varname, ' ~ s(yr) + s(dyear, k=3, bs="ts") + s(yr, dyear, k=36, bs="ts") ',
      ' + s(ca1, bs="ts") + s(ca2, bs="ts") ', 
      ' + s(t, bs="ts") + s(tmean, bs="ts") + s(tamplitude, bs="ts") + s(z, bs="ts")',
      ' + s(dZ, bs="ts") + s(ddZ, bs="ts")  + s(log.substrate.grainsize, bs="ts") ' )) 

    p$lbm_global_family = gaussian()
    p$lbm_local_family = gaussian()



    if (p$lbm_local_modelengine =="twostep") {

      # this is the time component (mostly) .. space enters as a rough constraint 
      p$lbm_local_modelformula = formula( paste(
        varname, '~ s(yr, k=5, bs="ts") + s(cos.w, k=3, bs="ts") + s(sin.w, k=3, bs="ts") ', 
          ' + s(cos.w, sin.w, yr, bs="ts", k=36) ',
          ' + s(plon, bs="ts") + s(plat, bs="ts") + s(plon, plat, k=25, bs="ts") ' ) )
      p$lbm_local_model_distanceweighted = TRUE

      # this is the spatial component
      # p$lbm_twostep_space = "spatial.process"
      # p$lbm_twostep_space = "fft"
      # p$lbm_twostep_space = "tps"
      p$lbm_twostep_space = "krige"

    }  else if (p$lbm_local_modelengine == "gam") {

      p$lbm_local_modelformula = formula( paste(
        varname, '~ s(yr, k=5, bs="ts") + s(cos.w, k=3, bs="ts") + s(sin.w, k=3, bs="ts") ', 
          ' + s(cos.w, sin.w, yr, bs="ts", k=36) ',
          ' + s(plon, bs="ts") + s(plat, bs="ts") + s(plon, plat, k=25, bs="ts") ' ) )    

      p$lbm_local_model_distanceweighted = TRUE
      p$lbm_gam_optimizer="perf"
      # p$lbm_gam_optimizer=c("outer", "bfgs") 
    
    }  else if (p$lbm_local_modelengine == "bayesx") {
 
      # bayesx families are specified as characters, this forces it to pass as is and 
      # then the next does the transformation internal to the "lbm__bayesx"
      p$lbm_local_family_bayesx = "gaussian" 

      # alternative models .. testing .. problem is that SE of fit is not accessible?
      p$lbm_local_modelformula = formula( paste( 
        varname, ' ~ sx(yr,   bs="ps") + sx(cos.w, bs="ps") + s(sin.w, bs="ps") +s(z, bs="ps") + sx(plon, bs="ps") + sx(plat,  bs="ps")', 
          ' + sx(plon, plat, cos.w, sin.w, yr, bs="te") ' )
          # te is tensor spline
      )
      p$lbm_local_model_bayesxmethod="MCMC"
      p$lbm_local_model_distanceweighted = FALSE
    
    } else {
    


      message( "The specified lbm_local_modelengine is not tested/supported ... you are on your own ;) ..." )

    }

  }

  return(p)




}
