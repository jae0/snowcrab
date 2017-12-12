
snowcrab.parameters = function( p=NULL, year.assessment=NULL, ... ) {

  # ---------------------
  # deal with additional passed parameters
  if ( is.null(p) ) p=list()
  p_add = list(...)
  if (length(p_add) > 0 ) p = c(p, p_add)
  i = which(duplicated(names(p), fromLast=TRUE))
  if ( length(i) > 0 ) p = p[-i] # give any passed parameters a higher priority, overwriting pre-existing variable

  # ---------------------
  # create/update library list
  p$libs = c( p$libs, RLibrary ( "colorspace",  "geosphere", "lattice",
    "maps", "mapdata", "maptools", "parallel",  "rgdal", "rgeos",  "sp", "splancs", "rgeos", "bigmemory", "numDeriv", "lubridate", "parallel", "fields", "mgcv" ) )
  p$libs = c( p$libs, project.library (
    "aegis.env", "bio.taxonomy", "stmv", "aegis",  "netmensuration", "bio.snowcrab" ) )
  p$libs = unique( p$libs )


  # ---------------------
  if (!exists("project.name", p)) p$project.name = "bio.snowcrab"
  if (!exists("data_root", p)) p$data_root = project.datadirectory( p$project.name ) 
  
  p$project.outputdir = project.datadirectory( p$project.name, "output" ) #required for interpolations and mapping
  p$transform_lookup = file.path( p$project.outputdir, "transform.lookup.rdata" ) # local storage of transforms for timeseries plots

  # ---------------------
  # define focal year. not required for pure spatial models but ignored by the spatial methods anyways
  if (!is.null(year.assessment) ) p$year.assessment = year.assessment
  if (!exists( "year.assessment", p)) p$year.assessment = lubridate::year(lubridate::now())

  
  if (!exists("yrs", p)) p$yrs = c(1999:p$year.assessment)

  p$seabird.yToload = intersect( p$yrs, 2012:p$year.assessment)
  p$minilog.yToload = intersect( p$yrs, 1999:p$year.assessment)
  p$netmind.yToload = intersect( p$yrs, 1999:p$year.assessment)
  p$esonar.yToload  = intersect( p$yrs, 2014:p$year.assessment)
  p$netmensuration.problems = c()

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

  p = spatial_parameters( p=p, spatial.domain="snowcrab" )  # data are from this domain .. so far

  p$spatial.domain.subareas = NULL # add cfa's as subareas .. TODO
  p$data.sources = c("groundfish", "snowcrab")

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

  return(p)

}
