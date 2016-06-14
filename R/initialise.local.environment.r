

  # ----------------------------------------------------------------------------------
  # NOTE to all: The year of "year.assessment must be changed every year before any other run
  #       It cannot be automatically loaded together with the "load.snowcrab.environment". This is because
  #       running in parallel mode requires overrding some parameters in "p" on occasion which cannot be done cleanly
  #       as "load.snowcrab.environment" is sourced with every initialisation of a  new CPU.
  #       Copying the following into each relevent file is not a solution as it is error prone and  repetitive.
  # ----------------------------------------------------------------------------------

initialise.local.environment = function( year.assessment=NULL, libs=NULL, p=NULL ) {

  Sys.setlocale("LC_COLLATE", "C")   # turn off locale-specific sorting,

  if (is.null(p)) p = list()

  p$project.name = "bio.snowcrab"
  p$project.outdir.root = project.datadirectory( p$project.name, "R" ) #required for interpolations and mapping

  p$libs = RLibrary ( c( "geosphere", "lubridate", "mgcv", "parallel", "DBI", "Cairo", "Hmisc",
      "vegan", "akima", "fields", "lattice", "gstat", "maptools",  "boot", "raster", "grid",
      "RColorBrewer", "rasterVis", "rgdal", "sp", "rgeos", "bigmemory" ) )

  p$libs = unique( c( p$libs, bioLibrary( "bio.snowcrab", "bio.spacetime", "bio.utilities",
      "bio.polygons", "bio.groundfish", "netmensuration", "bio.coastline",
      "bio.substrate", "bio.temperature", "bio.taxonomy", "bio.habitat", "bio.indicators",
      "bio.bathymetry" )) )

  if (!is.null(libs)) p$libs = unique( c(p$libs, RLibrary(libs) ) )

  if (is.null(year.assessment) ) year.assessment=as.numeric( substring(Sys.Date(),1,4))
  p$year.assessment = year.assessment
  p$seabird.yToload = 2012:p$year.assessment
  p$minilog.yToload = 1999:p$year.assessment
  p$netmind.yToload = 1999:p$year.assessment
  p$esonar.yToload  = 2014:p$year.assessment

  p = parameter.list.snowcrab ( p=p )
  p = spatial.parameters( p=p ) # region and lon/lats, projections
  p = gmt.parameters( p=p )

  return(p)
}

