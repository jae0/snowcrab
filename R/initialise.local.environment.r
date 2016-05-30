

  # ----------------------------------------------------------------------------------
  # NOTE to all: The year of "current.assessment.year must be changed every year before any other run
  #       It cannot be automatically loaded together with the "load.snowcrab.environment". This is because
  #       running in parallel mode requires overrding some parameters in "p" on occasion which cannot be done cleanly
  #       as "load.snowcrab.environment" is sourced with every initialisation of a  new CPU.
  #       Copying the following into each relevent file is not a solution as it is error prone and  repetitive.
  # ----------------------------------------------------------------------------------

initialize.local.environment = function( current.assessment.year=NULL, libs=NULL, p=NULL ) {

  if (is.null(p)) p = list()

  p$project.name = "snowcrab"
  p$project.outdir.root = project.datadirectory( p$project.name, "R" ) #required for interpolations and mapping

  p$libs = RLibrary ( c( "geosphere", "lubridate", "mgcv", "parallel", "DBI", "Cairo", "Hmisc", "chron",
      "vegan", "akima", "fields", "lattice", "gstat", "maptools",  "boot", "raster", "grid",
      "RColorBrewer", "rasterVis", "rgdal", "sp", "rgeos", "bigmemory" ) )

  p$libs = unique( c( p$libs, ecomodLibrary( "snowcrab", "ecomod_spacetime", "ecomod_utilities", "ecomod_parallel",
      "ecomod_polygons", "groundfish", "netmensuration", "ecomod_coastline",
      "ecomod_substrate", "ecomod_temperature", "ecomod_taxonomy", "ecomod_habitat", "ecomod_habitatsuitability",
      "ecomod_bathymetry", "ecomod_plottingmethods" )) )

  if (!is.null(libs)) p$libs = unique( c(p$libs, RLibrary(libs) ) )

  if (is.null(current.assessment.year) ) current.assessment.year=as.numeric( substring(Sys.Date(),1,4))
  p$current.assessment.year = current.assessment.year

  p = parameter.list.snowcrab ( p=p )
  p = spatial.parameters( p ) # region and lon/lats, projections
  p = gmt.parameters( p )

  return(p)
}

