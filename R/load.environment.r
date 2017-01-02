

  # ----------------------------------------------------------------------------------
  # NOTE to all: The year of "year.assessment must be changed every year before any other run
  #       It cannot be automatically loaded together with the "load.snowcrab.environment". This is because
  #       running in parallel mode requires overrding some parameters in "p" on occasion which cannot be done cleanly
  #       as "load.snowcrab.environment" is sourced with every initialisation of a  new CPU.
  #       Copying the following into each relevent file is not a solution as it is error prone and  repetitive.
  # ----------------------------------------------------------------------------------

load.environment = function( year.assessment=NULL, libs=NULL, p=NULL ) {

  Sys.setlocale("LC_COLLATE", "C")   # turn off locale-specific sorting,

  if (is.null(p)) p = list()

  p$project.name = "bio.snowcrab"
  p$project.outdir.root = project.datadirectory( p$project.name, "R" ) #required for interpolations and mapping


  rlibs = RLibrary ( "geosphere", "lubridate", "mgcv", "parallel", "DBI", "Cairo", "Hmisc",
      "vegan", "akima", "fields", "lattice", "gstat", "maptools",  "boot", "grid",
      "RColorBrewer", "raster","rgdal", "sp", "rgeos", "bigmemory" ,"numDeriv")

  blibs = bioLibrary( "bio.snowcrab", "bio.spacetime", "bio.utilities",
      "bio.polygons", "bio.groundfish", "netmensuration", "bio.coastline",
      "bio.substrate", "bio.temperature", "bio.taxonomy", "bio.habitat", "bio.indicators",
      "bio.bathymetry" )

  if (!is.null(libs)) RLibrary(libs)
  if ( exists("libs", p) ) libs = c(libs, p$libs)

  p$libs = unique( c( libs, rlibs, blibs ) )

  if (is.null(year.assessment) ) year.assessment=as.numeric( substring(Sys.Date(),1,4))
  p$year.assessment = year.assessment
  p$seabird.yToload = 2012:p$year.assessment
  p$minilog.yToload = 1999:p$year.assessment
  p$netmind.yToload = 1999:p$year.assessment
  p$esonar.yToload  = 2014:p$year.assessment
  p$netmensuration.problems = c()

  p$annual.results = file.path( project.datadirectory("bio.snowcrab"), "assessments", p$year.assessment ) # output location for year-specific results

  p$spatial.domain = "snowcrab"
  p$ext2 = extent(matrix(c(-66.4, 42.2, -57.2, 47.4), nrow=2, ncol=2)) #MG extent of mapping frame
  p$extUTM = extent(matrix(c(219287.2, 4677581, 937584, 5265946), nrow=2, ncol=2)) #MG UTM extent of mapping frame
  p$geog.proj = "+proj=longlat +ellps=WGS84"
  p$annot.cex=2
  p$do.parallel = TRUE
  p$clusters = c( rep("localhost", 1) )

  p$fisheries.grid.resolution = 2
  ## these are kriging related parameters:: the method is deprecated
  p$ofname = file.path(p$annual.results, paste("TSresults", p$year.assessment, "rdata", sep=".") )
  p$regions.to.model = c( "cfanorth", "cfasouth", "cfa4x", "cfaall" )

  p$vars.to.model = variable.list.expand("all.to.model")
  p$years.to.model = c(1998:p$year.assessment)

  p$yearswithTdata = c(1950:p$year.assessment)
  p$recode.data = TRUE
  p$map.results=TRUE

  p$prediction.dyear = 9/12  # time of year as fractional year to predict 1 Sept

  p$nw = 10  # from temperature.r, number of intervals in a year
  p$default.spatial.domain = "canada.east"  # for temperature/habitat lookups

  p$kformula = as.formula( "kv ~ z + t + tamp + wmin + dZ + ddZ + log.substrate.grainsize" )  # model in 2006-2008
  p$klocs = as.formula ( "~plon+plat" )
  p$vgm.dist = unique(sort(c( seq(10, 60, 4), seq(50, 100, 10), seq( 80, 160, 20) )))
  p$knmax=100  # must be greater than 30 for convergence
  p$krige.incremental = FALSE
  p$plot.variogram = FALSE
  p$transgaussian.kriging = TRUE
  p$n.conditional.sims = 100
  p$threshold.distance = 5  # in km for merging fisheries data into the trawl data for external drift kriging
  p$optimizers = c(  "bfgs", "nlm", "perf", "newton", "Nelder-Mead" )  # used by GAM

  p$plottimes=c("annual", "globalaverage")
  p$conversions=c("ps2png")

  p = spatial_parameters( p=p ) # region and lon/lats, projections
  p = gmt.parameters( p=p )

  return(p)
}

