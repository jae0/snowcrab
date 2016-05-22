
parameter.list.snowcrab = function ( p=NULL, set="default" ) {
   
  if (is.null( p) ) p = list()

  if ( set =="default" ) {

    if (!exists("annual.results", p ) ) p$annual.results = file.path( project.datadirectory("snowcrab"), "assessments", current.assessment.year ) # output location for year-specific results
    
    if (!exists("spatial.domain", p ) ) p$spatial.domain = "snowcrab"

    if (!exists("ext2", p ) ) p$ext2 = extent(matrix(c(-66.4, 42.2, -57.2, 47.4), nrow=2, ncol=2)) #MG extent of mapping frame
    if (!exists("extUTM", p ) ) p$extUTM = extent(matrix(c(219287.2, 4677581, 937584, 5265946), nrow=2, ncol=2)) #MG UTM extent of mapping frame
    if (!exists("geog.proj", p ) ) p$geog.proj = "+proj=longlat +ellps=WGS84"

    if (!exists("annot.cex", p ) ) p$annot.cex=2

    # ---------------------------------------------------
    # clusters definition
    if (!exists("do.parallel", p ) )     p$do.parallel = TRUE
   #    p$clusters = c( rep("io",8), rep("tethys",2) )  # put the fast cpu's first
   #    p$clusters = c( rep("localhost",8) )  # use local cluster to make file usage easier without refreshing each filesystem contantly
    if (!exists("clusters", p ) )     p$clusters = c( rep("localhost",1) )  # use local cluster to make file usage easier without refreshing each filesystem contantly

    if (!exists("fisheries.grid.resolution", p ) )   p$fisheries.grid.resolution = 2 

    ## these are kriging related parameters:: the method is deprecated
    if (!exists("ofname", p ) )   p$ofname = file.path(p$annual.results, paste("TSresults", p$current.assessment.year, "rdata", sep=".") )
    if (!exists("regions.to.model", p ) )   p$regions.to.model = c( "cfanorth", "cfasouth", "cfa4x", "cfaall" )
    
    if (!exists("vars.to.model", p ) )   p$vars.to.model = variable.list.expand("all.to.model")
    if (!exists("years.to.model", p ) )   p$years.to.model = c(1998:p$current.assessment.year)
      
    if (!exists("yearswithTdata", p ) )   p$yearswithTdata = c(1950:p$current.assessment.year)  
    if (!exists("recode.data", p ) )  p$recode.data = TRUE
    if (!exists("map.results", p ) )   p$map.results=TRUE

    if (!exists("prediction.dyear", p ) )  p$prediction.dyear = 9/12  # time of year as fractional year to predict 1 Sept

    if (!exists("nw", p ) )  p$nw = 10  # from temperature.r, number of intervals in a year
    if (!exists("default.spatial.domain", p ) )     p$default.spatial.domain = "canada.east"  # for temperature/habitat lookups
    
      # p$kformula = as.formula( "kv ~ z + t + total.cpue" )  # model in 2006-2008
      # p$kformula = as.formula( "kv ~ z + t + tamp + wmin + dZ + ddZ + substrate.mean + total.cpue" )  # model in 2006-2008
     if (!exists("kformula", p ) )  p$kformula = as.formula( "kv ~ z + t + tamp + wmin + dZ + ddZ + substrate.mean" )  # model in 2006-2008
#    p$kformula = as.formula( "kv ~ z + t + tamp + wmin + dZ + ddZ + substrate.mean " )  # model in 2006-2008
     if (!exists("klocs", p ) )  p$klocs = as.formula ( "~plon+plat" )
     if (!exists("vgm.dist", p ) )  p$vgm.dist = unique(sort(c( seq(10, 60, 4), seq(50, 100, 10), seq( 80, 160, 20) )))
     if (!exists("knmax", p ) )  p$knmax=100  # must be greater than 30 for convergence
     if (!exists("krige.incremental", p ) ) p$krige.incremental = FALSE
     if (!exists("plot.variogram", p ) )  p$plot.variogram = FALSE
     if (!exists("transgaussian.kriging", p ) )  p$transgaussian.kriging = TRUE
     if (!exists("n.conditional.sims", p ) )  p$n.conditional.sims = 100
      # p$inflection.t = 2
      # p$inflection.z = 250
      # p$linearise.drift.terms = T
     if (!exists("threshold.distance", p ) )  p$threshold.distance = 5  # in km for merging fisheries data into the trawl data for external drift kriging
     if (!exists("optimizers", p ) )  p$optimizers = c(  "bfgs", "nlm", "perf", "newton", "Nelder-Mead" )  # used by GAM

      # default plotting time format
      if (!exists("plottimes", p ) ) p$plottimes=c("annual", "globalaverage")

      # default figure generation (from maps)
      if (!exists("conversions", p ) ) p$conversions=c("ps2png")
    }

  return (p)
}
