
  # Figures obtained after completion of data assimilation and processing up to the end of "01.snowcrab.r"
  require(bio.base)
  p = bio.snowcrab::load.environment( year.assessment=2016)
  #loadfunctions('bio.snowcrab')


   # Think this is fixed now .. ?
   p$do.parallel = FALSE  # mapping in parallel is broken .. must fix ::TODO
  #p$do.parallel = TRUE  # mapping in parallel is broken .. must fix ::TODO

   p$clusters = rep("localhost", 24 )
   p$clusters = rep("localhost", 3 )

  # ------------------------------------------
   # Time-series: Fisheries landings
   figure.landings.timeseries( yearmax=p$year.assessment, outdir=file.path( p$annual.results,  "timeseries","fishery"), outfile="landings.ts", outfile2="landings.ts.sm" )

  # ------------------------------------------
  # Time-series: Fisheries effort
   figure.effort.timeseries( yearmax=p$year.assessment, outdir=file.path( p$annual.results,"timeseries", "fishery"), outfile="effort.ts", outfile2="effort.ts.sm" )

  # ------------------------------------------
  # Time-series: Fisheries CPUE
   figure.cpue.timeseries( yearmax=p$year.assessment, outdir=file.path( p$annual.results,"timeseries", "fishery"), outfile="cpue.ts", outfile2="cpue.sm.ts" )

  # ------------------------------------------
  # Size frequency distributions, broken down by moult category from at-sea observed data
    figure.observed.size.freq( regions = c("cfanorth", "cfasouth", "cfa4x"), years="all", outdir=file.path( p$annual.results, "figures", "size.freq", "observer")  )

  # ------------------------------------------
  # Size-frequency distributions of snow crab cw from trawl data, broken down by maturity classes
    histograms.size.maturity.update( outdir=file.path( p$annual.results, "figures", "size.freq", "survey"),  redo.data=T )
    histograms.size.maturity.single.area( outdir=file.path( p$annual.results, "figures", "size.freq", "survey"),  area='cfa4x',redo.data=T )

  # ------------------------------------------
  # Timeseries of all survey variables
    figure.timeseries.males( outdir=file.path(p$annual.results, "timeseries", "survey") )
    figure.timeseries.females( outdir=file.path(p$annual.results, "timeseries", "survey") )

  # ------------------------------------------
  # Timeseries of all survey variables
    #MG fix this one
    figure.timeseries.survey(p, areas=c("cfanorth", "cfasouth", "cfa4x"), from.file=F ) # goes to file.path( project.datadirectory("bio.snowcrab"), "R", "timeseries", "survey"

  #-----------------------------------------------
  # Time series of survey temperature
    figure.timeseries.raw.survey.temperature (outdir=file.path(p$annual.results, "timeseries", "survey"),plotyears=2001:2016)

  # ------------------------------------------
  # Timeseries: geometric mean density of R-1 recruits
    #figure.timeseries.recruits( outdir=file.path(p$annual.results, "timeseries", "survey") )
    #MG Fix R2
    figure.timeseries.R2( outdir=file.path(p$annual.results, "timeseries", "survey") )
    #MG Try to get a better measure of R1, not being used at this point because it seems to poorly represent recruits
    figure.timeseries.R1( outdir=file.path(p$annual.results, "timeseries", "survey") )

  #Timeseries: geometric mean density of R0 (fishable Biomass)
    figure.timeseries.R0( outdir=file.path(p$annual.results, "timeseries", "survey") )

  #To represent the reduced set of stations in 2014
  #figure.timeseries.R0.reduced.stations( outdir=file.path(p$annual.results, "timeseries", "survey") )

  # ------------------------------------------
  # Timeseries: geometric mean density of CC5 crab
    figure.timeseries.CC5( outdir=file.path(p$annual.results,  "timeseries", "survey") )

  # ------------------------------------------
  # Timeseries: sex ratios of immature crab
    figure.timeseries.sexratios( outdir=file.path(p$annual.results,  "timeseries", "survey" ), type="immature" )

  # ------------------------------------------
  # Timeseries: sex ratios of mature crab
    figure.timeseries.sexratios( outdir=file.path(p$annual.results,  "timeseries", "survey"), type="mature" )

  # ------------------------------------------
  # Timeseries: geometric mean carapace width from trawl surveys
    figure.timeseries.CW( outdir=file.path(p$annual.results,  "timeseries", "survey"), type="trawl" )

  # ------------------------------------------
  # Timeseries: geometric mean carapace width from at-sea observers
    figure.timeseries.CW( outdir=file.path(p$annual.results,  "timeseries", "survey"), type="observer" )

  #Timeseries: geometric mean biomass of by-catch from snow crab survey
    #cod, haddock, halibut, plaice, wolfish, thornyskate, smoothskate, winterskate, northernshrimp, jonahcrab, lessertoadcrab
    species = c(10, 11, 30, 40, 201, 50, 2521, 2511, 202, 204, 2211)
    figure.timeseries.bycatch(species, plotyears=2004:2016,outdir=file.path(p$annual.results, "timeseries", "survey")) 


   ##Predators
   #figure.timeseries.bycatch.halibut(outdir=file.path(p$annual.results, "timeseries", "survey")) #30
   #figure.timeseries.bycatch.cod(outdir=file.path(p$annual.results, "timeseries", "survey")) #10
   #figure.timeseries.bycatch.wolfish(outdir=file.path(p$annual.results, "timeseries", "survey")) #50
   #figure.timeseries.bycatch.thornyskate(outdir=file.path(p$annual.results, "timeseries", "survey")) #201
   #figure.timeseries.bycatch.smoothskate(outdir=file.path(p$annual.results, "timeseries", "survey")) #202
   #figure.timeseries.bycatch.winterskate(outdir=file.path(p$annual.results, "timeseries", "survey")) #204

   ##Competitors
   #figure.timeseries.bycatch.lessertoadcrab(outdir=file.path(p$annual.results, "timeseries", "survey")) #2512
   #figure.timeseries.bycatch.jonahcrab(outdir=file.path(p$annual.results, "timeseries", "survey")) #2511

   ##Prey
   #figure.timeseries.bycatch.northernshrimp(outdir=file.path(p$annual.results, "timeseries", "survey")) #2211

  # ------------------------------------------
  # Map: Basemap of the Scotian Shelf used by all other mapping routines
  #   creating a partial postscript file via GMT
  #   .. only required if changing resolution or spatial extent
  #  gmt.basemap (p)

    # ------------------------------------------
  # Map: Scotian Shelf with CFA lines and labels  .. using gmt
  # this is the basemap from map.r which is then post-labelled in sodipodi
  #  p$outdir = file.path(p$annual.results,"figures")
  #  p$outfile.basename = file.path(p$outdir, "map.CFAs")
    # p$basemap = file.path( project.datadirectory("bio.snowcrab"), "R", p$basemap)
  #  map.basemap.with.cfa.lines( p, conversions=c("ps2png")  )

  # ------------------------------------------
  # Map:  Interpolated mean/geometric mean of various variables in the set data table
  p$do.parallel=F
  outdir = file.path( project.datadirectory("bio.snowcrab"), "R", "maps", "survey","snowcrab","annual" ) 

    # just for the roadshow
    map.set.information( p, variables=c('totmass.male.com', 'totmass.female.mat'),mapyears=2014:2016,outdir=outdir) #,plot.method='gmt')
    map.set.information( p, variables='t',mapyears=2014:2016,outdir=outdir,log.variable=F,add.zeros=F,theta=100)   
    
    # all variables (geometric means)
    map.set.information( p, outdir=outdir)

    # Means
    # variables that shouldn't be logged 
    set = snowcrab.db( DS="set.biologicals")
    variables = variable.list.expand("all.data")
    variables = intersect( variables, names(set) )
 
    nolog.variables = c("t","z","sexratio.all","sexratio.mat","sexratio.imm","julian",variables[grep("cw",variables)])
    map.set.information( p, variables=nolog.variables,outdir=outdir,log.variable=F,add.zeros=F,theta=100)   

    # Geometric Means
    # all except variables that shouldn't be logged
    map.set.information( p, variables= variables[!variables%in%nolog.variables],outdir=outdir)

    #map.set.information.diff( p, outdir=file.path( project.datadirectory("bio.snowcrab"), "R", "maps", "survey.diff" )  )


  # ------------------------------------------
  # Map: Numerical density of by-catch species
   p$do.parallel=F
    map.cat.information( p, outdir=file.path( project.datadirectory("bio.snowcrab"), "R", "maps", "species" ) )

  # ------------------------------------------
  p=bio.snowcrab::load.environment()
  # Map:Fisheries logbook data (Effort, CPUE, Landings)
  # MG: This is not working properly, you must open the script and run through once by hand???
  # MG: To fix
    #map.fisheries.data( p)
    raster.map.variables( p, grid.fun=mean, variables = c("effort", "landings", "cpue"))
  # Map: fisheries logbook data (Effort, CPUE, Landings)  with a colour scale stretched only over the past two years
    raster.map.variables.year.assessment( p, grid.fun=mean, variables = c("effort", "landings", "cpue"), years=c('2014', '2015'))

  # ------------------------------------------
  # Map: Survey locations

    map.survey.locations( p, basedir=file.path(project.datadirectory("bio.snowcrab"), "R", "maps", "survey.locations"),  newyear=F, map.method="lattice"  )
    map.survey.locations( p, basedir=file.path(project.datadirectory("bio.snowcrab"), "R", "maps", "survey.locations"),  newyear=F, map.method="gmt"  )
    map.survey.locations( p, basedir=file.path(project.datadirectory("bio.snowcrab"), "R", "maps", "survey.locations"),  newyear=F, map.method="googleearth"  )

  # ------------------------------------------
  # Map: Observer locations
    map.observer.locations( p, basedir=file.path(project.datadirectory("bio.snowcrab"), "R", "maps","observer.locations" ), newyear=F , map.method="lattice"  )

  # ------------------------------------------
  # Map: Logbook recorded locations
    map.logbook.locations( p, basedir=file.path(project.datadirectory("bio.snowcrab"), "R", "maps","logbook.locations" ), newyear=F , map.method="lattice"  )

  # ------------------------------------------
  # Map: Crab movement from mark-recapture data
    #MG I think Brent is primarily mapping this stuff now. Not sure the data has been updated in a while
    map.movement( p, outdir=file.path(project.datadirectory("bio.snowcrab"), "R", "maps", "mark.recapture") )

  # ------------------------------------------
  # Map: Spatial representation of maturity patterns of snow crab
    #MG Not sure we use these maps either, check with Adam and Jae
    map.maturity( p, outdir=file.path(project.datadirectory("bio.snowcrab"), "R", "maps", "maturity"), newyear=T )

  ##########################################################################
  ###############################  Retired figures #########################
  # ------------------------------------------
  # Habitat usage comparisons (bivariate) ... requires the full "set.rdata" database and "logbook.dZ.rdata" database
    habitat.usage( usevar="totno.all", covariate="depth", outdir = file.path(p$annual.results, "habitat.templates") )

  # ------------------------------------------
  # Habitat usage comparisons (bivariate) ... requires the full "set.rdata" database and "logbook.dZ.rdata" database
    habitat.usage( usevar="totno.all", covariate="temperature", outdir = file.path(p$annual.results, "habitat.templates") )

  # ------------------------------------------
  # Habitat usage comparisons (bivariate) ... requires the full "set.rdata" database and "logbook.dZ.rdata" database
    habitat.usage( usevar="totno.all", covariate="bottom.slope", outdir = file.path(p$annual.results, "habitat.templates") )

  # ------------------------------------------
  # Habitat usage comparisons (bivariate) ... requires the full "set.rdata" database and "logbook.dZ.rdata" database
    habitat.usage( usevar="totno.all", covariate="bottom.curvature", outdir = file.path(p$annual.results, "habitat.templates") )

  # ------------------------------------------
  # Habitat usage comparisons (bivariate) ... requires the full "set.rdata" database and "logbook.dZ.rdata" database
    habitat.usage( usevar="totno.all", covariate="substrate", outdir = file.path(p$annual.results, "habitat.templates") )


  # ------------------------------------------
  # Timeseries: Larval brachyura from the SSIP data
    figure.timeseries.larvae( outdir=file.path(project.datadirectory("bio.snowcrab"), "R", "timeseries", "larvae") )

  # ------------------------------------------
  # Growth as a a function of instar for Scotian Shelf snow crab
    figure.growth.instar( outdir=file.path(project.datadirectory("bio.snowcrab"), "R", "growth") )


  # ------------------------------------------
  # Map: Larval distributions from the Scotian Shelf Ichtyoplankton Program data
    map.larvae( p, outdir=file.path(project.datadirectory("bio.snowcrab"), "R", "maps", "larvae"), conversions=conversions )


