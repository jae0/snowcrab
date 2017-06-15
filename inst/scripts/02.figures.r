
  # Figures obtained after completion of data assimilation and processing up to the end of "01.snowcrab.r"
  p = bio.snowcrab::load.environment( year.assessment=2016)
  #loadfunctions('bio.snowcrab')

#below lines set parallel mapping / figures
   p$do.parallel = FALSE  # mapping in parallel is broken .. must fix ::TODO
  #p$do.parallel = TRUE  # mapping in parallel is broken .. must fix ::TODO

# If parallel=T, below lines dtermine number of clusters to dedicate  
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
  #BZ- TODO add n= to the individual histograms
    figure.observed.size.freq( regions = c("cfanorth", "cfasouth", "cfa4x"), years="all", outdir=file.path( p$annual.results, "figures", "size.freq", "observer")  )

  # ------------------------------------------
  # Size-frequency distributions of snow crab cw from trawl data, broken down by maturity classes
    histograms.size.maturity.update( outdir=file.path( p$annual.results, "figures", "size.freq", "survey"),  redo.data=T )
    histograms.size.maturity.single.area( outdir=file.path( p$annual.results, "figures", "size.freq", "survey"),  area='cfa4x',redo.data=T )

  # ------------------------------------------
  # Timeseries of all survey variables
  figure.timeseries.survey(outdir=file.path(p$annual.results, "timeseries", "survey"),variables="R0.mass",plotyears=2001:p$year.assessment) # just R0 to see
  #figure.timeseries.survey(outdir=file.path(p$annual.results, "timeseries", "survey"),variables=c("sexratio.all","sexratio.mat","sexratio.imm")) 
  figure.timeseries.survey(outdir=file.path(p$annual.results, "timeseries", "survey"),plotyears=2001:p$year.assessment) # all variables
  figure.timeseries.survey(outdir=file.path(p$annual.results, "timeseries", "observer"),plotyears=2001:p$year.assessment,type='observer') 
  figure.timeseries.survey(outdir=file.path(p$annual.results, "timeseries", "survey"),type='groundfish.t') # groundfish survey temperature
  #-----------------------------------------------
  
  #Timeseries: geometric mean biomass of by-catch from snow crab survey

  # predators and competitors
    #cod, haddock, halibut, plaice, wolfish, thornyskate, smoothskate, winterskate, northernshrimp, jonahcrab, lessertoadcrab
    species = c(10, 11, 30, 40, 201, 50, 2521, 2511, 202, 204, 2211)
    figure.timeseries.bycatch(species, plotyears=2004:p$year.assessment,outdir=file.path(p$annual.results, "timeseries", "survey")) 

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
  p$corners = data.frame(plon=c(220, 990), plat=c(4750, 5270) )

  outdir = file.path( project.datadirectory("bio.snowcrab"), "R", "maps", "survey","snowcrab","annual" ) 

  #BZ TODO add a variable to p for mapyears  
  # just for the roadshow
    map.set.information( p, variables=c('totmass.male.com', 'totmass.female.mat'),mapyears=2014:p$year.assessment,outdir=outdir) #,plot.method='gmt')
    map.set.information( p, variables='t',mapyears=2014:p$year.assessment,outdir=outdir,log.variable=F,add.zeros=F,theta=100)   
    
    # bycatch (geometric means)
    bc.vars = c(paste("ms.mass",species,sep='.'),paste("ms.no",species,sep='.'))
    map.set.information( p, variables=bc.vars, mapyears=2014:p$year.assessment, outdir=outdir,probs=c(0,0.975)) #


    # all variables (geometric means)
    #map.set.information( p, outdir=outdir) # takes a long time

    # Means
    # variables that shouldn't be logged 
    set = snowcrab.db( DS="set.biologicals")
    variables = variable.list.expand("all.data")
    variables = intersect( variables, names(set) )
 
    nolog.variables = c("t","z","sexratio.all","sexratio.mat","sexratio.imm","julian",variables[grep("cw",variables)])
    map.set.information( p, variables=nolog.variables,outdir=outdir,log.variable=F,add.zeros=F,theta=100)   
    # logit transform for ratios
    map.set.information( p, variables=c("sexratio.all","sexratio.mat","sexratio.imm"),outdir=outdir,log.variable=F,add.zeros=F,theta=100) 

    # Geometric Means
    # all except variables that shouldn't be logged
    mass.vars = variables[!variables%in%nolog.variables][grep('mass',variables[!variables%in%nolog.variables])]
    no.vars = variables[!variables%in%nolog.variables][grep('no',variables[!variables%in%nolog.variables])]
    map.set.information( p, variables= mass.vars,outdir=outdir)
    map.set.information( p, variables= no.vars,outdir=outdir,probs=c(0,0.975))

    #map.set.information.diff( p, outdir=file.path( project.datadirectory("bio.snowcrab"), "R", "maps", "survey.diff" )  )


  # ------------------------------------------
  # Map: Survey locations

    map.survey.locations( p, basedir=file.path(project.datadirectory("bio.snowcrab"), "R", "maps", "survey.locations"),  newyear=F, map.method="lattice"  )
   # map.survey.locations( p, basedir=file.path(project.datadirectory("bio.snowcrab"), "R", "maps", "survey.locations"),  newyear=F, map.method="gmt"  )
  #  map.survey.locations( p, basedir=file.path(project.datadirectory("bio.snowcrab"), "R", "maps", "survey.locations"),  newyear=F, map.method="googleearth"  )

  # ------------------------------------------
  # Map: Observer locations
    map.observer.locations( p, basedir=file.path(project.datadirectory("bio.snowcrab"), "R", "maps","observer.locations" ), newyear=F , map.method="lattice"  )

  # ------------------------------------------
  # Map: Logbook recorded locations
    map.logbook.locations( p, basedir=file.path(project.datadirectory("bio.snowcrab"), "R", "maps","logbook.locations" ), newyear=F , map.method="lattice"  )




  # Map: Logbook data
  outdir = file.path( project.datadirectory("bio.snowcrab"), "R", "maps", "logbook","snowcrab","annual" ) 
  p$corners = data.frame(plon=c(220, 990), plat=c(4750, 5270) )
    map.fisheries.data( p, variable= 'effort',outdir=outdir,FUN=sum,probs=c(0,0.975))
    map.fisheries.data( p, variable= 'cpue',outdir=outdir,FUN=mean,probs=c(0,0.975))
    map.fisheries.data( p, variable= 'landings',outdir=outdir,FUN=sum,probs=c(0,0.975))










  # ------------------------------------------
  # Map: Numerical density of by-catch species
   p$do.parallel=F
    map.cat.information( p, outdir=file.path( project.datadirectory("bio.snowcrab"), "R", "maps", "species" ) )

  # ------------------------------------------
  p=bio.snowcrab::load.environment()
  # Map:Fisheries logbook data (Effort, CPUE, Landings)
  #BZ 2017. raster approach deprecated. consider removing in 2018
    # MG: This is not working properly, you must open the script and run through once by hand???
  # MG: To fix
    #map.fisheries.data( p)
    #raster.map.variables( p, grid.fun=mean, variables = c("effort", "landings", "cpue"))
  # Map: fisheries logbook data (Effort, CPUE, Landings)  with a colour scale stretched only over the past two years
    #raster.map.variables.year.assessment( p, grid.fun=mean, variables = c("effort", "landings", "cpue"), years=c('2014', '2015'))

  # ------------------------------------------
  # Map: Survey locations

    map.survey.locations( p, basedir=file.path(project.datadirectory("bio.snowcrab"), "R", "maps", "survey.locations"),  newyear=F, map.method="lattice"  )
    #map.survey.locations( p, basedir=file.path(project.datadirectory("bio.snowcrab"), "R", "maps", "survey.locations"),  newyear=F, map.method="gmt"  )
   # map.survey.locations( p, basedir=file.path(project.datadirectory("bio.snowcrab"), "R", "maps", "survey.locations"),  newyear=F, map.method="googleearth"  )

  # ------------------------------------------
  # Map: Observer locations
    map.observer.locations( p, basedir=file.path(project.datadirectory("bio.snowcrab"), "R", "maps","observer.locations" ), newyear=F , map.method="lattice"  )

  # ------------------------------------------
  # Map: Logbook recorded locations
    map.logbook.locations( p, basedir=file.path(project.datadirectory("bio.snowcrab"), "R", "maps","logbook.locations" ), newyear=F , map.method="lattice"  )

  # ------------------------------------------
  
  ##########################################################################
  ###############################  Retired figures #########################
  # ------------------------------------------
  # Habitat usage comparisons (bivariate) ... requires the full "set.rdata" database and "logbook.dZ.rdata" database
   # habitat.usage( usevar="totno.all", covariate="depth", outdir = file.path(p$annual.results, "habitat.templates") )

  # ------------------------------------------
  # Habitat usage comparisons (bivariate) ... requires the full "set.rdata" database and "logbook.dZ.rdata" database
    #habitat.usage( usevar="totno.all", covariate="temperature", outdir = file.path(p$annual.results, "habitat.templates") )

  # ------------------------------------------
  # Habitat usage comparisons (bivariate) ... requires the full "set.rdata" database and "logbook.dZ.rdata" database
   # habitat.usage( usevar="totno.all", covariate="bottom.slope", outdir = file.path(p$annual.results, "habitat.templates") )

  # ------------------------------------------
  # Habitat usage comparisons (bivariate) ... requires the full "set.rdata" database and "logbook.dZ.rdata" database
    #habitat.usage( usevar="totno.all", covariate="bottom.curvature", outdir = file.path(p$annual.results, "habitat.templates") )

  # ------------------------------------------
  # Habitat usage comparisons (bivariate) ... requires the full "set.rdata" database and "logbook.dZ.rdata" database
    #habitat.usage( usevar="totno.all", covariate="substrate", outdir = file.path(p$annual.results, "habitat.templates") )


  # ------------------------------------------
  # Timeseries: Larval brachyura from the SSIP data
    ##figure.timeseries.larvae( outdir=file.path(project.datadirectory("bio.snowcrab"), "R", "timeseries", "larvae") )

  # ------------------------------------------
  # Growth as a a function of instar for Scotian Shelf snow crab
    figure.growth.instar( outdir=file.path(project.datadirectory("bio.snowcrab"), "R", "growth") )


  # ------------------------------------------
  # Map: Larval distributions from the Scotian Shelf Ichtyoplankton Program data
    map.larvae( p, outdir=file.path(project.datadirectory("bio.snowcrab"), "R", "maps", "larvae"), conversions=conversions )

    # Map: Crab movement from mark-recapture data
    #MG I think Brent is primarily mapping this stuff now. Not sure the data has been updated in a while
    # map.movement( p, outdir=file.path(project.datadirectory("bio.snowcrab"), "R", "maps", "mark.recapture") )
    
    # ------------------------------------------
    # Map: Spatial representation of maturity patterns of snow crab
    #MG Not sure we use these maps either, check with Adam and Jae
    # map.maturity( p, outdir=file.path(project.datadirectory("bio.snowcrab"), "R", "maps", "maturity"), newyear=T )
    
