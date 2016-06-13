# ----------------------------------------------------------------------------------
# NOTE :
#
#  1. The year in the file "year.assessment.r" must be changed every year.
#     It must be kept separate from "initialise.local.environment.r" as running in parallel mode
#     requires occasionally overrding some parameters in "p". This override cannot be completed as
#     "initialise.local.environment.r" is sourced with every initialisation of a a new CPU.
  # ----------------------------------------------------------------------------------

  # load required functions and parameters

  p = bio.snowcrab::initialise.local.environment( year.assessment=2016 )

  debug = FALSE
  if (debug) {
    p$do.parallel =F
    p$clusters= c("localhost")
  }

    # get data tables from Oracle server and store local copies
  # !!!!!! --------- these should be run on a windows machine: !!!!!!!!! <--------- READ THIS
  if (obtain.database.snapshot) {
    snowcrab.db( DS="set.odbc.redo", yrs=1996:p$year.assessment ) # Copy over datadirectory ("bio.snowcrab"), "data", "trawl", "SNCRABSETS"
    snowcrab.db( DS="det.odbc.redo", yrs=1996:p$year.assessment ) # Copy over datadirectory ("bio.snowcrab"), "data", "trawl", "SNCRABDETAILS"
    snowcrab.db( DS="cat.odbc.redo", yrs=1996:p$year.assessment ) # Copy over datadirectory ("bio.snowcrab"), "data", "trawl", "SNTRAWLBYCATCH"
    logbook.db(  DS="odbc.logbook.redo", yrs=1996:p$year.assessment ) #Copy over datadirectory ("bio.snowcrab"), "data", "logbook", "datadump"
    logbook.db(  DS="odbc.licence.redo" ) #Copy over datadirectory ("bio.snowcrab"), "data", "logbook", "lic.datadump.rdata"
    logbook.db(  DS="odbc.areas.redo" ) #Copy over datadirectory ("bio.snowcrab"), "data", "observer", "datadump"
    observer.db( DS="odbc.redo", yrs=1996:p$year.assessment )
  }

# -------------------------------------------------------------------------------------
# produce base data files from bio.snowcrab logbook database (marfis) and historical data
# and the observer databases:
# this needs to be done after the above datadumps to refresh locally stored databases

  if (make.fisheries.data) {
    observer.db( DS="odb.redo", p=p ) # 3 minutes
    # fishing ground are used for determination of contraints for interpolation
    logbook.db( DS="logbook.redo", p=p )
    logbook.db( DS="logbook.filtered.positions.redo", p=p )
    logbook.db( DS="fishing.grounds.redo",  p=p )
    logbook.db( DS="logbook.gridded.redo", p=p )
  }


# -------------------------------------------------------------------------------------
# create base set data and add all historical data fixes

  if (get.base.data) {
    # sequence is important ... do not change
    # creates initial rdata and sqlite db
    snowcrab.db( DS="setInitial.redo", p=p ) # this is required by the seabird.db (but not minilog and netmind)
      # few sanity checks on the initial data pulled from the raw tables
      problems = data.quality.check( type="stations")  # duplicates
      problems = data.quality.check( type="count.stations")
      problems = data.quality.check( type="position") #MG try checking the end position of the tow, if there is an error

      # was in above. Split off as a separate function it is not essential
      # and can break without the right ESRI drivers. JC
      export.to.shapefile(
        inp=snowcrab.db( DS="setInitial", p=p ),
        out=file.path(project.datadirectory("bio.snowcrab"), "maps", "shapefiles", "survey"),
        fn="SurveyDataUpdate"
      )

    # check/choose the years
    (p$year.assessment)

    # if just updating a single year, run the following, else all years will be run by default
    p$seabird.yToload = p$year.assessment
    p$minilog.yToload = p$year.assessment
    p$netmind.yToload = p$year.assessment
    p$esonar.yToload  = p$year.assessment


    seabird.db( DS="load", Y=p$seabird.yToload ) # this begins 2012;
    minilog.db( DS="load", Y=p$minilog.yToload ) # minilog data series "begins" in 1999 -- 60 min?

    netmind.db( DS='esonar2netmind.conversion',Y=p$esonar.yToload )
    netmind.db( DS="load", Y=p$netmind.yToload) # netmind data series "begins" in 1998 -- 60 min?
      #JC note: 1998:2002 have about 60 files with no data, just a short header

    #MG I'm not sure why these stats are not being written automatically, neet to set it in the code above to run these after data is loaded
    seabird.db (DS="stats.redo", Y=p$seabird.yToload )
    minilog.db (DS="stats.redo", Y=p$minilog.yToload )
    netmind.db (DS="stats.redo", Y=p$netmind.yToload )


    # merge in netmind, minilog, seabird, esonar data and do some sanity checks
    snowcrab.db( DS="set.clean.redo", p=p )  # sanity checks
      problems = data.quality.check( type="minilog.mismatches", p=p )
      problems = data.quality.check( type="minilog.load", p=p)
      problems = data.quality.check( type="minilog.dateproblems", p=p)
      problems = data.quality.check( type="minilog", p=p) # Check for duplicate timestamps
      problems = data.quality.check( type="netmind.load", p=p)
      problems = data.quality.check( type="netmind.mismatches", p=p )
      problems = data.quality.check( type="tow.duration", p=p)
      problems = data.quality.check( type="tow.distance", p=p)
      problems = data.quality.check( type="seabird.mismatches", p=p )
      problems = data.quality.check( type="seabird.load", p=p)
      problems = data.quality.check( type="netmind.timestamp" , p=p)


    #MG det.initial.redo updates and processes morphology. This code now identifies morphology errors, which must be
    #checked with written logs, then sent to database and put in debugging here and re-run
    snowcrab.db( DS="det.initial.redo", p=p )
    snowcrab.db( DS="det.georeferenced.redo" )

    snowcrab.db( DS="cat.initial.redo", p=p )
    snowcrab.db( DS="cat.georeferenced.redo" )

    snowcrab.db( DS="set.merge.det.redo" )
    snowcrab.db( DS="set.merge.cat.redo" )

  }  # end base data



# -------------------------------------------------------------------------------------
# External Dependencies: (must be completed before the final lookup/mathcing phase)

  p=initialise.local.environment()

# Bathymetry data ::
  loadfunctions("bathymetry", functionname="bathymetry.r" ) # if necessary
# Substrate type  ::
  loadfunctions("substrate", functionname="substrate.r" ) # if necessary
# Groundfish data ::
# NOTE  groundfish.db( DS="odbc.redo" ) must first be done manually
# on a windows machine and data snapshots moved to local system
  loadfunctions( "groundfish", functionname="1.groundfish.r" ) #MG took 15-25 minutes to run
# Taxonomy ::
  loadfunctions("taxonomy", functionname="taxonomy.r" ) # if necessary #MG takes about 30 minutes to run
## The following are very SLOW:
# Temperatures ::
  loadfunctions ( "temperature", functionname="temperature.r" )  # days to run
# Habitat data ... environmentals only as it is used by bio.db etc
  loadfunctions ( "habitat", functionname="habitat.r") #MG fairly quick to run
  #loadfunctions ( "habitat", functionname="habitat.temperatures.r" )
# BIO db update ::
# must come after temperature interpolations to permit temperature lookups
  loadfunctions ( "bio", functionname="bio.r" )  #MG took about 20 minutes to run
# the folllowing depends upon bio.db and temperature
  #MG species area took 2 days to run in parallel, run some things on server if possible. It's quicker to run some things in serial though, ask Jae
  loadfunctions ( "speciesarea", functionname="speciesarea.r" )
  loadfunctions ( "speciescomposition", functionname="speciescomposition.r" )
  loadfunctions ( "sizespectrum", functionname="sizespectrum.r" )
  loadfunctions ( "metabolism", functionname="metabolism.r" )
  loadfunctions ( "condition", functionname="condition.r" ) #MG this one took 8 days to run on the laptop, in serial...
# Habitat data :: NOTE:: This glues all the above together in planar coord system
# to allow fast lookup of data for matching with set, logbook data
  loadfunctions ( "habitat", functionname="habitat.complete.r" )



# -------------------------------------------------------------------------------------
# Final data lookup/matching .. AFTER refreshing all above tables (where relevent/possible)

  p=initialise.local.environment()

  logbook.db( DS  ="fisheries.complete.redo", p=p )
  snowcrab.db( DS ="set.complete.redo", p=p )
  snowcrab.db( DS ="set.logbook.redo", yrs=1996:p$year.assessment ) # add gridded fisheries data
  # snowcrab.db( DS ="set.logbook", yrs=1996:p$year.assessment )

  #make.timeseries.data(p=p, areas=p$regions )  #  timeseries of means of all survey data
  #in 2014 as there were reduced stations for comparison
  #make.timeseries.data(p=p, areas=p$regions,reduced.stations=F, vars='R0.mass' ) #  timeseries of means of all survey data
  make.timeseries.data(p=p, areas=NULL,reduced.stations=F, vars=NULL) #  timeseries of means of all survey data
  #make.timeseries.data(p=p, areas=NULL,reduced.stations=F, vars=c('ms.mass.10', 'ms.mass.30', 'ms.mass.201', 'ms.mass.50', 'ms.mass.2521', 'ms.mass.2511', 'ms.mass.202', 'ms.mass.204', 'ms.mass.2211'), outfile = file.path( project.datadirectory("bio.snowcrab"), "R", "tsbycatch.rdata" )) #  timeseries of means of all survey data

  #  tsdata = snowcrab.db("set.timerseries")

# create a new lookuptable for data transformations after refreshing set data/ranges
  REPOS = recode.variable.initiate.db ( db="snowcrab" )


# -------------------------------------------------------------------------------------
# snow crab found in external databases tapped into for habitat determination
  #for ( vs in c( "R0.mass", "male.large", "male.small", "female.large", "female.small" ) ) {
    ### -------- not yet finished this one ...  TODO
    vs="R0.mass"
    snowcrab.external.db(p=p, DS="set.snowcrab.in.groundfish.survey.redo", vname=vs, year.current=p$year.assessment )

    # ---- TODO !!! must replace this with bio.db processing step


# simple geometric means of raw data:  used by indicators ordination and some figures
  # takes many hours ... need to make parallel  TODO
  tsdata =  get.time.series ( x=snowcrab.db( DS="set.logbook"),
  regions=p$regions, vars=variable.list.expand("all.data"), from.file=F, trim=0 )



#  ----- experimental / tests / ideas
  s =  snowcrab.db( DS ="set.complete" )
  d =   snowcrab.db( DS ="det.georeferenced" )
  l = merge( d, s[, c("trip", "set", "t")], by=c("trip", "set"), all.x=T, all.y=F)
  rm(s,d); gc()
  l = l[ which( as.numeric(as.character(l$mat)) %in% c(mature, immature)  &
                l$sex %in% c(male, female) ) , ]

  l$sex= factor( as.character(l$sex) )
  l$mat = factor( as.character( l$mat))

  m = glm( t~ as.factor(mat) * as.factor(sex), data= l, family=gaussian())
  require(car)

  Anova( m)
  require(effects)
  k=all.effects(m, l)
  plot(k)


  # -------------------------------------------------------------------------------------
  # test: make size at maturity estimates in a spatial context

  if( make.maturity.db ) {
    maturity = make.maturity.spatial( distance=50 )
    save(maturity, file="maturity.rdata", compress=T)
    # load(file.path( project.datadirectory("bio.snowcrab"), "R", "maturity.rdata"))
  }


