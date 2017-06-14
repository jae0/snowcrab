
require(bio.base)
require(raster)  ## TODO :: remove raster-based routines

p = bio.snowcrab::load.environment( year.assessment=2016 )
#loadfunctions('bio.snowcrab')

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

  # sequence is important ... do not change
  # creates initial rdata and sqlite db
  snowcrab.db( DS="setInitial.redo", p=p ) # this is required by the seabird.db (but not minilog and netmind)
    # few sanity checks on the initial data pulled from the raw tables
    problems = data.quality.check( type="stations", p=p)  # duplicates
    problems = data.quality.check( type="count.stations", p=p)
    problems = data.quality.check( type="position", p=p) #MG try checking the end position of the tow, if there is an error

    # was in above. Split off as a separate function it is not essential
    # and can break without the right ESRI drivers. JC
    export.to.shapefile(
      inp=snowcrab.db( DS="setInitial", p=p ),
      out=file.path(project.datadirectory("bio.snowcrab"), "maps", "shapefiles", "survey"),
      fn="SurveyDataUpdate"
    )

  # check/choose the years
  # if just updating a single year, run the following, else all years will be run by default
  (p$year.assessment)
  p$seabird.yToload = p$year.assessment
  p$minilog.yToload = p$year.assessment
  p$netmind.yToload = p$year.assessment
  p$esonar.yToload  = p$year.assessment

  #_________________________________________________________________#
  #                                                                 #
  #  BH: These functions need documentation
  #_________________________________________________________________#
  #                                                                 #

  seabird.db( DS="load", Y=p$seabird.yToload ) # this begins 2012;
  minilog.db( DS="load", Y=p$minilog.yToload ) # minilog data series "begins" in 1999 -- 60 min?
  netmind.db( DS='esonar2netmind.conversion',Y=p$esonar.yToload )
  netmind.db( DS="load", Y=p$netmind.yToload) # netmind data series "begins" in 1998 -- 60 min?

  #JC note: 1998:2002 have about 60 files with no data, just a short header

  #MG I'm not sure why these stats are not being written automatically, neet to set it in the code above to run these after data is loaded -- JC: mostly as "stats" can fail and need to be re-run. No need to re-run "load" steps.
  p$netmensuration.problems = c( "add trouble id's here") # add troublesome id's here .. eventually move into their respective functions

  seabird.db (DS="stats.redo", Y=p$seabird.yToload )
  minilog.db (DS="stats.redo", Y=p$minilog.yToload )  # note no depth in minilog any more (since 2014 ..  useful for temperature only)
  netmind.db (DS="stats.redo", Y=p$netmind.yToload )


  # merge in netmind, minilog, seabird, esonar data and do some sanity checks
  snowcrab.db( DS="set.clean.redo", p=p )  # sanity checks
    problems = data.quality.check( type="minilog.mismatches", p=p )
    problems = data.quality.check( type="minilog.load", p=p)
    problems = data.quality.check( type="minilog.dateproblems", p=p)
    problems = data.quality.check( type="minilog", p=p)   # Check for duplicate timestamps
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
  snowcrab.db( DS="det.georeferenced.redo", p=p )
  snowcrab.db( DS="cat.initial.redo", p=p )
  snowcrab.db( DS="cat.georeferenced.redo", p=p )

  snowcrab.db( DS="set.biologicals.redo", p=p )
  snowcrab.db( DS="set.complete.redo", p=p )


  # update a database of simple transformation ranges, etc.. for plotting range, etc.
  REPOS = bio.indicators::recode.variable.initiate.db ( db="snowcrab" )

  # create simple timeseries
  snowcrab.timeseries.db( DS="observer.redo" )
  snowcrab.timeseries.db( DS="biologicals.redo" )  # approx 30 min in 2015, JC
  snowcrab.timeseries.db( DS="biologicals.2014.redo" )  # reduced subset that matches 2014 station id's ..

  # example: or to get a quick one for a few vars of interest, region of interest ... no saving to file
  snowcrab.timeseries.db( DS="biologicals.direct", regions='cfa4x', vn=c('R0.mass'), trim=0 )  # returns the data


  # update data summaries for use with biomass estimation using simple indices:
  # this can be used should the habitat based approach fail or to get a sense of the patterns
  # while you wait for it to complete .. ;)
  biomass.summary.survey.nosa.db("complete.redo", p=p) #Uses the geometric mean biomass from the survey

