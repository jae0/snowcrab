
require(aegis)

if (!exists("year.assessment")) {
  year.assessment=lubridate::year(Sys.Date())      # year.assessment
  year.assessment=lubridate::year(Sys.Date()) -1   # or year previous to current
}
p = bio.snowcrab::load.environment( year.assessment=year.assessment )
#loadfunctions('bio.snowcrab')

# get data tables from Oracle server and store local copies
#BZ 2017- Can now be run on a Linux machine. rawdata connection changed to ROracle

if (obtain.database.snapshot) {
  # yrs = 1996:p$year.assessment  # to redo all years
  yrs = p$year.assessment # to update only the current year
  snowcrab.db( DS="set.rawdata.redo", yrs=yrs ) #  datadirectory ("bio.snowcrab"), "data", "trawl", "SNCRABSETS"
  snowcrab.db( DS="det.rawdata.redo", yrs=yrs ) #  datadirectory ("bio.snowcrab"), "data", "trawl", "SNCRABDETAILS"
  snowcrab.db( DS="cat.rawdata.redo", yrs=yrs ) #  datadirectory ("bio.snowcrab"), "data", "trawl", "SNTRAWLBYCATCH"
  logbook.db(  DS="rawdata.logbook.redo", yrs=yrs ) #  datadirectory ("bio.snowcrab"), "data", "logbook", "datadump"
  logbook.db(  DS="rawdata.licence.redo" ) # datadirectory ("bio.snowcrab"), "data", "logbook", "lic.datadump.rdata"
  logbook.db(  DS="rawdata.areas.redo" ) # datadirectory ("bio.snowcrab"), "data", "observer", "datadump"
  observer.db( DS="rawdata.redo", yrs=yrs )
}

# -------------------------------------------------------------------------------------
# produce base data files from bio.snowcrab logbook database (marfis) and historical data
# and the observer databases:
# this needs to be done after the above datadumps to refresh locally stored databases

  if (make.fisheries.data) {
    observer.db( DS="odb.redo", p=p ) # 3 minutes
    logbook.db( DS="logbook.redo", p=p )
    logbook.db( DS="logbook.filtered.positions.redo", p=p )
    # fishing ground are used for determination of contraints for interpolation
    logbook.db( DS="fishing.grounds.redo",  p=p )
    logbook.db( DS="logbook.gridded.redo", p=p )
  }


# -------------------------------------------------------------------------------------
# create base set data and add all historical data fixes

  # sequence is important ... do not change
  # creates initial rdata and sqlite db
  snowcrab.db( DS="setInitial.redo", p=p ) # this is required by the seabird.db (but not minilog and netmind)
  snowcrab.db( DS="set.clean.redo", p=p )

    # few sanity checks on the initial data pulled from the raw tables
    problems = data.quality.check( type="stations", p=p)  # duplicates
    problems = data.quality.check( type="count.stations", p=p)
    problems = data.quality.check( type="position", p=p) #MG try checking the end position of the tow, if there is an error
    #BZ ToDo 2019- getting an error with line 52 but runing the function directly works

    # was in above. Split off as a separate function it is not essential
    # and can break without the right ESRI drivers. JC
              #export.to.shapefile(
                #inp=snowcrab.db( DS="setInitial", p=p ),
                #out=file.path(project.datadirectory("bio.snowcrab"), "maps", "shapefiles", "survey"),
               # fn="SurveyDataUpdate"
              #)


# -------------------------------------------------------------------------------------
# process the net configuration and the temperatures from seabird, netmind, etc ..
# if just updating a single year, run the following, else all years will be run by default
  if (updating.year.assessment) {
    p$seabird.yToload = p$year.assessment
    p$minilog.yToload = p$year.assessment
    p$netmind.yToload = p$year.assessment
    p$esonar.yToload  = p$year.assessment
  }

  seabird.db( DS="load", Y=p$seabird.yToload ) # this begins 2012;duplicates are often due to seabird files not being restarted each morning

  minilog.db( DS="load", Y=p$minilog.yToload ) # minilog data series "begins" in 1999 -- 60 min?

  #netmind.db( DS='esonar2netmind.conversion',Y=p$esonar.yToload ) #may no longer need to be run BZ
  netmind.db( DS="load", Y=p$netmind.yToload) # netmind data series "begins" in 1998 -- 60 min?
  #JC note: 1998:2002 have about 60 files with no data, just a short header

  p$netmensuration.problems = c( "add trouble id's here") # add troublesome id's here .. eventually move into their respective functions
  seabird.db (DS="stats.redo", Y=p$seabird.yToload )
  minilog.db (DS="stats.redo", Y=p$minilog.yToload )  # note no depth in minilog any more (since 2014 ..  useful for temperature only)
  netmind.db (DS="stats.redo", Y=p$netmind.yToload )



# -------------------------------------------------------------------------------------
# merge in netmind, minilog, seabird, esonar data and do some sanity checks
# Can add any datachecks that might improve overall data quality
# BZ for 2017- If errors repeat and are not actually a problem, create an override

    #problems = data.quality.check( type="minilog.mismatches", p=p )
    problems = data.quality.check( type="minilog.load", p=p)
    problems = data.quality.check( type="minilog.dateproblems", p=p) #track down why ~all sets are giving mismatches
    problems = data.quality.check( type="minilog", p=p)   # Check for duplicate timestamps
    problems = data.quality.check( type="netmind.load", p=p)
    problems = data.quality.check( type="netmind.mismatches", p=p )
    problems = data.quality.check( type="tow.duration", p=p)
    problems = data.quality.check( type="tow.distance", p=p)
    problems = data.quality.check( type="seabird.mismatches", p=p )
    problems = data.quality.check( type="seabird.load", p=p)
    problems = data.quality.check( type="netmind.timestamp" , p=p)


# -------------------------------------------------------------------------------------
# local lookup trables are required for the following snowcrab.db steps

library(aegis.bathymetry)
    p = bio.snowcrab::snowcrab_carstm( DS="parameters", assessment.years=1999:year.assessment ) 
    pB = bathymetry_carstm( p=p, DS="parameters_override" )
    M = bathymetry.db( p=pB, DS="aggregated_data", redo=TRUE ) 
    
library(aegis.substrate)
    pS = substrate_carstm(p=p, DS="parameters_override" )
    M = substrate.db( p=pS, DS="aggregated_data", redo=TRUE ) 
    
library(aegis.temperature)
    pT = temperature_carstm(p=p, DS="parameters_override" )
    M = temperature.db( p=pT, DS="aggregated_data", redo=TRUE )  

#need to create a symlink through windows to the bathymetry file create through aegis.bathymetry
#run the following line in windows command prompt
    
#mklink "C:/bio.data/bio.snowcrab/data/bathymetry.canada.east.lonlat.rawdata.rdata" "C:/bio.data/aegis/bathymetry/data/bathymetry.canada.east.lonlat.rawdata.rdata"

#reinitialize parameters specific to the snow crab project    
p = bio.snowcrab::load.environment( year.assessment=year.assessment ) #reestablish parameters for snowcrab processes
        
# -------------------------------------------------------------------------------------
# Finalize the data sets
# MG det.initial.redo updates and processes morphology. This code now identifies morphology errors, which must be
# checked with written logs, then sent to database and put in debugging here and re-run

    
  snowcrab.db( DS="det.initial.redo", p=p )
  snowcrab.db( DS="det.georeferenced.redo", p=p )

  snowcrab.db( DS="cat.initial.redo", p=p )
  snowcrab.db( DS="cat.georeferenced.redo", p=p )

  snowcrab.db( DS="set.biologicals.redo", p=p )
  snowcrab.db( DS="set.complete.redo", p=p )


# -------------------------------------------------------------------------------------
# update a database of simple transformation ranges, etc.. for plotting range, etc.
  snowcrab.db( DS="data.transforms.redo", p=p)

# -------------------------------------------------------------------------------------
# create some simple/crude timeseries by each CFA
  snowcrab.timeseries.db( DS="observer.redo", p=p )
  snowcrab.timeseries.db( DS="biologicals.redo", p=p )
  snowcrab.timeseries.db(DS="groundfish.t.redo", p=p )
  # snowcrab.timeseries.db( DS="biologicals.2014.redo" )  # reduced subset that matches 2014 station id's ..

  # example: or to get a quick one for a few vars of interest, region of interest ... no saving to file
  snowcrab.timeseries.db( DS="biologicals.direct", p=p, regions='cfa4x', vn=c('R0.mass'), trim=0 )  # returns the data
