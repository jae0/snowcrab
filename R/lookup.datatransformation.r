
  lookup.datatransformation = function( ) {
    # determine data transformations based upon category of data and data source

    log.transform = bio.snowcrab::snowcrab.variablelist("log.transform")
    scaled.centered = bio.snowcrab::snowcrab.variablelist("scaled.centered")
    sn = bio.snowcrab::snowcrab.variablelist("all.data")
    set =  bio.snowcrab::snowcrab.db(DS="set.complete") # base transform characteristics
    logs =  bio.snowcrab::logbook.db(DS='logbook') 
    repository = file.path( project.datadirectory("bio.snowcrab"), "output", "transform.lookup.rdata" )
    return( list( log.transform=log.transform, sn=sn, set=set,logs = logs, repository=repository ) )

  }


