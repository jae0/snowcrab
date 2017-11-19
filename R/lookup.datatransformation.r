
  lookup.datatransformation = function( ) {
    # determine data transformations based upon category of data and data source

    log.transform = bio.snowcrab::variable.list.expand("log.transform")
    scaled.centered = bio.snowcrab::variable.list.expand("scaled.centered")
    sn = bio.snowcrab::variable.list.expand("all.data")
    set = bio.snowcrab::snowcrab.db(DS="set.merge.cat") # base transform characteristics
    logs = bio.snowcrab::logbook.db(DS='logbook')
    repository = file.path( project.datadirectory("bio.snowcrab"), "output", "transform.lookup.rdata" )
    return( list( log.transform=log.transform, sn=sn, set=set,logs = logs, repository=repository ) )

  }


