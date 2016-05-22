
  map.fisheries.data = function(p, gridfun, variables) {
    x = logbook.db( DS="logbook" )
    #x$landings = x$landings/1000  # convert kg to ton #MG code this out so that it calculates quantiles properly
    x$sa = 1  # this a dummy variable required by the mapping routine
    x = x [filter.region.polygon( x, region="isobath1000m"),]
    x = x[ which(x$effort <= 300) ,]
    x = x[ which(x$cpue < 500),]

    polydir = file.path(project.datadirectory("polygons"), "data", "Basemaps", "Marine", "Coastline")
    shpdir = file.path(project.datadirectory("snowcrab"), "maps", "shapefiles", "logbook")
    rasdir = file.path(project.datadirectory("snowcrab"), "maps", "rasters", "logbook")
    mapdir = file.path(project.datadirectory("snowcrab"), "maps", "images", "logbook")

    variables = c("effort", "landings", "cpue")
    grid.fun=mean
    raster.map.variables(x, p, variables, p$plottimes, polydir, shpdir, rasdir, mapdir, grid.fun=mean)
    return("Mapping Completed")
  }


