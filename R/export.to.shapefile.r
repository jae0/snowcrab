export.to.shapefile = function( inp, out=path.expand( file.path("~", "tmp") ), fn="export", proj="+proj=longlat +ellps=WGS84" ) {
  # moved from setInitial
  ## michelle:: please do not place hard-links into the code as this will force a fail for others ..
  ## this is probably better created as a functions and you can send to the data for a save into OGR format
  ## to a location of your choice? ..
  #MG Save the trawl file to a shapefile to bring into ArcGIS
  shape.set <- inp
  shape.set$lon <- -shape.set$lon
  shape.set$timestamp <- as.character(shape.set$timestamp)

  set.cords <- shape.set[, c("lon", "lat")]
  sdf.set <- SpatialPointsDataFrame(set.cords, data=shape.set)
  proj4string(sdf.set) <- CRS( proj )
  shpdir = out
  dir.create ( shpdir, recursive=TRUE, showWarnings=FALSE)
  setwd(shpdir)
  writeOGR(sdf.set, ".", fn, driver="ESRI Shapefile", overwrite=T)
  return ( paste(fn, "shapefile created at", shpdir, sep=" "))
}


