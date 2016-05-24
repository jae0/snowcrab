set = snowcrab.db( DS ="set.complete", p=p )   

set = set [ , c("trip", "set", "station", "yr", "lon", "lat", "t", "R0.mass", "R0.no" ) ]

shape.set <- set
#shape.set$lon <- -shape.set$lon
shape.set$chron <- as.character(shape.set$chron)

set.cords <- shape.set[, c("lon", "lat")]
sdf.set <- SpatialPointsDataFrame(set.cords, data=shape.set)
proj4string(sdf.set) <- CRS(p$geog.proj)
shpdir = file.path(project.datadirectory("snowcrab"), "maps", "shapefiles", "survey")
setwd(shpdir)

writeOGR(sdf.set, ".", "SurveyDataR0.mass", driver="ESRI Shapefile", overwrite=T)
setwd("/home/michelle/tmp")  
shp.path <- paste("SurveyDataUpdate shapefile created at", shpdir, sep=" ")
print(shp.path)
