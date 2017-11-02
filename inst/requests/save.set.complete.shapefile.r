
p = bio.snowcrab::load.environment()

set = snowcrab.db( DS ="set.biologicals", p=p )

set = set [ , c("trip", "set", "station", "yr", "lon", "lat", "t", "R0.mass", "R0.no" ) ]

shape.set <- set
#shape.set$lon <- -shape.set$lon
shape.set$timestamp <- as.character(shape.set$timestamp)

set.cords <- shape.set[, c("lon", "lat")]
sdf.set <- SpatialPointsDataFrame(set.cords, data=shape.set)
proj4string(sdf.set) <- CRS("+proj=longlat +ellps=WGS84")

shpdir = file.path(project.datadirectory("bio.snowcrab"), "maps", "shapefiles", "survey")
setwd(shpdir)

writeOGR(sdf.set, ".", "SurveyDataR0.mass", driver="ESRI Shapefile", overwrite=T)
shp.path <- paste("SurveyDataUpdate shapefile created at", shpdir, sep=" ")
print(shp.path)
