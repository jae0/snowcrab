

#p = bio.snowcrab::initialise.local.environment()

 setwd( file.path( project.datadirectory("bio.snowcrab"), "R" )
 load("det.georef.rdata")
 load("set.complete.rdata")

 set = set[, c("trip", "set", "timestamp", "julian", "z", "t" )]
 det = merge( x=det, y=set, by=c("trip", "set"), all.x=T, all.y=F )
 save(det, file="det_ben.rdata", compress=T)


