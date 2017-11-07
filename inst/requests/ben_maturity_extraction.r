

#p = bio.snowcrab::load.environment()

 setwd( file.path( project.datadirectory("bio.snowcrab"), "output" )
 load("det.georef.rdata")
 load("set.complete.rdata")

 set = set[, c("trip", "set", "timestamp", "julian", "z", "t" )]
 det = merge( x=det, y=set, by=c("trip", "set"), all.x=T, all.y=F )
 save(det, file="det_ben.rdata", compress=T)


