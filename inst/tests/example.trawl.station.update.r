
  # used to identify points to modify and obtain coords for new points
   	
	loadfunctions( "snowcrab", functionname="initialise.local.environment.r") 
  
  set = snowcrab.db("set")
  set = set[ which(set$yr==2007), ]
  set = set[, c("station", "lon", "lat")]
  set = set[ is.finite( rowSums(set) ) ,]

  plot(set$lon, set$lat, pch=20)
  
  todrop = NULL
  todrop = identify(set$lon, set$lat)
  # todrop = c( 56  65  66  79 160 162 168 169 184 185 234 294 295 301)
  points( set$lon[todrop], set$lat[todrop], col="red")

  tomove = identify(set$lon, set$lat)
  for (i in tomove) {
    points(set$lon[i], set$lat[i], col="green",pch=20)
  }
   
  newlocations = NULL
  newlocations = locator()
  newlocations = as.data.frame(sapply(  newlocations, unlist))
  for (i in 1:nrow(newlocations) ) {
    points(  newlocations$x[i], newlocations$y[i], col="blue",pch=22)

  }
   

  print (newlocations) 
  
  ll = as.data.frame( array(unlist(newlocations), dim=c(35,2)) )

  l = cbind(set[tomove,], ll)
  names(l) = c("station", "lon.2005", "lat.2005", "lon.2006", "lat.2006")
  l = as.data.frame(l)

  write.table(l, file=file.path( project.datadirectory("snowcrab"), "R", "trawl.stations.toadd.txt") ,sep=";")
  write.table(set[todrop,c("station", "lon", "lat")], file=file.path( project.datadirectory("snowcrab"), "R", "trawl.stations.toremove.txt"), sep=";")
   


