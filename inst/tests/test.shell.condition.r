  
  require(lattice)
 
	loadfunctions( "snowcrab", functionname="initialise.local.environment.r") 
 
  k = snowcrab.db( DS ="det.georeferenced")
  
  i12 = which(k$shell %in% c(1,2) )
  i34 = which(k$shell %in% c(3,4) )
  i5 = which(k$shell %in% c(5) )
  
  k$sc = NA
  k$sc[i12] = "sc12"
  k$sc[i34] = "sc34"
  k$sc[i5] = "sc5"

  xyplot( durometer ~ cw | sex, data=k )

  histogram( ~ durometer | sc+sex, data=k )
  densityplot( ~ durometer | sc+sex , data=k )

