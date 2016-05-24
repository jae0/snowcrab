	
	loadfunctions( "snowcrab", functionname="initialise.local.environment.r") 
   
  w = snowcrab.db( DS ="set.complete", p=p )   
  w = w[ filter.region.polygon( w, "cfanorth"),]

  out = data.frame( cbind( 
    no.male = tapply( w$totno.male, w$yr, sum, na.rm=T) ,  
    no.female = tapply( w$totno.female, w$yr, sum, na.rm=T)  
  ) )

  out$year = rownames( out )

  write.csv( out, file="cfanorth.counts.csv")



