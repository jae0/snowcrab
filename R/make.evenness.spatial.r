
make.evenness.spatial = function( distance=50 ) {

  require(MASS)
  require(doBy)

  #  distance = seq(10,60,10)
  #  distance = 50

  #  subsample size-at-maturity estimates in a spatio-temporal context and then map it
  # mapping is continued in functions,figures.r and "map.maturity" in 2.figures.r
  
  det = snowcrab.db("det.georeferenced")
  to.keep = c("trip", "set", "sex", "cw", "mass", "abdomen", "chela", "mat","lon", "lat", "yr", "sa","plon","plat")
  det = det[,to.keep]
  gc()

  det$X 	<- det$lon
  det$Y 	<- det$lat
  det$lon 	<- det$plon
  det$lat 	<- det$plat
  loc = c("lon", "lat")

  years = sort( unique( det$yr ) )
  years = years[years>=1999]
  nyears = length( years )
  output = NULL

  for (y in 1:nyears) {
    	j = which( det$yr == years[y] )
    	dj = det[j,]
    	dj$uid <- paste(dj$trip,dj$set,sep="")
    	sets = sort( unique(dj$uid ) )
    	nset = length(sets)
	    print (years[y])
    for (s in 1:nset) {
      	i = which( dj$uid == sets[s] )
      dji = dj[i,]
      coord0 = dji[1,loc]  # focal pt
      browser()
      d.ss = which( geodist(
        point=coord0, locations=dj[,loc], method="vincenty") <= distance
      )
      dji = dj[d.ss,]
   for (sex in c(male, female)) {
        sexi = which(dji$sex == sex )
        out = inf = NULL
        
        ##to here December 11, 2013 
        out = est.size.at.maturity (dji[sexi,], sex)

        if (!is.null(out) ) {
          out = data.frame( t(out ))
          inf = cbind( sex=sex, distance=distance, uniqueid=I(sets[s]),
                       yr=years[y], lon=coord0[1], lat=coord0[2], out
                )
          output = rbind( output , inf )
  }}}}

  good = which(
    output$converged==1 &
    output$cw50 >  20 &
    output$cw50 < 150 &
    output$deviance < 400 &
    output$aic < 400
  )

  maturity = output[good,]
  return (maturity )
}
                  

