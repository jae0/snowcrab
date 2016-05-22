

  filter.prediction.locations = function( DS="", PS=NULL, y=NULL, p=NULL  ) {
    
    ips = NULL
    require (fields)

    if( DS=="default" ) {
      # choose good habitat and 
      # remove poor prediction locations
      # remove the western most parts of 4X 
      ips = which(  
        ( PS$habitat > 0.5 ) & 
        ( PS$habitat.se < PS$habitat ) & 
        ( PS$plon > 295) 
      )  
      return (ips)
    }
    
    if ( DS=="limit.to.near.survey.stations" ) {
      ## add points within X km of a survey station
      S = snowcrab.db( DS="set.logbook") 
      if ( ! is.null(S) ) { 
				S = S[ , c("plon", "plat") ]
				distances =  rdist( PS[,c("plon", "plat")], S)
				distances[ which(distances < p$threshold.distance) ] = NA
				ips = which( !is.finite( rowSums(distances) ) ) 
			}
			return (ips)
     }


    return ( "Wrong data source chosen" )
  }



