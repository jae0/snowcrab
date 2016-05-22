
  data.quality.check = function( set, type ) {
    
 
     
    if (type=="seabird.load") {
        SS = snowcrab.db( DS="setInitial")
        sb = seabird.db( DS="set.seabird.lookuptable" )
        SS = merge( SS, sb, by=c("trip", "set"), all.x=TRUE, all.y=FALSE, sort=TRUE ) 
        isb = which( is.na( SS$seabird_uid) & SS$yr %in% seabird.yToload & SS$yr >= 2012 )
        print( "Missing seabird matches: ")
        print( SS[ isb,] )
        return( SS[isb,])
    }
  
 
    if (type=="minilog.load") {
        SS = snowcrab.db( DS="setInitial")
        ml = minilog.db( DS="set.minilog.lookuptable" )
        SS = merge( SS, ml, by=c("trip", "set"), all.x=TRUE, all.y=FALSE, sort=TRUE ) 
        iml = which( is.na( SS$minilog_uid) & SS$yr %in% minilog.yToload & SS$yr >= 2004 )
        print( "Missing minilog matches: ")
        print( SS[ iml,] )
        return( SS[iml,])
    }
  
    if (type=="netmind.load") {
        SS = snowcrab.db( DS="setInitial")
        nm = netmind.db( DS="set.netmind.lookuptable" )
        SS = merge( SS, nm, by=c("trip", "set"), all.x=TRUE, all.y=FALSE, sort=TRUE ) 
        inm = which( is.na( SS$netmind_uid) & SS$yr %in% netmind.yToload & SS$yr >= 2004 )
        print( "Missing netmind matches: ")
        print( SS[ inm,] )
        return( SS[inm,])
    }




    if (type=="tow.duration") {
      e0 = which( ( set$dt > 9  | set$dt < 3.5 )  & set$yr >=2004 )
      if  (length(e0)>0 ) {
        print( "The following have rather short/long tow times (dt)" )
        print( set[e0, c("trip", "set", "station", "dt", "chron")] )
        return (set[e0,] ) 
      }
    }


    if (type=="tow.distance") {
      # expected = 2 knots * 5 min = 2 * 1.852 * 5/60 = 0.309 km ( so a good range is {-25%, +75%} = (0.232, 0.5408)
      e0 = which( ( set$distance > 0.541  | set$distance < 0.232 )  & set$yr >=2004 )
      if  (length(e0)>0 ) {
        print( "The following have rather short/long tow distances" )
        print( set[e0, c("trip", "set", "station", "distance", "chron")] )
        return (set[e0,] ) 
      }
    }

  
    if (type=="netmind.timestamp") {
      #  check times/data and merge remaining data using datestamps and {station, set} 
      time.diff = set$netmind.timestamp - set$chron
      time.thresh = 30/60/24  # in days
      i = which( abs( as.numeric( time.diff)) > time.thresh )  
      if (length(i)>0) {
        print ("Potential date/time mismatches::")
        print( set[i, ] )
        return(set[i,])
      } 
    }
  
    # netmind mismatches
    if(type=="netmind.mismatches") {
      q = which( set$yr > 2004 & (set$netmind_uid==""| is.na(set$netmind_uid) )  )
      if ( length (q) > 0 ) {
        print( "No netmind matches for the following sets:")
        print ( set[q,c("trip", "set", "station", "t0", "chron") ] )
        return( set[q,c("trip", "set", "station", "t0", "chron") ] )
      }
    }


    # duplicated stations
    if (type=="stations") {
      must.be.unique = paste( set$yr, set$station, sep="~" )
      dups = duplicates.toremove( must.be.unique )
      x = set[ dups, c( "trip", "set", "station" ) ]
      x = x[ is.finite(x$station), ]
      dup.stations = x[ order(x$station) ,] 
      print( "Duplicated stations:" )
      print( dup.stations  )
      return (dup.stations)
    }    
       
    # poor minilog matches
    if(type=="minilog") {
      must.be.unique = set$t0
      dups = duplicates.toremove( must.be.unique )
      x = set[ dups, c( "trip", "set", "station", "t0" ) ]
      dup.t0 = x[ is.finite(x$t0), ]
      print( "Duplicated minilog times:" )
      print( dup.t0 )
      return (dup.t0)
    }
    
    # minilog mismatches
    if(type=="minilog.mismatches") {
      q = which( set$yr > 2004 & (set$minilog_uid==""| is.na(set$minilog_uid) )  )
      if ( length (q) > 0 ) {
        print( "No minilog matches for the following sets:")
        print ( set[q,c("trip", "set", "station", "t0", "chron") ] )
        return( set[q,c("trip", "set", "station", "t0", "chron") ] )
      }
    }

    if(type=="minilog.dateproblems") {
      time.thresh = 1/24  # in days
      ii = which( ( abs( as.numeric( set$t0 - set$chron)) > time.thresh ) )   
      if ( length (ii) > 0 ) {
        print( "Minilog date mismatches with Trip ID, using Trip id as default pre 2005:" )
        print ( set[ii, c("trip", "set", "station", "t0", "chron", "minilog_uid")] )
        return ( set[ii, c("trip", "set", "station", "t0", "chron", "minilog_uid")] )
      }
    }

    # positional information
    if(type=="position") {
      plot(set$lon, set$lat) 
      inside = filter.region.polygon( set[, c("lon", "lat") ], "cfaall")
      if (length (inside) == nrow(set) ) {
        print("All data are within positional bounds")
        return (NULL)
      } else {  
        outside = setdiff( 1:nrow(set), inside )
        points( set$lon[outside], set$lat[outside], col="red")
        print( "------------- The following are out of the cfa bounds: "  )
        print( set[ outside, ] )
        return(  set[ outside, ]  )
      }
    }

    # counts of stations by area
    if(type=="count.stations") {
      years = sort( unique( set$yr ) )
      nyears = length(years)
      nregions = length(p$regions)
      res = matrix( NA, nrow=nyears, ncol=nregions)
      for (r in 1:nregions) {
        nr = filter.region.polygon(x=set, region=recode.areas(p$regions[r]), planar=F)
        for (y in 1:nyears) {
          ni = which( set$yr==years[y] )  
          res[y,r] = length( unique( intersect (nr, ni) ) )
      }}
      x = as.data.frame(res)
      names(x) = c(p$regions)
      x$yr = years
      x = x[ , c("yr", p$regions)]
      print( "Number of stations: ")
      print (x)
      return(x)
    }
  if(type=='na.spread') {
    	ii <- which(is.na(set$spread))
    	x = set[ii,c('trip','set','station','netmind')]
    	print(x)
    	return(x)   
    }
    if(type=='na.distance') {
    	ii <- which(is.na(set$distance))
    	x = set[ii,c('trip','set','station','netmind')]
    	print(x)
    	return(x)   
    }
  }
  
