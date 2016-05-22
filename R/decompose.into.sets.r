
  decompose.into.sets = function( x, data.series="seabird" ) {
         
    # determine an empirical threshold for fishing vs not fishing 
    # this algorithm uses empircal CDF's and quantiles: most of the time is spent near the surface not fishing
    # set a threshold of the 0.05 quantile as the cutoff.. and make sure there is no drift either too
  

    ndata = length(x)
    X = data.frame( dataorder=1:ndata, X=x, fishing=1, id=NA )

    # based upon empirical observation from 2012 data, near-surface values should be ~ 0.014 decibar
    # detect any drift in the initial pressure values
    if (data.series=="seabird") {
      # error checking values to detect drift in pressure sensors at surface
      expectedValue=0.014 
      PrThreshold = 0.05  # median value ~ 0.055 decibars
      fishing.threshold = 10 # 5 decibars will remove most non-fishing events
      lscale =  60  # ==60 / 5 *5, the number of pings expected for a 5 min tow at 5 pings per second == 60/5 *5)
      # for the seabird there is a ping every 5 sec ; 2min minimum = 2*60/5 = 24 pings required 
    }

    threshold = quantile( X$X, probs=PrThreshold ) 
    if (abs (threshold - expectedValue ) > (2 * sd( X$X[ X$X <= threshold] )) ) {
      print( "----------")
      print( "A significant deviation was detected in the data stream for the data file, above" )
      print( "Assuming it was drift, but the data should be checked." )
      print( paste( "observed:", threshold) ) 
      print( "expected: ~ 0.1 decibar" )
      print( paste( "mean: ", mean( x[x<=threshold] ) ) )
      print( paste( "sd: ", sd( x[x<=threshold] ) ) )
      print( "----------")
      
      if (abs(threshold) > 1) {
        print( "In fact, the deviation was > 1 decibar" )
        print( "Stopping now as there may be something wrong with the data" )
        stop()
      }
      
      X$X = X$X - threshold 
    
    }

    
    # flag to indicate is fishing is occurring via bottom contact
    not.fishing = which( X$X < fishing.threshold)   # 5 decibars is a safe limit given the survey
    if (length(not.fishing) > 0) X$fishing[ not.fishing ] = 0
  
    # by taking the sum of fishing with fishing offset by 1, we can 
    # get the boundaries where fishing starts and stop by looking for 1's
    X$fishing.boundary = X$fishing + c(0, X$fishing[1:(ndata-1)] ) 
    intervals = which( X$fishing.boundary == 1 ) # these are the boundaries between fishing and non-fishing events

  
    count = 0
    for (o in 1:(length((intervals))-1)) { 
      i0 = intervals[o]
      i1 = intervals[o+1]
      if ( (i1 - i0) < lscale*0.75 ) next()  # accept only if the number of pings > 75% of expected number to be safe 
      # number of pings is sufficient ... now check if fishing or not
      if (mean( X$X[i0:i1], na.rm=T) < 2*fishing.threshold ) next()
      count = count + 1
      X$id [ i0:i1 ] = count
    }
    return( X$id[ order(X$dataorder) ])
  }


