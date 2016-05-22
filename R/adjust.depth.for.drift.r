
  adjust.depth.for.drift = function( x, data.series="seabird" ) {
         
    # most of the time is spent near the surface not fishing

    # based upon empirical observation from 2012 data, near-surface values should be ~ 0.014 decibar
    # detect any drift in the initial pressure values

    ndata = length(x)
    X = data.frame( dataorder=1:ndata, X=x, fishing=1, id=NA )

    if (data.series=="seabird") {
      # error checking values to detect drift in pressure sensors at surface
      expectedValue=0.014 
      PrThreshold = 0.05  
    }

    threshold = quantile( X$X, probs=PrThreshold ) 
    if (abs (threshold - expectedValue ) > (2 * sd( X$X[ X$X <= threshold] )) ) {
      print( paste( "Adjusting depth sensors for drift: observed:", threshold, "expected: ~ 0.1 decibar" ) )
      print( paste( "mean: ", mean( x[x<=threshold] ), "sd: ", sd( x[x<=threshold] ) ) )
      print( "----------")
      
      if (abs(threshold) > 1) {
        print( "In fact, the deviation was > 1 decibar" )
        print( "Stopping now as there may be something wrong with the data" )
        
        stop()
      }
      
      X$X = X$X - threshold 
    
    }

    return(X$X)

}
