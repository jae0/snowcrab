  

  jags.extract = function( jo, var, func, cnames=c( "cfanorth", "cfasouth", "cfa4x" ) ) {
    o =  summary( jo[[var]], func)[["stat"]] 
    if ( is.vector(o) ) o = t(o)
    colnames( o) = cnames
    return(o)
  }
 

