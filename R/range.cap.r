
  range.cap = function( x, r ) {
    x[ which( x < r[1] ) ] = r[1]
    x[ which( x > r[2] ) ] = r[2]
    return(x)
  } 


