  filter.durometer = function( x, durometer ) {
    i = which( is.finite(x) & (x>= durometer[1]) & (x < durometer[2] ) )
    return(i)
  }


