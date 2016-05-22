  filter.size = function( x, size, index=T ) {
    i = which( is.finite(x) & (x >= size[1] ) & (x < size[2]) )
    if (index) return(i) else return(x[i,])
  }


