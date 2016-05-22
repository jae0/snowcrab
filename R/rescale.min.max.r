 
    rescale.min.max = function( x, xmax=1 ) {
      y = x - min(x, na.rm=T)
      y = y / max(y, na.rm=T)
    }


