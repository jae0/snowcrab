
  filter.shell = function( x, shell ) {
    i = which( is.finite(x) & (x %in% shell) )
    return(i)
  }


