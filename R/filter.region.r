  filter.region = function( x, region, index=T ) {
    area = region.decompose(region)
    i = which( as.character(x) %in% area )
    if (index) return(i) else return(x[i,])
  }


