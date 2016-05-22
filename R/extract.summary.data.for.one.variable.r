  
  #  extract summary information for ONE variable at the time 

   extract.summary.data.for.one.variable = function( x, y ) {
    # extract various (unadjusted) summary stats: medians, quantiles, etc
      d0 = tapply(X=x, INDEX=y, FUN=mean, na.rm=T, simplify=T)
      d1 = tapply(X=x, INDEX=y, FUN=sd, na.rm=T, simplify=T)
      d2 = tapply(X=x, INDEX=y, FUN=function(x) length(x[is.finite(x)]), simplify=T)
      d3 = tapply(X=x, INDEX=y, FUN=median, na.rm=T, simplify=T)
      d4 = tapply(X=x, INDEX=y, FUN=function(x) quantile(x, probs=0.25, na.rm=T), simplify=T)
      d5 = tapply(X=x, INDEX=y, FUN=function(x) quantile(x, probs=0.75, na.rm=T), simplify=T)
      out = rbind( mean=exp(d0)-1, sd=exp(d1), n=d2, median=exp(d3)-1, q25=exp(d4)-1, q75=exp(d5)-1 )
      print (out)
    return (out)
  }


