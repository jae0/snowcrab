
  observer.get.counts.by.class = function( p, x ) {
    
    #good.data = which(  x$cw>=95 & x$sex==male & ( x$shell>=2 | x$durometer >= 68 ) &  x$prodcd_id==0 & is.finite(x$mat) & is.finite(x$cw)  )
    good.data = which(  x$cw>=95 & x$sex==male & ( x$shell>=2 | x$durometer >= 68 ) &  x$prodcd_id==0 & is.finite(x$mat) & is.finite(x$cw)  )
    
    x = x[ good.data ,]
    x$dummy =1
    
    classnames = p$nodes
    classes = p$nodelabels
    
    nclasses = length(classes)
    x$class = NA
    for (c in 1:nclasses) x$class[ filter.class(x, classes[c]) ] = classnames[c]
     
    s = tapply2( x=x, indices=c("fishyr", "class"), var="dummy", newvars=c("count"), func=function(x){ sum(x) } )

    names(s) = c("yr", "vars", "n")
    s = factor2number(s, c("yr", "n" ) )
    return(s)
  }


