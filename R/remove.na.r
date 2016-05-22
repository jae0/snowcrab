
  remove.na = function( x, .model, not.factor=NA) {
  
    all.variables = all.vars(.model)
    dependent.var = all.variables[1] 
    main.effects = setdiff( all.variables, dependent.var )

    x = y = x[, all.variables]
 
    x = ifelse(!is.na(x) | !is.finite(x), 1, NA)
    good.data = which( is.finite( rowSums( x ) ) ) 
    y = y[good.data,]
    
    factor.variables = all.variables
    if ( is.null(not.factor) | length(not.factor >0) ) factor.variables = setdiff( all.variables, not.factor )
    
#    factor.variables = ifelse( is.na(not.factor), all.variables, setdiff( all.variables, not.factor ) )
#    factor.variables = factor.variables[ which( is.finite(factor.variables) ) ]
    if (length(   factor.variables ) >0 )  {
      for ( v in factor.variables ) {
        y[,v] = as.factor(y[,v])
        y[,v] = drop.levels(y[,v], reorder=F)
      }
    }

    return( y )
  }



