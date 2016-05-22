
LDx = function( model, p, error.method, covariates) {
  # julian day at 50% transition of X-value
  
  tr = NULL
  ncov = length(covariates)

  if (  ncov==1) { 
    ncov1 = length( covariates[[1]] )
    for ( i in 1:ncov1 ) {
      covars = c(covariates[[1]][i]  )
      tr = rbind( tr, LDx.core( obj=model, p=p, method=error.method, covars=covars ) )  # ######## must choose correct method depending on model specification
    }
  }  
    

  if (  ncov==2) { 
    ncov1 = length( covariates[[1]] )
    ncov2 = length( covariates[[2]] )
    
    for ( i in 1:ncov1 ) {
    for ( j in 1:ncov2 ) {
       covars = c(covariates[[1]][i], covariates[[2]][j] )
      tr = rbind( tr, LDx.core( obj=model, p=p, method=error.method, covars=covars ) )  # ######## must choose correct method depending on model specification
    }}
  }  
    
  
  if (  ncov==3) { 
    ncov1 = length( covariates[[1]] )
    ncov2 = length( covariates[[2]] )
    ncov3 = length( covariates[[3]] )
    
    for ( i in 1:ncov1 ) {
    for ( j in 1:ncov2 ) {
    for ( k in 1:ncov3 ) {
      covars = c( covariates[[1]][i], covariates[[2]][j],  covariates[[3]][k] )
      tr = rbind( tr, LDx.core( obj=model, p=p, method=error.method, covars=covars ) )  # ######## must choose correct method depending on model specification
    }}}
  }  
       
   names( tr) = c( "Prob", "Dose", "SD", names(covariates) ) 

   print (" -------------------- " )
   print( "Expecting: y ~ aX + bT + cU + dXT + eXU ..." )
   print( "Order of terms is important: main effects first, then first order effects")
   print( "Where X is the variable of interest ... dose, size, date, etc)" )
   print( "When there are interaction terms, the first interaction term must be the second main effect" )
   print (" -------------------- " )
   print( " ")
   print( tr )

  return(tr)

}




