
  get.structured.data = function(ageclass, yrs, cats, varmap, K, outvar="total") {

    M=NULL
    n.age = length( ageclass )
    n.yrs = length( yrs )
    n.cats = length( cats )

    ageclass.indices =  c( 1:n.age )

    M = array( data=0, dim=c(n.age, n.cats, n.yrs), dimnames=list(ageclass, cats, yrs) )

    for ( yr in 1:n.yrs ) {
    for ( ac in 1:n.age ) {
    for ( vn in 1:n.cats ) {
      vname = i = NULL
      vname = unlist ( varmap[ cats[vn] ] ) [ac]
      i = which( K$vars==vname & K$yr==yrs[yr] )
      if (length(i) == 0) next
      if (length(i) > 1) {
          print (paste ( "Problem: ", yrs[yr], ageclass[ac], cats[vn] ))
          print (K[i,])
          break
      }
      M[ac, vn, yr] = K[i, outvar]
      }}}

  return (M)
  }

