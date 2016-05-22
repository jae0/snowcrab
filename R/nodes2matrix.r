
  nodes2matrix = function(x, cl) {

    nx = names(x)

    ageclass = cl$yclass
    cats = cl$cats
    varmap = cl$varmap

    n.age = length( ageclass )
    n.cats = length( cats )

    ageclass.indices =  c( 1:n.age )

    M = array( data=0, dim=c(n.age, n.cats), dimnames=list(ageclass, cats) )

    for ( ac in 1:n.age ) {
    for ( vn in 1:n.cats ) {
      vname = i = NULL
      vname = paste( cats[vn], ageclass[ac], sep=".")
      i = which( nx==vname )
      if (length(i) == 0) next
      if (length(i) > 1) {
          print (paste ( "Problem (duplication): ", ageclass[ac], cats[vn] ))
          print (x[i,])
          break
      }
      M[ac, vn] = x[i]
      }}

      M[!is.finite(M)] = 0

  return (M)

  }


