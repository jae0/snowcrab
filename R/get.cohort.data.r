
  get.cohort.data = function(ageclass, yrs, cats, varmap, M, outvar="total") {

    yrs = as.numeric(yrs)
    n.age = length( ageclass )
    n.yrs = length( yrs )

    YRS = c( (max(yrs) - (2*length(yrs))+1 ): max(yrs) )
    n.YRS = length( YRS )
    n.cats = length( cats )

  #  ageclass.indices =  c( 1:n.age )

    Q = array( data=0, dim=c(n.age, n.cats, n.YRS), dimnames=list(ageclass, cats, YRS) )

    for ( iyc in n.age:1) {
      for ( jyc in iyc:n.age ) {
        ac = ageclass[jyc]

        for ( iyr in 1:n.yrs ) {  # the year class (starting with yc5)
          yr = yrs[iyr]
          yr2 = yr - jyc + 1

          Q[ac,, as.character(yr2)] = M[ac,, as.character(yr) ]
        }
      }
   }

  return (Q)
  }


