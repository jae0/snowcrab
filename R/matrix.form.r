
  matrix.form = function( x, type="number" ) {

   # rescale variables .. log transforms and express per km2 before ordinations
    # ...  large magnitude numbers are needed for ordinations

    # form datasets
    # numbers and weights have already been converted to per km2 and with vessel corrections

      k = 1e4         # a large constant number to make xtabs work
      if (type == "number") {
        o = which(is.finite(x$totno))
        m = xtabs( as.numeric(as.integer(totno*k)) ~ as.factor(id) + as.factor(spec), data=x[o,] ) / k
      }
      if (type == "biomass") {
        o = which(is.finite(x$totmass))
        m = xtabs(as.numeric( as.integer(totwgt*k)) ~ as.factor(id) + as.factor(spec), data=x[o,] ) / k
      }

      threshold= min( x$totno[ which(x$totno>0) ], na.rm=T ) / nrow( m )
      finished.j = finished.i = F
      while( !(finished.j & finished.i) ) {
        nr = nrow(m)
        nc = ncol(m)
        rowsm = rowSums(m)
        colsm = colSums(m)
        i = unique( c( which( rowsm/nr <= threshold ), which(rowsm==0 ) ))
        j = unique( c( which( colsm/nc <= threshold ), which(colsm==0 ) ))
        if (length(i) > 0 ) {
          m = m[ -i , ]
        } else {
          finished.i = T 
        }
        if (length(j) > 0 ) {
          m = m[ , -j ]
        } else {
          finished.j = T
        }
      }

      n = array(m, dim=dim(m), dimnames(m))

      return(n)
    }



