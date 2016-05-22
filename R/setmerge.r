
   setmerge = function(X, Y, varname, filter, variable, index=NULL) {
    factors = c("trip", "set")
    #browser()
    print(varname)
    if (!is.null(filter)) {
      i = filter.class(Y, filter )
    } else {
      i = index  # "index" is the override mechanism
    }
    if (length(i)>0) {
      y = sum.data(Y[i,], factors, variable)
      names(y) = c(factors, varname)
      X = merge(x=X, y=y, by=factors, all.x=T )
      bi = X[which(X$trip == 'S01122015' & X$set==1 & X$station ==101),]
      #print('before')
      #print(bi$R0.mass)
      #print(bi$sa)
      X[,varname] = X[,varname] / X$sa   # express as x / km2
      #ci = X[which(X$trip == 'S01122015' & X$set==1 & X$station ==101),]
      #print(ci$R0.mass)
      X[!is.finite(X[,varname]),varname] = 0
      #if (varname=='R0.mass'){
      #  stop("......")
      #}
    } else {
      dummy = rep(0, dim(X)[1])
      oldnames = names(X)
      X = cbind(X, dummy)
      names(X) = c(oldnames, varname)
    }
    return(X)
  }


