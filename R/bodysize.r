
  # ------------------------
  # mean body size

  bodysize = function(x, factors=c("trip", "set"), variable, logtransform=T) {
    if (logtransform) x[,variable]=log10(x[,variable])

    x = x[is.finite(x[,variable]),]

    m = as.data.frame.table(tapply(X=x[,variable], INDEX=x[,factors], FUN=mean, simplify=T, na.rm=T))
    n = as.data.frame.table(tapply(X=x[,variable], INDEX=x[,factors], FUN=length, simplify=T))
    v = as.data.frame.table(tapply(X=x[,variable], INDEX=x[,factors], FUN=var, simplify=T, na.rm=T))

    names(m) = c(factors, paste(variable,"mean",sep="."))
    names(n) = c(factors, paste(variable,"n",sep="."))
    names(v) = c(factors, paste(variable,"var",sep="."))

    m = factor2character(m, factors)
    v = factor2character(v, factors)
    n = factor2character(n, factors)

    all = merge(x=m, y=v, by=factors, all=T)
    all = merge(x=all, y=n, by=factors, all=T)

    return(all)
  }


