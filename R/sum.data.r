
  sum.data = function(x, factors=c("trip", "set"), variable) {

    if (variable=="number") x[,variable] = 1  # create a dummy variable to allow counting of numbers
    m = as.data.frame.table(
          tapply( X=x[,variable],
                  INDEX=x[,factors],
                  FUN=sum,
                  simplify=T,
                  na.rm=T ) )
    names(m) = c(factors, paste(variable,"tot",sep="."))
    m = factor2character(m, factors)

    return(m)
  }


