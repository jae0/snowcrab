  # ------------------------------
  # extract sex ratios

  sex.ratios = function(x, factors=c("trip", "set")) {

    mi = which(x$sex==male)
    nx = rep(1,dim(x[mi,])[1])
    m = as.data.frame(xtabs( nx ~ as.factor(trip) + as.factor(set), data=x[mi,] ) )
    names(m) = c(factors, "no.male")
    m = factor2character(m, factors)

    fi = which(x$sex==female)
    nx = rep(1,dim(x[fi,])[1])
    f = as.data.frame(xtabs( nx ~ as.factor(trip) + as.factor(set), data=x[fi,] ) )
    names(f) = c(factors, "no.female")
    f = factor2character(f, factors)

    all = merge(x=m, y=f, by=factors, all=T)
    all$no.female[!is.finite(all$no.female)] = 0
    all$no.male[!is.finite(all$no.male)] = 0
    all$sexratio = all$no.female / (all$no.female + all$no.male)
    all$sexratio[!is.finite(all$sexratio)] = NA

    return (all)
  }


