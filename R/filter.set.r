  filter.set = function(x, f) {
    i = c(1:dim(x)[1])
    if (!is.null(f$towquality)) i = intersect(i, filter.towquality(x$towquality, f$towquality))
    if (!is.null(f$region))     i = intersect(i, filter.region(x$cfa, f$region))
    if (!is.null(f$season))     i = intersect(i, filter.season(x$julian, f$season))
    return (i)
  }

