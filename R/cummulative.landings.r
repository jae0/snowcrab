
  cummulative.landings = function(x) {

    x$uniqueid = paste(x$gridid, x$yr, sep="~")
    ids = sort(unique(x$uniqueid))
    x$week = weeks(x$chrono) # need to to add chron objebct or week to grid before implementation

    v = "landings"
    ti = "week"
    for (i in ids) {
      j = x[ which(x$uniqueid==i),]
      y = as.data.frame(xtabs( as.integer(x[,v]) ~ as.factor(x[,
                ti]), exclude="" ))
      names(y) = c(ti, v)
      y$Csum = cumsum(y[, v])
# --- must add a flexible model that can do convex/concave relationships (log relationshpi may be sufficient)

#      extract params (t0, t50, t100, slope parameter)


    }
#    redo for total effort as well

  }


