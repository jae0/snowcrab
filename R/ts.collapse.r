
  ts.collapse = function ( x, vars ) {
    out = NULL
    for (v in vars) {
      u = x[which(x$variable==v),]
      wm = NULL
      for (w in sort(unique(as.numeric(as.character(u$yr))))) {
        q = u[ which(u$yr==w), ]
        mean = wtd.mean(q$mean, q$sumwgt, normwt=T, na.rm=T)
        variance = wtd.var(q$mean, q$sumwgt, normwt=T, na.rm=T) # from Hmisc
        nsets = length(is.finite(q$mean))
        wm = data.frame(mean=mean, variance=variance, nsets=nsets, yr=w, variable=v )
        out = rbind(out, wm)
      }
    }
  return (out)
  }


