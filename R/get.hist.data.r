    get.hist.data = function(set, det, save=F) {
      dd = det[is.finite(det$cw),]
      hdata = merge(dd, set, by=c("trip", "set"), all.x=T, all.y=F, sort=F)
      if (save) save(hdata, file="hist.Rdata", compress=T)
      return (hdata)
    }


