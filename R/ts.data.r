

  ts.data = function ( vars, times="annual", db, delta=1) {
    V = dbGetQuery(gs, paste('select * from', db))
    V = merge.stratum.sa(V)
    out = NULL

    # lookup.table = snowcrab.db( p=p, DS="data.transforms" )

    for (v in vars) {
      j = which(colnames(V)==v)
      u = recode.time( V, times, delta=delta )
#      u = bio.snowcrab::variable.recode( V, v, direction="forward", lookuptable=lookuptable )
      u = u[is.finite(u$yr) ,]
      for (i in sort(unique( as.numeric(as.character(u$yr ))))) {
        q = u[ which(u$yr==i), ]
        wm =  means.strata(v=q[,j], strata=q$strat, w=q$area)
        wm$yr = i
        wm$variable = v
        out = rbind(out, wm)
      }
    }
  out = out[is.finite(out$mean),]
  return (out)
  }


