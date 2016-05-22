
  grid.logbook = function(x, res) {

    logbook.gridded = change.resolution (x, res=res)
    logbook.gridded$gridid = paste(logbook.gridded$lat, logbook.gridded$lon, logbook.gridded$year, sep="~")

    v = "pro_rated_slip_wt_lbs"
      x = logbook.gridded[is.finite(logbook.gridded[,v]),]
      grid = as.data.frame(xtabs( as.integer(x[,v]) ~ as.factor(x[,"gridid"]), exclude="" ))
      names(grid) = c("gridid", "landings")
      grid$landings = grid$landings * 0.45359237  # convert to kg

    v = "num_of_traps"
      x = logbook.gridded[is.finite(logbook.gridded[,v]),]
      effort = as.data.frame(xtabs( as.integer(x[,v]) ~ as.factor(x[,"gridid"]), exclude="" ))
      names(effort) = c("gridid", "notraps")

    grid = merge(grid, effort, by="gridid", all=T, sort=F)
    grid$cpue = grid$landings / grid$notraps

    tmp = matrix(unlist(strsplit(as.character(grid$gridid), "~")), ncol=3, byrow=T)
    grid$lat = as.numeric(tmp[,1])
    grid$lon = as.numeric(tmp[,2])
    grid$yr = as.numeric(tmp[,3])

#    grid = grid[is.finite(grid$yr),]

    return(grid)
  }


