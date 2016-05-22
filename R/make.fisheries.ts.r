
  make.fisheries.ts = function(landings) {

      l =NULL
      l = aggregate(landings$landings, list(yr=landings$yr), function(x) sum(x, na.rm=T))
      names(l) = c("yr", "landings")
      l = factor2number(l, c("yr", "landings"))
      out = l

      l =NULL
      l = aggregate(landings$effort, list(yr=landings$yr), function(x) sum(x, na.rm=T))
      names(l) = c("yr", "effort")
      l = factor2number(l, c("yr", "effort") )
      out = merge(out, l, by="yr", all=T)

#      l =NULL
#      l = aggregate(landings$tac.tons, list(yr=landings$yr), function(x) sum(x, na.rm=T))
#      names(l) = c("yr", "tac.tons")
#      l = factor2number(l, c("yr", "tac.tons") )
#      out = merge(out, l, by="yr", all=T)

      out$landings.kg = out$landings
      out$cpue = out$landings.kg / out$effort
      out$landings.kt = out$landings / 1000 / 1000
#      out$tac.kt = out$tac.tons / 1000

      return (out)
    }


