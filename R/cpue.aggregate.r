
  cpue.aggregate = function( format="default"  ) {

    cpue = snowcrab_landings_db() 
    # year is year of capture
    # yr is "fishing year" relative to the assessment cycle 

    # if (format=="default") {
    #   landings.north = landings[ which(landings$cfa =="cfanorth"), ]
    #   ln = aggregate(landings.north$landings, list(yr=landings.north$yr), function(x) sum(x, na.rm=T))
    #   names(ln) = c("yr", "landings")
    #   ln = factor2number(ln, c("yr", "landings"))
    #   ln$id = paste(ln$yr, "cfanorth", sep=".")
    # 
    #   landings.south = landings[ which(landings$cfa =="cfasouth") , ]
    #   ls = aggregate(landings.south$landings, list(yr=landings.south$yr), function(x) sum(x, na.rm=T))
    #   names(ls) = c("yr", "landings")
    #   ls = factor2number(ls, c("yr", "landings"))
    #   ls$id = paste(ls$yr, "cfasouth", sep=".")
    # 
    #   landings.4x = landings[ which(landings$cfa =="cfa4x" ) , ]
    #   lx = aggregate(landings.4x$landings, list(yr=landings.4x$yr), function(x) sum(x, na.rm=T))
    #   names(lx) = c("yr", "landings")
    #   lx = factor2number(lx, c("yr", "landings"))
    #   lx$id = paste(lx$yr, "cfa4x", sep=".")
    # 
    #   # lx$landings[2:(nrow(lx))] = lx$landings[1:(nrow(lx)-1)]  # offset due to the fishing season being after the Survey
    # 
    #   la = rbind(ln, ls, lx)
    #   la$yr = NULL
    #   la$landings.kt = as.numeric(as.character(la$landings))/1000/1000
    #   return(la)
    # }

    if (format=="bugs") {
      regs =  c( "cfanorth", "cfasouth", "cfa4x" )
      cpue = cpue[ which (cpue$cfa %in% regs) , ]
      yrs = sort(unique( cpue$yr ) ) 
      L = tapply( cpue$landings, INDEX=cpue[,c("yr", "cfa")], FUN=sum, na.rm=T )
      L2 = tapply( cpue$effort, INDEX=cpue[,c("yr", "cfa")], FUN=sum, na.rm=T )
      L=L/L2
      L1 = tapply( cpue$landings, INDEX=cpue[,"yr"], FUN=sum, na.rm=T )
      L2 = tapply( cpue$effort, INDEX=cpue[,"yr"], FUN=sum, na.rm=T )
      cfaall = L1/L2
      L = cbind( L, cfaall )
      return (L)
    }

  }


