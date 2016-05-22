
  logbook.fisheries.stats.merge = function( Z ) {
 
    # merge snow crab fishery stats to incoming data frame Z
    # using plon/plat/yr discretization
    # data exists for 1996 and more recent 
    
    nZ = nrow(Z)
    Z0 = Z[ which( Z$yr < 1996 ) ,]
    
    yrs = sort( unique( Z$yr ), decreasing=T )  # most recent will have data and is used to init var lists
    yrs = yrs[ which(yrs >= 1996) ]
    out = NULL
    for ( yr in yrs )  {
    
      X = Z[ which( Z$yr == yr ) ,]
      X$plon = floor(X$plon)
      X$plat = floor(X$plat)
      X$gridid = paste( 
        X$plon%/%p$fisheries.grid.resolution * p$fisheries.grid.resolution, 
        X$plat%/%p$fisheries.grid.resolution * p$fisheries.grid.resolution, 
        sep="." 
      )
  
      # fisheries data regridded
      fg = logbook.db( DS="logbook.gridded", p=p, yrs=yr )
      fg$gridid = as.character( fg$gridid )
      X = merge( X, fg, by="gridid", all.x=T, all.y=F, sort=F)
      X$gridid = NULL # no longer needed
      
      for (vl in setdiff(names(fg), c("plon", "plat","gridid") )) {
        X[!is.finite(X[,vl]),vl] = 0 # make sure NA's created by merge statement are set to 0
      }
      rm (fg); gc()

      out = rbind( out, X )  
    }

    # fill years without appropriate data with NAs
    if (nrow(Z0) > 0 ) { 
      out.names = names ( out ) 
      newvars = setdiff( out.names, names(Z0) )
      for (nv in newvars) Z0[,nv] = NA
      out = rbind(out, Z0[, names(out)] ) 
    }

    if (nrow(out) != nZ ) stop ("Merge error")
    
    return (out)

  }
