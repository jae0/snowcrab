

  clean.surface.area = function (set,  qreject = c( 0.025, 0.975 )) {
     # data extremes to accept/reject
  
    # final edits / sanity checks
    # many good tows where surface area was not recorded in the historical data
    # filter : set missing values and any large extremes to year-specific medians
    # historical data not touched except if missing
   

    # ----------------------------------------------------------
    # tow distance
    # quantile truncation of net mensuration-based estimates
    distance.thresholds = quantile( set$distance, probs=qreject, na.rm=T )  # recalc after removing unrealistic extremes
    dist.throwout = which ( set$distance < distance.thresholds[1] | set$distance > distance.thresholds[2] )
    if (length( dist.throwout) > 0 ) set$distance[ dist.throwout ] = NA

    set$dist  = set$dist / 1000 # m->km; these are historical estimates from manually-determined examination (Moncton)
    set$dist0 = set$dist0 / 1000 # distanced based upon logged records from touchdown/lift off records

    # historical data assumed correct (Moncton source)
    j = which( is.finite( set$dist) )
    # plot( dist0 ~ dist, set[j,] )
    # plot( distance ~ dist, set[j,] )
    if (length(j) > 0 ) set$distance[j] = set$dist[j]

    # where distance is missing, assume ship log start-end positions correct
    j = which( is.finite( set$dist0 ) & !is.finite( set$distance) ) 
    # plot( dist0 ~ dist, set[j,] )
    if (length(j) > 0 ) set$distance[j] = set$dist0[j]

    # final pass for distances.. where missing, fill with average for each year
    j = which( !is.finite( set$distance) ) 
    if (length(j) > 0 ) {
      distyear = tapply( set$distance, set$yr, mean, na.rm=TRUE )
      set$distance[j] = distyear[ as.character(set$yr[j]) ]
    }
    
    j = which( !is.finite( set$distance) ) # should be zero 
    if (length(j) > 0 ) {
      set$distance[j] = mean( set$distance, na.rm=TRUE )
    }
 

    #---------------------------------------------------------------
    # net spread 
    spread.thresholds = c(0.004, 0.016) # initially fixed at realisitc values based upon net configuration 
    j = which( set$spread < spread.thresholds[1] | set$spread > spread.thresholds[2] )
    if (length( j >0 )) set$spread [ j] = NA
    
    # quantile truncation
    spread.thresholds = quantile( set$spread, probs=qreject, na.rm=T ) # recalc after removing unrealistic extremes 
    j = which( set$spread < spread.thresholds[1]  | set$spread > spread.thresholds[2] ) 
    if (length(j)>0) set$spread[ j ] = NA 
 
    # .. where missing, fill with average for each year
    j = which( !is.finite( set$spread) ) 
    if (length(j) > 0 ) {
      spreadyear = tapply( set$spread, set$yr, mean, na.rm=TRUE )
      set$spread[j] = spreadyear[ as.character(set$yr[j]) ]
    }
   
    j = which( !is.finite( set$spread) ) 
    if (length(j) > 0 ) {
      set$spread[j] = mean( set$spread, na.rm=TRUE )
    }
     
 
    # ------------------------------------------
    # time 
    time.thresholds = c(4,10) # hard limits
    ith = which( set$dt < time.thresholds[1]  | set$dt > time.thresholds[2] )  
    if (length( ith >0 )) set$dt [ ith ] = NA
    
    time.quants = quantile( set$dt, probs=qreject, na.rm=T ) 
    ith = which( set$dt < time.quants[1]  | set$dt > time.quants[2] )  
    if (length( ith >0 )) set$dt [ ith ] = NA
  
    j = which( !is.finite(set$dt) ) 
    if (length(j) > 0) set$sa[ j ] = NA  ## if time is unreliable then so is surface area
   
    # .. where missing, fill with average for each year
    j = which( !is.finite( set$dt) ) 
    if (length(j) > 0 ) {
      dtyear = tapply( set$dt, set$yr, mean, na.rm=TRUE )
      set$dt[j] = dtyear[ as.character(set$yr[j]) ]
    }
   
    j = which( !is.finite( set$dt) ) 
    if (length(j) > 0 ) {
      set$dt[j] = mean( set$dt, na.rm=TRUE )
    }

 
    #-----------------------------------------
    # surface area
    # up to this piint set$sa contains historical data, determined by "hand", now merge with modern data streams
    
    # hard limits (gating)
    qsah = c( 0.001, 0.01 ) # more extreme than these are really not possible -- force as errors
    iqsah = which( set$sa < qsah[1] | set$sa > qsah[2] ) # historical data
    if (length( iqsah)>0 )  set$sa[iqsah] = NA
    
    iqsah2 = which( set$surfacearea < qsah[1] | set$surfacearea > qsah[2] ) # modern data
    if (length( iqsah2)>0 )  set$surfacearea[iqsah2] = NA
    
    # quantile truncation
    qsa = quantile( set$sa, qreject, na.rm=T )    # historical data
    iqsa = which( set$sa < qsa[1] | set$sa > qsa[2] ) 
    if (length( iqsa)>0 )  set$sa[iqsa] = NA
    
    qsa2 = quantile( set$surfacearea, qreject, na.rm=T ) # ~ 0.0028 0.0057    # modern data streams ...
    iqsa2 = which( set$surfacearea < qsa2[1] | set$surfacearea > qsa2[2] ) 
    if (length( iqsa2)>0 )  set$surfacearea[iqsa2] = NA

    # surfacearea is determined by automatic mechanisms .. fill where no historical data
    j = which( is.finite(set$surfacearea) & !is.finite( set$sa) )  
    set$sa[j] = set$surfacearea[j]   

    j = which( !is.finite( set$sa) )  
    set$sa[j] = set$distance[j] * set$spread[j]
   
    j = which( !is.finite( set$sa ) )
    if (length(j) > 0 ) {
      sayear = tapply( set$sa, set$yr, mean, na.rm=TRUE )
      set$sa[j] = sayear[ as.character(set$yr[j]) ]
    }

    j = which( !is.finite( set$sa ) )
    if (length(j) > 0 ) set$sa[ j ] = mean( set$sa, na.rm=TRUE )

    return (set)
  }



