 
  get.fishery.stats.by.region = function( Reg="cfaall", y=NULL ) {
    
    landings = landings.db()
    
    if (is.null(y)) y = sort(unique(landings$yr) )
    out = data.frame( yr=y )

    if (Reg=="cfaall")  region = sort( unique(landings$cfa) ) # all data
    if (Reg=="cfanorth") region = c("cfa20", "cfa21", "cfa22", "cfanorth", "north")
    if (Reg=="cfasouth") region = c("cfa23", "cfa24", "cfasouth", "cfaslope")
    if (Reg=="cfa4x") region = "cfa4x"
    
    lnd = landings[ which(landings$cfa %in% region) ,]
     
    l = aggregate( lnd$landings, list(yr=lnd$yr), function(x) sum(x, na.rm=T))
    names(l) = c("yr", "landings")
    
    out = merge(out, l, by="yr", all.x=T, all.y=F, sort=T)

    lnd$cpue_direct = lnd$landings / lnd$effort
    lnd$cpue_direct[ which( lnd$cpue_direct > (650*0.454)) ] = NA  # same rule as in landings.db -- 650lbs/trap is a reasonable upper limit

    cpue = aggregate( lnd$cpue, list(yr=lnd$yr), function(x) mean(x, na.rm=T))
    names(cpue) = c("yr", "cpue")
  
    out = merge (out, cpue, by="yr", all.x=T, all.y=F, sort=T)
    
    out$effort = out$landings / out$cpue  ## estimate effort level as direct estimates are underestimates (due to improper logbook records)
    rownames(out) = out$yr
    
    return(out)
    
  }


