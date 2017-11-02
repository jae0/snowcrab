 
  # local function used below
  aggregated.fish.stats = function (sm ) { 
    # w is a time metric
    # lon, lat, w , catch , effort and cpue are expected in sm
    l = aggregate( sm$catch, list(w=sm$w, lon=sm$lon, lat=sm$lat), function(x) sum(x, na.rm=T))
    names(l) = c("w", "lon", "lat", "catch")
    l = factor2number(l, c("w", "lon", "lat", "catch"))
    rownames(l) = NULL
    l$id = paste( l$w, l$lon, l$lat, sep="~" )
    
    e = aggregate( sm$effort, list(w=sm$w, lon=sm$lon, lat=sm$lat), function(x) sum(x, na.rm=T))
    names(e) = c("w", "lon", "lat", "effort")
    e = factor2number(e, c("w", "lon", "lat", "effort"))
    rownames(e) = NULL
    e$id = paste( e$w, e$lon, e$lat, sep="~" )
      
    cc = aggregate( sm$cpue, list(w=sm$w, lon=sm$lon, lat=sm$lat), function(x) mean(x, na.rm=T))
    names(cc) = c("w", "lon", "lat", "cpue")
    cc = factor2number(cc, c("w", "lon", "lat", "cpue"))
    rownames(cc) = NULL
    cc$id = paste( cc$w, cc$lon, cc$lat, sep="~" )
  
    x = merge( l, e, by="id", sort=F, suffixes=c("",".e") )
    x = merge( x, cc,  by="id", sort=F, suffixes=c("",".cc") )
    x = x[ , c("w", "lon", "lat", "catch", "effort", "cpue") ]
    return (x)
  }


