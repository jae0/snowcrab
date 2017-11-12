
snowcrab.timeseries.db = function( DS="default", p=NULL, regions=c( "cfa4x", "cfanorth", "cfasouth", "cfaall" ), trim=0, vn=NULL, sdci=F ) {

  tsoutdir = file.path( p$project.outputdir, "timeseries" )
  dir.create(tsoutdir, showWarnings=FALSE, recursive=TRUE)

  if (DS == "default") return( snowcrab.timeseries.db( DS="biologicals") )

  if (DS == "biologicals.direct" ) {
    # \\ no saving .. just a direct one-off
    dat = snowcrab.db( DS ="set.biologicals", p=p )
    dat$year = as.character(dat$yr)

    if (is.null(vn)) vn = c( "R0.mass", "t", "R1.mass" )
    yrs = sort(unique(dat$yr))

    #area designations
    for (a in regions) {
      dat[,a] = NA
      ai = NULL
      ai = emaf::polygon_inside(dat, a)
      if (length(ai) > 0) dat[ai,a] = a
    }
    tsdata = expand.grid( region=regions, year=yrs, variable=vn, stringsAsFactors=FALSE, KEEP.OUT.ATTRS=FALSE )
    tsdata$year = as.character( tsdata$year)
    tsdata$mean = NA
    tsdata$se = NA
    tsdata$sd = NA
    tsdata$n = NA
    tsdata$ub = NA
    tsdata$lb = NA

    for (vi in 1:length(vn) ) {
      v = vn[vi]
      if ( !is.numeric( dat[,v] ) ) next()
      print( paste( vi, v) )
      XX = emaf::variable.recode( dat[,v], v, direction="forward", db="snowcrab" ) # transform variables where necessary
      for (r in regions) {
        ri = which( dat[,r] == r)
        if (length(ri)==0) next()
        XXmean = tapply( XX[ri], INDEX=dat$year[ri], FUN=mean, na.rm=TRUE )
        XXn =  tapply( XX[ri], INDEX=dat$year[ri], FUN=function(x) length(which(is.finite(x))) )
        XXse = tapply( XX[ri], INDEX=dat$year[ri], FUN=sd, na.rm=TRUE ) / XXn
        XXsd = tapply( XX[ri], INDEX=dat$year[ri], FUN=sd, na.rm=TRUE ) 
        tsi = which(tsdata$variable==v & tsdata$region==r)

        tsdata[ tsi,"mean"] = emaf::variable.recode (XXmean[ tsdata[ tsi, "year"] ], v, direction="backward", db="snowcrab" )
        tsdata[ tsi,"n"] = XXn[ tsdata[ tsi, "year"] ]
        tsdata[ tsi,"se"] = emaf::variable.recode (XXse[ tsdata[ tsi, "year"] ], v, direction="backward", db="snowcrab" )
        tsdata[ tsi,"sd"] = emaf::variable.recode (XXsd[ tsdata[ tsi, "year"] ], v, direction="backward", db="snowcrab" )
        XXlb = XXmean - XXse* 1.96
        XXub = XXmean + XXse* 1.96
        if(sdci){
          XXlb = XXmean - XXsd* 1.96
          XXub = XXmean + XXsd* 1.96
        }
        tsdata[ tsi,"lb"] = emaf::variable.recode (XXlb[ tsdata[ tsi, "year"] ], v, direction="backward", db="snowcrab" )
        tsdata[ tsi,"ub"] = emaf::variable.recode (XXub[ tsdata[ tsi, "year"] ], v, direction="backward", db="snowcrab" )
      }
    }
    tsdata$year = as.numeric( tsdata$year)
    return(tsdata)
  }


  # -------------------


  if (DS %in% c( "biologicals", "biologicals.redo" ) ) {
    fn = file.path( tsoutdir, "snowcrab.timeseries.rdata" )
    if (DS=="biologicals") {
      tsdata = NULL
      if (file.exists( fn) ) load(fn)
      return(tsdata)
    }

    dat = snowcrab.db( DS ="set.biologicals", p=p )
    dat$year = as.character(dat$yr)

    if (is.null(vn)) vn = setdiff( colnames(dat), c("trip", "set", "set_type", "station", "lon1", "lat1", "towquality", "lon", "lat", "plon", "plat", "seabird_uid", "minilog_uid", "netmind_uid" ) )
    yrs = sort(unique(dat$yr))

    print( "This might take a bit of time... " )
    print( paste( "There are", length(vn), "variables" ) )

    #area designations
    for (a in regions) {
      dat[,a] = NA
      ai = NULL
      ai = emaf::polygon_inside(dat, a) 
      if (length(ai) > 0) dat[ai,a] = a
    }

    tsdata = expand.grid( region=regions, year=yrs, variable=vn, stringsAsFactors=FALSE, KEEP.OUT.ATTRS=FALSE )
    tsdata$year = as.character( tsdata$year)
    tsdata$mean = NA
    tsdata$se = NA
    tsdata$sd = NA
    tsdata$n = NA
    tsdata$ub = NA
    tsdata$lb = NA

    for (vi in 1:length(vn) ) {
      v = vn[vi]
      if ( !is.numeric( dat[,v] ) ) next()
      print( paste( vi, v) )
      XX = emaf::variable.recode( dat[,v], v, direction="forward", db="snowcrab" ) # transform variables where necessary
      for (r in regions) {
        ri = which( dat[,r] == r)
        if (length(ri)==0) next()
        XXmean = tapply( XX[ri], INDEX=dat$year[ri], FUN=mean, na.rm=TRUE )
        XXn =  tapply( XX[ri], INDEX=dat$year[ri], FUN=function(x) length(which(is.finite(x))) )
        XXse = tapply( XX[ri], INDEX=dat$year[ri], FUN=sd, na.rm=TRUE ) / XXn
        XXsd = tapply( XX[ri], INDEX=dat$year[ri], FUN=sd, na.rm=TRUE ) 
        tsi = which(tsdata$variable==v & tsdata$region==r)

        tsdata[ tsi,"mean"] = emaf::variable.recode (XXmean[ tsdata[ tsi, "year"] ], v, direction="backward", db="snowcrab" )
        tsdata[ tsi,"n"] = XXn[ tsdata[ tsi, "year"] ]
        tsdata[ tsi,"se"] = emaf::variable.recode (XXse[ tsdata[ tsi, "year"] ], v, direction="backward", db="snowcrab" )
        tsdata[ tsi,"sd"] = emaf::variable.recode (XXsd[ tsdata[ tsi, "year"] ], v, direction="backward", db="snowcrab" )
        XXlb = XXmean - XXse* 1.96
        XXub = XXmean + XXse* 1.96
        if(sdci){
          XXlb = XXmean - XXsd* 1.96 # confidence intervals for population instead of mean
          XXub = XXmean + XXsd* 1.96
        }
        tsdata[ tsi,"lb"] = emaf::variable.recode (XXlb[ tsdata[ tsi, "year"] ], v, direction="backward", db="snowcrab" )
        tsdata[ tsi,"ub"] = emaf::variable.recode (XXub[ tsdata[ tsi, "year"] ], v, direction="backward", db="snowcrab" )
      }
    }
    tsdata$year = as.numeric( tsdata$year)
    save( tsdata, file=fn, compress=TRUE )
    return( fn)
  }


  # -------------------


  if (DS %in% c( "biologicals.2014", "biologicals.2014.redo" ) ) {
    #\\ "reduced" subset of stations found in 2014 ... to be comparable with smaller survey
    fn = file.path( tsoutdir, "snowcrab.timeseries.2014.rdata" )
    if (DS=="biologicals.2014") {
      tsdata = NULL
      if (file.exists( fn) ) load(fn)
      return(tsdata)
    }

    dat = snowcrab.db( DS ="set.biologicals", p=p )
    dat$year = as.character(dat$yr)
    stations.in.2014 = unique( dat$station[ which(dat$yr==2014) ] )
    dat = dat[ which(dat$station %in% stations.in.2014 ),]

    vn = setdiff( colnames(dat), c("trip", "set", "set_type", "station", "lon1", "lat1", "towquality",
                                   "lon", "lat", "plon", "plat", "seabird_uid", "minilog_uid", "netmind_uid" ) )
    yrs = sort(unique(dat$yr))

    print( "This will take a bit of time... " )
    print( paste( "There are", length(vn), "variables" ) )

    #area designations
    for (a in regions) {
      dat[,a] = NA
      ai = NULL
      ai = emaf::polygon_inside(dat, a)
      if (length(ai) > 0) dat[ai,a] = a
    }

    tsdata = expand.grid( region=regions, year=yrs, variable=vn, stringsAsFactors=FALSE, KEEP.OUT.ATTRS=FALSE )
    tsdata$year = as.character( tsdata$year)
    tsdata$mean = NA
    tsdata$se = NA
    tsdata$n = NA
    tsdata$ub = NA
    tsdata$lb = NA

    for (vi in 1:length(vn) ) {
      v = vn[vi]
      if ( !is.numeric( dat[,v] ) ) next()
      print( paste( vi, v) )
      XX = emaf::variable.recode( dat[,v], v, direction="forward", db="snowcrab" ) # transform variables where necessary
      for (r in regions) {
        ri = which( dat[,r] == r)
        if (length(ri)==0) next()
        XXmean = tapply( XX[ri], INDEX=dat$year[ri], FUN=mean, na.rm=TRUE )
        XXn =  tapply( XX[ri], INDEX=dat$year[ri], FUN=function(x) length(which(is.finite(x))) )
        XXse = tapply( XX[ri], INDEX=dat$year[ri], FUN=sd, na.rm=TRUE ) / XXn
        tsi = which(tsdata$variable==v & tsdata$region==r)

        tsdata[ tsi,"mean"] = emaf::variable.recode (XXmean[ tsdata[ tsi, "year"] ], v, direction="backward", db="snowcrab" )
        tsdata[ tsi,"n"] = XXn[ tsdata[ tsi, "year"] ]
        tsdata[ tsi,"se"] = emaf::variable.recode (XXse[ tsdata[ tsi, "year"] ], v, direction="backward", db="snowcrab" )
        XXlb = XXmean - XXse* 1.96
        XXub = XXmean + XXse* 1.96
        tsdata[ tsi,"lb"] = emaf::variable.recode (XXlb[ tsdata[ tsi, "year"] ], v, direction="backward", db="snowcrab" )
        tsdata[ tsi,"ub"] = emaf::variable.recode (XXub[ tsdata[ tsi, "year"] ], v, direction="backward", db="snowcrab" )
      }
    }
    tsdata$year = as.numeric( tsdata$year)
    save( tsdata, file=fn, compress=TRUE )
    return( fn)
  }


  # -------------------


  if (DS %in% c( "observer", "observer.redo" ) ) {
    #\\ "reduced" subset of stations found in 2014 ... to be comparable with smaller survey
    fn = file.path( tsoutdir, "snowcrab.observer.timeseries.rdata" )
    if (DS=="observer") {
      tsdata = NULL
      if (file.exists( fn) ) load(fn)
      return(tsdata)
    }

    dat = observer.db( DS="odb" )
    dat$year = as.character(dat$yr)
    dat = dat[ which( dat$cw >= 95),]
    vn = c( "cw", "totmass", "abdomen", "chela", "shell", "durometer",  "cpue.kg.trap", "mass", "mat" )
    yrs = sort(unique(dat$yr))

    print( "This will take a bit of time... " )
    print( paste( "There are", length(vn), "variables" ) )

    #area designations
    for (a in regions) {
      dat[,a] = NA
      ai = NULL
      ai = emaf::polygon_inside(dat, a)
      if (length(ai) > 0) dat[ai,a] = a
    }

    tsdata = expand.grid( region=regions, year=yrs, variable=vn, stringsAsFactors=FALSE, KEEP.OUT.ATTRS=FALSE )
    tsdata$mean = NA
    tsdata$se = NA
    tsdata$n = NA
    tsdata$ub = NA
    tsdata$lb = NA
    tsdata$year = as.character( tsdata$year)

    for (vi in 1:length(vn) ) {
      v = vn[vi]
      if ( !is.numeric( dat[,v] ) ) next()
      print( paste( vi, v) )
      XX = emaf::variable.recode( dat[,v], v, direction="forward", db="snowcrab" ) # transform variables where necessary
      for (r in regions) {
        ri = which( dat[,r] == r)
        if (length(ri)==0) next()
        XXmean = tapply( XX[ri], INDEX=dat$year[ri], FUN=mean, na.rm=TRUE )
        XXn =  tapply( XX[ri], INDEX=dat$year[ri], FUN=function(x) length(which(is.finite(x))) )
        XXse = tapply( XX[ri], INDEX=dat$year[ri], FUN=sd, na.rm=TRUE ) / XXn
        tsi = which(tsdata$variable==v & tsdata$region==r)

        tsdata[ tsi,"mean"] = emaf::variable.recode (XXmean[ tsdata[ tsi, "year"] ], v, direction="backward", db="snowcrab" )
        tsdata[ tsi,"n"] = XXn[ tsdata[ tsi, "year"] ]
        tsdata[ tsi,"se"] = emaf::variable.recode (XXse[ tsdata[ tsi, "year"] ], v, direction="backward", db="snowcrab" )
        XXlb = XXmean - XXse* 1.96
        XXub = XXmean + XXse* 1.96
        tsdata[ tsi,"lb"] = emaf::variable.recode (XXlb[ tsdata[ tsi, "year"] ], v, direction="backward", db="snowcrab" )
        tsdata[ tsi,"ub"] = emaf::variable.recode (XXub[ tsdata[ tsi, "year"] ], v, direction="backward", db="snowcrab" )
      }
    }
    tsdata$year = as.numeric( tsdata$year)
    save( tsdata, file=fn, compress=TRUE )
    return( fn)
  }

  # -------------------


  if (DS %in% c( "groundfish.t", "groundfish.t.redo" ) ) {
    #\\ "reduced" subset of stations found in 2014 ... to be comparable with smaller survey
    fn = file.path( tsoutdir, "groundfish.t.rdata" )
    if (DS=="groundfish.t") {
      tsdata = NULL
      if (file.exists( fn) ) load(fn)
      return(tsdata)
    }
    tsdata = data.frame(r=NA,yrs=NA,V3='t',meanval=NA,se=NA,n=NA,ub=NA,lb=NA)
    h = groundfish.db('gshyd')
    g = groundfish.db('gsinf')
    g = g[,c('id','sdate','lon','lat','bottom_temperature')]
    h = h[,c('id','temp')]
    names(h)[2] <- 'bottom_temperature'
    f <- merge(g,h,by='id',all.x=T)
    i <- which(is.na(f$bottom_temperature.x) & !is.na(f$bottom_temperature.y))
    f[i,'bottom_temperature.x'] <- f[i,'bottom_temperature.y']
    f$yr <- as.numeric(format(f$sdate,'%Y'))
    f <- fishing.area.designations(f)
    ar <- unique(f$cfa)
    yy <- unique(f$yr)
    yy <- yy[order(yy)]
    for (r in ar) {
      for (yrs in yy) {
        y = f[which(f$yr == yrs & f$cfa ==r & !is.na(f$bottom_temperature.x)),]
        if(nrow(y)>3) {
          ym <- min(y$bottom_temperature.x[y$bottom_temperature.x>0])
          q = log(y$bottom_temperature.x+ym)
          m =  mean (q, na.rm=T)
          n = length(q)
          se = sd(q, na.rm=T)/ sqrt(n-1)
          meanval = exp(m)-ym
          ub = exp(m+se*1.96)-ym
          lb = exp(m-se*1.96)-ym
          j = as.data.frame(cbind(r, yrs, 't',meanval, se, n, ub, lb))
          tsdata <- rbind(tsdata,j)
        }
      }
    }
    #browser()
    colnames(tsdata) = c("region", "year", "variable","mean", "se","n", "ub", "lb")
    numbers = c("year", "mean", "se", "n", "ub", "lb")
    tsdata = factor2number(tsdata, numbers)
    tsdata <- tsdata[!is.na(tsdata$year),]

    save(tsdata, file=fn, compress=T)
    return(fn)
  }

}



