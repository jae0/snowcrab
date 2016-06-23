
snowcrab.timeseries.db = function( DS="default", p=NULL, regions=c( "cfa4x", "cfanorth", "cfasouth", "cfaall" ), trim=0, vn=NULL ) {

  tsoutdir = project.datadirectory( "bio.snowcrab", "R" )

  if (DS == "default") return( snowcrab.timeseries.db( DS="biologicals") )

  if (DS == "biologicals.direct" ) {
    # \\ no saving .. just a direct one-off
    dat = snowcrab.db( DS ="set.biologicals", p=p )
    if (is.null(vn)) vn = c( "R0.mass", "t", "R1.mass" )
    yrs = sort(unique(dat$yr))
    #area designations
    for (a in regions) {
      dat[,a] = NA
      ai = NULL
      ai = filter.region.polygon(dat, a)
      if (length(ai) > 0) dat[ai,a] = a
    }
    tsdata = expand.grid( region=regions, year=yrs, variable=vn, stringsAsFactors=FALSE, KEEP.OUT.ATTRS=FALSE )
    tsdata$mean = NA
    tsdata$se = NA
    tsdata$n = NA
    tsdata$ub = NA
    tsdata$lb = NA

    for (v in vn) {
      if ( !is.numeric( dat[,v] ) ) next()
      XX = bio.indicators::variable.recode( dat[,v], v, direction="forward", db="snowcrab" ) # transform variables where necessary
      for (r in regions) {
        ri = which( dat[,r] == r)
        if (length(ri)==0) next()
        XXmean = tapply( XX[ri], INDEX=dat$year[ri], FUN=mean, na.rm=TRUE )
        XXn =  tapply( XX[ri], INDEX=dat$year[ri], FUN=function(x) length(which(is.finite(x))) )
        XXse = tapply( XX[ri], INDEX=dat$year[ri], FUN=sd, na.rm=TRUE ) / XXn
        tsi = which(tsdata$variable==v & tsdata$region==r)

        tsdata[ tsi,"mean"] = bio.indicators::variable.recode (XXmean[ tsdata[ tsi, "year"] ], v, direction="backward", db="snowcrab" )
        tsdata[ tsi,"n"] = XXn[ tsdata[ tsi, "year"] ]
        tsdata[ tsi,"se"] = bio.indicators::variable.recode (XXse[ tsdata[ tsi, "year"] ], v, direction="backward", db="snowcrab" )
        XXlb = XXmean - XXse* 1.96
        XXub = XXmean + XXse* 1.96
        tsdata[ tsi,"lb"] = bio.indicators::variable.recode (XXlb[ tsdata[ tsi, "year"] ], v, direction="backward", db="snowcrab" )
        tsdata[ tsi,"ub"] = bio.indicators::variable.recode (XXub[ tsdata[ tsi, "year"] ], v, direction="backward", db="snowcrab" )
      }
    }
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
    print( "This will take a bit of time... make this parallel ..." )

    dat = snowcrab.db( DS ="set.biologicals", p=p )
    vn = setdiff( colnames(dat), c("trip", "set", "set_type", "station", "lon1", "lat1", "towquality",
                                   "lon", "lat", "plon", "plat", "seabird_uid", "minilog_uid", "netmind_uid" ) )
    yrs = sort(unique(dat$yr))

    dat$year = as.character(dat$yr)

    #area designations
    for (a in regions) {
      dat[,a] = NA
      ai = NULL
      ai = filter.region.polygon(dat, a)
      if (length(ai) > 0) dat[ai,a] = a
    }

    tsdata = expand.grid( region=regions, year=as.character(yrs), variable=vn, stringsAsFactors=FALSE, KEEP.OUT.ATTRS=FALSE )
    tsdata$mean = NA
    tsdata$se = NA
    tsdata$n = NA
    tsdata$ub = NA
    tsdata$lb = NA

    for (v in vn) {
      if ( !is.numeric( dat[,v] ) ) next()
      XX = bio.indicators::variable.recode( dat[,v], v, direction="forward", db="snowcrab" ) # transform variables where necessary
      for (r in regions) {
        ri = which( dat[,r] == r)
        if (length(ri)==0) next()
        XXmean = tapply( XX[ri], INDEX=dat$year[ri], FUN=mean, na.rm=TRUE )
        XXn =  tapply( XX[ri], INDEX=dat$year[ri], FUN=function(x) length(which(is.finite(x))) )
        XXse = tapply( XX[ri], INDEX=dat$year[ri], FUN=sd, na.rm=TRUE ) / XXn
        tsi = which(tsdata$variable==v & tsdata$region==r)

        tsdata[ tsi,"mean"] = bio.indicators::variable.recode (XXmean[ tsdata[ tsi, "year"] ], v, direction="backward", db="snowcrab" )
        tsdata[ tsi,"n"] = XXn[ tsdata[ tsi, "year"] ]
        tsdata[ tsi,"se"] = bio.indicators::variable.recode (XXse[ tsdata[ tsi, "year"] ], v, direction="backward", db="snowcrab" )
        XXlb = XXmean - XXse* 1.96
        XXub = XXmean + XXse* 1.96
        tsdata[ tsi,"lb"] = bio.indicators::variable.recode (XXlb[ tsdata[ tsi, "year"] ], v, direction="backward", db="snowcrab" )
        tsdata[ tsi,"ub"] = bio.indicators::variable.recode (XXub[ tsdata[ tsi, "year"] ], v, direction="backward", db="snowcrab" )
      }
    }

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
    print( "This will take a bit of time... " )

    dat = snowcrab.db( DS ="set.biologicals", p=p )
    stations.in.2014 = unique( dat$station[ which(dat$yr==2014) ] )
    dat = dat[ which(dat$station %in% stations.in.2014 ),]

    vn = setdiff( colnames(dat), c("trip", "set", "set_type", "station", "lon1", "lat1", "towquality",
                                   "lon", "lat", "plon", "plat", "seabird_uid", "minilog_uid", "netmind_uid" ) )
    yrs = sort(unique(dat$yr))

    #area designations
    for (a in regions) {
      dat[,a] = NA
      ai = NULL
      ai = filter.region.polygon(dat, a)
      if (length(ai) > 0) dat[ai,a] = a
    }

    tsdata = expand.grid( region=regions, year=yrs, variable=vn, stringsAsFactors=FALSE, KEEP.OUT.ATTRS=FALSE )
    tsdata$mean = NA
    tsdata$se = NA
    tsdata$n = NA
    tsdata$ub = NA
    tsdata$lb = NA

    for (v in vn) {
      if ( !is.numeric( dat[,v] ) ) next()
      XX = bio.indicators::variable.recode( dat[,v], v, direction="forward", db="snowcrab" ) # transform variables where necessary
      for (r in regions) {
        ri = which( dat[,r] == r)
        if (length(ri)==0) next()
        XXmean = tapply( XX[ri], INDEX=dat$year[ri], FUN=mean, na.rm=TRUE )
        XXn =  tapply( XX[ri], INDEX=dat$year[ri], FUN=function(x) length(which(is.finite(x))) )
        XXse = tapply( XX[ri], INDEX=dat$year[ri], FUN=sd, na.rm=TRUE ) / XXn
        tsi = which(tsdata$variable==v & tsdata$region==r)

        tsdata[ tsi,"mean"] = bio.indicators::variable.recode (XXmean[ tsdata[ tsi, "year"] ], v, direction="backward", db="snowcrab" )
        tsdata[ tsi,"n"] = XXn[ tsdata[ tsi, "year"] ]
        tsdata[ tsi,"se"] = bio.indicators::variable.recode (XXse[ tsdata[ tsi, "year"] ], v, direction="backward", db="snowcrab" )
        XXlb = XXmean - XXse* 1.96
        XXub = XXmean + XXse* 1.96
        tsdata[ tsi,"lb"] = bio.indicators::variable.recode (XXlb[ tsdata[ tsi, "year"] ], v, direction="backward", db="snowcrab" )
        tsdata[ tsi,"ub"] = bio.indicators::variable.recode (XXub[ tsdata[ tsi, "year"] ], v, direction="backward", db="snowcrab" )
      }
    }

    save( tsdata, file=fn, compress=TRUE )
    return( fn)
  }


  # -------------------


  if (DS %in% c( "observer", "observer.redo" ) ) {
    #\\ "reduced" subset of stations found in 2014 ... to be comparable with smaller survey
    fn = file.path( tsoutdir, "snowcrab.observer.timeseries.rdata" )
    if (DS=="odb") {
      tsdata = NULL
      if (file.exists( fn) ) load(fn)
      return(tsdata)
    }

    dat = observer.db( DS="odb" )
    dat = dat[ which( dat$cw >= 95),]
    vn = c( "cw", "totmass", "abdomen", "chela", "shell", "durometer",  "cpue.kg.trap", "mass", "mat" )
    yrs = sort(unique(dat$yr))

    #area designations
    for (a in regions) {
      dat[,a] = NA
      ai = NULL
      ai = filter.region.polygon(dat, a)
      if (length(ai) > 0) dat[ai,a] = a
    }

    tsdata = expand.grid( region=regions, year=yrs, variable=vn, stringsAsFactors=FALSE, KEEP.OUT.ATTRS=FALSE )
    tsdata$mean = NA
    tsdata$se = NA
    tsdata$n = NA
    tsdata$ub = NA
    tsdata$lb = NA


    for (v in vn) {
      if ( !is.numeric( dat[,v] ) ) next()
      XX = bio.indicators::variable.recode( dat[,v], v, direction="forward", db="snowcrab" ) # transform variables where necessary
      for (r in regions) {
        ri = which( dat[,r] == r)
        if (length(ri)==0) next()
        XXmean = tapply( XX[ri], INDEX=dat$year[ri], FUN=mean, na.rm=TRUE )
        XXn =  tapply( XX[ri], INDEX=dat$year[ri], FUN=function(x) length(which(is.finite(x))) )
        XXse = tapply( XX[ri], INDEX=dat$year[ri], FUN=sd, na.rm=TRUE ) / XXn
        tsi = which(tsdata$variable==v & tsdata$region==r)

        tsdata[ tsi,"mean"] = bio.indicators::variable.recode (XXmean[ tsdata[ tsi, "year"] ], v, direction="backward", db="snowcrab" )
        tsdata[ tsi,"n"] = XXn[ tsdata[ tsi, "year"] ]
        tsdata[ tsi,"se"] = bio.indicators::variable.recode (XXse[ tsdata[ tsi, "year"] ], v, direction="backward", db="snowcrab" )
        XXlb = XXmean - XXse* 1.96
        XXub = XXmean + XXse* 1.96
        tsdata[ tsi,"lb"] = bio.indicators::variable.recode (XXlb[ tsdata[ tsi, "year"] ], v, direction="backward", db="snowcrab" )
        tsdata[ tsi,"ub"] = bio.indicators::variable.recode (XXub[ tsdata[ tsi, "year"] ], v, direction="backward", db="snowcrab" )
      }
    }


    save( tsdata, file=fn, compress=TRUE )
    return( fn)
  }

}



