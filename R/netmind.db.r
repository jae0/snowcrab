netmind.db = function( DS, Y=NULL, plotdata=FALSE ) {

  netmind.dir = project.datadirectory("bio.snowcrab", "data", "netmind" )
  netmind.rawdata.location = file.path( netmind.dir, "archive" )

  if (!is.null(Y)) {
    iY = which( Y>=1999 )  # no historical data prior to 1998
    if (length(iY)==0) return ("No data for specified years")
    Y = Y[iY]
  }

  if(DS =='esonar2netmind.conversion') {

    if(is.null(Y) | any(Y < 2014)) stop('This only begins in 2014')
    for(y in Y) {
      #Changing to convert esonar directly to netmind. -Brent
      esonar.raw.location = file.path(project.datadirectory("bio.snowcrab", "data", "esonar", "archive" ), y)
      flist = list.files(path=esonar.raw.location, full.names=T, recursive=FALSE)
      for(fl in  flist){
        esonar2netmind(fl)
      }
    }
  }

  # -----------------------

  if ( DS %in% c("basedata", "metadata", "load") ) {

    if (DS=="basedata" ){
      flist = list.files(path=netmind.dir, pattern="basedata", full.names=T, recursive=FALSE)
      if (!is.null(Y)) {
        mm = NULL
        for (yy in Y ) {
          ll = grep( yy, flist)
          if (length(ll)>0 ) mm = c( mm, ll)
        }
        if (length(mm) > 0 ) flist= flist[mm]
      }
      out = NULL
      for ( i in flist ) {
        load( i )
        out= rbind( out, basedata )
      }
      return( out )
    }

    if (DS=="metadata" ){
      flist = list.files(path=netmind.dir, pattern="metadata", full.names=T, recursive=FALSE)
      if (!is.null(Y)) {
        mm = NULL
        for (yy in Y ) {
          ll = grep( yy, flist)
          if (length(ll)>0 ) mm = c( mm, ll)
        }
        if (length(mm) > 0 ) flist= flist[mm]
      }
      out = NULL
      for ( i in flist ) {
        load( i )
        out= rbind( out, metadata )
      }
      return( out )
    }

    # default is to "load"
    #

    if (any( Y < 2004) ) {
      print( "Net metrics and bottom contact stats (distance towed) were processed manually by Gulf Region until 2004 ")
      print( "and now stored in 'SNTOWS'. This is therefore redundant for historical data and only fills in time of ")
      print( "bottom contact, etc for the sake of completeness" )
    }

    dirlist = list.files(path=netmind.rawdata.location, full.names=T, recursive=T)
    # process every data file ... even bad tows .. marginal overhead in order to be complete (sometimes file names are wrong)

    nfiles = length(dirlist)
    filelist = matrix( NA, ncol=3, nrow=nfiles)
    for (f in 1:nfiles) {
      yr = netmindDate( fnNetmind=dirlist[f] )
      if ( is.null(yr) ) next()
      if ( yr %in% Y ) filelist[f,] = c( f, dirlist[f], yr )
    }
    filelist = filelist[ which( !is.na( filelist[,1] ) ) , ]

    set = snowcrab.db( DS="setInitial" )  # UTC

    for ( yr in Y ) {
      print(yr)
      fn.meta = file.path( netmind.dir, paste( "netmind", "metadata", yr, "rdata", sep="." ) )
      fn.raw = file.path( netmind.dir, paste( "netmind", "basedata", yr, "rdata", sep="." ) )
      fs = filelist[ which( as.numeric(filelist[,3])==yr ) , 2 ]

      if (length(fs)==0) next()

      basedata = NULL
      metadata = NULL

      for (f in 1:length(fs) ) {
        j = load.netmind.rawdata( fs[f], f=f, set=set )  # variable naming conventions in the past
        if (is.null(j)) next()
        metadata = rbind( metadata, j$metadata)
        basedata = rbind( basedata, j$basedata)
      }
      save( metadata, file=fn.meta, compress=TRUE )
      save( basedata, file=fn.raw, compress=TRUE )

    }
    # now that it is complete, refresh the set/uid lookup table
    netmind.db( DS="set.netmind.lookuptable.redo" )

    return ( netmind.dir  )
  }


   # ------------------


  if (DS %in% c("stats", "stats.redo" ) ) {

    if (DS %in% c("stats") ){

      flist = list.files(path=netmind.dir, pattern="stats", full.names=T, recursive=FALSE)
      if (!is.null(Y)) { # if Y is declared then subset ... default is to return all
        mm = NULL
        for (yy in Y ) {
          ll = grep( yy, flist)
          if (length(ll)>0 ) mm = c( mm, ll)
        }
        if (length(mm) > 0 ) flist= flist[mm]
      }
      netmind.stat = NULL
      for ( i in flist ) {
        load( i )
        netmind.stat = rbind( netmind.stat, Stats )
      }

      netmind.stat$yr = NULL
      nm = netmind.db( DS="set.netmind.lookuptable" )
      res = merge( nm, netmind.stat,  by="netmind_uid", all.x=TRUE, all.y=FALSE, sort=FALSE )
      # not really required but just in case missing values cause confusion with rbind
      #res$t0 = as.POSIXct( res$t0, origin=lubridate::origin, tz="UTC" )
      #res$t1 = as.POSIXct( res$t1, origin=lubridate::origin, tz="UTC" )
      #res$dt = difftime( res$t1, res$t0 )  # reset in case time info gets lost with rbind
      return (res)
    }

    # "stats.redo" is the default action
    # bring in stats from each data stream and then calculate netmind stats
    # bring in minilog and seabird data that has t0, t1 times for start and stop of bottom contact

    set = snowcrab.db( DS="setInitial")  # UTC

    sbStats =  seabird.db( DS="stats" )
    sbv = c('trip','set', "z", "zsd", "t", "tsd", "n", "t0", "t1", "dt" )
    set_sb = merge( set[, c("trip", "set") ], sbStats[,sbv], by=c("trip","set"), all.x=TRUE, all.y=FALSE, sort=FALSE )

    mlStats =  minilog.db( DS="stats" )
    mlv =  c('trip','set', "z",    "zsd",    "t",    "tsd",    "n",    "t0",    "t1",    "dt" )
    set_ml = merge( set[, c("trip", "set") ], mlStats[,mlv], by=c("trip","set"), all.x=TRUE, all.y=FALSE, sort=FALSE )

    set = merge( set, set_sb, by=c("trip", "set" ), all.x=TRUE, all.y=FALSE, sort=FALSE, suffixes=c("", ".sb" ) )
    set = merge( set, set_ml, by=c("trip", "set" ), all.x=TRUE, all.y=FALSE, sort=FALSE, suffixes=c("", ".ml" ))

    # use seabird data as the standard, replace with minilog data where missing
    ii = which(!is.finite( set$t0) )
    if (length(ii) > 0 )  set$t0[ ii] = set$t0.ml[ii]

    ii = which(!is.finite( set$t1) )
    if (length(ii) > 0 )  set$t1[ ii] = set$t1.ml[ii]

    ii = which(!is.finite( set$z) )
    if (length(ii) > 0 )  set$z[ ii] = set$z.ml[ii]

    ii = which(!is.finite( set$zsd) )
    if (length(ii) > 0 )  set$zsd[ ii] = set$zsd.ml[ii]

    ii = which(!is.finite( set$t) )
    if (length(ii) > 0 )  set$t[ ii] = set$t.ml[ii]

    ii = which(!is.finite( set$tsd) )
    if (length(ii) > 0 )  set$tsd[ ii] = set$tsd.ml[ii]

    ii = which(!is.finite( set$dt) )
    if (length(ii) > 0 )  set$dt[ ii] = set$dt.ml[ii]

    tokeep = grep( "\\.ml$", colnames(set), invert=TRUE )
    set = set[, tokeep]
    set$n = NULL

    nm = netmind.db( DS="set.netmind.lookuptable" )
    set = merge( set, nm, by=c("trip","set"), all.x=T, all.y=F, sort=F, suffixes=c("", ".netmind") )

    # add more data .. t0,t1, dt where missing and width and SA estimates where possible
    for ( yr in Y ) {
      print(yr)
      fn = file.path( netmind.dir, paste( "netmind.stats", yr, "rdata", sep=".") )
      Stats = NULL
      basedata = netmind.db( DS="basedata", Y=yr )
      ii = which( set$yr==yr & !is.na(set$netmind_uid) )
      nii =  length( ii )
      if ( nii== 0 ) next()
      rid = set[ ii,]
      Stats = NULL
      for ( i in 1:nii  ){
       print(i)
        id = rid$netmind_uid[i]
        print(rid[i,])
        bdi = which( basedata$netmind_uid==id )
        if (length(bdi) < 5 ) next()
        l = net.configuration( basedata[ bdi ,], t0=rid$t0[i], t1=rid$t1[i], set_timestamp=rid$timestamp[i], yr=yr, plotdata=plotdata )
        #if(is.na(l$surfacearea))browser()
        l$netmind_uid = id
        # not really required but just in case missing values cause confusion with rbind
        Stats = rbind( Stats, l )
      }
      if (is.null(Stats)) next()
      Stats$t0 = as.POSIXct(Stats$t0, origin=lubridate::origin, tz="UTC" )
      Stats$t1 = as.POSIXct(Stats$t1, origin=lubridate::origin, tz="UTC")
      Stats$dt = difftime( Stats$t1, Stats$t0 )
      save( Stats, file=fn, compress=TRUE )
    }
    return ( netmind.dir )
  }


  # -------------------


  if (DS %in% c("set.netmind.lookuptable", "set.netmind.lookuptable.redo") ) {

    fn = file.path( netmind.dir, "set.netmind.lookuptable.rdata" )

    if (DS=="set.netmind.lookuptable" ) {
      B = NULL
      if ( file.exists( fn) ) load (fn)
      return (B)
    }

    B = netmind.db( DS="metadata" )

    # double check .. should not be necessary .. but in case
    uuid = paste( B$trip, B$set, sep="." )
    dups = which( duplicated( uuid) )


    if (length(dups > 0 ) ) {
      toremove =NULL
      for (i in dups) {
        di = which( uuid == uuid[i] )
        w <- B[di,]
        tdiff = difftime( B$set_timestamp[di], B$netmind_timestamp[di])
        oo = which.min( abs( tdiff) )
        toremove = c(toremove, di[-oo] )
        print("----")
        print( "Matching based upon closest time stamps")
        print(B[di, ])
        print( "Choosing: ")
        print(B[di[oo], ])
        print("")

      }

      B = B[-toremove, ]

    }

    # double check .. should not be necessary .. but in case

    B = B[ , c("trip", "set", "netmind_uid" )]

    save(B, file=fn, compress=TRUE )
    return(fn)
  }

}
