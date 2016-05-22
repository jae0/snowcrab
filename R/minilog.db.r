
  minilog.db = function( DS="", Y=NULL ){

    minilog.dir = project.datadirectory("snowcrab", "data", "minilog" )
    minilog.rawdata.location = file.path( minilog.dir, "archive" )

    if (!is.null(Y)) {
      iY = which( Y>=1999 )  # no historical data prior to 1999
      if (length(iY)==0) return ("No data for specified years")
      Y = Y[iY]
    }

    if ( DS %in% c("basedata", "metadata", "load") ) {
      if (DS=="basedata" ){
        flist = list.files(path=minilog.dir, pattern="basedata", full.names=T, recursive=FALSE)
        if (!is.null(Y)) {
          mm = NULL
          for (yy in Y ) {
            ll = grep( yy, flist)
            if (length(ll)==0) return( NULL) # nothing to do
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
        flist = list.files(path=minilog.dir, pattern="metadata", full.names=T, recursive=FALSE)
        if (!is.null(Y)) {
          mm = NULL
          for (yy in Y ) {
            ll = grep( yy, flist)
            if (length(ll)==0) return( NULL ) # nothing to do
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
      dirlist = list.files(path=minilog.rawdata.location, full.names=T, recursive=T)
      oo = grep("backup", dirlist)
      if (length(oo) > 0) {
        backups = dirlist[ oo ]
        dirlist = dirlist[-oo]
      }

      nfiles = length(dirlist)
      filelist = matrix( NA, ncol=3, nrow=nfiles)

      for (f in 1:nfiles) {
        yr = minilogDate( fnMini=dirlist[f] )
        if (is.null(yr) ) next()
        if ( yr %in% Y ) filelist[f,] = c( f, dirlist[f], yr )
      }
      filelist = filelist[ which( !is.na( filelist[,1] ) ) , ]

      set = snowcrab.db( DS="setInitial" )  # set$chron is in local time America/Halifax

      for ( yr in Y ) {
        print(yr)
        fn.meta = file.path( minilog.dir, paste( "minilog", "metadata", yr, "rdata", sep="." ) )
        fn.raw = file.path( minilog.dir, paste( "minilog", "basedata", yr, "rdata", sep="." ) )
        fs = filelist[ which( as.numeric(filelist[,3])==yr ) , 2 ]

        if (length(fs)==0) next()

        basedata = NULL
        metadata = NULL
        for (f in 1:length(fs)) {
          if( yr >= 2014 ) {
            j = load.minilog.rawdata.one.file.per.day( fn=fs[f], f=f, set=set)  # variable naming conventions in the past
          } else {
            j = load.minilog.rawdata( fn=fs[f], f=f, set=set)  # variable naming conventions in the past
          }
          if (is.null(j)) next()
          metadata = rbind( metadata, j$metadata)
          basedata = rbind( basedata, j$basedata)
        }

        # now do a last pass for the "backups" ....
        # incomplete ....
        add.backup.minilogs=FALSE
        if (add.backup.minilogs) {
          fb = backups[ which( as.numeric(backups[,3])==yr ) , 2 ]
          for (f in 1:length(fb)) {
            j = load.minilog.rawdata.backups( fn=fb[f], f=f, set=set)  # variable naming conventions in the past
            if (is.null(j)) next()
            metadata = rbind( metadata, j$metadata)
            basedata = rbind( basedata, j$basedata)
          }
        }

        save( metadata, file=fn.meta, compress=TRUE )
        save( basedata, file=fn.raw, compress=TRUE )

      }

      minilog.db( DS="set.minilog.lookuptable.redo" )

      return ( minilog.dir )
    }

    # -----------------------------------------------

    if (DS %in% c("stats", "stats.redo" ) ) {

      if (DS %in% c("stats") ){
        flist = list.files(path=minilog.dir, pattern="stats", full.names=T, recursive=FALSE)
        if (!is.null(Y)) {
          mm = NULL
          for (yy in Y ) {
            ll = grep( yy, flist)
            if (length(ll)==0) return(NULL) # nothing to do
            if (length(ll)>0 ) mm = c( mm, ll)
          }
          if (length(mm) > 0 ) flist= flist[mm]
        }
        mini.stat = NULL
        for ( i in flist ) {
          load( i )
          mini.stat = rbind( mini.stat, miniStats )
        }
        mini.meta = minilog.db( DS="metadata", Y=Y )
        res = merge( mini.meta, mini.stat,  by="minilog_uid", all.x=TRUE, all.y=FALSE, sort=FALSE )
        if(any(duplicated(res[,c('trip','set')]))) {
            res = removeDuplicateswithNA(res,cols=c('trip','set'),idvar='dt')
          }

        # TODO:: move the following to the load.minilog funcition .. and remove chron dependence
        res$t0 = as.POSIXct( as.chron(res$t0), tz=tzone, origin=lubridate::origin )
        res$t1 = as.POSIXct( as.chron(res$t1), tz=tzone, origin=lubridate::origin )
        res$dt = as.numeric( res$t1 - res$t0  )
        res$timestamp = lubridate::ymd_hms( res$timestamp)

        return (res)
       }


      # "stats.redo" is the default action
      tzone = "America/Halifax"

      for ( yr in Y ) {
        print (yr )
        fn = file.path( minilog.dir, paste( "minilog.stats", yr, "rdata", sep=".") )
        miniStats = NULL
        miniRAW = minilog.db( DS="basedata", Y=yr )
        miniRAW$timestamp = lubridate::ymd_hms( miniRAW$chron)

        mta = minilog.db( DS="metadata", Y=yr )
        mta$timestamp = ymd_hms( mta$timestamp )

        rid = minilog.db( DS="set.minilog.lookuptable" )
        rid = data.frame( minilog_uid=rid$minilog_uid, stringsAsFactors=FALSE )
        rid = merge( rid, mta, by="minilog_uid", all.x=TRUE, all.y=FALSE )
        rid = rid[ which(rid$yr== yr) ,]
        #rid = rid[grepl('S19092004',rid$minilog_uid),]
        if (nrow(rid) == 0 ) next()

        for ( i in 1:nrow(rid)  ) {

          id = rid$minilog_uid[i]
          sso.trip = rid$trip[i]
          sso.set = rid$set[i]
          sso.station = rid$station[i]

          Mi = which( miniRAW$minilog_uid == id )
          if (length( Mi) == 0 ) next()
          M = miniRAW[ Mi, ]

          M$timestamp = as.POSIXct( M$chron, tz=tzone, origin=lubridate::origin )
          settimestamp= as.POSIXct( rid$setChron[i] , tz=tzone , origin=lubridate::origin )
          time.gate =  list( t0=settimestamp - dminutes(5), t1=settimestamp + dminutes(11) )

          print( paste( i, ":", id) )

          # default, empty container
          res = data.frame(z=NA, t=NA, zsd=NA, tsd=NA, n=NA, t0=NA, t1=NA, dt=NA)

          bad.list = c(
#'minilog.S20052000.10.NA.NA.NA.13',
#'minilog.S19092004.8.389.NA.NA.321',
#'minilog.S19062000.8.NA.NA.NA.165' ,
#"minilog.S07092002.12.NA.NA.NA.245",
#"minilog.S08092002.10.NA.NA.NA.254",
#'minilog.S12102002.8.NA.15.59.349',
#'minilog.S28052002.10.NA.19.30.445'
# "minilog.S24112009.4.370.NA.NA.276",
# "minilog.S08092010.3.178.NA.NA.170",
# "minilog.S21102010.9.341.14.51.252",
# "minilog.S25092010.8.36.NA.NA.33",
# "minilog.S27102010.3.918.8.11.423"
          )

          if (! ( id %in% bad.list ) ) {

            ndat = length(M$depth[!is.na(M$depth)])
            if( ndat > 15 ) {
              # defaults appropriate for more modern scanmar data have > 3500 pings
              # depth resolution is about 4-5 m
              bcp = list(
                id=id, nr=nrow(M), YR=yr,
                tdif.min=3, tdif.max=9, time.gate=time.gate, depth.min=20, depth.range=c(-20,30), eps.depth = 1
              )

              #if(id=="minilog.S18092004.6.392.13.9.326") browser()
              # if(id=="minilog.S22061999.8.NA.NA.NA.84") browser()
              # if(id=="minilog.S04102007.12.903.17.10.378") browser()

              bcp = bottom.contact.parameters( bcp ) # add other default parameters .. not specified above
              bc =  NULL
              bc = bottom.contact( x=M, bcp=bcp )
              ## bottom.contact.plot (bc)

              if ( is.null(bc) || ( exists( "res", bc) && ( ( !is.finite(bc$res$t0 ) || !is.finite(bc$res$t1 ) ) ) )) {
                 bc = bottom.contact( x=M, bcp=bcp )
              }

              if ( is.null(bc) || ( exists( "res", bc) && ( ( !is.finite(bc$res$t0 ) || !is.finite(bc$res$t1 ) ) ) )) {
                bcp$noisefilter.inla.h =0.1
                bc = bottom.contact( x=M, bcp=bcp )
              }

              if ( is.null(bc) || ( exists( "res", bc) && ( ( !is.finite(bc$res$t0 ) || !is.finite(bc$res$t1 ) ) ) )) {
                M$depth = jitter( M$depth, amount = bcp$eps.depth/10 )
                bcp$noisefilter.inla.h = 0.01
                bc = bottom.contact( x=M, bcp=bcp )
              }

              if ( is.null(bc) || ( exists( "res", bc) && ( ( !is.finite(bc$res$t0 ) || !is.finite(bc$res$t1 ) ) ) )) {
                M$depth = jitter( M$depth, amount = bcp$eps.depth/10 )
                bcp$noisefilter.inla.h = 0.1
                bc = bottom.contact( x=M, bcp=bcp )
              }

              if ( !is.null(bc$res) & exists( "res", bc) ) res = bc$res
            }

            if( ndat == 0) {
              # nothing to do now ...
              # headerall = rid[i,'headerall'] , z = NA, t =  NA, zsd =NA, tsd =NA, n = NA,t0 =NA , t1=NA, dt =NA)
              #todo tie in the seabird data bottom contact to get secondary temperature data}
            }

          }

          if (FALSE) {
            # to visualize/debug
            bottom.contact.plot(bc)
            ## --- NOTE modal seems to work best ... but
            # no single best method .. use the default which is the mean of all methods
            ##  likely due to greater precision and data density relative to minilog
            #            res$t0 = bc$smooth.method[1]
            #            res$t0 = bc$smooth.method[2]
            #            res$dt = bc$smooth.method[2] -  bc$smooth.method[1]
          }
          res$t0 = as.POSIXct(res$t0,origin=lubridate::origin, tz=tzone )
          res$t1 = as.POSIXct(res$t1,origin=lubridate::origin, tz=tzone )
          res$dt = as.numeric(res$dt)
          miniStats = rbind(miniStats, cbind( minilog_uid=id, res ) )
        }

        miniStats$minilog_uid =  as.character(miniStats$minilog_uid)
        minidt = miniStats$dt
        miniStats$dt = NA
        i = which(!is.na( minidt ) )
        if (length(i) >0 ) miniStats$dt[i] = times( minidt[i] )

        save( miniStats, file=fn, compress=TRUE )
      }

      return ( minilog.dir )
    }

    # --------------------------------

    if (DS %in% c("set.minilog.lookuptable", "set.minilog.lookuptable.redo") ) {

      fn = file.path( minilog.dir, "set.minilog.lookuptable.rdata" )

      if (DS=="set.minilog.lookuptable" ) {
        B = NULL
        if ( file.exists( fn) ) load (fn)
        return (B)
      }

      B = minilog.db( DS="metadata" )

      # double check .. should not be necessary .. but in case
      uuid = paste( B$trip, B$set, sep="." )
      dups = which( duplicated( uuid) )

      if (length(dups > 0 ) ) {
        toremove =NULL
        for (i in dups) {
          di = which( uuid == uuid[i] )
          tdiff = B$setChron[di] - B$timestamp[di]
          oo = which.min( abs( tdiff) )
          toremove = c(toremove, di[-oo] )
          print("----")
          print( "Matching based upon closest time stamps")
          print(B[di, ])
          print( "Choosing: ")
          print(B[di[oo], ])
          print("")
          toremove = c(toremove, di[-oo] )
        }
        B = B[ -toremove, ]
      }
      B = B[, c("trip", "set", "minilog_uid" )]
      save(B, file=fn, compress=TRUE )
      return(fn)
    }
	}


