
	# mostly a copy over of the MINILOG functions with variable nemaes being replaced

  seabird.db = function( DS="", Y=NULL, plotdata=TRUE ){

    sb.dir = project.datadirectory("bio.snowcrab", "data", "seabird" )
    seabird.rawdata.location = file.path( sb.dir, "archive" )

    if (!is.null(Y)) {
        iY = which( Y>=2012 )  # no historical data prior to 2012
        if (length(iY)==0) return ("No data for specified years")
        Y = Y[iY]
    }


    if ( DS %in% c("basedata", "metadata", "load") ) {

      if (DS=="basedata" ){
        flist = list.files(path=sb.dir, pattern="basedata", full.names=T, recursive=FALSE)
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
        flist = list.files(path=sb.dir, pattern="metadata", full.names=T, recursive=FALSE)
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
          out= rbind( out, metadata )
        }
        return( out )
      }

      # default is to "load"

      dirlist = list.files(path=seabird.rawdata.location, full.names=T, recursive=T)
      nfiles = length(dirlist)
      filelist = matrix( NA, ncol=3, nrow=nfiles)
      for (f in 1:nfiles) {
        yr = seabirdDate( fnSeaBird=dirlist[f] )
        if (is.null(yr) ) next()
        if ( yr %in% Y ) filelist[f,] = c( f, dirlist[f], yr )
      }
      fli = which( !is.na( filelist[,1] ) )
      if (length(fli) == 0) return( "No data found matching the criteria")

      filelist = filelist[ fli , ]

      set = snowcrab.db( DS="setInitial" )

      for ( yr in Y ) {
        print(yr)
        fn.meta = file.path( sb.dir, paste( "seabird", "metadata", yr, "rdata", sep="." ) )
        fn.raw = file.path( sb.dir, paste( "seabird", "basedata", yr, "rdata", sep="." ) )
        fs = filelist[ which( as.numeric(filelist[,3])==yr ) , 2 ]
        if (length(fs)==0) next()
        basedata = NULL
        metadata = NULL
        for (f in 1:length(fs)) {
          j = load.seabird.rawdata( fn=fs, f=f, set=set, plotdata=plotdata )
          if (is.null(j)) next()
          metadata = rbind( metadata, j$metadata)
          basedata = rbind( basedata, j$basedata)
        }

        save( metadata, file=fn.meta, compress=TRUE )
        save( basedata, file=fn.raw, compress=TRUE )
      }

      seabird.db( DS="set.seabird.lookuptable.redo" )

      return ( sb.dir )
    }


    if (DS %in% c("stats", "stats.redo" ) ) {

      dir.create( file.path( sb.dir, "figures" ), recursive=TRUE, showWarnings=FALSE )

      if (DS=="stats") {

        flist = list.files(path=sb.dir, pattern="stats", full.names=T, recursive=FALSE)
        if (!is.null(Y)) {
          mm = NULL
          for (yy in Y ) {
            ll = grep( yy, flist)
            if (length(ll)==0) return( NULL) # nothing to do
            if (length(ll)>0 ) mm = c( mm, ll)
          }
          if (length(mm) > 0 ) flist= flist[mm]
        }

        sb.stat = NULL
        for ( i in 1:length(flist )) {
          load( flist[i] )
          sb.stat = rbind( sb.stat, sbStats )
        }

        sb.meta = seabird.db( DS="metadata", Y=Y )

        res = merge( sb.meta, sb.stat,  by="seabird_uid", all.x=TRUE, all.y=FALSE, sort=FALSE )

        if(any(duplicated(res[,c('trip','set')]))) {
          res = removeDuplicateswithNA(res,cols=c('trip','set'),idvar='dt')
        }

        tzone = "UTC"
        res$t0 = as.POSIXct( lubridate::ymd_hms(res$t0), origin=lubridate::origin, tz=tzone )
        res$t1 = as.POSIXct( lubridate::ymd_hms(res$t1), origin=lubridate::origin, tz=tzone )
        # res$timestamp = lubridate::ymd_hms( res$timestamp)

        return(res)
      }

      # default action  is "stats.redo"
      for ( yr in Y ) {
        print (yr )

        fn = file.path( sb.dir, paste( "seabird.stats", yr, "rdata", sep=".") )
        sbRAW = mta = sbStats = NULL

        sbRAW = seabird.db( DS="basedata", Y=yr )
        mta = seabird.db( DS="metadata", Y=yr )

        if (is.null(sbRAW) | is.null(mta)) next()

        rid = seabird.db( DS="set.seabird.lookuptable" )
        rid = data.frame( seabird_uid=rid$seabird_uid, stringsAsFactors=FALSE )
        rid = merge( rid, mta, by="seabird_uid", all.x=TRUE, all.y=FALSE )
        rid = rid[ which(rid$yr== yr) ,]

        if (nrow(rid) == 0 ) next()

  #      if (length(fli) == 0) return( "No data found matching the criteria")

        #prune down the rids to only a subset
        #rid = rid[grepl('S30092013',rid$seabird_uid),]
        for ( i in 1:nrow(rid) ) {
          #browser()
          print(i)

          id = rid$seabird_uid[i]

          if (id=="seabird.S06112014.8.210.18.44.9") next()  # strange problems

          sso.trip = rid$trip[i]
          sso.set = rid$set[i]
          sso.station = rid$station[i]

          Mi = which( sbRAW$seabird_uid == id )
          if (length( Mi) == 0 ) next()
          M = sbRAW[ Mi, ]

          settimestamp= rid$timestamp[i]
          time.gate =  list( t0=settimestamp - dminutes(5), t1=settimestamp + dminutes(10) )

          bcp = list(id=id, nr=nrow(M), YR=yr, tdif.min=3, tdif.max=9, time.gate=time.gate, depth.min=20, depth.range=c(-40,20),
                     smooth.windowsize=5, modal.windowsize=5, noisefilter.trim=0.01, noisefilter.target.r2=0.9, noisefilter.quants=c(0.025, 0.975) )

          bcp = bottom.contact.parameters( bcp ) # add other default parameters

          print(id)
          bc = NULL
          bc = bottom.contact( x=M, bcp=bcp )
          #browser()

          if ( is.null(bc) || ( !is.null(bc$res)  && ( ( !is.finite(bc$res$t0 ) || !is.finite(bc$res$t1 ) ) ) )) {
             bc = bottom.contact( x=M, bcp=bcp )
          }
          if ( is.null(bc) || ( !is.null(bc$res) && ( ( !is.finite(bc$res$t0 ) || !is.finite(bc$res$t1 ) ) ) )) {
            # try once more with random settings
            bcp$noisefilter.inla.h = 0.01
            bc = bottom.contact( x=M, bcp=bcp )
          }
          if ( is.null(bc) || ( !is.null(bc$res) && ( ( !is.finite(bc$res$t0 ) || !is.finite(bc$res$t1 ) ) ) )) {
            M$depth = jitter( M$depth, amount = bcp$eps.depth/10 )
            bcp$noisefilter.inla.h = 0.1
            bc = bottom.contact( x=M, bcp=bcp )
          }
          if ( is.null(bc) || ( !is.null(bc$res) && ( ( !is.finite(bc$res$t0 ) || !is.finite(bc$res$t1 ) ) ) )) {
            M$depth = jitter( M$depth, amount = bcp$eps.depth/10 )
            bcp$noisefilter.inla.h = 0.01
            bc = bottom.contact( x=M, bcp=bcp )
          }

          # default, empty container
          res = data.frame(z=NA, t=NA, zsd=NA, tsd=NA, n=NA, t0=NA, t1=NA, dt=NA)
          if ( !is.null( bc$res ) ) res = bc$res

          if (!is.null(bc)) {
            if (plotdata) {
              bottom.contact.plot( bc )
              plotfn = file.path( sb.dir, "figures", paste(id, "pdf", sep="." ) )
              print (plotfn)
              dev.flush()
              dev.copy2pdf( file=plotfn )
            }
          }

          if ( !is.null(bc)) if ( !is.null(bc$res )) res = bc$res

          sbStats = rbind( sbStats, cbind( seabird_uid=id, res ) )
        }

        sbStats$seabird_uid =  as.character(sbStats$seabird_uid)

        save( sbStats, file=fn, compress=TRUE)

      }

      return ( sb.dir )
    }


    # ----------------------

    if ( DS %in% c("set.seabird.lookuptable", "set.seabird.lookuptable.redo" ) ) {

      fn = file.path( sb.dir, "set.seabird.lookuptable.rdata" )

      if ( DS=="set.seabird.lookuptable" ) {
        B = NULL
        if ( file.exists( fn) ) load (fn)

        return (B)
      }

      B = seabird.db( DS="metadata" )

      # duplicates are created when datafiles are not reset completely
      # and the same data is reloaded during the survey ... cleanup at this stage

      # seabird_uid cannot be used as it has a file number attached
      # to make each one unique, even with duplicated data streams
      uuid = paste( B$trip, B$set, sep="." )
      dups = which( duplicated( uuid) )

      if (length(dups > 0 ) ) {
        toremove =NULL
        for (i in dups) {
          di = which( uuid == uuid[i] )

          tdiff = difftime( B$set_timestamp[di],  B$timestamp[di] )
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

      B =  B[, c("trip", "set", "seabird_uid" )]
      save(B, file=fn, compress=TRUE )
      return(fn)

    }

	}

  # --------------------


