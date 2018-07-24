
  minilog.db = function( DS="", Y=NULL, plotdata=TRUE ){

    minilog.dir = project.datadirectory("bio.snowcrab", "data", "minilog" )
    minilog.rawdata.location = file.path( minilog.dir, "archive" )
    plotdir = project.datadirectory("bio.snowcrab", "data", "minilog", "figures" )

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
      fli = which( !is.na( filelist[,1] ) )
      if ( length(fli) == 0) return( "No files matching the criteria.")

      filelist = filelist[ fli , ]

      set = snowcrab.db( DS="setInitial" )  # set$timestamp is in UTC

      for ( yr in Y ) {
        print(yr)
        #browser()
        fn.meta = file.path( minilog.dir, paste( "minilog", "metadata", yr, "rdata", sep="." ) )
        fn.raw = file.path( minilog.dir, paste( "minilog", "basedata", yr, "rdata", sep="." ) )
        fs = filelist[ which( as.numeric(filelist[,3])==yr ) , 2 ]

        if (length(fs)==0) next()

        basedata = NULL
        metadata = NULL
        for (f in 1:length(fs)) {
          if( yr >= 2014 ) {
            j = load.minilog.rawdata.one.file.per.day( fn=fs[f], f=f, set=set)
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
          stop( "TODO")
          fb = backups[ which( as.numeric(backups[,3])==yr ) , 2 ]
          for (f in 1:length(fb)) {
            j = load.minilog.rawdata.backups( fn=fb[f], f=f, set=set)
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
            res = removeDuplicateswithNA(res,cols=c('trip','set'),idvar='t')
          }
        #res$t0 = as.POSIXct( res$t0, tz="UTC", origin=lubridate::origin )
        #res$t1 = as.POSIXct( res$t1, tz="UTC", origin=lubridate::origin )
        #res$dt = difftime( res$t1, res$t0 )
       ids = paste(res$trip, res$set, sep=".") 
       uu = which( duplicated( ids ) )
       if (length(uu)>0 ) {
          message( "Duplicated trip/set found .. please fix this at the data level:" )
          toshow = which( ids %in% ids[uu] )
          print( res[ toshow,] )
       }

       return (res)
      }

      # "stats.redo" is the default action

      #      bad.list = c(
      #"minilog.S02112006.9.151.22.14.142",
      #"minilog.S27042001.7.NA.18.7.17",
      #"minilog.S08112008.9.55.NA.NA.55",
      #"minilog.S12102011.12.129.NA.NA.145",
      #"minilog.S18102007.11.226.18.44.198",
      #"minilog.S23102007.6.308.13.28.232",
      #"minilog.S27092007.9.86.NA.NA.87"
      #'minilog.S12071999.1.NA.NA.NA.190',
      #'minilog.S20052000.10.NA.NA.NA.13',
      #'minilog.S19092004.8.389.NA.NA.321',
      #'minilog.S19062000.8.NA.NA.NA.165',
      #'minilog.S07092002.12.NA.NA.NA.245',
      #'minilog.S08092002.10.NA.NA.NA.254',
      #'minilog.S12102002.8.NA.15.59.349',
      #'minilog.S28052002.10.NA.19.30.445',
      #'minilog.S24112009.4.370.NA.NA.276',
      #'minilog.S08092010.3.178.NA.NA.170',
      #'minilog.S21102010.9.341.14.51.252',
      #'minilog.S25092010.8.36.NA.NA.33',
      #'minilog.S27102010.3.918.8.11.423' '
      #      )
      bad.list = NULL
      bad.list = unique( c(bad.list, p$netmensuration.problem ) )

      for ( yr in Y ) {
        print (yr )

        fn = file.path( minilog.dir, paste( "minilog.stats", yr, "rdata", sep=".") )
        mta = miniRAW = miniStats = NULL
        miniRAW = minilog.db( DS="basedata", Y=yr )
        mta = minilog.db( DS="metadata", Y=yr )
        if (is.null(mta) | is.null(miniRAW)) next()

        rid = minilog.db( DS="set.minilog.lookuptable" )
        rid = data.frame( minilog_uid=rid$minilog_uid, stringsAsFactors=FALSE )
        rid = merge( rid, mta, by="minilog_uid", all.x=TRUE, all.y=FALSE )
        rid = rid[ which(rid$yr== yr) ,]
        if (nrow(rid) == 0 ) next()

        for ( i in 1:nrow(rid)  ) {
          id = rid$minilog_uid[i]
          sso.trip = rid$trip[i]
          sso.set = rid$set[i]
          sso.station = rid$station[i]

          Mi = which( miniRAW$minilog_uid == id )
          if (length( Mi) == 0 ) next()
          M = miniRAW[ Mi, ]

          settimestamp= rid$set_timestamp[i]
          time.gate =  list( t0=settimestamp - dminutes(6), t1=settimestamp + dminutes(12) )

          print( paste( i, ":", id) )

          # default, empty container
          res = data.frame(z=NA, t=NA, zsd=NA, tsd=NA, n=NA, t0=NA, t1=NA, dt=NA)

          rii = which( M$timestamp > settimestamp &  (M$timestamp < settimestamp+dminutes(5)) )
          # first estimate in case the following does not work
          if (length(rii) > 30) {
            res$z = mean(M$depth[rii], na.rm=TRUE)
            res$t = mean(M$temperature[rii], na.rm=TRUE)
            res$zsd = sd(M$depth[rii], na.rm=TRUE)
            res$tsd = sd(M$temperature[rii], na.rm=TRUE)
          }
          if (! ( id %in% bad.list ) ) {

            ndat = length( which( !is.na(M$depth) ))
            if (ndat ==0 ) print ("No depth data in minilogs")
            if( ndat < 30 ) {
              miniStats = rbind(miniStats, cbind( minilog_uid=id, res ) )

              next()
            } else {

              
              
              bcp = list(id=id, nr=nrow(M), YR=yr, tdif.min=3, tdif.max=11, time.gate=time.gate,
                         depth.min=20, depth.range=c(-25,15), eps.depth = 2 ,
                         smooth.windowsize=5, modal.windowsize=5,
                         noisefilter.trim=0.025, noisefilter.target.r2=0.85, noisefilter.quants=c(0.025, 0.975) )

              if(yr<2007)bcp$from.manual.archive=FALSE # manual touchdown only done since 2007
            
              bcp = bottom.contact.parameters( bcp ) # add other default parameters .. not specified above
              
              #BC: Determine if this station was done yet, if not then we want user interaction.
              if(file.exists(file.path(bcp$from.manual.archive, "clicktouchdown_all.csv"))){
                manualclick = read.csv(file.path(bcp$from.manual.archive, "clicktouchdown_all.csv"), as.is=TRUE)
                station = unlist(strsplit(bcp$id, "\\."))[4]
                sta.ind = which(manualclick$station == station & manualclick$year == bcp$YR)
                if(length(sta.ind == 1)) bcp$user.interaction = FALSE
                else bcp$user.interaction = TRUE
              }
  
              bc =  NULL

              bc = bottom.contact( x=M, bcp=bcp )

              redo = FALSE
              if ( is.null(bc) ) redo =TRUE
              if ( !is.null(bc) && exists("res", bc)) {
                if ( !is.finite(bc$res$t0 ) || !is.finite(bc$res$t1 ) ) redo = TRUE
              }
              if (redo) {
                 bcp$noisefilter.target.r2=0.8
                 bc = bottom.contact( x=M, bcp=bcp )
                 redo = FALSE
              }

              if ( is.null(bc) ) redo =TRUE
              if ( !is.null(bc) && exists("res", bc)) {
                if ( !is.finite(bc$res$t0 ) || !is.finite(bc$res$t1 ) ) redo = TRUE
              }
              if (redo) {
                 bcp$noisefilter.target.r2=0.75
                 bcp$noisefilter.trim=0.05
                 bcp$noisefilter.quants=c(0.025, 0.975)
                 bc = bottom.contact( x=M, bcp=bcp )
                 redo = FALSE
              }

              if (!is.null(bc) ) {
                if (plotdata) {
                  bottom.contact.plot( bc )
                  plotfn = file.path( plotdir, paste(id, "pdf", sep="." ) )
                  print (plotfn)
                  dev.flush()
                  dev.copy2pdf( file=plotfn )
                }
              }
              if ( !is.null(bc) && !is.null(bc$res) ) {
                res = bc$res
                miniStats = rbind(miniStats, cbind( minilog_uid=id, res ) )
              }
            } #end if dat

          } # end if badlist

        } #end nrow id
        # time needs to be reset as posix as it gets lost with rbind/cbind
        miniStats$minilog_uid =  as.character(miniStats$minilog_uid)
        miniStats$t0 = as.POSIXct(miniStats$t0,origin=lubridate::origin, tz="UTC" )
        miniStats$t1 = as.POSIXct(miniStats$t1,origin=lubridate::origin, tz="UTC")
        miniStats$dt = difftime( miniStats$t1, miniStats$t0 )

        # minidt = miniStats$dt
        # miniStats$dt = NA
        # i = which(!is.na( minidt ) )
        # if (length(i) >0 ) miniStats$dt[i] = minidt[i]

        save( miniStats, file=fn, compress=TRUE )
      } # end for year

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
          tdiff = difftime( B$set_timestamp[di], B$timestamp[di])
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


