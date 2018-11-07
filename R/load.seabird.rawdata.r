
  # ---------------------------------------- low-level data access

	load.seabird.rawdata = function(fn, f, set, plotdata=F) {

    out = NULL

    filename= fn[f]
    header =  readLines(filename, n=75)
    headerall = paste( header[1:12], collapse="~")

    if (length(header) < 75 )  return( out )

    data.start = grep( "start sample number", header, perl=T )
    if (grepl( "test", header[ 7 ], ignore.case=T) ) return(out)

    fileinfo = tolower(unlist(strsplit(filename, split="/")))
    seabird = as.data.frame(read.table( file=filename, sep=",", as.is=T, colClasses="character", header=F, skip= data.start, strip.white=T))

    if ( nrow(seabird) < 10 ) return( NULL )
    if(ncol(seabird)==1) {
      #AMC ADDED 2013 FOR FUNCTIONALITY AND IF SEABIRD DATA NOT SEPARAETED OUT BY COLUMNS
      seabird = as.data.frame(matrix(apply(seabird,2,function(x) unlist(strsplit(x,","))),ncol=4,byrow=T))
    }

    tz.system = "America/Halifax"
    tz.seabird = "America/Halifax"
    tz.snowcrab = "UTC"

    # if time is not correct in seabird file and is due to DST fix by BMC
    mhead = readLines(filename, n=40)
    #local time
    lt = mhead[grep("System UpLoad Time",mhead)]
    lt = unlist(strsplit(unlist(strsplit(lt,"= "))[2], " "))
    lt  = lubridate::mdy_hms( paste( paste(lt[1],lt[2],lt[3],sep="-"), lt[4]), tz=tz.system )

    #seabird time
    st = mhead[grep("SERIAL NO", mhead)]
    st = unlist(strsplit(unlist(strsplit(st,"    "))[2], "  "))
    st  = lubridate::dmy_hms( paste( st[1], st[2]), tz=tz.seabird )

    #add difference to seabird times
    colnames(seabird) = c( "temperature", "pressure", "mdate", "mtime")

    seabird$timestamp = lubridate::dmy_hms( paste( trimWhiteSpace(seabird$mdate), trimWhiteSpace(seabird$mtime) ), tz=tz.seabird )

    seabird$timestamp = with_tz( seabird$timestamp, tz.snowcrab )  # from now on, seabird is UTC

  # BC 2018 - Seabird header data not populated in single instance so comparing system time with seabird
  #           time is not accurate to get time offset betweem system and seabird. Seabird is time is synced 
  #           to system time every morning and should not be more than a few seconds in difference. 
  # seabird$timestamp = seabird$timestamp + difftime(lt,st)

    numerics = c("temperature", "pressure")
    seabird = factor2number(seabird, numerics)

    # check time first
    seabird.date.range = range( seabird$timestamp )

    setxi = which( set$timestamp >= seabird.date.range[1] & set$timestamp <= seabird.date.range[2] )
    if ( length( setxi ) == 0 ) return(NULL)

    yr = lubridate::year( seabird$timestamp[1] )

    # First pass:
    # break down multi-set records into separate records using a simple depth rule
    print (filename)

    seabird$pressure = adjust.depth.for.drift( seabird$pressure )
    surface = quantile( seabird$pressure, probs=0.5 ) # pressure at which it is assumed the sensor in at surface (most data will be shallow)

    # initiate new variables
    seabird$seabird_uid = NA
    seabird$depth = NA

    metadata = NULL
    for ( ssid in 1:length(setxi) ) {

	    filename2 = tolower( fileinfo[length(fileinfo)] )

      iset = setxi[ssid]
      setx =  set[iset,] # matching trip/set/station

      # find sets by running until some threshold depth is found, stating with set$chron
      istart = which.min( abs( difftime( seabird$timestamp, set$timestamp[iset] )) )

      j0 = istart
      while ( j0 > 1 ) {
        if (seabird$pressure[j0] <= surface ) break()
        j0 = j0 - 1
      }
      tdiff0 = abs( difftime( seabird$timestamp[j0], set$timestamp[iset] ))
      if ( tdiff0 > lubridate::minutes(10) ) {
        # then something went wrong .. default to a few minutes before set$chron
        j0 = which.min( abs( difftime( seabird$timestamp, (set$timestamp[iset] - lubridate::minutes(5) ) )) )
      }


      j1 = istart
      while ( j1 < nrow(seabird) ) {
        if (seabird$pressure[j1] <= surface ) break()
        j1 = j1 + 1
      }
      tdiff1 = abs( difftime( seabird$timestamp[j1], set$timestamp[iset] ))
      if ( tdiff1 > lubridate::minutes(15) ) {
        # then something went wrong .. default to a few minutes before set$chron
        j1 = which.min( abs( difftime( seabird$timestamp, (set$timestamp[iset] + lubridate::minutes(12) ) )) )
      }

      if(j0==j1){ #skip some entries and then re calculate
        j1=j1+20
          while ( j1 < nrow(seabird) ) {
        if (seabird$pressure[j1] <= surface ) break()
        j1 = j1 + 1
      }
      }

      o = j0:j1
      error = "" #dummy value
      if (length(o) < 30) {
        # no matching time range found ... use closest match and raise a flag
        # data stream did not start appropriately .. use minilogs
        ot = which.min( abs( difftime( seabird$timestamp, setx$timestamp) ) )
        o =  which( seabird$timestamp >= (seabird$timestamp[ot]- lubridate::minutes(2) ) &
                    seabird$timestamp <= (seabird$timestamp[ot] + lubridate::minutes(10) )  )  # 5 min tow + 5 min in case
        print( "No set data matched for seabird data:")
        print( filename )
        print( head( seabird[ o,] ) )
        print( "The following is the closest matching in set:")
        print( setx)
          if(plotdata){
            setinfo = set[setxi,]
            fp =  file.path( project.datadirectory("bio.snowcrab"), "data", "seabird", "mismatches", yr)
            dir.create( fp, recursive = TRUE, showWarnings = FALSE )
            pdf(file.path(fp,paste(setinfo$trip[1],'pdf',sep='.') ),12,8)
            plot(seabird$timestamp, -seabird$pressure,type='l',main=setinfo$trip[1])
            points(setinfo$timestamp, rep(0,nrow(setinfo)),pch=16,col='red')
            dev.off()
          }
        next()
      }

      zmaxi = which.max( as.numeric( seabird$pressure[o] ) )
      if (length(zmaxi)==0) stop('need depth data use esonar') #zmaxi = which.min( as.numeric( seabird$temperature[o]) ) #should use esonar depth if no depth in seabird dont do this step ever
      if (length(zmaxi)==0) zmaxi = floor( length(o) / 2 )  # take midpoint

      tstamp = seabird$timestamp[o[zmaxi]]
      uid =  paste( "seabird", setx$trip, setx$set, setx$station, lubridate::hour(tstamp), lubridate::minute(tstamp), f, sep=".")
      seabird$seabird_uid[o] = uid
      seabird$depth[o] = decibar2depth ( P=seabird$pressure[o], lat=setx$lat )
      studyid = paste( setx$trip, setx$set, setx$station, sep="." )
      out = data.frame( uid, yr, seabird$timestamp[o[zmaxi]], setx$trip, setx$set, setx$station, studyid, setx$Zx, setx$timestamp, error, filename2, headerall, stringsAsFactors=FALSE )
      names( out ) = c( "seabird_uid", "yr", "timestamp", "trip", "set", "station", "studyid", "setZx", "set_timestamp", "error", "filename", "headerall" )
      metadata = rbind( metadata, out )
    }

    basedata = seabird[ which( !is.na( seabird$seabird_uid) ) ,]

    return( list( metadata=metadata, basedata=basedata ) )
  }


