	load.minilog.rawdata = function(fn, f, set ) {
   	
    out = NULL
    minilog=NULL
    
    filename = fn

    header =  readLines(filename, n=21)
    headerall = paste( header[1:6], collapse="~")

    if (length(header) < 20)  return( NULL )
      
    l.study = grep( "Study ID=", header, perl=T )
    studyid = tolower( gsub("^.*Study ID=", "", header[ l.study ] ) )
    if (grepl( "test", studyid, ignore.case=T) ) return( NULL )
    if (grepl( "testage", studyid, ignore.case=T) ) return( NULL )
    

    minute = 1 / 24 / 60 # 1 minute in chron terms

    minilog = as.data.frame(read.table( file=filename, sep=",", as.is=T, colClasses="character", header=F, skip=7))
   
    if ( nrow(minilog) < 10 ) return( NULL )
   
	  fileinfo = tolower(unlist(strsplit(filename, split="/")))
	 
    
    stationid = tolower(unlist(strsplit(basename(filename), "\\."))[1])
      if (grepl("^ep", stationid, ignore.case=T)) {
        stationid = gsub( "(^ep[[:space:]]*)([[:digit:]]*)([[:space:]]*.*$)", "\\2", stationid, ignore.case=T )
      } else if (grepl("asc", stationid, ignore.case=T)) {
        # older data series used different file naming conventions
        # sequence is slightly important: e.g., "post" must come before "pos"
        possible.prefixes = c( "postr", "posr", "pos t", "post", "pos", "ep", "station" )
        stationid = studyid
        for (pp in possible.prefixes ) stationid = gsub( pp, "QQQ", stationid, ignore.case=T )
        if (grepl("QQQ", stationid )) {
          stationid = gsub( "(^.*QQQ[[:space:]]*)([[:digit:]]*)([[:space:]]*.*$)", "\\2", stationid, ignore.case=T )
        }
        # manual fixes  .. still in the asc* series
        if ( gsub("[[:space:]]", "", studyid) == gsub("[[:space:]]", "", "ens 2001 zone 20 167") ) stationid = "167"
        if ( gsub("[[:space:]]", "", studyid) == gsub("[[:space:]]", "", "ens 2001 zone 20 162") ) stationid = "162"
      }
      
    
      if (nchar(stationid)>=1) {
        stationid = as.numeric( gsub("[[:alpha:]]", "", stationid) ) 
      } else {
        stationid = NA
      }


	  filename2 = fileinfo[length(fileinfo)] #Changed from postions in file location to the last entry since last entry is the filename
    filename2 = tolower(filename2)
	  
	  if (ncol(minilog)== 3)  minilog[,4] = NA  # fill with NA's when depth is not recorded (e.g. 1998)

    colnames(minilog) = c( "mdate", "mtime", "temperature", "depth")
    numerics = c("temperature", "depth")
    minilog = factor2number(minilog, numerics)

    # depth offets as this can be large sometimes (Esp 2009 ~ 50 m)
    surface = quantile( minilog$depth, probs=0.01, na.rm=TRUE )
    if ( is.finite( surface) ) minilog$depth = minilog$depth - surface
    
    # obtain date format from the minilog header
	  headerdateformat = minilogDate( header=header, outvalue="format"  ) 
    if (is.null(headerdateformat) ) return( NULL )
    
    minilog$mdate = gsub("^80-", "2000-", minilog$mdate )  # there are incorrect year entries
    date.format = c( dates=headerdateformat, times="h:m:s")
    minilog$chron = chron( dates.=minilog$mdate, times.=minilog$mtime, format=date.format, out.format=dateformat.snow )
    
    yr = as.numeric( as.character( years(minilog$chron[1]) ))
    if (!is.finite(yr) ) yr = minilogDate( header=header, outvalue="year"  ) 
    

    # check time first
    minilog.date.range = range( minilog$chron )
    setxi = which( set$chron >= minilog.date.range[1] & set$chron <= minilog.date.range[2] )
    if ( length(setxi ) ==0 ) return (NULL)
  
    minilog$minilog_uid = NA 
    metadata = NULL
    for ( ssid in 1:length(setxi) ) {
    
      iset = setxi[ssid]
      setx =  set[iset,] # matching trip/set/station  

      # reduce size of minilog data stream
      istart = which.min( abs(minilog$chron - setx$chron ))
      
      j0 = istart
      while ( j0 > 1 ) {
        if (minilog$depth[j0] <= surface ) break() 
        j0 = j0 - 1
      }
      tdiff0 = abs(minilog$chron[j0] - setx$chron )
      if ( tdiff0 > (10 * minute) ) { 
        # then something went wrong .. default to a few minutes before set$chron
        j0 = which.min( abs( minilog$chron - (setx$chron - 5 * minute) )) 
      }

    
      j1 = istart
      while ( j1 < nrow(minilog) ) {
        if (minilog$depth[j1] <= surface ) break() 
        j1 = j1 + 1
      }
      tdiff1 = abs(minilog$chron[j1] - setx$chron )
      if ( tdiff1 > (15 * minute) ) { 
        # then something went wrong .. default to a few minutes before set$chron
        j1 = which.min( abs( minilog$chron - (setx$chron + 12 * minute) )) 
      }

      o = j0:j1

      if (length(o) < 30) {
        # no matching time range found ... use closest match and raise a flag
        # data stream did not start appropriately .. use minilogs
        ot = which.min( abs( minilog$chron - setx$chron) )
        o =  which( minilog$chron >= (minilog$chron[ot]- 2*minute) & 
                 minilog$chron <= (minilog$chron[ot] + 10*minute)  )  # 5 min tow + 5 min in case
        print( "No matching set entry found for:")
        print( filename )
        print( head( minilog[ o,] ) )
        return( NULL )
      }

      error = ""
      if (any(is.na(minilog$chron))) error = paste(error, "Ambiguous time format" )
      if (years(minilog$chron[1]) != yr)  error = paste(error, "Years inconsistent" )
      strangedatacheck = try( lm(temperature ~ depth, data=minilog,  na.action="na.omit"), silent=T )
      if ( "try-error" %in% class( strangedatacheck ) ) error=paste(error, "no depth data?")
      
      # error corrections:  This one is hard to fix without changing the raw data file
      if  (yr == 2006) {
        if ( setx$station==124 ) {
          d0 = chron( dates.="2006-10-11", times.="16:37:16", format=c(dates="y-m-d", times="h:m:s"))
          d1 = chron( dates.="2006-10-16", times.="05:40:50", format=c(dates="y-m-d", times="h:m:s"))
          offset = times(d1) - times(d0)
          minilog$chron = chron( times( minilog$chron ) + offset, out.format=dateformat.snow )
      }}
     
      
      zmaxi = which.max( as.numeric( minilog$depth) )
      if (length(zmaxi)==0) zmaxi = which.min( as.numeric( minilog$temperature) )
      if (length(zmaxi)==0) zmaxi = floor( nrow(minilog) / 2 )  # take midpoint
      if ( !(length(zmaxi)==1) ) stop( filename )
      tstamp = minilog$chron[o[zmaxi]]
      minilog_uid = paste( "minilog",  setx$trip, setx$set, setx$station, hours(tstamp), minutes(tstamp), f, sep=".")
      minilog$minilog_uid[o] = minilog_uid

      out = data.frame( minilog_uid, yr, minilog$chron[zmaxi], setx$trip, setx$set, setx$station, studyid, setx$Zx, setx$chron, error, filename2, headerall, stringsAsFactors=FALSE)

      names( out ) = c( "minilog_uid", "yr", "timestamp", "trip", "set", "station", "studyid", "setZx", "setChron",  "error", "filename", "headerall" )
      
      metadata = rbind( metadata, out )

    }
    
    basedata = minilog[ which( !is.na( minilog$minilog_uid) ) ,]

    return( list( metadata=out, basedata=basedata ) )
  }
 

