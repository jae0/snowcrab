	load.minilog.rawdata.one.file.per.day = function(fn, f, set ) {

    tz.minilog = "America/Halifax"
    tz.snowcrab = "UTC"
  
    out = NULL
    minilog=NULL

    filename = fn
    ind_skip_header = 7
    header =  readLines(filename, n=21)
    headerall = paste( header[1:6], collapse="~")

    if (length(header) < 20)  return( NULL )

    
    if(any(grepl("Study ID=", header, perl=T))){
      l.study = grep( "Study ID=", header, perl=T )
      studyid = tolower( gsub("^.*Study ID=", "", header[ l.study ] ) )
    }
    if(any(grepl("Study Description:", header, perl=T))){
      l.study = grep( "Study Description:", header, perl=T )
      studyid = trimws(tolower( gsub("Study Description:", "", header[ l.study ] ) ))
      ind_skip_header = 8 
    }
    
    if (grepl( "test", studyid, ignore.case=T) ) return( NULL )
    if (grepl( "testage", studyid, ignore.case=T) ) return( NULL )

    minilog = as.data.frame(read.table( file=filename, sep=",", as.is=T, colClasses="character", header=F, skip=ind_skip_header))

    if ( nrow(minilog) < 10 ) return( NULL )

	  fileinfo = tolower(unlist(strsplit(filename, split="/")))
    tripid =   tolower(unlist(strsplit(basename(filename), "\\."))[1])

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
    minilog$timestamp = lubridate::parse_date_time ( paste(minilog$mdate, minilog$mtime), orders=paste(headerdateformat, "H:M:S" ), tz=tz.minilog )
    minilog$timestamp = with_tz(  minilog$timestamp, tz.snowcrab )  # set/setInitial is now in UTC

    yr = lubridate::year( minilog$timestamp[1])
    if (!is.finite(yr) ) yr = minilogDate( header=header, outvalue="year"  )
    print(filename2)
    #break up minilog by station
    setS = set[which(tolower(set$trip)==tripid),]
    # minilog$timestamp1 = trunc(minilog$timestamp, 'mins')

    metadata = NULL
    basedata = NULL

    for(pp in 1:nrow(setS)){
        xS = setS[pp,]
        xSii = which( minilog$timestamp >= (xS$timestamp - lubridate::minutes(10)) &
                      minilog$timestamp <= (xS$timestamp + lubridate::minutes(20)))
        if (length(xSii) < 30 ) next()
        mi = minilog[ xSii , ] #three minutes before the start of the tow and 15 min after from setinfo
        mi$minilog_uid = xS$minilog_uid = paste('minilog',tripid,xS$station,xS$set,sep=".")
        meta = data.frame(
          minilog_uid=xS$minilog_uid,
          yr=xS$yr,
          timestamp=xS$timestamp,
          trip=xS$trip,
          set=xS$set,
          station=xS$station,
          studyid=sprintf("ep%03d",xS$station),
          setZx=xS$Zx,
          set_timestamp=xS$timestamp,
          error= NA,
          filename=filename2,
          headerall=headerall,
          stringsAsFactors=FALSE)
        base = data.frame(
          mdate=mi$mdate,
          mtime=mi$mtime,
          temperature=mi$temperature,
          depth=mi$depth,
          timestamp=mi$timestamp,
          minilog_uid = mi$minilog_uid,
          stringsAsFactors=FALSE)
        metadata = rbind(metadata, meta )
        basedata = rbind(basedata, base )
    }

    return( list( metadata=metadata, basedata=basedata ) )
  }


