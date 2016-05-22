	load.minilog.rawdata.one.file.per.day = function(fn, f, set ) {
   	
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
    date.format = c( dates=headerdateformat, times="h:m:s")
    minilog$chron = chron( dates.=minilog$mdate, times.=minilog$mtime, format=date.format, out.format=dateformat.snow )
    
    yr = as.numeric( as.character( years(minilog$chron[1]) ))
    if (!is.finite(yr) ) yr = minilogDate( header=header, outvalue="year"  ) 
    print(filename2)
    #break up minilog by station
    setS = set[which(tolower(set$trip)==tripid),]
    minilog$chron1 = trunc(minilog$chron,'minutes')
    metadata = NULL
    basedata = NULL
    for(pp in 1:nrow(setS)){
        xS = setS[pp,]
        mi = minilog[which(minilog$chron1>=xS$chron - (minute*10) & minilog$chron1<=xS$chron + (minute*20) ),] #three minutes before the start of the tow and 15 min after from setinfo
        mi$minilog_uid = xS$minilog_uid = paste('minilog',tripid,xS$station,xS$set,sep=".")
        xS = data.frame(minilog_uid = xS$minilog_uid,yr=xS$yr,timestamp = xS$chron,trip=xS$trip,set=xS$set,station=xS$station,studyid=sprintf("ep%03d",xS$station),setZx=xS$Zx,setChron=xS$chron,error= NA,filename=filename2,headerall=headerall,stringsAsFactors=FALSE)
        mi = data.frame(mdate=mi$mdate,mtime=mi$mtime,temperature=mi$temperature,depth=mi$depth,chron=mi$chron,minilog_uid = mi$minilog_uid,stringsAsFactors=FALSE)
        metadata = rbind(metadata,xS )
        basedata = rbind(basedata,mi)       
    }

    return( list( metadata=metadata, basedata=basedata ) )
  }
 

