
    netmindDate = function( fnNetmind=NULL, header=NULL, outvalue="year", linenumber=20 ) {
      # input can be file name or the file header
      out = NULL
      if (!is.null( fnNetmind) && file.exists(fnNetmind) ) {
        header = readLines( fnNetmind, n=20 )
        if ( !( (any(grepl("FileName", header))) & (any(grepl("Local", header))) & (any(grepl("Ship", header))) & (length(header) > 15) ) ) {
		print(paste(fnNetmind,' error in file header;'))
          return( out )
        }
      }

      # local time vs actual in the same file can be wrong: mostly in the year 2000,
      # the Local Time (computer) entry in incorrect vs the GPS date/time (stored in GMT/UTC)
      # eg:  /home/jae/ecomod/snowcrab/data/netmind/archive/2000/pos065.txt
      # just in case this occurs elsewhere, use the time stamp in the actual data series rather than "Local Time"
      # determine time from "Local Time:"
      # for offsets, assume that at least the time of day component is correct in the local computer time

      if (outvalue=="localtime") {
        lineloc = grep( "Local Time:", header, ignore.case=T  )
        ndt = stringsplit( header[lineloc], "[[:space:]]+" )
        if ( length(ndt) != 7) stop ( paste("Local Time error:", header[lineloc], fnNetmind ) )
        dateobject = dates( gsub( "[[:space:]]", "", paste( ndt[7], ndt[4], ndt[5], sep="-" ) ), format="y-mon-d")
        timeobject = times( gsub( "[[:space:]]", "", ndt[6] ) )
        loctime = chron( dateobject, timeobject, out.format=c("y-m-d", "hh:mm:ss") )
        return (loctime)
      }

      if (outvalue=="linetime") {
        rec = stringsplit( header[linenumber], "[[:space:]]+" )
        recdate = paste(substring(rec[1],1,2), substring(rec[1],3,4), substring(rec[1],5,6), sep="-")
        recyr = (substring(rec[2],1,2))
        rectime = paste(recyr, substring(rec[2],3,4), substring(rec[2],5,6), sep=":")
        recchron = chron( dates(recdate, format="y-m-d" ), times(rectime), out.format=c("y-m-d", "hh:mm:ss") )
        return( recchron )
      }

      if (outvalue=="year") {
        out = netmindDate( header=header, outvalue="date", linenumber=linenumber )
        out = as.numeric( as.character( years( out  ) ) )
        return (out )
      }

      if (outvalue=="date") {
        out = netmindDate( header=header, outvalue="linetime", linenumber=linenumber )
        return (out)
      }

      if (outvalue=="timeoffset") {
        locdate = netmindDate( header=header, outvalue="localtime" )
        recdate = netmindDate( header=header, outvalue="linetime", linenumber=linenumber )
        toffset = as.numeric(locdate - recdate)  # in "days" .. keep only up to hour differences
        toffset = round( toffset *24) / 24 
        return (toffset )
      }

      return (out)
    }


