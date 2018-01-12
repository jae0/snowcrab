
    netmindDate = function( fnNetmind=NULL, header=NULL, outvalue="year", linenumber=20 ) {
      # input can be file name or the file header
      out = NULL
      if (!is.null( fnNetmind) && file.exists(fnNetmind) ) {
        header = readLines( fnNetmind, n=20, encoding="UTF-8", skipNul=TRUE )
        if ( !( (any(grepl("FileName", header, useBytes=TRUE)))
              & (any(grepl("Local", header, useBytes=TRUE)))
              & (any(grepl("Ship", header, useBytes=TRUE)))
              & (length(header) > 15) ) ) {
		      print(paste(fnNetmind,' There is an error in file header.'))
          return( out )
        }
      }

      # local time vs actual in the same file can be wrong: mostly in the year 2000,
      # the Local Time (computer) entry in incorrect vs the GPS date/time (stored in GMT/UTC)
      # eg:  /home/jae/bio/bio.snowcrab/data/netmind/archive/2000/pos065.txt
      # just in case this occurs elsewhere, use the time stamp in the actual data series rather than "Local Time"
      # determine time from "Local Time:"
      # for offsets, assume that at least the time of day component is correct in the local computer time

      if (outvalue=="localtime") {
        lineloc = grep( "Local Time:", header, ignore.case=T  )
        ndt = strsplit( header[lineloc], "[[:space:]]+" )
        if ( length(ndt) != 7) stop ( paste("Local Time error:", header[lineloc], fnNetmind ) )
        dateobject = gsub( "[[:space:]]", "", paste( ndt[7], ndt[4], ndt[5], sep="-" ) )
        timeobject = gsub( "[[:space:]]", "", ndt[6] )
        loctime = lubridate::ymd_hms( paste( dateobject, timeobject), tz="America/Halifax" )
        return (loctime)
      }


      if (outvalue=="linetime") {
        rec = strsplit( header[linenumber], "[[:space:]]+" )
        recdate = paste(substring(rec[1],1,2), substring(rec[1],3,4), substring(rec[1],5,6), sep="-")
        recyr = (substring(rec[2],1,2))
        rectime = paste(recyr, substring(rec[2],3,4), substring(rec[2],5,6), sep=":")
        rects = lubridate::ymd_hms( paste( recdate, rectime), tz="UTC" )
        return( rects )
      }

      if (outvalue=="year") {
        out = netmindDate( header=header, outvalue="date", linenumber=linenumber )
        out = lubridate::year( out  )
        return (out )
      }

      if (outvalue=="date") {
        out = netmindDate( header=header, outvalue="linetime", linenumber=linenumber )
        return (out)
      }

      if (outvalue=="timeoffset") {
        locdate = netmindDate( header=header, outvalue="localtime" )
        recdate = netmindDate( header=header, outvalue="linetime", linenumber=linenumber )
        toffset = difftime(  with_tz( locdate, "UTC"),  with_tz(recdate, "UTC"), units="hours" )
        toffset = round( toffset)
        return (toffset )
      }

      return (out)
    }


