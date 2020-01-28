

    minilogDate = function( fnMini=NULL, header=NULL, outvalue="year" ) {
      # input can be file name or the file header
      out = NULL
      if (!is.null( fnMini) && file.exists(fnMini) ) {
        header = readLines( fnMini, n=8)
        if ( ! any( grepl("Minilog", header) )) return( out )
      }
      
      lineformat = grep( "^\\*.Date\\(", header, perl=T )
      lineno = grep( "Start Time=", header, perl=T )
      
      if(length(lineformat) == 0)lineformat = grep( "Date\\(", header, perl=T )
      if(length(lineno) == 0)lineno = grep( "Start Time:", header, perl=T )
      
      if ( !length(lineformat)==1 | !length(lineno)==1 ) return(out)

      # minilog date format is variable must extract and convert to an identifiable format
      date.format = gsub( "^\\*.Date\\(", "", header[lineformat])
      if(any(grepl("Date", date.format))) date.format = gsub( "Date\\(", "", date.format)
      date.format = gsub( "\\).*$", "", date.format)
      date.format = gsub( "dd", "d", date.format )
      date.format = gsub( "mm", "m", date.format )
      date.format = gsub( "yyyy", "y", date.format )
      date.format = gsub( "yy", "y", date.format )

      Mts = gsub( "^.*Time=", "", header[lineno])
      if(any(grepl("Time", Mts))) Mts = trimws(gsub( "Study Start Time:", "", Mts))
      Mdate = gsub( "[[:space:],]+.*$", "", Mts ) # break on space or comma (convention changes in 2000-2001)
      Mtime = gsub( "^.*[[:space:],]+", "", Mts ) # break on space or comma (convention changes in 2000-2001)

      Mtimestamp = paste( Mdate, Mtime )

      out = lubridate::parse_date_time ( Mtimestamp, orders=paste( date.format, "H:M:S" ), tz="America/Halifax"  )
      out = with_tz( out, "UTC" )
      if(months(out) %in% c("January")) out = out-months(1)
      if (outvalue=="year") out = lubridate::year( out )
      if (outvalue=="format") out = date.format
      return (out)
    }



