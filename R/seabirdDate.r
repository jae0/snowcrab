

    seabirdDate = function( fnSeaBird=NULL, header=NULL, outvalue="year" ) {
      # read dates from file name or the file header and then return the dates
      # or the date format for seabirds
      out = NULL
      if (!is.null( fnSeaBird) && file.exists(fnSeaBird) ) {
        header = readLines( fnSeaBird, n=52)
      }

      if ( any( grepl("Sea-Bird", header, useBytes=TRUE ) )) {
        lineno = grep( "^start\\ time\\ =", header, perl=TRUE, useBytes=TRUE  )
        if ( length(lineno)==1 ) {
          Mdate = gsub( "^start\\ time\\ =", "", header[lineno])
          Mdate = gsub( "^[[:space:]]{1,}", "", Mdate) # empty space at beginning
          Mdate = gsub( "[[:space:]]{1,}", " ", Mdate) # multiple spaces into one space
          y =  matrix(unlist(strsplit(Mdate, " ")), nco=4, byrow=T)
          dstring = paste(y[1], y[2], y[3], sep=" ")
          tstring = y[4]
          out = lubridate::parse_date_time( paste(dstring, tstring ), orders="d b y H:M:S" ) # b = mon, in POSIXct
          if (outvalue=="year"){
            if(months(out) %in% c("January")) out = out-months(1)
             return( lubridate::year( out ) )
          }
          if (outvalue=="posix") return(out)
        }
      }
      return (out)
    }


