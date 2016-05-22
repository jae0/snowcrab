

    minilogDate = function( fnMini=NULL, header=NULL, outvalue="year" ) {
      # input can be file name or the file header
      out = NULL
      if (!is.null( fnMini) && file.exists(fnMini) ) {
        header = readLines( fnMini, n=8)
        if ( ! any( grepl("Minilog", header) )) return( out )
      }  
      lineformat = grep( "^\\*.Date\\(", header, perl=T )
      lineno = grep( "Start Time=", header, perl=T )
      
      if ( !length(lineformat)==1 | !length(lineno)==1 ) return(out)
      # minilog date format is variable must extract and convert to chron format
      date.format = gsub( "^\\*.Date\\(", "", header[lineformat])
      date.format = gsub( "\\).*$", "", date.format)
      date.format = gsub( "dd", "d", date.format )
      date.format = gsub( "mm", "m", date.format )
      date.format = gsub( "yyyy", "y", date.format )
      date.format = gsub( "yy", "y", date.format )

      Mdate = gsub( "^.*Time=", "", header[lineno])
      Mdate = gsub( "[[:space:],]+.*$", "", Mdate ) # break on space or comma (convention changes in 2000-2001)
      
      out = dates( Mdate, format=date.format)  

      if (outvalue=="year") out = as.numeric( as.character( years( out ) ) ) 
      if (outvalue=="date") out = chron( out, out.format="y-m-d" )
      if (outvalue=="format") out = date.format
      return (out)
    }



