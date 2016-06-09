
# (Brent)
# 2013 +. A boolean function to test if a file is an eSonar data file.
# If it is a eSonar file then conversion to netmind is needed.

  isEsonar = function( fn=NULL) {
    # input can be file name or the file header
    out = FALSE
    if (!is.null( fn) && file.exists(fn) ) {
      header = readLines( fn, n=5)
      if(grepl("Validity", header[5]))
        if(grepl("CPU Date and Time", header[5]))
          if(grepl("Comments", header[3]))
            if(grepl("CPU Date and Time", header[5]))
              if(grepl("Tow Number", header[1]))
                out = TRUE
    }
    return(out)
  }

