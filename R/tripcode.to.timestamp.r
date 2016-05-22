
  tripcode.to.timestamp = function(x, y, tzone="America/Halifax" ) {
    #require( lubridate) 
    # take a snow crab trip code and time to create a POSIXct object using lubridate
    z = nchar(x[1])
    yr = as.numeric(substring( x, z-3, z))
    mth = as.numeric(substring( x, z-5, z-4))
    day = as.numeric(substring( x, 2, z-6))
    datecode = paste( yr, mth, day, sep="-")
    datetime = ymd_hms( paste( datecode, y) )
    tz( datetime ) = tzone
    return( datetime )
  }


