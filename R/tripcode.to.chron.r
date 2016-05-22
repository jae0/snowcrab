  
  tripcode.to.chron = function(x, y) {
    require( chron) 
    # take a snowcrab trip code and time to create a chron object 
    z = nchar(x[1])
    yr = as.numeric(substring( x, z-3, z))
    mth = as.numeric(substring( x, z-5, z-4))
    day = as.numeric(substring( x, 2, z-6))
    datecode = paste( yr, mth, day, sep="-")
    out = chron( dates.=datecode, times.=y, format=c(dates="y-m-d", times="h:m:s"), out.format=dateformat.snow )
    return(out)
  }


