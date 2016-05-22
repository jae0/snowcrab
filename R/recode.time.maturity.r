
recode.time.maturity = function( x, time.resolution="daily" ) {
  if (time.resolution == "daily" ) {
    zz = 1
  } else if( time.resolution == "weekly" ) {
    zz = 52/365
  } else if( time.resolution == "monthly" ) {
    zz = 12/365
  }
  x = round( x * zz )
  return(x)
}



