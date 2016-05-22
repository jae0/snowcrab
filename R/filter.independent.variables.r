
  filter.independent.variables = function( x, vtype="R0.mass")  {
    if (vtype %in% c("R0.mass", "R1.no") ) {
     
      if (exists( "t", x)) {
        ii =  which( x$t > 14  )
        if (length(ii) > 0)  x$t[ii] = 14
      }

      if (exists( "tmean", x)) {
        ii =  which( x$tmean > 10  )
        if (length(ii) > 0)  x$tmean[ ii ] = 10
      }

      if (exists( "z", x)) {
        ii =  which( x$z > log(600) ) 
        if (length(ii) > 0)  x$z[ii] = log(600)
        ii =  which( x$z < log(25) ) 
        if (length(ii) > 0) x$z[ii] = log(25)
      }
      
      if (exists( "dZ", x)) {
        ii =  which(  x$dZ < -8 ) 
        if (length(ii) > 0) x$dZ[ii] = -8
        ii =  which(  x$dZ > -2 ) 
        if (length(ii) > 0) x$dZ[ii] = -2
      }
      
      if (exists( "tamp", x)) {
        ii =  which( x$tamp> 12   ) 
        if (length(ii) > 0) x$tamp[ ii ] = 12
      }
      
      if (exists( "dt.annual", x)) {
        ii =  which( x$dt.annual < -3  ) 
        if (length(ii) > 0) x$dt.annual[ ii ] = -3
        ii =  which( x$dt.annual > 3  ) 
        if (length(ii) > 0) x$dt.annual[ ii ] = 3
      }
      
      if (exists( "dt.seasonal", x)) {
        ii =  which( x$dt.seasonal < -3  ) 
        if (length(ii) > 0) x$dt.seasonal[ii ] = -3
        ii =  which( x$dt.seasonal > 3   ) 
        if (length(ii) > 0) x$dt.seasonal[ ii ] = 3
      }
      
      if (exists( "wmin", x)) {
        ii =  which( x$wmin< -7   ) 
        if (length(ii) > 0) x$wmin[ ii ] = -7
      }
      
      if (exists( "smr", x)) {
        ii =  which( x$smr < 0.0045  ) 
        if (length(ii) > 0) x$smr[ ii ] = 0.0045 
        ii =  which( x$smr > 0.0065  ) 
        if (length(ii) > 0) x$smr[ii ] = 0.0065
      }
      
      if (exists( "mr", x)) {
        ii =  which( x$mr < 0  ) 
        if (length(ii) > 0) x$mr[ ii ] = 0
        ii =  which( x$mr > 50  ) 
        if (length(ii) > 0) x$mr[ ii ] = 50
      }

      if (exists( "ca1", x)) {
        ii =  which( x$ca1 < -2.5  ) 
        if (length(ii) > 0) x$ca1[ ii] = -2.5
        ii =  which( x$ca1 >  2  ) 
        if (length(ii) > 0) x$ca1[ ii ] = 2
      }

      if (exists( "ca2", x)) {
        ii =  which( x$ca2 < -2.0  ) 
        if (length(ii) > 0) x$ca2[ ii ] = -2.0
        ii =  which( x$ca2 >  2.0  ) 
        if (length(ii) > 0) x$ca2[ ii ] =  2.0
      }
      
      if (exists( "Npred", x)) {
        ii =  which( x$Npred > 130  ) 
        if (length(ii) > 0) x$Npred[ ii] = 130
      }
      
      if (exists( "Z", x)) {
        ii =  which(  x$Z > 0.41  ) 
        if (length(ii) > 0) x$Z[ ii ] = 0.41
      }
  
    }
    return (x)
  }


