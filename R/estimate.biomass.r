
  estimate.biomass = function(Q, R, S, y, r, use.global.average=F) {
  
    Y = as.character(y)
      
#    S[,"mass.mean",,] = S[,"mass.mean",, ] 
#    S[,"mass.sd",,]   = S[,"mass.sd",, ]
    
    m.q = S[,"mass.mean",Y,r]
    error.q = S[,"mass.sd",Y,r]
    j.q = which ( S[,"n",Y,r] <= 30 )
    m.q[ j.q ] = NA
    error.q[ j.q ] = NA

    m.r = S[,"mass.mean",Y,"cfaall"]
    error.r = S[,"mass.sd",Y,"cfaall"]

#    j.r = which ( S[,"n",Y,"cfaall"] <= 30 )
#    m.r[ j.r ] = NA
#    error.r[ j.r ] = NA

    m.s =  apply( S[,"mass.mean",,"cfaall"], MARGIN=1, mean, na.rm=T)
    error.s = apply( S[,"mass.sd",,"cfaall"], MARGIN=1, mean, na.rm=T)

    m.o = m.q
    error.o = error.q
   
    if (use.global.average) {
      m.o = m.r
      error.o = error.r
   } else {
      i = which (!is.finite(m.o) )
      m.o[ i ] = m.r[ i ]
      error.o[ i ] = error.r[ i ]
     
      j = which (!is.finite(m.o) )
      m.o[ j ] = m.s[ j ]
      error.o[ j ] = error.s[ j ]
   }

    if ( !is.finite( m.o["CC5.13"] ) ) {
      m.o[ "CC5.13" ] =  m.o ["CC3to4.13"]
      error.o[ "CC5.13" ] =  error.o ["CC3to4.13"]
    }

    o = m.o * Q
    oerror =  o * sqrt( (error.o/m.o) ^2 + (R/Q)^2)
      
    names(o) = names(Q)
    o[ o < 0.0001] = 0
    o[ !is.finite(o) ] = 0

    names(oerror) = names(Q)
    oerror [ oerror < 0.0001] = 0
    oerror [ !is.finite(oerror) ] = 0

    out = list( x=o / 1000 , error=oerror / 1000 ) # convert g to kg
    
    return(out)
  }



