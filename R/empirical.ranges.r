
  empirical.ranges = function( db="snowcrab", var, remove.zeros=F, probs=c(0.025, 0.975) ) {
    if (db =="snowcrab") {
      x = snowcrab.db("set.complete")[, var]  
      if (remove.zeros) x = x[ which( x>0 ) ]
      return ( quantile( x, probs=probs, na.rm=T ) )
    }    
  }


