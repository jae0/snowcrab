
  error.propagate = function(x=NULL, error.x=NULL, y=NULL, error.y=NULL, type="") {
    if (type=="product") { # expects y to be the transition matrix
      error.x.sq = (error.x/x)^2
      error.x.sq[ !is.finite(error.x.sq) ] = 0
      error.y.sq =  (error.y/y)^2
      error.y.sq[ !is.finite(error.y.sq) ] = 0
   
      if (is.array(y) | is.matrix(y) ) { 
        error.z = sqrt( as.vector(error.y.sq %*% error.x.sq) )
      } else {
        error.z = sqrt( error.x.sq + error.y.sq )
      }
      
    } else if (type=="sum") {
      error.x.sq = error.x^2
      error.x.sq[ !is.finite(error.x.sq) ] = 0
      error.y.sq = error.y^2
      error.y.sq[ !is.finite(error.y.sq) ] = 0
      
      if (is.array(y) | is.matrix(y) ) { 
        error.z = sqrt( as.vector( error.x.sq %*% t(error.y.sq) ))
      } else {
        error.z = sqrt( error.x.sq + error.y.sq )
      }
       
    } else if ( type=="power" ) {
      error.z = abs(y) * error.x/x
      error.z[ !is.finite(error.z) ] = 0
    }
    return(error.z)
  }



