
  interpolate.missing.xyz = function( xyz, method="default",  nmax=30, maxdist=50, idp=0.5 ) {
    
    out = NULL
    names( xyz) =  c("x", "y", "z")

    missing = which( !is.finite( xyz$z ) )
    if (any( (!is.finite( xyz$x)) | (!is.finite( xyz$y)) )) {
      print( "Missing values in coords not allowed" )  
      return (NULL)
    }

    if ( length( missing ) ==0 ) {
      print( "Nothing to interpolate" )
      return (NULL)
    }

    if ( (nrow( xyz ) - length(missing)) < 30 ) {
      print ("Not enough data to interpolate")  
      return ( NULL  )
    }

    data.range = range( xyz$z, na.rm =T ) 

    if ( method %in% c("default", "linear.akima" ) ) { # interpolation is simple linear
      require( akima)
      gs = try( interp( x=xyz$x[-missing], y=xyz$y[-missing], z=xyz$z[-missing], xo=xyz$x[missing], yo=xyz$y[missing], linear=T, duplicate="mean" ) )
      if (class(gs) =="try-error" )  return (NULL )
      k =  as.matrix(gs[[3]])
      dimnames(k) = list( gs[[1]], gs[[2]] ) 
      pred = as.data.frame.table(k)
      names( pred ) =  names( xyz)
      pred = pred [ which(is.finite( pred$z ) ), ]
      xyz$order = c( 1:nrow(xyz) )
      out = merge( xyz, pred, by=c( "x", "y" ), all.x=T, all.y=F, sort=F, suffixes=c("",".new") )
      if (nrow(xyz) != nrow(out) ) {
        print( "Merge failure at interpolation" )
        stop()  # error
      }
      out$z[missing] = out$z.new[missing]
      out = out[ order( out$order) , c("x", "y", "z" ) ]
      out = range.cap( out, data.range )
      return (out )
    }
    
    if ( method =="bicubic.spline.akima" ) { 
      require( akima)
      gs = try( interp( x=xyz$x[-missing], y=xyz$y[-missing], z=xyz$z[-missing], xo=xyz$x[missing], yo=xyz$y[missing], linear=F, extrap=T, duplicate="mean" ) )
      if (class(gs) =="try-error" )  return (NULL )
      k =  as.matrix(gs[[3]])
      dimnames(k) = list( gs[[1]], gs[[2]] ) 
      pred = as.data.frame.table(k)
      names( pred ) =  names( xyz)
      pred = pred [ which(is.finite( pred$z ) ), ]
      xyz$order = c( 1:nrow(xyz) )
      xyz = merge( xyz, pred, by=c( "x", "y" ), all.x=T, all.y=F, sort=F, suffixes=c("",".new") )
      xyz$z[missing] = xyz$z.new[missing]
      xyz = xyz[ order( xyz$order) , c("x", "y", "z" ) ]
      out = range.cap( xyz, data.range )
      return (out )
    }

    if ( method =="inverse.distance.gstat" ) {
      require( gstat )
      m = gstat( id="interpolation", formula=z~1, locations=~x+y, data=xyz[-missing,], nmax=nmax, maxdist=maxdist, set=list(idp=idp)) 
      pred = predict( object=m, newdata=xyz[missing,] ) [,3]
      if (class(pred) =="try-error" )  return (NULL )
      pred = range.cap( pred, data.range )
      xyz$z[missing] = pred 
      return (xyz$z )
    }

    return ( "Method not found" )

  }
  

