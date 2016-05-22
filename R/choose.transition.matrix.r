
choose.transition.matrix = function(tm, region, threshold.to.delete=5, use.global.average=T) {
    
    eps = 1e-2
    trim = 0

      tm[ !is.finite(tm) ] = NA
      tm[ tm > threshold.to.delete ] = NA  # large estimates in transfer function is not believable
      tm[ tm < eps ] = NA  # very small estimates in transfer function are temporarily not trusted
   
      xm = tm[,,,region]
       
      # region specific values
      XMm = apply(xm, c(1,2), function(b) {  mean( b, na.rm=T)  } )
      XMmn = apply(xm, c(1,2), function(b) { length( b[is.finite(b)] ) } )
      XMmsd = apply(xm, c(1,2), function(b) {  sd (b, na.rm=T) } ) 
      
      XM = XMm
      XMsd = XMmsd
     
      # shelf-wide averages
      XMg = apply(tm, c(1,2), function(b) { mean( b, na.rm=T) } )
      XMng = apply(tm, c(1,2), function(b) { length( b[is.finite(b)] ) } )
      XMsdg = apply(tm, c(1,2), function(b) { sd( as.vector(b), na.rm=T)  } )  
      
      if (use.global.average)  {
        XM = XMg 
        XM[ !is.finite( XM ) ] = 0
        XMsd = XMsdg 
        XMsd[ !is.finite( XMsd ) ] = 0
      } else {
        iii = which( !is.finite( XM )| (XM <= eps) ) 
        XM[ iii ] = XMg [iii]  # try using shelf-wide means when no recent data exists
        XM[ !is.finite( XM ) ] = 0
     
        jjj = which( !is.finite( XMsd )| (XM < eps) ) 
        XMsd[ iii ] = XMsdg [iii]  # try using shelf-wide se when no recent data exists
        XMsd[ !is.finite( XMsd ) ] = 0
       }

  return (list( XM=XM, XMsd=XMsd) )
  
  }


