 
  # surplus production analysis
  # Maxmimum Likelihood optimization approach (aka, "timeseries modelling") 

 	
  loadfunctions( "snowcrab", functionname="initialise.local.environment.r") 
  
   
  ###  all data follow this sequence: c("cfanorth", "cfasouth", "cfa4x")
  redo.data=F
  if (redo.data) { 
    biomass.summary.db("redo", p=p )
  }
  
  res = biomass.summary.db(p=p)
  
  sc = res$BI[ , which(colnames(res$BI)=="cfasouth")]
  cc = res$L[,  which(colnames(res$L)=="cfasouth")]
  yrs = res$BI[ , which(colnames(res$BI)=="yr")]
  params = c(r=1, K=60, q=1 , B0=40 )
 
  # lognormal errors are used:   
  res = optim( params, fn=biomass.logistic, O=sc, C=cc, errorType="lognormal" ) 

  # to obtain accessory estimates using solution parameters obtained above
  o = biomass.logistic.recursion( res$par, O=sc, C=cc, returnErrorOnly=F, TAC=rep(8, 3), yrs=yrs ) 
   lognormal timeseries method
  #

  res0 = optim( params, fn=biomass.logistic.recursion, O=sc, C=cc, errorType="lognormal" ) 
  res1 = nlminb( start=params, objective=biomass.logistic.recursion, O=sc, C=cc,
      errorType="lognormal" , control=list(eval.max=500, iter.max=500 ))
  res2 = nlm( p=params, f=biomass.logistic.recursion, O=sc, C=cc, errorType="lognormal"  )
  
  # Schnute's method
  res3 = biomass.logistic.schnute.noequilibrium( O=sc, C=cc ) 
   

  res = res0$par
  # eg --
  # normal errors
  # res = optim( params, fn=biomass.logistic, O=sc, C=cc, errorType="normal" ) 
  #          r            K            q           B0 
  #   2.386061 26398.577857     0.223553  7000.074956 
  
  # in this run, lognormal errors are used:   

  # to obtain accessory estimates using solution parameters obtained above
  o = biomass.logistic( r=res["r"], K=res["K"], q=res["q"], B0=res["B0"],
    O=sc, C=cc, TAC=rep(8, 5), yrs=yrs ) 
  
  # a few useful plots
  ymax = max( c( o$res$O, o$res$Op ), na.rm=T )
  plot( Op ~ yrs, o$res, type="b", ylim=c(0, ymax), pch=20 )  # predicted catch rate
  points( O ~ yrs, o$res, pch=22, col="red" )  # observed catch rate
  lines( O ~ yrs, o$res, lty="dashed", col="red" )  # observed catch rate
  abline ( v=yrs[ length(yrs)]+0.5, lty="dotted" )
  

  # catch
  plot( C ~ yrs, o$res, type="b" )
  abline ( v=yrs[ length(yrs)]+0.5, lty="dotted" )

  # biomass
  plot( B ~ yrs, o$res, type="b" ) 
  abline ( v=yrs[ length(yrs)]+0.5, lty="dotted" )


  # fishing mortality
  plot ( F ~ yrs, o$res, type="b" )
  abline ( v=yrs[ length(yrs)]+0.5, lty="dotted" )
  
  
  # surplus production 
  plot( P ~ yrs, o$res, type="b", pch=20 )
  abline ( v=yrs[ length(yrs)]+0.5, lty="dotted" )
  abline ( h=0, lty="dotted" )

  # surplus production 
  plot( P ~ B, o$res, type="p", pch=20 )
  tmp = na.omit( o$res[, c("P","B")] )
  tmp = tmp[ order( tmp$B ), ]
  lines( predict( loess(  P ~ B, tmp, span=0.6 ) ) ~ B, tmp )
  abline ( h=0, lty="dotted" )

  # surplus production 
  plot( P ~ C, o$res )
  tmp = na.omit( o$res[, c("P","C")] )
  tmp = tmp[ order( tmp$C ), ]
  lines( predict( loess(  P ~ C, tmp, span=0.6 ) ) ~ C, tmp )
  abline ( h=0, lty="dotted" )




  # a few useful plots
  ymax = max( c( o$res$O, o$res$Op ), na.rm=T )
  plot( Op ~ yrs, o$res, type="b", ylim=c(0, ymax), pch=20 )  # predicted catch rate
  points( O ~ yrs, o$res, pch=22, col="red" )  # observed catch rate
  lines( O ~ yrs, o$res, lty="dashed", col="red" )  # observed catch rate
  abline ( v=sc$yr[ length(sc$yr)]+0.5, lty="dotted" )
  
  plot( C ~ yrs, o$res, type="b" )
  abline ( v=sc$yr[ length(sc$yr)]+0.5, lty="dotted" )

  plot( B ~ yrs, o$res, type="b" ) 
  abline ( v=sc$yr[ length(sc$yr)]+0.5, lty="dotted" )

  plot ( F ~ yrs, o$res, type="b" )
  abline ( v=sc$yr[ length(sc$yr)]+0.5, lty="dotted" )

  plot( P ~ yrs, o$res, type="b", pch=20 )
  abline ( v=sc$yr[ length(sc$yr)]+0.5, lty="dotted" )
  abline ( h=0, lty="dotted" )

  plot( P ~ B, o$res, type="p", pch=20 )
  tmp = na.omit( o$res[, c("P","B")] )
  tmp = tmp[ order( tmp$B ), ]
  lines( predict( loess(  P ~ B, tmp) ) ~ B, tmp )
  abline ( h=0, lty="dotted" )

  plot( P ~ C, o$res )
  tmp = na.omit( o$res[, c("P","C")] )
  tmp = tmp[ order( tmp$C ), ]
  lines( predict( loess(  P ~ C, tmp ) ) ~ C, tmp )
  abline ( h=0, lty="dotted" )



