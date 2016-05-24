
  # tests using geoR and geoRglm

  require(geoR)
  require(geoRglm)

	
	loadfunctions( "spacetime", "utility", "parallel" )
  loadfunctions( "snowcrab", functionname="initialise.local.environment.r") 
 
	
	set = snowcrab.db("set")
  ll = set[ which(set$yr==2000), c("plon","plat","R0.mass")]

  rownames(ll) = c(1:nrow(ll))
  ad = dup.coords(ll[,c(1,2)])
  for ( i in 1:ncol(ad) ) {
    dups.i = as.numeric(ad[,i])
    mm = NULL
    mm = colMeans( ll[ dups.i, ], na.rm=T )
    ll = ll[ - dups.i, ]
    ll = rbind( ll, mm )
  }

  a = as.geodata(ll)

  plot(a)
  av = variog(a, estimator.type="modulus", max.dist=150)
  av.env = variog.mc.env(a, obj.variog=av, nsim=1000) 
 
  x11()
  par(mfrow=c(1,2), mar=c(3,3,1.5,0.5))
  plot(av,type="b", envelope=av.env )

  av4 = variog4(a, estimator.type ="modulus", max.dist=150 )
  direction=pi/3
  
  title(main=expression(paste("directional,angle =",60*degree)))
  plot(av4)
  

  a.model = likfit(a, 
       trend = "cte", ini.cov.pars=c(0.5, 0.5), fix.nugget = FALSE, nugget = 0,
       fix.kappa = TRUE, kappa = 0.5, fix.lambda = TRUE, lambda = 1,
       fix.psiA = TRUE, psiA = 0, fix.psiR = TRUE, psiR = 1, 
       cov.model = "matern", method.lik = "REML", components = FALSE,
       nospatial = TRUE, limits = pars.limits(),
       print.pars = FALSE)

  str(a.model)




