
  	
	loadfunctions( "snowcrab", functionname="initialise.local.environment.r") 
 
  require(sp)
  
  
  o = snowcrab.db( DS ="set.complete", p=p )   
  
  p = list()
  p$interval = 5
  p$drange=120
  p$nmax=0
  p$variogram.cressie=F
  p$model="Sph"
  p$vgm.dist = c( seq( 10, p$drange, p$interval) )

  v = "totmass.male.com"

  ow = o[ , c(v, "plon", "plat", "z", "t", "yr" ) ]
  names (ow)[ which( names(ow)==v ) ] = "Y"
 


  require(MASS)
  nb.start = list(size=mean(ow$Y)^2/(var(ow$Y) - mean(ow$Y) ), mu=mean(ow$Y) )
  nb.start = list(size=1, mu=1.5 )
  of = fitdistr( ceiling(ow$Y), "Negative Binomial", start=nb.start )
  
  # mean to variance
  Ym = tapply( ceiling(ow$Y), ow$yr, fitdistr, 
      densfun="Negative Binomial", start=list( mu=of$estimate[["mu"]], size=of$estimate[["size"]] ) )
  Ymm = t( matrix(  unlist( Ym ), ncol=length(Ym ) ))
  colnames(Ymm) = names( unlist(Ym[[1]]) )
  rownames(Ymm) = names(Ym)
  Ymm = as.data.frame(Ymm)
  Ymm$prob = Ymm$estimate.size / (Ymm$estimate.size + Ymm$estimate.mu)
  Ymm$variance =  Ymm$estimate.mu +  Ymm$estimate.mu^2/Ymm$estimate.size

  plot( estimate.size ~ estimate.mu, Ymm, type="n" )
  text( Ymm$estimate.mu , Ymm$estimate.size, rownames(Ymm) )




  set.seed(123)
  ow$plon = jitter( ow$plon)
  ow$plat = jitter( ow$plat)

  omin = min( ow$Y[ ow$Y>0] )/10
  ow$Y[ ow$Y<=0] = omin
  hist(ow$Y,"fd")
  hist(log(ow$Y),"fd")
  ow$Y = log(ow$Y )


  params = NULL
  out = NULL
  yrs = sort(unique( ow$yr ))
  yrs = yrs[ yrs>=1999 ]

  regions= c("cfanorth", "cfasouth", "cfa4x" )

  for (y in yrs ) { 
  for (r in regions) {  
    
    ir = filter.region.polygon(x=ow[ , c("plon", "plat")], region=r, planar=T)
    it = which(ow$yr==y )
    iy = intersect( ir, it )
    if (length(iy) < 10)  next()
    oy = ow[ iy , c("Y", "plon", "plat" ) ]
    
    oz = empirical.variogram( oy, p )
    ozp = oz[[1]]
    ozp$yr = y
    ozp$region = r
    out = rbind( out, ozp )
    
    ozq = c( unlist( c(oz[[2]][["psill"]], oz[[2]][["range"]] )) , r, y, min(ozp$gamma), max(ozp$gamma) )
    names(ozq) = c( "nugget", "psill","dummy", "range", "region", "yr", "gamma.min", "gamma.max" )
    
    params = rbind( params, ozq)
  }}

  params = as.data.frame( params, stringsAsFactors=F )
  rownames( params ) = NULL
  for (i in setdiff(names(params), "region") ) params[,i] = as.numeric( params[,i]  )
  
  params$ratio  = params$gamma.min / params$gamma.max
  xyplot(ratio ~ yr | region, data=params, type="b" )

  
  params$nugget2  = params$nugget / params$gamma.max
  plot(params$nugget2 ~yrs, type="b" )



  out = as.data.frame(out)
  out$dist =  cut( out$dist, breaks=c(0,p$vgm.dist), labels=p$vgm.dist )

  require(lattice)
  xyplot( gamma~dist | yr, out )
  vg.mat = xtabs( gamma ~ yr+dist, out )
  cbrks = seq( min( vg.mat), max(vg.mat), length.out=100 )
  image(vg.mat, col=terrain.colors(length(cbrks)+1) )


  empirical.variogram = function( oy, p ) {
    g  = gstat(id=v, formula=Y~1, loc=~plon+plat, data=oy, maxdist=p$drange, nmax=p$nmax, omax=floor(p$nmax/4) )
    vg.dat = variogram( g, cutoff=p$drange, cressie=p$variogram.cressie, boundaries=p$vgm.dist)
    var.dat = var( oy$Y )
    psill.start = 0.5*var.dat
    nugget.start = 0.5*var.dat
    vg.mod = vgm( psill=psill.start, model=p$model, range=p$drange, nugget=nugget.start)variogram
    vg.fit = fit.variogram(vg.dat, vg.mod )
    # plot(vg.dat, model=vg.fit)
    return ( list(vg.dat, vg.fit) )
  }


  require(nlme)

  mrange = mean( params$range[ which(params$range < p$drange )] )
  mpsill = mean( params$psill[ which(params$psill < p$drange )]  )
  sp.cor = corSpher(c( vg.fit$range[2], vg.fit$psill[1] /vg.fit$psill[2] ), form= ~plon+plat | yr, nugget=TRUE )

  sp.cor = Initialize( sp.cor, oy )
  
  plot( Variogram( sp.cor, distance=p$vgm.dist  ))

  sp.cor = corSpher( form= ~plon+plat | yr )
  oo <- lme( Y ~ yr, ow, random= ~ 1|yr, correlation=sp.cor )



