
 
  require(rjags)
  rjags::load.module("dic")
  rjags::load.module("glm")

  
	loadfunctions( "bayesian" )
	loadfunctions( "snowcrab", functionname="initialise.local.environment.r") 


  ###  all data follow this sequence: c("cfanorth", "cfasouth", "cfa4x")
  
  # curve( dnorm(x, mean=1, sd=1*0.5), from=0.1, to=4  )
  # curve( dlnorm(x, meanlog=log(1), sdlog=1),  from=0.1, to=100  )
  # curve( dbeta(x, 10, 10),  from=0.01, to=1  )
  # curve( dgamma(x, shape=0.001, rate=0.001),  from=0.01, to=2  )

  
  redo.data=F
  if (redo.data) { 
    biomass.summary.db("complete.redo", p=p)
  }
  
  res = biomass.summary.db(p=p)
     
  sb = list( 
    b.min = 0.001, # scaled to 1 but allow overshooting
    b.max = 1.1, # scaled to 1 but allow overshooting
    q.min = 0.1,  # max value of q , anything larger is not believable
    q.max = 2.0,  # max value of q , anything larger is not believable
    r.min = 0.1,
    r.max = 2.0,  
    rec.max= c( 10^3, 10^4, 10^2 ), 
    K.min = c(1, 10, 0.1 ),  # max carrying capacity estimate: 
    K.max = c(10, 100, 5 ),  # max carrying capacity estimate: 
    b0.min = c(0.5, 0.5, 0.2),  # prior: mean value possible in  N,S,4X 
    b0.max = c(0.8, 0.8, 0.6),  # prior: mean value possible in  N,S,4X 
    cv.normal.min = 0.05, # upper limit of CV for normally distributed variables ~ 0.5 covers a reasonably large range, try:   curve( dnorm(x, mean=1, sd=0.5), from=0.1, to=4  )
    cv.normal.max = 0.4, # upper limit of CV for normally distributed variables ~ 0.5 covers a reasonably large range, try:   curve( dnorm(x, mean=1, sd=0.5), from=0.1, to=4  )
    cv.lognormal.min = 0.05, #  curve( dlnorm(x, meanlog=log(1), sdlog=0.25), from=0.01, to=2 )
    cv.lognormal.max = 0.4, #  curve( dlnorm(x, meanlog=log(1), sdlog=0.25), from=0.01, to=2 )
  # for lognormal: cv = sqrt(exp(sigma^2) - 1); or sigma = sqrt(log(cv^2+ 1) ) ==> sigma = sqrt( log(0.25^2 + 1)) = 0.246 ~ cv -- i.e. cv ~ sd
    IOA = as.matrix(res$B), # observed index of abundance
    IOAcv = as.matrix(res$B.sd ), # observed index of log abundance SD estimates ~ CV
    IREC = as.matrix(res$R), # observed index of abundance
    IRECcv = as.matrix(res$R.sd ), # observed index of log abundance SD estimates ~CV
    CAT = as.matrix(res$L) , # catches  , assume 20% handling mortality and illegal landings
    CAT.min = apply( res$L, 2, min, na.rm=T), 
    CAT.max = apply( res$L, 2, max, na.rm=T), 
    # p.g = 2,  # dgamma parameter = 2 is an exponential family with mean = p.g * p.h; p.h=mean/p.g
    # cv = 0.4, #default 1/cv  -- muliplier for dgamma on inverse scale
    # logsdmax = 1,  # upper bound of sd of lognormal distributions 1== assume < 1 order of magnitude variation and modal 
    er = 0.2,  # target exploitation rate
    U = ncol( res$B),  # number of regions
    N = nrow( res$B) , # no years with data  
    M = 3, # no years for projections  
    ty = 7,  # index of the transition year (2004) between spring and fall surveys 
    cfa4x = 3, # index of cfa4x
    eps = 1e-4  # small non-zero number
  )
 
# set poor data to NA's and maximize associated CV's
  
  cfa.north =  1 # column index
  cfa.north.bad.data = which( as.numeric(rownames(sb$IOA)) <= 1997 )
  sb$IOA[ cfa.north.bad.data, cfa.north ] = NA 
  sb$IOAcv[ cfa.north.bad.data, cfa.north ] = mean( sb$IOAcv[,cfa.north ] ,na.rm=T )
 
  cfa.south =  2 # column index
  cfa.south.bad.data = which( as.numeric(rownames(sb$IOA)) <= 1998 )
  sb$IOA[ cfa.south.bad.data, cfa.south ] = NA 
  sb$IOAcv[ cfa.south.bad.data, cfa.south ] = mean( sb$IOAcv[,cfa.south ] ,na.rm=T )
 
  cfa.nodata =   which( as.numeric(rownames(sb$IOA)) <= 2003 )
  sb$IOA[ cfa.nodata , sb$cfa4x ] = NA 
  sb$IOAcv[ cfa.nodata , sb$cfa4x ] = mean( sb$IOAcv[,sb$cfa4x] ,na.rm=T )
 
  # recruitment index has no process model right now ... use min value
  sb$IREC[ cfa.nodata , sb$cfa4x ] = min( sb$IREC[,sb$cfa4x] ,na.rm=T )
  sb$IRECcv[ cfa.nodata , sb$cfa4x ] = mean( sb$IRECcv[,sb$cfa4x] ,na.rm=T )

  sb$IREC = log( sb$IREC )
  
 

  # MCMC/Gibbs parameters
  n.adapt = 4000 # burn-in  .. 4000 is enough for the full model but in case ...
  n.iter = 3000 
  n.chains = 3
  n.thin = 100 # use of uniform distributions causes high autocorrelations ? 
  n.iter.final = n.iter * n.thin
  fnres = file.path( project.datadirectory("snowcrab"), "R", paste( "surplus.prod.mcmc", p$current.assessment.year,"rdata", sep=".") )
 

  debug =F
  if (debug) {
    n.adapt = 1000 # burn-in
    n.iter = 4000 
    n.chains = 3
    n.thin = 10
    n.iter.final = n.iter 
    fnres = file.path( project.datadirectory("snowcrab"), "R", "surplus.prod.mcmc.debug.rdata" )
  }


  # -------------------
  ##  simple surplus production with observation and process error
  # m = jags.model( file=fishery.model.jags ( DS="biomassdynamic_basic_2013.bugs" ), data=sb, n.chains=n.chains, n.adapt=n.adapt )  
  # m = jags.model( file=fishery.model.jags ( DS="biomassdynamic_recruitment_2013.bugs" ), data=sb, n.chains=n.chains, n.adapt=n.adapt )  # "recruitment" + spring/summer q's

  m = jags.model( file=fishery.model.jags ( DS="biomassdynamic_full_2013.bugs" ), data=sb, n.chains=n.chains, n.adapt=n.adapt ) # recruitment + spring/summer q's + all observed CVs
 
 
  tomonitor =  c( 
    "r", "K", "q", "qs",
    "r.mu", "r.sd",
    "K.mu", "K.sd",
    "q.mu", "q.sd", 
#    "qs.mu", "qs.sd",
    "b", 
    "bp.sd", "bo.sd", 
    "b0", "b0.sd",
#    "bm", 
    "rem", "rem.sd", "rem.mu", 
#    "ill",
    "REM", 
    "MSY", "BMSY", "FMSY", "Fcrash", "Bdrop", "BX2MSY",
    "F", "TAC",  "C", "P", "B" )
  
  tomonitor = intersect( variable.names (m), tomonitor )
  coef(m)
  

  # ----------------

  dic.samples(m, n.iter=n.iter ) # pDIC

  
  # ----------------

  y = jags.samples(m, variable.names=tomonitor, n.iter=n.iter.final, thin=n.thin) # sample from posterior
  
  figure.bugs( type="timeseries", vname="biomass", y=y, sb=sb, fn=file.path(dir.output, "biomass.timeseries.png" ) ) 
  
  figure.bugs( type="timeseries", vname="fishingmortality", y=y, sb=sb, fn=file.path(dir.output, "fishingmortality.timeseries.png" ) ) 
   
  graphics.off() ; x11()
  layout( matrix(c(1,2,3), 3, 1 )); par(mar = c(5, 4, 0, 2))
  for( i in 1:3) hist(y$cv.r[i,,], "fd")
           
  # ----------------
  # convergence testing -- by 1000 to 1500 convergence observed by Gelman shrink factor diagnostic
    y = jags.samples(m, variable.names=tomonitor, n.iter=6000, thin=1 )
    
    gelman.plot(y[["r"]])
    gelman.plot(y[["K"]])
    gelman.plot(y[["q"]])  # about 6-8000 runs required to converge
    gelman.plot(y[["r.sd"]])
    gelman.plot(y[["K.sd"]]) 
    gelman.plot(y[["bo.sd"]])
    gelman.plot(y[["bp.p"]])
    geweke.plot(y[["r"]])


  # update if not yet converged 
  #  update(m, n.iter=n.iter ) # above seems enough for convergence but a few more to be sure


  # ------------------
  # determine autocorrelation thinning
    y = coda.samples(m, variable.names=c("K", "r", "q"), n.iter=20000, thin=10) # sample from posterior
    autocorr.plot(y)   # about 10 to 20 required
    # plot(y, ask=T)
    # autocorr(y, lags = c(0, 1, 5, 10, 50), relative=TRUE)


  # final sampling from the posteriors
  #  y = jags.samples(m, variable.names=tomonitor, n.iter=10000, thin=20) # sample from posterior
    y = jags.samples(m, variable.names=tomonitor, n.iter=n.iter.final, thin=n.thin) # sample from posterior
		
    
    fnres =  file.path( project.datadirectory("snowcrab"), "R", "surplus.prod.mcmc.2013_final.rdata" )
		# fnres =  file.path( project.datadirectory("snowcrab"), "R", "surplus.prod.mcmc.2012_final.rdata" )
    # fnres =  file.path( project.datadirectory("snowcrab"), "R", "surplus.prod.mcmc.2012a.rdata" )
    save(y, file=fnres, compress=T)
    # load( fnres )


  # Figures
    graphics.off()
    
    dir.output = file.path( dirname(p$ofname), "figures", "bugs")   
    dir.create( dir.output, recursive=T, showWarnings=F )

    # frequency density of key parameters
    figure.bugs( "K", y=y, sb=sb, fn=file.path(dir.output, "K.density.png" ) ) 
    # figure.bugs( "r", y=y, sb=sb, fn=file.path(dir.output, "r.density.png" ) ) 
    figure.bugs( "r.ts", y=y, sb=sb, fn=file.path(dir.output, "r.ts.density.png" ) ) 
    figure.bugs( "q", y=y, sb=sb, fn=file.path(dir.output, "q.density.png" ) ) 
    figure.bugs( "qs", y=y, sb=sb, fn=file.path(dir.output, "qs.density.png" ) ) 
    figure.bugs( "FMSY", y=y, sb=sb, fn=file.path(dir.output, "FMSY.density.png" ) ) 
  
    # timeseries
    figure.bugs( type="timeseries", vname="biomass", y=y, sb=sb, fn=file.path(dir.output, "biomass.timeseries.png" ) ) 
    figure.bugs( type="timeseries", vname="fishingmortality", y=y, sb=sb, fn=file.path(dir.output, "fishingmortality.timeseries.png" ) ) 
      
    # Harvest control rules
    figure.bugs( type="hcr", vname="default", y=y, sb=sb, fn=file.path(dir.output, "hcr.default.png" ) ) 
    figure.bugs( type="hcr", vname="simple", y=y, sb=sb, fn=file.path(dir.output, "hcr.simple.png" ) ) 
     
    # diagnostics
    figure.bugs( type="diagnostic.production", y=y, sb=sb, fn=file.path(dir.output, "diagnostic.production.png" ) ) 
    figure.bugs( type="diagnostic.errors", y=y, sb=sb, fn=file.path(dir.output, "diagnostic.errors.png" ) ) 
    figure.bugs( type="diagnostic.phase", y=y, sb=sb, fn=file.path(dir.output, "diagnostic.phase.png" ) ) 


    ndata = sb$N

    # densities of biomass estimates for the current year
      for (i in 1:3) plot(density(y$B[ndata,i,,] ), main="")
      qs = apply( y$B[ndata,,,], 1, quantile, probs=c(0.025, 0.5, 0.975) )
      qs

      # densities of biomass estimates for the previous year
      for (i in 1:3) plot(density(y$B[ndata-1,i,,] ), main="")
      qs = apply( y$B[ndata-1,,,], 1, quantile, probs=c(0.025, 0.5, 0.975) )
      qs

      # densities of F in assessment year
      for (i in 1:3) plot(density( y$F[ndata,i,,] ), xlim=c(0.05, 0.5), main="")
      qs = apply( y$F[ndata,,,], 1, quantile, probs=c(0.025, 0.5, 0.975) )
      qs
      qs = apply( y$F[ndata,,,], 1, mean )

      # densities of F in previous year
      for (i in 1:3) plot(density( y$F[ndata-1,i,,] ), main="")
      qs = apply( y$F[ndata-1,,,], 1, quantile, probs=c(0.025, 0.5, 0.975) )
      qs

      # F for table ---
      summary(y$F, median)


    debug = F
    if (debug) {
      
      graphics.off()
      
      x11()
      layout( matrix(c(1,2,3), 3, 1 ))
      par(mar = c(5, 4, 0, 2))

      for( i in 1:3) hist(y$r[i,,], "fd")
      for( i in 1:3) hist(y$q.sd[i,,], "fd")
      for( i in 1:3) hist(y$K.sd[i,,], "fd")
      for( i in 1:3) hist(y$b0.sd[i,,], "fd")
      for( i in 1:3) hist(y$r.sd[i,,], "fd")
      for( i in 1:3) hist(y$b0[i,,], "fd")
   
    }


