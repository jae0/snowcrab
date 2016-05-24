
 
  require(rjags)
  rjags::load.module("dic")
  rjags::load.module("glm")
	
	loadfunctions( "bayesian") 
	loadfunctions( "snowcrab", functionname="initialise.local.environment.r") 


  ###  all data follow this sequence: c("cfanorth", "cfasouth", "cfa4x")
  redo.data=F
  if (redo.data) { 
    biomass.summary.db("complete.redo", p=p ) 
  }
  res = biomass.summary.db( p=p )  


  sb = list( 
    FB0x = c(0.8, 0.6, 0.1),  # priors of the starting fishable biomass : N,S,4X 
    q0x = c(1, 1, 1),  # priors of "catchability" 
    K0x = c(4, 50, 1 ),  # mean carrying capacity estimate 
    IOA = as.matrix(res$B), # observed index of fishable biomass 
    CAT = as.matrix(res$L) , # observed landings
    REC = as.matrix(res$R) , # observed recruitment
    qREC0 = 1, # q correction for recruitment
    sigma = 0.3, # default upper limit in error (CV) of data 
    sigmaB = 0.3, # default upper limit in error (CV) of FB data
    cv = 0.5, # ... 
    er = 0.2,  # target exploitation rate
    MFB = 0.8, # prior estimate of natural mortality for FB
    MFB.sd = 0.8, # prior estimate of SD in natural mortality for FB
    MREC = 0.8, # prior estimate of natural mortality for FB
    MREC.sd = 0.8, # prior estimate of SD in natural mortality for FB
    brodycoef = 2,  # proportional increase in weight with each age
    brodysd0 = 0.2,  # proportional increase in weight with each age
    N = nrow( res$B) , # no years with data  
    M = 5, # no years for projections  
    R = ncol( res$B),  # number of regions
    ty=7,  # index of the transition year (2004) between spring and fall surveys 
    cfa4x=3, # index of cfa4x
    eps = 1e-4  # small non-zero number
  )

  # MCMC/Gibbs parameters
  n.adapt = 5000 # burn-in
  n.chains = 3
  n.thin = 100
  n.iter = 10000


  if (modelrun == "simple" ) {
    ## model 1 --- simple delay difference with observation and process error
    m = jags.model( file=fishery.model.jags ( DS="delay.difference" ), data=sb, n.chains=n.chains, n.adapt=n.adapt )  
    coef(m)
    tomonitor = c("FB","K", "REC", "Znat", "Ztot", "q.biomass", "qREC", "Catch", "Cr" )
    dic.samples(m, n.iter=n.iter ) # pDIC
    fnres = file.path( project.datadirectory("snowcrab"), "R", "delaydifference.mcmc.simple.rdata" )

  }

  
  if (modelrun == "simple.illegal" ) {
    ## model 2 --- simple delay difference with observation and process error and illegal landings
    m = jags.model( file=fishery.model.jags ( DS="delay.difference.illegal" ), data=sb, n.chains=n.chains, n.adapt=n.adapt )  
    coef(m)
    tomonitor =  c("FB","K", "REC", "Znat", "Ztot", "q.biomass", "qREC", "Catch", "Cr" )
    dic.samples(m, n.iter=n.iter ) # pDIC
    fnres = file.path( project.datadirectory("snowcrab"), "R", "delaydifference.mcmc.simple.illegal.rdata" )
  
  }
 

  tomonitor = intersect( variable.names (m), tomonitor )


  # convergence testing -- by 1000 to 1500 convergence observed by Gelman shrink factor diagnostic
  convergence.test = F
  if (convergence.test) {
    y = jags.samples(m, variable.names=tomonitor, n.iter=10000, thin=10)
    gelman.plot(y[["Ztot"]])
    gelman.plot(y[["K"]])
    gelman.plot(y[["q.biomass"]])
    gelman.plot(y[["qREC"]])
    gelman.plot(y[["sd.K"]])
    gelman.plot(y[["sd.o"]])
    gelman.plot(y[["sd.p"]])
    geweke.plot(y[["r"]])
  }
  

  # autocorrelation thinning
    y = coda.samples(m, variable.names=c("K", "r", "q"), n.iter=10000, thin=10) # sample from posterior
    autocorr.plot(y)
    # plot(y, ask=T)
    # autocorr(y, lags = c(0, 1, 5, 10, 50), relative=TRUE)


  # update if not yet converged 
    update(m, n.iter=n.iter ) # above seems enough for convergence but a few more to be sure


  # final sampling from the posteriors
    n.iter.final = n.iter * n.thin
    # n.iter.final = n.iter 
    y = jags.samples(m, variable.names=tomonitor, n.iter=n.iter.final, thin=n.thin) # sample from posterior


  # save(y, file=fnres, compress=T)
  # load( fnres )


  # Figures
    dir.output = file.path( dirname(p$ofname), "figures", "bugs")   
    dir.create( dir.output, recursive=T, showWarnings=F )


    # frequency density of key parameters
    figure.bugs( "K", y=y, fn=file.path(dir.output, "K.density.png" ) ) 
    figure.bugs( "r", y=y, fn=file.path(dir.output, "r.density.png" ) ) 
    figure.bugs( "q", y=y, fn=file.path(dir.output, "q.density.png" ) ) 
    figure.bugs( "BMSY", y=y, fn=file.path(dir.output, "BMSY.density.png" ) ) 
  
    # timeseries
    figure.bugs( type="timeseries", var="biomass", y=y, fn=file.path(dir.output, "biomass.timeseries.png" ) ) 
    figure.bugs( type="timeseries", var="fishingmortality", y=y, fn=file.path(dir.output, "fishingmortality.timeseries.png" ) ) 
    
    # Harvest control rules
    figure.bugs( type="hcr", var="default", y=y, fn=file.path(dir.output, "hcr.default.png" ) ) 
    figure.bugs( type="hcr", var="simple", y=y, fn=file.path(dir.output, "hcr.simple.png" ) ) 
     
    # diagnostics
    figure.bugs( type="diagnostics.production", y=y, fn=file.path(dir.output, "diagnostic.production.png" ) ) 
    figure.bugs( type="diagnostics....", y=y, fn=file.path(dir.output, "diagnostic. ... .png" ) ) 

   


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

      # densities of F in previous year
      for (i in 1:3) plot(density( y$F[ndata-1,i,,] ), main="")
      qs = apply( y$F[ndata-1,,,], 1, quantile, probs=c(0.025, 0.5, 0.975) )
      qs

      # F for table ---
      summary(y$F, median)



  
  # B timeseries with error
  graphics.off()
  layout( matrix(c(1,2,3), 3, 1 ))
  par(mar = c(5, 4, 0, 2))

  vv = "B"
  Xm = jags.extract( y, vv, mean )
  X = y[[vv]]
  nR = dim(X)[2]
  for (r in 1:nR) { 
    nn = create.histograms.yearly( X[,r,,] )
    nn = nn[ 1:(nrow(nn)-1) ,]  # truncate larges and smallest data points as they are not linear increments
    nnr = as.numeric(rownames(nn) )
    if (r==2) nn = nn[ which( nnr < 100 ) ,]
    if (r==3) nn = nn[ which( nnr < 2 ) ,]
    colnames(nn) = yrs
    plot.histograms.temporal( nn, overlaydata=Xm[,r], xlab="Year", ylab="Biomass; kt", barscale=.9, barcol="gray" ) 
  }

  
  
  # F timeseries with error
  graphics.off()
  layout( matrix(c(1,2,3), 3, 1 ))
  par(mar = c(5, 4, 0, 2))

  vv = "F"
  Xm = jags.extract( y, vv, mean )
  X = y[[vv]]
  X[ X>1] = 1
  X = X[ 1:length(yrs0), ,, ]

  nR = dim(X)[2]
  for (r in 1:nR) { 
    nn = create.histograms.yearly( X[,r,,] )
    nn = nn[ 2:(nrow(nn)-1) ,]  # truncate larges and smallest data points as they are not linear increments
    uu = scale( nn, center=FALSE, scale=colSums(nn) ) 
    colnames(nn) = yrs0
    plot.histograms.temporal( nn, overlaydata=Xm[1:length(yrs0),r], xlab="Year", ylab="Fishing mortality", barscale=0.9, barcol="gray" ) 
  }

  



