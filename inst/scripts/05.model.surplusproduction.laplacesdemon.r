
p = bio.snowcrab::load.environment( year.assessment=2016)


# choose data
# res =  biomass.summary.survey.db(p=p)
res = biomass.summary.survey.nosa.db(p=p)


# install_github("lazycrazyowl/LaplacesDemon")
# install_github("lazycrazyowl/LaplacesDemonCpp")
require(LaplacesDemon)
# require(LaplacesDemonCpp)

# data object
region = "cfasouth"
biomassindex = res$B[[region]]
catch = res$L[[region]]
yrs = as.numeric(rownames( res$B))


sb = list(
  O = biomassindex, # observed index of abundance
  Omissing0 = mean(biomassindex, na.rm=TRUE ),
  removals = catch , # removalsches  , assume 20% handling mortality and illegal landings
  removalsmissing0 = mean(catch, na.rm=TRUE ),
  er = 0.2,  # target exploitation rate
  N = length( biomassindex ) , # no years with data
  M = 5, # no years for projections
  MN = length( biomassindex ) + 5,
  ty = which(yrs==2004) ,  # index of the transition year (2004) between spring and fall surveys
  r0= 1,
  K0= 100*1.25,
  q0= mean(biomassindex, na.rm=TRUE)/100,
  S0= 0.6, # normalised
  cv = 0.5,
  smax =1.25,
  ii = 2:length(biomassindex),
  jj = 1:(length(biomassindex)-1),
  eps = 1e-6
)

loadfunctions( "bio.snowcrab", "R", "surplus.production.simple.laplacesdemon.setup.r" )


# set up the model
sb = surplus.production.simple.laplacesdemon.setup( sb )  

sb$Model( parm=sb$PGF(sb), Data=sb )



n.iter = 500 
tol = 1.0E-6

f.ml = optim( par=sb$PGF(sb), fn=sb$Model.ML, Data=sb, hessian=TRUE, control=list(maxit=5000)  )
names(f.ml$par ) = sb$parm.names
f.ml$par # check convergence
print(sqrt( diag( solve(f.ml$hessian) )) ) # asymptotic standard errors


# PML .. better ("penalised ML")
# penalized maximum likelihood .. better but still a little unstable depending on algorithm
#f.pml = optim( par=sb$PGF(sb), fn=sb$Model.PML, Data=sb,  method="L-BFGS-B", control=list(maxit=5000), hessian=TRUE )
#f.pml = optim( par=sb$PGF(sb), fn=sb$Model.PML, Data=sb,  method="BFGS", control=list(maxit=5000), hessian=TRUE )
f.pml = optim( par=sb$PGF(sb), fn=sb$Model.PML, Data=sb,  control=list(maxit=5000), hessian=TRUE )
names(f.pml$par ) = sb$parm.names
f.pml$par

str(f.pml) # check convergence
print(sqrt( diag( solve(f.pml$hessian) )) ) # assymptotic standard errors


f = LaplacesDemon(sb$Model, Data=sb, Initial.Values=sb$PGF(sb), Iterations=1000000, Status=2000, Thinning=1000)

# quick solution: "burn-in"
f = LaplaceApproximation(sb$Model, Data=sb, parm=sb$PGF(sb) ) 

f = LaplaceApproximation(sb$Model, Data=sb, parm=as.initial.values(f), Method="optim", Stop.Tolerance=1e-8, Iterations = 5000, method="Nelder-Mead", control=list(maxit=5, reltol=1e-9 ) )
f

f = LaplaceApproximation(sb$Model, Data=sb, parm=as.initial.values(f), Method="optim", Stop.Tolerance=1e-8, Iterations = 1000, method="L-BFGS-B", control=list(maxit=10, reltol=1e-8 ) )


f = LaplaceApproximation(sb$Model, Data=sb, parm=as.initial.values(f), Method="optim", optim.params=list(method="BFGS") , Stop.Tolerance=1e-9, Iterations = 1000  )


f = LaplaceApproximation(sb$Model, Data=sb, parm=as.initial.values(f)  ) # fast spin up of paramters
f = LaplaceApproximation(sb$Model, Data=sb, parm=as.initial.values(f), Method="BFGS", Stop.Tolerance=1e-8  ) 

Deviance  1.530190e+01  0    0 1000  1.530190e+01  1.530190e+01  1.530190e+01
LP       -1.778748e+03  0    0 1000 -1.778748e+03 -1.778748e+03 -1.778748e+03
r         7.970775e-01  0    0 1000  7.970775e-01  7.970775e-01  7.970775e-01
K         6.277875e+01  0    0 1000  6.277875e+01  6.277875e+01  6.277875e+01
q         8.223097e-03  0    0 1000  8.223097e-03  8.223097e-03  8.223097e-03

f = LaplaceApproximation(sb$Model, Data=sb, parm=as.initial.values(f), Method="LBFGS"  ) 


# NOTES: NM is slow
# algorithms = c( "TR", "LM", "CG", "NR", "NM", "SPG", "LBFGS", "BFGS", "PSO", "SR1", "HAR", "DFP", "BHHH", "HJ", "Rprop" ) # available
# algorithms = c( "CG", "NM", "SPG", "LBFGS", "BFGS", "PSO", "SR1", "HAR", "DFP", "Rprop" ) # reliably working
algorithms = c( "CG", "LBFGS", "BFGS", "PSO", "SR1", "HAR", "DFP", "Rprop" ) # reliably working
for (a in algorithms) {
  print(a)
  ft = try( LaplaceApproximation(sb$Model, Data=sb, parm=as.initial.values(f), Method=a, Stop.Tolerance=tol, Iterations=n.iter ) )
  if (! class(ft) %in% "try-error" ) if (ft$Converged) if( ft$LP.Final > f$LP.Final )  {
    f = ft
    f$Method = a 
  }
}
str(f)
plot(f, Data=sb)

# MCMC ... look at ACF
f = LaplacesDemon(sb$Model, Data=sb, Initial.Values=as.initial.values(f), Covar=f$Covar, 
  Iterations=1000, Status=100, Thinning=1, Algorithm="CHARM" )
Consort(f)
plot(f, Data=sb)

 PosteriorChecks(f)
     #caterpillar.plot(f, Parms="beta")
     #plot(f, Data, PDF=FALSE)
     #Pred <- predict(f, Model, Data, CPUs=1)

  #summary(Pred, Discrep="Chi-Square")
     #plot(Pred, Style="Covariates", Data=Data)
     #plot(Pred, Style="Density", Rows=1:9)
     #plot(Pred, Style="Fitted")
     #plot(Pred, Style="Jarque-Bera")
     #plot(Pred, Style="Predictive Quantiles")
     #plot(Pred, Style="Residual Density")
     #plot(Pred, Style="Residuals")
     #Levene.Test(Pred)
     #Importance(f, Model, Data, Discrep="Chi-Square")

# MCMC: run with appropriate thinning and other options:
f = LaplacesDemon(sb$Model, Data=sb, as.initial.values(f),
  Covar=NULL, Iterations=10000, Status=1000, Thinning=1000, Algorithm="CHARM", Specs=NULL)

f = LaplacesDemon(sb$Model, Data=sb, as.initial.values(f),
  Covar=f$Covar , Iterations=50000, Status=10204, Thinning=1000, Algorithm="CHARM", Specs=list(alpha.star=0.44))


f = LaplacesDemon(sb$Model, Data=sb, 
  Initial.Values=as.initial.values(f), Covar=f$Covar, Iterations=10000, Status=1000, Thinning=100)

f = LaplacesDemon(sb$Model, Data=sb, 
  Initial.Values=as.initial.values(f), Covar=f$Covar, Iterations=30000, Status=1000, Thinning=35, 
  Algorithm="AMWG", Specs=list(B=NULL, n=1000, Periodicity=35 ) )

f <- LaplacesDemon(sb$Model, Data=sb,
  Initial.Values=as.initial.values(f), Covar=f$Covar, Iterations=10000, Status=1000, Thinning=105, 
  Algorithm="AFSS", Specs=list(A=Inf, B=NULL, m=100, n=0, w=1))

# medium speed .. increase Iterations til convergence
f = VariationalBayes(sb$Model, Data=sb,  parm=as.initial.values(f), Iterations=100,  Samples=10, CPUs=5 )
f = VariationalBayes(sb$Model, Data=sb,  parm=as.initial.values(f), Iterations=1000,  Samples=1000, CPUs=5 Stop.Tolerance=1e-9)

# slow
f = IterativeQuadrature(sb$Model, Data=sb, parm=as.initial.values(f), Iterations=10, Algorithm="AGH",
 Specs=list(N=5, Nmax=7, Packages=NULL, Dyn.libs=NULL) )



# construct data structured for LaplacesDemon:
  sb = list(
    b.min = 0.001, # scaled to 1 but allow overshooting
    b.max = 1.1, # scaled to 1 but allow overshooting
    q.min = rep(0.001,3),
    q.max = rep(1,3),
    #q.min = rep(0.1,3),  # max value of q , anything larger is not believable
    #q.max = rep(2,3),  # max value of q , anything larger is not believable
   # r.min = 0.3, #for hyperprior
   # r.max = 2.0,  #for hyperprior
    bo.mup=rep(-2.5579,3),
    bo.sdp=rep(0.47726,3),
    bp.mup=rep(-2.5579,3),
    bp.sdp=rep(0.47726,3),
    rec.max= c( 10^3, 10^4, 10^2 ),
    K.mu = c( 1.831139,4.170013,0.784308), #for ln
    K.sd = c(0.04595339,0.04350642,0.02571607), #for ln
    r.mu = rep(0.96,3),
    r.sd = rep(0.01041271,3),
   # K.min = c(1, 10, 0.1 ),  # max carrying capacity estimate: #for hyperprior
   # K.max = c(10, 100, 5 ),  # max carrying capacity estimate: #for hyperprior
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
  sb$IRECcv[ cfa.north.bad.data, cfa.north ] = mean( sb$IRECcv[,cfa.north ] ,na.rm=T )
  sb$IREC[ cfa.north.bad.data, cfa.north ] = mean( sb$IREC[,cfa.north ] ,na.rm=T )


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
  fnres = file.path( project.datadirectory("bio.snowcrab"), "R", paste( "surplus.prod.mcmc", p$year.assessment,"rdata", sep=".") )


  debug =T
  if (debug) {
    n.adapt = 1000 # burn-in
    n.iter = 4000
    n.chains = 3
    n.thin = 10
    n.iter.final = n.iter
    fnres = file.path( project.datadirectory("bio.snowcrab"), "R", "surplus.prod.mcmc.debug.rdata" )
  }


  # -------------------

  m = jags.model( file=fishery.model.jags ( DS="biomassdynamic_nonhyper_2014.bugs" ), data=sb, n.chains=n.chains, n.adapt=n.adapt ) # recruitment + spring/summer q's + all observed CVs


  tomonitor =  c(
    "r", "K", "q", "qs",
    "r.mu", "r.sd",
   # "K.mu", "K.sd",
    #"q.mu", "q.sd",
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
  dir.output = file.path(project.datadirectory('bio.snowcrab'),"assessments","2015")
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


    fnres =  file.path( project.datadirectory("bio.snowcrab"), "R", "surplus.prod.mcmc.2016.survey_final.rdata" )
    # fnres =  file.path( project.datadirectory("bio.snowcrab"), "R", "surplus.prod.mcmc.2012_final.rdata" )
    # fnres =  file.path( project.datadirectory("bio.snowcrab"), "R", "surplus.prod.mcmc.2012a.rdata" )
    save(y, file=fnres, compress=T)
   load( fnres )


  # Figures
    graphics.off()

    dir.output = file.path( dirname(p$ofname), "figures", "bugs","survey")
    dir.create( dir.output, recursive=T, showWarnings=F )

    # frequency density of key parameters
    figure.bugs( "K", y=y, sb=sb, fn=file.path(dir.output, "K.density.png" ) )
    figure.bugs( "r", y=y, sb=sb, fn=file.path(dir.output, "r.density.png" ) )
    #figure.bugs( "r.ts", y=y, sb=sb, fn=file.path(dir.output, "r.ts.density.png" ) )
    figure.bugs( "q", y=y, sb=sb, fn=file.path(dir.output, "q.density.png" ) )
    #figure.bugs( "qs", y=y, sb=sb, fn=file.path(dir.output, "qs.density.png" ) )
    figure.bugs( "FMSY", y=y, sb=sb, fn=file.path(dir.output, "FMSY.density.png" ) )
    figure.bugs( "bo.sd", y=y, sb=sb, fn=file.path(dir.output, "bo.sd.density.png" ) )
    figure.bugs( "bp.sd", y=y, sb=sb, fn=file.path(dir.output, "bp.sd.density.png" ) )

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

    # densities of biomass estimates for the year.assessment
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

