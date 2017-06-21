
  
if (!exists("current.year")) current.year=lubridate::year(Sys.Date())
p = bio.snowcrab::load.environment( year.assessment=current.year)

require(rjags)
rjags::load.module("dic")
rjags::load.module("glm")

#bioLibrary( "bio.models" )

# MCMC/Gibbs parameters
n.adapt = 4000 # burn-in  .. 4000 is enough for the full model but in case ...
n.iter = 30000
n.chains = 3
n.thin = 100 # use of uniform distributions causes high autocorrelations ?
n.iter.final = n.iter * n.thin


sb = switch( as.character(p$year.assessment),
  "2013" = surplusproduction.db( DS="jags.2013", sourcedata="default" ) ,
  "2014" = surplusproduction.db( DS="jags.2014", sourcedata="nosa" ) ,
  "2015" = surplusproduction.db( DS="jags.2015", sourcedata="nosa" ) ,
  "2016" = surplusproduction.db( DS="jags.2016", sourcedata="default" ) 
  )



dir.output = file.path(project.datadirectory('bio.snowcrab'),"assessments",p$year.assessment)
dir.create( dir.output, recursive=T, showWarnings=F )
fnres = file.path( project.datadirectory("bio.snowcrab"), "R", paste( "surplus.prod.mcmc", p$year.assessment,"rdata", sep=".") )
#fnres = file.path( project.datadirectory("bio.snowcrab"), "R", paste( "surplus.prod.mcmc", p$year.assessment,"survey_final.rdata", sep=".") )


m = jags.model( file=fishery.model.jags ( DS=sb$jagsmodelname ), data=sb, n.chains=n.chains, n.adapt=n.adapt ) # recruitment + spring/summer q's + all observed CVs
tomonitor = intersect( variable.names (m), sb$tomonitor )
y = jags.samples(m, variable.names=tomonitor, n.iter=n.iter.final, thin=n.thin) # sample from posterior
dic.samples(m, n.iter=n.iter ) # pDIC

graphics.off() ; x11(); layout( matrix(c(1,2,3), 3, 1 )); par(mar = c(5, 4, 0, 2))
for( i in 1:3) hist(y$cv.r[i,,], "fd")

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

# determine autocorrelation thinning
y = coda.samples(m, variable.names=c("K", "r", "q"), n.iter=20000, thin=10) # sample from posterior
autocorr.plot(y)   # about 10 to 20 required
# plot(y, ask=T)
# autocorr(y, lags = c(0, 1, 5, 10, 50), relative=TRUE)

# final sampling from the posteriors
#  y = jags.samples(m, variable.names=tomonitor, n.iter=10000, thin=20) # sample from posterior
y = jags.samples(m, variable.names=tomonitor, n.iter=n.iter.final, thin=n.thin) # sample from posterior
save(y, file=fnres, compress=T)
# load( fnres )



## 
b1=apply( y$B[,1,,], 1,quantile , probs=c(0.025,0.5,0.975), na.rm=T  )
f1=apply( y$F[,1,,], 1,quantile , probs=c(0.5), na.rm=T  )
NENS=data.frame(t(b1),F=f1,row.names=1999:2019)
NENS$U=1-exp(-NENS$F)
b2=apply( y$B[,2,,], 1,quantile , probs=c(0.025,0.5,0.975), na.rm=T  )
f2=apply( y$F[,2,,], 1,quantile , probs=c(0.5), na.rm=T  )
SENS=data.frame(t(b2),F=f2,row.names=1999:2019)
SENS$U=1-exp(-SENS$F)

b3=apply( y$B[,3,,], 1,quantile , probs=c(0.025,0.5,0.975), na.rm=T  )
f3=apply( y$F[,3,,], 1,quantile , probs=c(0.5), na.rm=T  )
CFA4X=data.frame(t(b3),F=f3,row.names=1999:2019)
CFA4X$U=1-exp(-CFA4X$F)

NENS
SENS
CFA4X
# Figures
graphics.off()


# frequency density of key parameters
figure.bugs( "K", y=y, sb=sb, fn=file.path(dir.output, "K.density.png" ) )
figure.bugs( "r", y=y, sb=sb, fn=file.path(dir.output, "r.density.png" ) )
figure.bugs( "q", y=y, sb=sb, fn=file.path(dir.output, "q.density.png" ) ,xrange=c(0,2))
figure.bugs( "FMSY", y=y, sb=sb, fn=file.path(dir.output, "FMSY.density.png" ) )
figure.bugs( "bo.sd", y=y, sb=sb, fn=file.path(dir.output, "bo.sd.density.png" ) )
figure.bugs( "bp.sd", y=y, sb=sb, fn=file.path(dir.output, "bp.sd.density.png" ) )

# timeseries
figure.bugs( type="timeseries", vname="biomass", y=y, sb=sb, fn=file.path(dir.output, "biomass.timeseries.png" ), save.plot=T )
figure.bugs( type="timeseries", vname="fishingmortality", y=y, sb=sb, fn=file.path(dir.output, "fishingmortality.timeseries.png" ) )

# Harvest control rules
figure.bugs( type="hcr", vname="default", y=y, sb=sb, fn=file.path(dir.output, "hcr.default.png" ), save.plot=T  )
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


