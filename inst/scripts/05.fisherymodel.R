
## ---------
#### final estimation of biomass via fishery models and associated figures and tables:

#Pick whichever year reference below is correct (most often year.assessment...-1)
if (!exists("year.assessment")) {
  # year.assessment=lubridate::year(Sys.Date())
  # year.assessment=lubridate::year(Sys.Date()) - 1
}
p = bio.snowcrab::load.environment( year.assessment=year.assessment )


# update data summaries
p$vars.tomodel="R0.mass"
biomass.summary.db("complete.redo", p=p) #Uses the model results to create a habitat area expanded survey index


p$fishery_model = list()
p$fishery_model$method = "stan"  # "jags", etc.
p$fishery_model$outdir = file.path(project.datadirectory('bio.snowcrab'), "assessments", p$year.assessment )
p$fishery_model$fnres  = file.path(p$fishery_model$outdir, paste( "surplus.prod.mcmc", p$year.assessment, p$fishery_model$method, "rdata", sep=".") )
p$fishery_model$standata = fishery_model( p=p, DS="stan_data" )
p$fishery_model$stancode = fishery_model( p=p, DS="stan_surplus_production" )
p$fishery_model$stancode_compiled = rstan::stan_model( model_code=p$fishery_model$stancode )

# later:::ensureInitialized()  # solve mode error

res = fishery_model( p=p, DS="stan",
  chains=4,
  iter=10000,
  warmup=4000,
  refresh = 1000,
  control = list(adapt_delta = 0.96, max_treedepth=15)
)
  # warmup = 200,          # number of warmup iterations per chain
  # control = list(adapt_delta = 0.9),
  # # refresh = 500,          # show progress every 'refresh' iterations
  # iter = 1000,            # total number of iterations per chain
  # chains = 5,             # number of Markov chains
  # cores = 5              # number of cores (using 2 just for the vignette)

# load( p$fishery_model$fnres )


# frequency density of key parameters
figure.mcmc( "K", res=res, fn=file.path(p$fishery_model$outdir, "K.density.png" ) )
figure.mcmc( "r", res=res, fn=file.path(p$fishery_model$outdir, "r.density.png" ) )
figure.mcmc( "q", res=res, fn=file.path(p$fishery_model$outdir, "q.density.png" ) ,xrange=c(0,2))
figure.mcmc( "FMSY", res=res, fn=file.path(p$fishery_model$outdir, "FMSY.density.png" ) )
figure.mcmc( "bosd", res=res, fn=file.path(p$fishery_model$outdir, "bosd.density.png" ) )
figure.mcmc( "bpsd", res=res, fn=file.path(p$fishery_model$outdir, "bpsd.density.png" ) )

# timeseries
figure.mcmc( type="timeseries", vname="biomass", res=res, fn=file.path(p$fishery_model$outdir, "biomass.timeseries.png" ), save.plot=T )
figure.mcmc( type="timeseries", vname="fishingmortality", res=res, fn=file.path(p$fishery_model$outdir, "fishingmortality.timeseries.png" ) )

# Harvest control rules
figure.mcmc( type="hcr", vname="default", res=res, fn=file.path(p$fishery_model$outdir, "hcr.default.png" ), save.plot=T  )
figure.mcmc( type="hcr", vname="simple", res=res, fn=file.path(p$fishery_model$outdir, "hcr.simple.png" ) )

# diagnostics
figure.mcmc( type="diagnostic.production", res=res, fn=file.path(p$fishery_model$outdir, "diagnostic.production.png" ) )
figure.mcmc( type="diagnostic.errors", res=res, fn=file.path(p$fishery_model$outdir, "diagnostic.errors.png" ) )
figure.mcmc( type="diagnostic.phase", res=res, fn=file.path(p$fishery_model$outdir, "diagnostic.phase.png" ) )

# K
plot.new()
layout( matrix(c(1,2,3), 3, 1 ))
par(mar = c(4.4, 4.4, 0.65, 0.75))
for (i in 1:3) plot(density(res$mcmc$K[,i] ), main="")
( qs = apply(  res$mcmc$K[,], 2, quantile, probs=c(0.025, 0.5, 0.975) ) )

# R
plot.new()
layout( matrix(c(1,2,3), 3, 1 ))
par(mar = c(4.4, 4.4, 0.65, 0.75))
for (i in 1:3) plot(density(res$mcmc$r[,i] ), main="")
( qs = apply(  res$mcmc$r[,], 2, quantile, probs=c(0.025, 0.5, 0.975) ) )

# q
plot.new()
layout( matrix(c(1,2,3), 3, 1 ))
par(mar = c(4.4, 4.4, 0.65, 0.75))
for (i in 1:3) plot(density(res$mcmc$q[,i] ), main="")
( qs = apply(  res$mcmc$q[,], 2, quantile, probs=c(0.025, 0.5, 0.975) ) )

# FMSY
plot.new()
layout( matrix(c(1,2,3), 3, 1 ))
par(mar = c(4.4, 4.4, 0.65, 0.75))
for (i in 1:3) plot(density(res$mcmc$FMSY[,i] ), main="")
( qs = apply(  res$mcmc$FMSY[,], 2, quantile, probs=c(0.025, 0.5, 0.975) ) )


# densities of biomass estimates for the year.assessment
plot.new()
layout( matrix(c(1,2,3), 3, 1 ))
par(mar = c(4.4, 4.4, 0.65, 0.75))
for (i in 1:3) plot(density(res$mcmc$B[,res$sb$N,i] ), main="")
( qs = apply(  res$mcmc$B[,res$sb$N,], 2, quantile, probs=c(0.025, 0.5, 0.975) ) )

# densities of biomass estimates for the previous year
plot.new()
layout( matrix(c(1,2,3), 3, 1 ))
par(mar = c(4.4, 4.4, 0.65, 0.75))
for (i in 1:3) plot(density( res$mcmc$B[,res$sb$N-1,i] ), main="")
( qs = apply(  res$mcmc$B[,res$sb$N-1,], 2, quantile, probs=c(0.025, 0.5, 0.975) ) )

# densities of F in assessment year
plot.new()
layout( matrix(c(1,2,3), 3, 1 ))
par(mar = c(4.4, 4.4, 0.65, 0.75))
for (i in 1:3) plot(density(  res$mcmc$F[,res$sb$N,i] ), xlim=c(0.01, 0.6), main="")
( qs = apply(  res$mcmc$F[,res$sb$N,], 2, quantile, probs=c(0.025, 0.5, 0.975) ) )
( qs = apply(  res$mcmc$F[,res$sb$N,], 2, mean ) )

# densities of F in previous year
plot.new()
layout( matrix(c(1,2,3), 3, 1 ))
par(mar = c(4.4, 4.4, 0.65, 0.75))
for (i in 1:3) plot(density(  res$mcmc$F[,res$sb$N-1,i] ), xlim=c(0.01, 0.6), main="")
( qs = apply(  res$mcmc$F[,res$sb$N-1,], 2, quantile, probs=c(0.025, 0.5, 0.975) ) )
( qs = apply(  res$mcmc$F[,res$sb$N-1,], 2, mean ) )

# F for table ---
summary( res$mcmc$F, median)
