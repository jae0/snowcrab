
## ---------
#### final estimation of biomass via fishery models and associated figures and tables:

year.assessment = 2020   # NOTE: for 4X, the season 2019-2020 -> 2019


p = bio.snowcrab::snowcrab_parameters(
  project_class="carstm",
  assessment.years=2000:year.assessment,
 # carstm_model_label="nonseparable_simple",  # to choose alt carstm models
  tag="default"
)

p$tag = "tesselation"

p$fishery_model = fishery_model( DS = "logistic_parameters", p=p, tag=p$areal_units_type )
p$fishery_model$stancode = fishery_model( p=p, DS="stan_surplus_production" )
p$fishery_model$stancode_compiled = rstan::stan_model( model_code=p$fishery_model$stancode )


if (0) {
  # testing:
  # later:::ensureInitialized()  # solve mode error
  str( p$fishery_model$standata )
  p$fishery_model$standata$Kmu = c(4, 40, 1)
  p$fishery_model$standata$rmu = c(1, 1, 1)
  p$fishery_model$standata$qmu = c(1, 1, 1)
  p$fishery_model$standata$Ksd = c(0.5, 0.5, 0.5) * p$fishery_model$standata$Kmu  # c( 2, 20, 0.5)
  p$fishery_model$standata$rsd = c(0.5, 0.5, 0.5) * p$fishery_model$standata$rmu  # rep( 0.3, 3)
  p$fishery_model$standata$qsd = c(0.5, 0.5, 0.5) * p$fishery_model$standata$qmu  # rep( 0.3, 3)

  # to recompile
  p$fishery_model$stancode_compiled = rstan::stan_model( model_code=p$fishery_model$stancode )
}


str( p$fishery_model)

res = fishery_model( p=p, DS="logistic_model", tag=p$tag,
  chains=3, cores=3, iter=20000, warmup=12000, refresh=1000,
  control=list(adapt_delta=0.99, max_treedepth=18)
)
# res = fishery_model( p=p, DS="logistic_samples", tag=tag )  # to get samples



p$fishery_model$outdir = file.path( p$modeldir, p$carstm_model_label, "fishery_model_results" )


# frequency density of key parameters
fishery_model( DS="plot", type="K", res=res, fn=file.path(p$fishery_model$outdir, "K.density.png" ) )
fishery_model( DS="plot", type="r", res=res, fn=file.path(p$fishery_model$outdir, "r.density.png" ) )
fishery_model( DS="plot", type="q", res=res, fn=file.path(p$fishery_model$outdir, "q.density.png" ) ,xrange=c(0,5.5))
fishery_model( DS="plot", type="FMSY", res=res, fn=file.path(p$fishery_model$outdir, "FMSY.density.png" ) )
# fishery_model( DS="plot", type="bosd", res=res, fn=file.path(p$fishery_model$outdir, "bosd.density.png" ) )
# fishery_model( DS="plot", type="bpsd", res=res, fn=file.path(p$fishery_model$outdir, "bpsd.density.png" ) )

# timeseries
fishery_model( DS="plot", type="timeseries", vname="biomass", res=res, fn=file.path(p$fishery_model$outdir, "biomass.timeseries.png" ), save.plot=T )
fishery_model( DS="plot", type="timeseries", vname="fishingmortality", res=res, fn=file.path(p$fishery_model$outdir, "fishingmortality.timeseries.png" ) )

# Summary table of mean values for inclusion in document
biomass.summary.table(x)

# Harvest control rules
fishery_model( DS="plot", type="hcr", vname="default", res=res, fn=file.path(p$fishery_model$outdir, "hcr.default.png" ), save.plot=T  )
fishery_model( DS="plot", type="hcr", vname="simple", res=res, fn=file.path(p$fishery_model$outdir, "hcr.simple.png" ) )

# diagnostics
# fishery_model( DS="plot", type="diagnostic.errors", res=res, fn=file.path(p$fishery_model$outdir, "diagnostic.errors.png" ) )
# fishery_model( DS="plot", type="diagnostic.phase", res=res, fn=file.path(p$fishery_model$outdir, "diagnostic.phase.png" ) )


NN = res$p$fishery_model$standata$N

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
for (i in 1:3) plot(density(res$mcmc$B[,NN,i] ), main="")
( qs = apply(  res$mcmc$B[,NN,], 2, quantile, probs=c(0.025, 0.5, 0.975) ) )

# densities of biomass estimates for the previous year
plot.new()
layout( matrix(c(1,2,3), 3, 1 ))
par(mar = c(4.4, 4.4, 0.65, 0.75))
for (i in 1:3) plot(density( res$mcmc$B[,NN-1,i] ), main="")
( qs = apply(  res$mcmc$B[,NN-1,], 2, quantile, probs=c(0.025, 0.5, 0.975) ) )

# densities of F in assessment year
plot.new()
layout( matrix(c(1,2,3), 3, 1 ))
par(mar = c(4.4, 4.4, 0.65, 0.75))
for (i in 1:3) plot(density(  res$mcmc$F[,NN,i] ), xlim=c(0.01, 0.6), main="")
( qs = apply(  res$mcmc$F[,NN,], 2, quantile, probs=c(0.025, 0.5, 0.975) ) )
( qs = apply(  res$mcmc$F[,NN,], 2, mean ) )

# densities of F in previous year
plot.new()
layout( matrix(c(1,2,3), 3, 1 ))
par(mar = c(4.4, 4.4, 0.65, 0.75))
for (i in 1:3) plot(density(  res$mcmc$F[,NN-1,i] ), xlim=c(0.01, 0.6), main="")
( qs = apply(  res$mcmc$F[,NN-1,], 2, quantile, probs=c(0.025, 0.5, 0.975) ) )
( qs = apply(  res$mcmc$F[,NN-1,], 2, mean ) )

# F for table ---
summary( res$mcmc$F, median)
