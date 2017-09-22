
## --------- 
#### final estimation of biomass via fishery models and associated figures and tables:


if (!exists("year.assessment")) year.assessment=lubridate::year(Sys.Date())
p = bio.snowcrab::load.environment( year.assessment=year.assessment )



p$fishery_model = list()
p$fishery_model$method = "stan"  # "jags", etc.
p$fishery_model$outdir = file.path(project.datadirectory('bio.snowcrab'), "assessments", p$year.assessment )
p$fishery_model$fnres  = file.path(p$fishery_model$outdir, paste( "surplus.prod.mcmc", p$year.assessment, p$fishery_model$method, "rdata", sep=".") )


res = fishery_model( p=p, DS=p$fishery_model$method )
# load( p$fishery_model$fnres )


# frequency density of key parameters
figure.mcmc( "K", res=res, fn=file.path(p$fishery_model$outdir, "K.density.png" ) )
figure.mcmc( "r", res=res, fn=file.path(p$fishery_model$outdir, "r.density.png" ) )
figure.mcmc( "q", res=res, fn=file.path(p$fishery_model$outdir, "q.density.png" ) ,xrange=c(0,2))
figure.mcmc( "FMSY", res=res, fn=file.path(p$fishery_model$outdir, "FMSY.density.png" ) )
figure.mcmc( "bo.sd", res=res, fn=file.path(p$fishery_model$outdir, "bo.sd.density.png" ) )
figure.mcmc( "bp.sd", res=res, fn=file.path(p$fishery_model$outdir, "bp.sd.density.png" ) )

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

# densities of biomass estimates for the year.assessment
for (i in 1:3) plot(density(res$mcmc$B[res$sb$N,i,,] ), main="")
( qs = apply(  res$mcmc$B[res$sb$N,,,], 1, quantile, probs=c(0.025, 0.5, 0.975) ) )

# densities of biomass estimates for the previous year
for (i in 1:3) plot(density( res$mcmc$B[res$sb$N-1,i,,] ), main="")
( qs = apply(  res$mcmc$B[res$sb$N-1,,,], 1, quantile, probs=c(0.025, 0.5, 0.975) ) )

# densities of F in assessment year
for (i in 1:3) plot(density(  res$mcmc$F[res$sb$N,i,,] ), xlim=c(0.05, 0.5), main="")
( qs = apply(  res$mcmc$F[res$sb$N,,,], 1, quantile, probs=c(0.025, 0.5, 0.975) ) )
( qs = apply(  res$mcmc$F[res$sb$N,,,], 1, mean ) )

# densities of F in previous year
for (i in 1:3) plot(density(  res$mcmc$F[res$sb$N-1,i,,] ), main="")
( qs = apply(  res$mcmc$F[res$sb$N-1,,,], 1, quantile, probs=c(0.025, 0.5, 0.975) ) )

# F for table ---
summary( res$mcmc$F, median)




