### --------- 
#### final estimation of biomass via fishery model:

  
if (!exists("current.year")) current.year=lubridate::year(Sys.Date())
p = bio.snowcrab::load.environment( year.assessment=current.year)

p$surplusproduction_model = list()
p$surplusproduction_model$method = "stan"  # "jags", etc.
p$surplusproduction_model$outdir = file.path(project.datadirectory('bio.snowcrab'), "assessments", p$year.assessment ),
p$surplusproduction_model$fnres  = file.path(p$surplusproduction_model$outdir, paste( "surplus.prod.mcmc", p$year.assessment, p$surplusproduction_model$method, "rdata", sep=".") )



if (p$surplusproduction_model$method=="stan"){
  
  res = surplusproduction_model( p=p, DS="stan" )
  # load( p$surplusproduction_model$fnres )

  outdir = p$surplusproduction_model$outdir
    # frequency density of key parameters
    figure.stan( "K", res=res, fn=file.path(outdir, "K.density.png" ) )
    figure.stan( "r", res=res, fn=file.path(outdir, "r.density.png" ) )
    figure.stan( "q", res=res, fn=file.path(outdir, "q.density.png" ) ,xrange=c(0,2))
    figure.stan( "FMSY", res=res, fn=file.path(outdir, "FMSY.density.png" ) )
    figure.stan( "bo.sd", res=res, fn=file.path(outdir, "bo.sd.density.png" ) )
    figure.stan( "bp.sd", res=res, fn=file.path(outdir, "bp.sd.density.png" ) )

    # timeseries
    figure.stan( type="timeseries", vname="biomass", res=res, fn=file.path(outdir, "biomass.timeseries.png" ), save.plot=T )
    figure.stan( type="timeseries", vname="fishingmortality", res=res, fn=file.path(outdir, "fishingmortality.timeseries.png" ) )

    # Harvest control rules
    figure.stan( type="hcr", vname="default", res=res, fn=file.path(outdir, "hcr.default.png" ), save.plot=T  )
    figure.stan( type="hcr", vname="simple", res=res, fn=file.path(outdir, "hcr.simple.png" ) )

    # diagnostics
    figure.stan( type="diagnostic.production", res=res, fn=file.path(outdir, "diagnostic.production.png" ) )
    figure.stan( type="diagnostic.errors", res=res, fn=file.path(outdir, "diagnostic.errors.png" ) )
    figure.stan( type="diagnostic.phase", res=res, fn=file.path(outdir, "diagnostic.phase.png" ) )

    # densities of biomass estimates for the year.assessment
    for (i in 1:3) plot(density(y$B[res$sb$N,i,,] ), main="")
    ( qs = apply(  res$stan$B[res$sb$N,,,], 1, quantile, probs=c(0.025, 0.5, 0.975) ) )

    # densities of biomass estimates for the previous year
    for (i in 1:3) plot(density( res$stan$B[res$sb$N-1,i,,] ), main="")
    ( qs = apply(  res$stan$B[res$sb$N-1,,,], 1, quantile, probs=c(0.025, 0.5, 0.975) ) )

    # densities of F in assessment year
    for (i in 1:3) plot(density(  res$stan$F[res$sb$N,i,,] ), xlim=c(0.05, 0.5), main="")
    ( qs = apply(  res$stan$F[res$sb$N,,,], 1, quantile, probs=c(0.025, 0.5, 0.975) ) )
    ( qs = apply(  res$stan$F[res$sb$N,,,], 1, mean ) )

    # densities of F in previous year
    for (i in 1:3) plot(density(  res$stan$F[res$sb$N-1,i,,] ), main="")
    ( qs = apply(  res$stan$F[res$sb$N-1,,,], 1, quantile, probs=c(0.025, 0.5, 0.975) ) )

    # F for table ---
    summary( res$stan$F, median)


}


if (p$surplusproduction_model$method=="jags"){
  
  res = surplusproduction_model( p=p, DS="jags" )  # jags method is deprecated
  # load( p$surplusproduction_model$fnres )

  outdir = p$surplusproduction_model$outdir

    # frequency density of key parameters
    figure.bugs( "K", res=res, fn=file.path(outdir, "K.density.png" ) )
    figure.bugs( "r", res=res, fn=file.path(outdir, "r.density.png" ) )
    figure.bugs( "q", res=res, fn=file.path(outdir, "q.density.png" ) ,xrange=c(0,2))
    figure.bugs( "FMSY", res=res, fn=file.path(outdir, "FMSY.density.png" ) )
    figure.bugs( "bo.sd", res=res, fn=file.path(outdir, "bo.sd.density.png" ) )
    figure.bugs( "bp.sd", res=res, fn=file.path(outdir, "bp.sd.density.png" ) )

    # timeseries
    figure.bugs( type="timeseries", vname="biomass", res=res, fn=file.path(outdir, "biomass.timeseries.png" ), save.plot=T )
    figure.bugs( type="timeseries", vname="fishingmortality", res=res, fn=file.path(outdir, "fishingmortality.timeseries.png" ) )

    # Harvest control rules
    figure.bugs( type="hcr", vname="default", res=res, fn=file.path(outdir, "hcr.default.png" ), save.plot=T  )
    figure.bugs( type="hcr", vname="simple", res=res, fn=file.path(outdir, "hcr.simple.png" ) )

    # diagnostics
    figure.bugs( type="diagnostic.production", res=res, fn=file.path(outdir, "diagnostic.production.png" ) )
    figure.bugs( type="diagnostic.errors", res=res, fn=file.path(outdir, "diagnostic.errors.png" ) )
    figure.bugs( type="diagnostic.phase", res=res, fn=file.path(outdir, "diagnostic.phase.png" ) )

    # densities of biomass estimates for the year.assessment
    for (i in 1:3) plot(density(y$B[res$sb$N,i,,] ), main="")
    ( qs = apply(  res$jags$B[res$sb$N,,,], 1, quantile, probs=c(0.025, 0.5, 0.975) ) )

    # densities of biomass estimates for the previous year
    for (i in 1:3) plot(density( res$jags$B[res$sb$N-1,i,,] ), main="")
    ( qs = apply(  res$jags$B[res$sb$N-1,,,], 1, quantile, probs=c(0.025, 0.5, 0.975) ) )

    # densities of F in assessment year
    for (i in 1:3) plot(density(  res$jags$F[res$sb$N,i,,] ), xlim=c(0.05, 0.5), main="")
    ( qs = apply(  res$jags$F[res$sb$N,,,], 1, quantile, probs=c(0.025, 0.5, 0.975) ) )
    ( qs = apply(  res$jags$F[res$sb$N,,,], 1, mean ) )

    # densities of F in previous year
    for (i in 1:3) plot(density(  res$jags$F[res$sb$N-1,i,,] ), main="")
    ( qs = apply(  res$jags$F[res$sb$N-1,,,], 1, quantile, probs=c(0.025, 0.5, 0.975) ) )

    # F for table ---
    summary( res$jags$F, median)

}



