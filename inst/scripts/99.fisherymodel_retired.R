
## ---------
#### final estimation of biomass via fishery models and associated figures and tables:

year.assessment = 2020   # NOTE: for 4X, the season 2019-2020 -> 2019


p = bio.snowcrab::snowcrab_parameters(
  project_class="carstm",
  assessment.years=2000:year.assessment,
 # carstm_model_label="nonseparable_simple",  # to choose alt carstm models
  tag="default"
)

p$fishery_model = fishery_model( DS = "logistic_parameters", p=p, tag=p$areal_units_type )
p$fishery_model$stancode = stan_initialize( stan_code=fishery_model( p=p, DS="stan_surplus_production" ) )


res = fishery_model( 
  DS="logistic_model", 
  p=p, 
  tag=p$areal_units_type,
  # from here down are params for cmdstanr::sample()
  data=p$fishery_model$standata, 
  iter_warmup = 2000,
  iter_sampling = 1000,
  seed = 123,
  chains = 4,
  parallel_chains = 4,  # The maximum number of MCMC chains to run in parallel.
  max_treedepth = 18,
  adapt_delta = 0.975,
  refresh = 500
)

# res = fishery_model( p=p, DS="samples", tag=p$areal_units_type )  # to get samples

  if (0) {
    fit = fishery_model( p=p, DS="fit", tag=p$areal_units_type )  # to get samples
  
    print( fit, max_rows=30 )
    # fit$summary("K", "r", "q")
    
    fit$cmdstan_diagnose()
    fit$cmdstan_summary()
  }


      # frequency density of key parameters
      fishery_model( DS="plot", vname="K", res=res )
      fishery_model( DS="plot", vname="r", res=res )
      fishery_model( DS="plot", vname="q", res=res, xrange=c(0.5, 2))
      fishery_model( DS="plot", vname="FMSY", res=res  )
      # fishery_model( DS="plot", vname="bosd", res=res  )
      # fishery_model( DS="plot", vname="bpsd", res=res  )

      # timeseries
      fishery_model( DS="plot", type="timeseries", vname="biomass", res=res  )
      fishery_model( DS="plot", type="timeseries", vname="fishingmortality", res=res)

      # Summary table of mean values for inclusion in document
      biomass.summary.table(x)

      # Harvest control rules
      fishery_model( DS="plot", type="hcr", vname="default", res=res  )
      fishery_model( DS="plot", type="hcr", vname="simple", res=res  )

      # diagnostics
      # fishery_model( DS="plot", type="diagnostic.errors", res=res )
      # fishery_model( DS="plot", type="diagnostic.phase", res=res  )



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
