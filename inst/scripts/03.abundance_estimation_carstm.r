

# Snow crab --- Areal unit modelling of habitat  -- no reliance upon stmv fields


# -------------------------------------------------
# Part 1 -- construct basic parameter list defining the main characteristics of the study
  p = bio.snowcrab::snowcrab_carstm( DS="parameters", assessment.years=1999:2018 )

  # misc run params adjustments here:
  p$inla_num.threads = 4
  p$inla_blas.num.threads = 4


# -------------------------------------------------
# Part 2 -- polygon structure
  sppoly = areal_units( p=p )  # to reload
  plot(sppoly)
  spplot( sppoly, "au_sa_km2", main="AUID", sp.layout=p$coastLayout )
  if (0) {
    # create if not yet made
    MS = snowcrab.db( p=p, DS="biological_data"  )  # domain is  sse     # the underlying observations/data
    # ensure if polys exist and create if required
    for (au in c("cfanorth", "cfasouth", "cfa4x", "cfaall" )) plot(polygons_managementarea( species="snowcrab", au))
    sppoly = areal_units( p=p, areal_units_constraint=MS[, c("lon", "lat")], redo=TRUE )  # create constrained polygons with neighbourhood as an attribute
    coastLayout = aegis.coastline::coastline_layout( p=p, redo=TRUE )
    MS = NULL
  }


# -------------------------------------------------
# Part 3 -- create covariate field for bathymetry
# bathymetry -- ensure the data assimilation in bathymetry is first completed :: 01.bathymetry_data.R
# about 4.4 hrs to redo; 15 configs @ 0.5 hrs each

  pB = bathymetry_carstm( p=p, DS="parameters_override" )
  M = bathymetry.db( p=pB, DS="aggregated_data", redo=TRUE )
  M = bathymetry_carstm( p=pB, DS="carstm_inputs", redo=TRUE  ) # will redo if not found
  M = NULL; gc()
  res = carstm_model( p=pB, M='bathymetry_carstm( p=pB, DS="carstm_inputs" )', DS="redo", carstm_model_label="production"  ) # run model and obtain predictions

  if (0) {
    # to use a saved instance
    res = carstm_model( p=pB, DS="carstm_modelled", carstm_model_label="production" ) # run model and obtain predictions
    fit = carstm_model( p=pB, DS="carstm_modelled_fit", carstm_model_label="production" )  # extract currently saved model fit
    # maps of some of the results
    vn = paste(pB$variabletomodel, "predicted", sep=".")
    carstm_plot( p=pB, res=res, vn=vn )

    vn = paste(pB$variabletomodel, "predicted_se", sep=".")
    carstm_plot( p=pB, res=res, vn=vn )

    vn = paste(pB$variabletomodel, "random_auid_nonspatial", sep=".")
    carstm_plot( p=pB, res=res, vn=vn )

    vn = paste(pB$variabletomodel, "random_auid_spatial", sep=".")
    carstm_plot( p=pB, res=res, vn=vn )

    plot(fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE, single=TRUE )

    # Time used:
    #     Pre = 4.71, Running = 10928, Post = 4.51, Total = 10938
    # Fixed effects:
    #              mean    sd 0.025quant 0.5quant 0.975quant  mode kld
    # (Intercept) 7.883 0.001      7.882    7.883      7.885 7.883   0

    # Random effects:
    #   Name	  Model
    #     auid BYM2 model

    # Model hyperparameters:
    #                                             mean     sd 0.025quant 0.5quant 0.975quant    mode
    # Precision for the lognormal observations 752.050  3.198    745.879  752.003    758.450 751.849
    # Precision for auid                     286.200 23.004    239.367  287.143    329.336 291.721
    # Phi for auid                             0.979  0.021      0.922    0.985      0.999   0.996

    # Expected number of effective parameters(stdev): 698.31(1.78)
    # Number of equivalent replicates : 161.11

    # Deviance Information Criterion (DIC) ...............: 1348750.89
    # Deviance Information Criterion (DIC, saturated) ....: 113171.52
    # Effective number of parameters .....................: 698.85

    # Marginal log-Likelihood:  -674832.61
    # Posterior marginals for the linear predictor and
    #  the fitted values are computed

  }




# -------------------------------------------------
# Part 4 -- create covariate field for  substrate
# ensure the data assimilation in substrate is first completed :: 01.substrate_data.R
# 25 configs @ 2 hrs each, total time 32 hrs
  pS = substrate_carstm(p=p, DS="parameters_override" )
  M = substrate.db( p=pS, DS="aggregated_data", redo=TRUE )  # used for data matching/lookup in other aegis projects that use substrate
  M = substrate_carstm( p=pS, DS="carstm_inputs", redo=TRUE )  # will redo if not found
  M = NULL; gc()
  res = carstm_model( p=pS, M='substrate_carstm( p=pS, DS="carstm_inputs")', DS="redo", carstm_model_label="production"  )  # run model and obtain predictions

  if(0) {
    res = carstm_model( p=pS, DS="carstm_modelled", carstm_model_label="production"   ) # run model and obtain predictions
    fit = carstm_model( p=pS, DS="carstm_modelled_fit", carstm_model_label="production" )  # extract currently saved model fit
    summary(fit)
    # Model hyperparameters:
    #                                           mean    sd 0.025quant 0.5quant 0.975quant  mode
    # Precision for the lognormal observations 1.405 0.006      1.392    1.405      1.417 1.405
    # Precision for zi                         4.318 2.467      1.105    3.822     10.495 2.771
    # Precision for auid                       0.840 0.124      0.652    0.820      1.134 0.770
    # Phi for auid                             0.959 0.034      0.867    0.969      0.995 0.986

    # Expected number of effective parameters(stdev): 191.31(0.21)
    # Number of equivalent replicates : 502.93

    # Deviance Information Criterion (DIC) ...............: 41922.78
    # Deviance Information Criterion (DIC, saturated) ....: 96422.26
    # Effective number of parameters .....................: 192.34

    # Marginal log-Likelihood:  -21357.76

    vn = paste(pS$variabletomodel, "predicted", sep=".")
    carstm_plot( p=pS, res=res, vn=vn ) # maps of some of the results

    vn = paste(pS$variabletomodel, "predicted_se", sep=".")
    carstm_plot( p=pS, res=res, vn=vn )

    vn = paste(pS$variabletomodel, "random_auid_nonspatial", sep=".")
    carstm_plot( p=pS, res=res, vn=vn )

    vn = paste(pS$variabletomodel, "random_auid_spatial", sep=".")
    carstm_plot( p=pS, res=res, vn=vn )



    plot(fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE, single=TRUE )

    #------------
    # 6567.727 sec = 1.8hrs

    pS$carstm_model_label = "production_inla.group_quantile_25"
    # DIC:
    #   Mean of Deviance ................. 33840.1
    #   Deviance at Mean ................. 33636.8
    #   Effective number of parameters ... 203.318
    #   DIC .............................. 34043.4
    # DIC (Saturated):
    #   Mean of Deviance ................. 95897.9
    #   Deviance at Mean ................. 95694.6
    #   Effective number of parameters ... 203.318
    #   DIC .............................. 96101.3




    pS$carstm_model_label = "production_inla.group_quantile_20"
    # DIC:
    # 	Mean of Deviance ................. 33995.4
    # 	Deviance at Mean ................. 33794.6
    # 	Effective number of parameters ... 200.881
    # 	DIC .............................. 34196.3
    # DIC (Saturated):
    # 	Mean of Deviance ................. 96187.3
    # 	Deviance at Mean ................. 95986.4
    # 	Effective number of parameters ... 200.881
    # 	DIC .............................. 96388.2


    pS$carstm_modelcall = paste('
      inla(
        formula =', pS$variabletomodel, ' ~ 1
          + f( inla.group(z, method="quantile", n=20) ,  model="rw2", scale.model=TRUE, hyper=H$rw2)
          + f(auid, model="bym2", graph=sppoly@nb, scale.model=TRUE, constr=TRUE, hyper=H$bym2),
        family = "lognormal",
        data= M,
        control.compute=list(dic=TRUE, config=TRUE),
        control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
        control.predictor=list(compute=FALSE, link=1 ),
        control.fixed=H$fixed,  # priors for fixed effects, generic is ok
        # control.inla=list( strategy="laplace", cutoff=1e-6, correct=TRUE, correct.verbose=FALSE ),  # extra work to get tails
        # control.inla = list( h=1e-6, tolerance=1e-12), # increase in case values are too close to zero
        # control.mode = list( restart=TRUE, result=RES ), # restart from previous estimates
        # control.inla = list(cmin = 0 ),
        # control.inla=list( strategy="laplace", cutoff=1e-6, correct=TRUE, correct.verbose=FALSE ),
        # control.inla = list( h=1e-6, tolerance=1e-12), # increase in case values are too close to zero
        # control.inla = list(h=1e-3, tolerance=1e-9, cmin=0), # restart=3), # restart a few times in case posteriors are poorly defined
        # control.mode = list( restart=TRUE, result=RES ), # restart from previous estimates
        verbose=TRUE
      ) ' )
    res = carstm_model( p=pS, M='substrate_carstm( p=pS, DS="carstm_inputs")', DS="redo", carstm_model_label=pS$carstm_model_label  )  # run model and obtain predictions
    fit = carstm_model( p=pS, DS="carstm_modelled_fit", carstm_model_label=pS$carstm_model_label )  # extract currently saved model fit

  }




# -------------------------------------------------
# Part 5 -- create covariate field for temperature
# ensure the data assimilation in temperature is first completed :: 01.temperature_data.R
# total: 30 min, 80 configs .. fast
  pT = temperature_carstm(p=p, DS="parameters_override" )
  M = temperature.db( p=pT, DS="aggregated_data", redo=TRUE )  #  used for data matching/lookup in other aegis projects that use temperature
  M = temperature_carstm( p=pT, DS="carstm_inputs", redo=TRUE )  # will redo if not found
  M = NULL; gc()
  res = carstm_model( p=pT, M='temperature_carstm( p=pT, DS="carstm_inputs")', DS="redo", carstm_model_label="production"  ) # run model and obtain predictions
  if (0) {
    res = carstm_model( p=pT, DS="carstm_modelled", carstm_model_label="production" ) # run model and obtain predictions
    fit = carstm_model(  p=pT, DS="carstm_modelled_fit", carstm_model_label="production" )  # extract currently saved model fit
    summary(fit)

# Time used:
#     Pre = 3.23, Running = 3680, Post = 16.4, Total = 3700
# Fixed effects:
#              mean    sd 0.025quant 0.5quant 0.975quant  mode kld
# (Intercept) 4.454 0.388      3.682    4.453      5.228 4.451   0

# Random effects:
#   Name	  Model
#     year_factor AR1 model
#    dyri AR1 model
#    inla.group(z, method = "quantile", n = 25) RW2 model
#    auid BYM2 model

# Model hyperparameters:
#                                                           mean    sd 0.025quant 0.5quant 0.975quant  mode
# Precision for the Gaussian observations                  0.429 0.003      0.422    0.429      0.436 0.429
# Precision for year_factor                                2.878 0.841      1.549    2.772      4.821 2.571
# Rho for year_factor                                      0.386 0.155      0.059    0.395      0.663 0.411
# Precision for dyri                                       3.092 1.126      1.383    2.934      5.749 2.622
# Rho for dyri                                             0.601 0.147      0.268    0.617      0.838 0.653
# Precision for inla.group(z, method = "quantile", n = 25) 0.128 0.028      0.079    0.126      0.190 0.122
# Precision for auid                                       0.410 0.033      0.352    0.408      0.482 0.401
# Phi for auid                                             0.996 0.005      0.983    0.998      1.000 1.000
# GroupRho for auid                                        0.723 0.025      0.670    0.725      0.767 0.731

# Expected number of effective parameters(stdev): 1563.56(31.92)
# Number of equivalent replicates : 22.46

# Deviance Information Criterion (DIC) ...............: 130925.13
# Deviance Information Criterion (DIC, saturated) ....: 36676.07
# Effective number of parameters .....................: 1565.21

# Marginal log-Likelihood:  -64612.59
# Posterior marginals for the linear predictor and
#  the fitted values are computed


  # Time used:
  #     Pre = 2.67, Running = 559, Post = 2.69, Total = 564
  # Fixed effects:
  #             mean    sd 0.025quant 0.5quant 0.975quant  mode kld
  # (Intercept) 5.127 0.309       4.52    5.127      5.734 5.127   0

  # Random effects:
  #   Name	  Model
  #     year_factor AR1 model
  #   dyri RW2 model
  #   zi RW2 model
  #   auid BYM2 model

  # Model hyperparameters:
  #                                         mean    sd 0.025quant 0.5quant 0.975quant  mode
  # Precision for the Gaussian observations 0.421 0.000      0.421    0.421   4.52e-01 0.421
  # Precision for year_factor               1.325 0.057      1.181    1.308   2.29e+05 1.344
  # Rho for year_factor                     0.458 0.022      0.398    0.452   1.00e+00 0.466
  # Precision for dyri                      1.131 0.061      1.005    1.121   1.30e+00 1.085
  # Precision for zi                        0.001 0.000      0.001    0.001   1.06e-01 0.001
  # Precision for auid                    0.220 0.001      0.214    0.220   2.23e-01 0.220
  # Phi for auid                          0.993 0.000      0.992    0.993   1.00e+00 0.994

  # Expected number of effective parameters(stdev): 215.49(0.00)
  # Number of equivalent replicates : 146.64

  # Deviance Information Criterion (DIC) ...............: 128672.37
  # Deviance Information Criterion (DIC, saturated) ....: 43266.47
  # Effective number of parameters .....................: 215.48

  # Marginal log-Likelihood:  -64703.07
  # Posterior marginals for the linear predictor and
  # the fitted values are computed



  # Time used:
  #     Pre = 1.5, Running = 848, Post = 2.71, Total = 852
  # Fixed effects:
  #             mean    sd 0.025quant 0.5quant 0.975quant  mode kld
  # (Intercept) 5.101 0.206      4.698    5.101      5.504 5.101   0

  # Random effects:
  #   Name	  Model
  #     year_factor AR1 model
  #   dyri RW2 model
  #   zi RW2 model
  #   auid BYM2 model

  # Model hyperparameters:
  #                                         mean    sd 0.025quant 0.5quant 0.975quant  mode
  # Precision for the Gaussian observations 0.355 0.003      0.349    0.355      0.361 0.355
  # Precision for year_factor               2.417 0.850      1.036    2.335      4.306 2.137
  # Rho for year_factor                     0.440 0.176      0.090    0.444      0.761 0.442
  # Precision for dyri                      0.556 0.305      0.167    0.492      1.331 0.375
  # Precision for zi                        0.001 0.000      0.000    0.001      0.001 0.001
  # Precision for auid                    0.328 0.027      0.280    0.325      0.387 0.319
  # Phi for auid                          0.992 0.010      0.966    0.995      1.000 1.000
  # GroupRho for auid                     0.843 0.015      0.810    0.844      0.868 0.848

  # Expected number of effective parameters(stdev): 1216.66(0.00)
  # Number of equivalent replicates : 25.97

  # Deviance Information Criterion (DIC) ...............: 123610.43
  # Deviance Information Criterion (DIC, saturated) ....: 32822.21
  # Effective number of parameters .....................: 1216.63

  # Marginal log-Likelihood:  -60756.52



    vn = paste(pT$variabletomodel, "predicted", sep=".")
    carstm_plot( p=pT, res=res, vn=vn, time_match=list(year="2000", dyear="0.85" ) )       # maps of some of the results

    vn = paste(pT$variabletomodel, "predicted_se", sep=".")
    carstm_plot( p=pT, res=res, vn=vn, time_match=list(year="2000", dyear="0.85" ) )       # maps of some of the results

    vn = paste(pT$variabletomodel, "random_auid_nonspatial", sep=".")
    carstm_plot( p=pT, res=res, vn=vn, time_match=list(year="2000"  ) )       # maps of some of the results
    vn = paste(pT$variabletomodel, "random_auid_spatial", sep=".")
    carstm_plot( p=pT, res=res, vn=vn, time_match=list(year="2000"  ) )       # maps of some of the results



    plot(fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE, single=TRUE )


        pT$carstm_model_label = "production"
        pT$carstm_modelcall = paste('
          inla(
            formula = ', pT$variabletomodel, ' ~ 1
              + f( year_factor, model="ar1", hyper=H$ar1 )
              + f( dyri, model="ar1", scale.model=TRUE, hyper=H$ar1 )
              + f( inla.group( z, method="quantile", n=25 ), model="rw2", scale.model=TRUE, hyper=H$rw2)
              + f( auid, model="bym2", graph=sppoly@nb, group=year_factor, scale.model=TRUE, constr=TRUE, hyper=H$bym2),
            family = "normal",
            data= M,
            control.compute=list(dic=TRUE, config=TRUE),
            control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
            control.predictor=list(compute=FALSE, link=1 ),
            control.fixed=H$fixed,  # priors for fixed effects, generic is ok
            # control.inla=list( strategy="laplace", cutoff=1e-6, correct=TRUE, correct.verbose=FALSE ),
            # control.inla = list( h=1e-6, tolerance=1e-12), # increase in case values are too close to zero
            # control.mode = list( restart=TRUE, result=RES ), # restart from previous estimates
            # control.inla = list(cmin = 0 ),
            # control.inla=list( strategy="laplace", cutoff=1e-6, correct=TRUE, correct.verbose=FALSE ),
            # control.inla = list( h=1e-6, tolerance=1e-12), # increase in case values are too close to zero
            # control.inla = list(h=1e-3, tolerance=1e-9, cmin=0), # restart=3), # restart a few times in case posteriors are poorly defined
            # control.mode = list( restart=TRUE, result=RES ), # restart from previous estimates
            verbose=TRUE
          ) ' )

        #  + f(tiyr2, model="seasonal", season.length=10 )
        #  + f(dyear, model="ar1", hyper=H$ar1 )
        #  + f(seasonal, model="seasonal", season.length=', pT$n.season, ', scale.model=TRUE )  # using seasonal effect is not recommended as it is not smoothed well .. rw2 is better



  }





# -------------------------------------------------
# Part 6 -- create covariate field for species composition 1
# ensure that survey data is assimilated : bio.snowcrab::01snowcb_data.R, aegis.survey::01.surveys.data.R , etc.
# 30 min, 150 configs
  pPC1 = speciescomposition_carstm(p=p, DS="parameters_override", varnametomodel="pca1" )
  M = speciescomposition_carstm( p=pPC1, DS="carstm_inputs", varnametomodel="pca1", redo=TRUE )  # will redo if not found
  M = NULL; gc()
  res = carstm_model( p=pPC1, M='speciescomposition_carstm( p=p, DS="carstm_inputs", varnametomodel="pca1" )', DS="redo", carstm_model_label="production"   ) # run model and obtain predictions
  if (0) {
    res = carstm_model( p=pPC1, DS="carstm_modelled", carstm_model_label="production"  ) # run model and obtain predictions
    fit = carstm_model( p=pPC1, DS="carstm_modelled_fit", carstm_model_label="production" )  # extract currently saved model fit
    summary(fit)
    vn = paste(pPC1$variabletomodel, "predicted", sep=".")
    carstm_plot( p=pPC1, res=res, vn=vn, time_match=list(year="2000" ), dyear="0.85" ) # maps of some of the results
    vn = paste(pPC1$variabletomodel, "random_auid_nonspatial", sep=".")
    carstm_plot( p=pPC1, res=res, vn=vn, time_match=list(year="2000" ), dyear="0.85" )       # maps of some of the results , dyear="0.85"
    vn = paste(pPC1$variabletomodel, "random_auid_spatial", sep=".")
    carstm_plot( p=pPC1, res=res, vn=vn, time_match=list(year="2000" ), dyear="0.85" )       # maps of some of the results , dyear="0.85"
  }





# -------------------------------------------------
# Part 7 -- create covariate field for species composition 2
# ensure that survey data is assimilated : bio.snowcrab::01snowcb_data.R, aegis.survey::01.surveys.data.R ,
# etc.30 min, 30 min
  pPC2 = speciescomposition_carstm(p=p, DS="parameters_override", varnametomodel="pca2" )
  M = speciescomposition_carstm( p=pPC2, DS="carstm_inputs", redo=TRUE )  # will redo if not found
  M = NULL; gc()
  res = carstm_model( p=pPC2, M='speciescomposition_carstm( p=p, DS="carstm_inputs", varnametomodel="pca2" )', DS="redo" , carstm_model_label="production"  ) # run model and obtain predictions
  if (0) {
    res = carstm_model( p=pPC2, DS="carstm_modelled", carstm_model_label="production"   ) # run model and obtain predictions
    fit = carstm_model( p=pPC2, DS="carstm_modelled_fit", carstm_model_label="production"  )  # extract currently saved model fit
    summary(fit)
    vn = paste(pPC2$variabletomodel, "predicted", sep=".")
    carstm_plot( p=pPC2, res=res, vn=vn, time_match=list(year="2000" ), dyear="0.85" )       # maps of some of the results
    vn = paste(pPC2$variabletomodel, "random_auid_nonspatial", sep=".")
    carstm_plot( p=pPC2, res=res, vn=vn, time_match=list(year="2000" ), dyear="0.85" )       # maps of some of the results , dyear="0.85"
    vn = paste(pPC2$variabletomodel, "random_auid_spatial", sep=".")
    carstm_plot( p=pPC2, res=res, vn=vn, time_match=list(year="2000" ), dyear="0.85" )       # maps of some of the results , dyear="0.85"
  }



# finished covariates ... move onto abundance index estimation



# -------------------------------------------------
# Part 8 -- Snow crab anbundance -- main mode used for production
  M = snowcrab_carstm( p=p, DS="carstm_inputs", redo=TRUE )  # will redo if not found
  M = NULL; gc()
  res = carstm_model( p=p, M='snowcrab_carstm( p=p, DS="carstm_inputs" )' ) # 151 configs and long optim .. 19 hrs

  # m = get("inla.models", INLA:::inla.get.inlaEnv())
	# m$latent$rw2$min.diff = NULL
	# assign("inla.models", m, INLA:::inla.get.inlaEnv())


    # extract results
    res = carstm_model( p=p, DS="carstm_modelled" ) # to load currently saved res
    fit = carstm_model( p=p, DS="carstm_modelled_fit" )  # extract currently saved model fit
    summary(fit)
    vn = paste(p$variabletomodel, "predicted", sep=".")
    carstm_plot( p=p, res=res, vn=vn, time_match=list(year="2000" ) )     # maps of some of the results

    plot(fit)
    plot(fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE, single=TRUE )
    s = summary(fit)
    s$dic$dic
    s$dic$p.eff

    # maps of some of the results
    vn = paste(p$variabletomodel, "random_sample_iid", sep=".")
    if (exists(vn, res)) carstm_plot( p=p, res=res, vn=vn, time_match=list(year="1950", dyear="0.05") )

    vn = paste(p$variabletomodel, "random_auid_nonspatial", sep=".")
    if (exists(vn, res)) {
      res_dim = dim( res[[vn]] )
      if (res_dim == 1 ) time_match = NULL
      if (res_dim == 2 ) time_match = list(year="2000")
      if (res_dim == 3 ) time_match = list(year="2000", dyear="0.85" )
      carstm_plot( p=p, res=res, vn=vn, time_match=time_match )
    }

    vn = paste(p$variabletomodel, "random_auid_spatial", sep=".")
    if (exists(vn, res)) {
      res_dim = dim( res[[vn]] )
      if (res_dim == 1 ) time_match = NULL
      if (res_dim == 2 ) time_match = list(year="2000")
      if (res_dim == 3 ) time_match = list(year="2000", dyear="0.85" )
      carstm_plot( p=p, res=res, vn=vn, time_match=time_match )
    }



    snowcrab_abundance_index( p=p, operation="compute", carstm_model_label=p$carstm_model_label )

    RES = snowcrab_abundance_index(p=p, operation="load_timeseries", carstm_model_label=p$carstm_model_label  )

    bio = snowcrab_abundance_index(p=p, operation="load_spacetime_biomass", carstm_model_label=p$carstm_model_label  )
    num = snowcrab_abundance_index(p=p, operation="load_spacetime_number", carstm_model_label=p$carstm_model_label  )
    wt = snowcrab_abundance_index(p=p, operation="load_spacetime_weights", carstm_model_label=p$carstm_model_label  )


    plot( cfaall ~ yrs, data=RES, lty=1, lwd=2.5, col="red", type="b")
    plot( cfasouth ~ yrs, data=RES, lty=1, lwd=2.5, col="red", type="b")
    plot( cfanorth ~ yrs, data=RES, lty=1, lwd=2.5, col="red", type="b")
    plot( cfa4x ~ yrs, data=RES, lty=1, lwd=2.5, col="red", type="b")

    p$boundingbox = list( xlim=p$corners$lon, ylim=p$corners$lat) # bounding box for plots using spplot

    p$coastLayout = aegis.coastline::coastline_layout(p=p)
    p$mypalette=RColorBrewer::brewer.pal(9, "YlOrRd")
      sppoly = areal_units( p=p )  # to reload


    # map it ..mean density
    vn = "pred"
    sppoly@data[,vn] = bio[,"2017"]
    brks = interval_break(X= sppoly[[vn]], n=length(p$mypalette), style="quantile")
    spplot( sppoly, vn, col.regions=p$mypalette, main=vn, at=brks, sp.layout=p$coastLayout, col="transparent" )


    plot( fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE )
    plot( fit$marginals.hyperpar$"Phi for auid", type="l")  # posterior distribution of phi nonspatial dominates
    plot( fit$marginals.hyperpar$"Precision for auid", type="l")
    plot( fit$marginals.hyperpar$"Precision for setno", type="l")





# end
