

# Snow crab --- Areal unit modelling of habitat  -- no reliance upon stmv fields


# -------------------------------------------------
# Part 1 -- construct basic parameter list defining the main characteristics of the study
require(aegis) 

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
# 25 configs @ 5 min each, total time 2 hrs
  pS = substrate_carstm(p=p, DS="parameters_override" )
  M = substrate.db( p=pS, DS="aggregated_data", redo=TRUE )  # used for data matching/lookup in other aegis projects that use substrate
  M = substrate_carstm( p=pS, DS="carstm_inputs", redo=TRUE )  # will redo if not found
  M = NULL; gc()
  res = carstm_model( p=pS, M='substrate_carstm( p=pS, DS="carstm_inputs")', DS="redo", carstm_model_label="production"  )  # run model and obtain predictions

  if(0) {
    res = carstm_model( p=pS, DS="carstm_modelled", carstm_model_label="production"   ) # run model and obtain predictions
    fit = carstm_model( p=pS, DS="carstm_modelled_fit", carstm_model_label="production" )  # extract currently saved model fit
    summary(fit)

# Fixed effects:
#               mean    sd 0.025quant 0.5quant 0.975quant   mode kld
# (Intercept) -1.038 0.024     -1.086   -1.038     -0.991 -1.038   0

# Random effects:
#   Name	  Model
#     inla.group(z, method = "quantile", n = 13) RW2 model
#    auid BYM2 model

# Model hyperparameters:
#                                                           mean    sd 0.025quant 0.5quant 0.975quant  mode
# Precision for the lognormal observations                 1.519 0.007      1.506    1.519      1.533 1.519
# Precision for inla.group(z, method = "quantile", n = 13) 5.833 2.643      2.047    5.401     12.175 4.490
# Precision for auid                                       0.903 0.110      0.714    0.893      1.144 0.871
# Phi for auid                                             0.962 0.034      0.871    0.972      0.996 0.989

# Expected number of effective parameters(stdev): 195.14(0.203)
# Number of equivalent replicates : 493.07

# Deviance Information Criterion (DIC) ...............: 34384.88
# Deviance Information Criterion (DIC, saturated) ....: 96394.93
# Effective number of parameters .....................: 196.18

# Marginal log-Likelihood:  -17603.13

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
          + f( inla.group(z, method="quantile", n=13) ,  model="rw2", scale.model=TRUE, hyper=H$rw2)
          + f(auid, model="bym2", graph=sppoly@nb, scale.model=TRUE, constr=TRUE, hyper=H$bym2),
        family = "lognormal",
        data= M,
        control.compute=list(dic=TRUE, config=TRUE),
        control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
        control.predictor=list(compute=FALSE, link=1 ),
        control.fixed=H$fixed,  # priors for fixed effects, generic is ok
        verbose=TRUE
      ) ' )
    res = carstm_model( p=pS, M='substrate_carstm( p=pS, DS="carstm_inputs")', DS="redo", carstm_model_label=pS$carstm_model_label  )  # run model and obtain predictions
    fit = carstm_model( p=pS, DS="carstm_modelled_fit", carstm_model_label=pS$carstm_model_label )  # extract currently saved model fit

  }




# -------------------------------------------------
# Part 5 -- create covariate field for temperature
# ensure the data assimilation in temperature is first completed :: 01.temperature_data.R
# total: 30 min, 80 configs .. fast
# -------------------------------------------------
# Part 1 -- construct basic parameter list defining the main characteristics of the study

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
#     Pre = 3.31, Running = 2959, Post = 15.4, Total = 2978
# Fixed effects:
#              mean    sd 0.025quant 0.5quant 0.975quant  mode kld
# (Intercept) 4.575 0.342      3.889    4.574      5.267 4.573   0

# Random effects:
#   Name	  Model
#     dyri AR1 model
#    inla.group(z, method = "quantile", n = 13) RW2 model
#    auid BYM2 model

# Model hyperparameters:
#                                                           mean    sd 0.025quant 0.5quant 0.975quant  mode
# Precision for the Gaussian observations                  0.420 0.003      0.413    0.420      0.427 0.419
# Precision for dyri                                       3.119 1.143      1.417    2.943      5.853 2.612
# Rho for dyri                                             0.589 0.154      0.243    0.605      0.838 0.642
# Precision for inla.group(z, method = "quantile", n = 13) 1.137 0.340      0.602    1.094      1.923 1.012
# Precision for auid                                       0.407 0.021      0.365    0.407      0.449 0.409
# Phi for auid                                             0.449 0.032      0.386    0.449      0.514 0.448
# GroupRho for auid                                        0.828 0.011      0.806    0.828      0.851 0.826

# Expected number of effective parameters(stdev): 1976.69(32.72)
# Number of equivalent replicates : 17.76

# Deviance Information Criterion (DIC) ...............: 132097.21
# Deviance Information Criterion (DIC, saturated) ....: 37062.41
# Effective number of parameters .....................: 1977.19

# Marginal log-Likelihood:  -65482.69



    vn = paste(pT$variabletomodel, "predicted", sep=".")
    carstm_plot( p=pT, res=res, vn=vn, time_match=list(year="2000", dyear="0.85" ) )       # maps of some of the results

    vn = paste(pT$variabletomodel, "predicted_se", sep=".")
    carstm_plot( p=pT, res=res, vn=vn, time_match=list(year="2000", dyear="0.85" ) )       # maps of some of the results

    vn = paste(pT$variabletomodel, "random_auid_nonspatial", sep=".")
    carstm_plot( p=pT, res=res, vn=vn, time_match=list(year="2000"  ) )       # maps of some of the results
    vn = paste(pT$variabletomodel, "random_auid_spatial", sep=".")
    carstm_plot( p=pT, res=res, vn=vn, time_match=list(year="2000"  ) )       # maps of some of the results



    plot(fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE, single=TRUE )


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

    plot( fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE )

    summary(fit)
    vn = paste(pPC1$variabletomodel, "predicted", sep=".")
    carstm_plot( p=pPC1, res=res, vn=vn, time_match=list(year="2017" ), dyear="0.85" ) # maps of some of the results
    vn = paste(pPC1$variabletomodel, "predicted_se", sep=".")
    carstm_plot( p=pPC1, res=res, vn=vn, time_match=list(year="2017" ), dyear="0.85" ) # maps of some of the results

    vn = paste(pPC1$variabletomodel, "random_auid_nonspatial", sep=".")
    carstm_plot( p=pPC1, res=res, vn=vn, time_match=list(year="2017" ), dyear="0.85" )       # maps of some of the results , dyear="0.85"
    vn = paste(pPC1$variabletomodel, "random_auid_spatial", sep=".")
    carstm_plot( p=pPC1, res=res, vn=vn, time_match=list(year="2017" ), dyear="0.85" )       # maps of some of the results , dyear="0.85"

# Time used:
#     Pre = 2.67, Running = 667, Post = 4.98, Total = 675
# Fixed effects:
#              mean    sd 0.025quant 0.5quant 0.975quant  mode kld
# (Intercept) 0.143 0.015       0.11    0.143      0.173 0.144   0

# Random effects:
#   Name	  Model
#     dyri AR1 model
#    inla.group(t, method = "quantile", n = 13) RW2 model
#    inla.group(z, method = "quantile", n = 13) RW2 model
#    inla.group(substrate.grainsize, method = "quantile", n = 13) RW2 model
#    auid BYM2 model

# Model hyperparameters:
#                                                                                mean       sd 0.025quant 0.5quant
# Precision for the Gaussian observations                                    1.91e+02 2.88e+00    185.644 1.91e+02
# Precision for dyri                                                         1.03e+03 6.60e+02    273.266 8.62e+02
# Rho for dyri                                                               5.50e-02 2.47e-01     -0.430 5.90e-02
# Precision for inla.group(t, method = "quantile", n = 13)                   1.70e+04 1.45e+04   2813.038 1.30e+04
# Precision for inla.group(z, method = "quantile", n = 13)                   1.61e+02 7.50e+01     53.252 1.49e+02
# Precision for inla.group(substrate.grainsize, method = "quantile", n = 13) 7.29e+02 2.72e+03      9.795 1.95e+02
# Precision for auid                                                         1.71e+02 1.21e+01    147.609 1.71e+02
# Phi for auid                                                               5.19e-01 3.80e-02      0.442 5.20e-01
# GroupRho for auid                                                          9.19e-01 7.00e-03      0.905 9.19e-01
#                                                                            0.975quant     mode
# Precision for the Gaussian observations                                      1.97e+02  191.036
# Precision for dyri                                                           2.75e+03  612.450
# Rho for dyri                                                                 5.19e-01    0.069
# Precision for inla.group(t, method = "quantile", n = 13)                     5.54e+04 7245.080
# Precision for inla.group(z, method = "quantile", n = 13)                     3.41e+02  121.730
# Precision for inla.group(substrate.grainsize, method = "quantile", n = 13)   4.73e+03   19.393
# Precision for auid                                                           1.95e+02  170.771
# Phi for auid                                                                 5.92e-01    0.523
# GroupRho for auid                                                            9.33e-01    0.919

# Expected number of effective parameters(stdev): 1129.49(30.74)
# Number of equivalent replicates : 9.58

# Deviance Information Criterion (DIC) ...............: -24998.14
# Deviance Information Criterion (DIC, saturated) ....: 11953.32
# Effective number of parameters .....................: 1132.66

# Marginal log-Likelihood:  13670.07


  }





# -------------------------------------------------
# Part 7 -- create covariate field for species composition 2
# ensure that survey data is assimilated : bio.snowcrab::01snowcb_data.R, aegis.survey::01.surveys.data.R ,
# etc.30 min, 30 min
 p = bio.snowcrab::snowcrab_carstm( DS="parameters", assessment.years=1999:2018 )

  # misc run params adjustments here:
  p$inla_num.threads = 4
  p$inla_blas.num.threads = 4


  pPC2 = speciescomposition_carstm(p=p, DS="parameters_override", varnametomodel="pca2" )
  M = speciescomposition_carstm( p=pPC2, DS="carstm_inputs", redo=TRUE )  # will redo if not found
  M = NULL; gc()
  res = carstm_model( p=pPC2, M='speciescomposition_carstm( p=p, DS="carstm_inputs", varnametomodel="pca2" )', DS="redo" , carstm_model_label="production"  ) # run model and obtain predictions
  if (0) {
    res = carstm_model( p=pPC2, DS="carstm_modelled", carstm_model_label="production"   ) # run model and obtain predictions
    fit = carstm_model( p=pPC2, DS="carstm_modelled_fit", carstm_model_label="production"  )  # extract currently saved model fit
    summary(fit)
    plot( fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE )

    vn = paste(pPC2$variabletomodel, "predicted", sep=".")
    carstm_plot( p=pPC2, res=res, vn=vn, time_match=list(year="2017" ), dyear="0.85" )       # maps of some of the results
    vn = paste(pPC2$variabletomodel, "predicted_se", sep=".")
    carstm_plot( p=pPC2, res=res, vn=vn, time_match=list(year="2017" ), dyear="0.85" ) # maps of some of the results

    vn = paste(pPC2$variabletomodel, "random_auid_nonspatial", sep=".")
    carstm_plot( p=pPC2, res=res, vn=vn, time_match=list(year="2017" ), dyear="0.85" )       # maps of some of the results , dyear="0.85"
    vn = paste(pPC2$variabletomodel, "random_auid_spatial", sep=".")
    carstm_plot( p=pPC2, res=res, vn=vn, time_match=list(year="2017" ), dyear="0.85" )       # maps of some of the results , dyear="0.85"

  #   Time used:
  #     Pre = 4.86, Running = 1768, Post = 11.8, Total = 1785
  # Fixed effects:
    #             mean    sd 0.025quant 0.5quant 0.975quant mode kld
    # (Intercept) 0.02 0.021     -0.024     0.02      0.065 0.02   0

    # Random effects:
    #   Name	  Model
    #     dyri AR1 model
    #   inla.group(t, method = "quantile", n = 13) RW2 model
    #   inla.group(z, method = "quantile", n = 13) RW2 model
    #   inla.group(substrate.grainsize, method = "quantile", n = 13) RW2 model
    #   auid BYM2 model

    # Model hyperparameters:
    #                                                                               mean       sd 0.025quant 0.5quant
    # Precision for the Gaussian observations                                     255.139 3.86e+00    247.336  255.239
    # Precision for dyri                                                          972.573 5.21e+02    267.742  875.460
    # Rho for dyri                                                                  0.601 2.19e-01      0.123    0.627
    # Precision for inla.group(t, method = "quantile", n = 13)                   7071.547 1.04e+04    757.266 4032.628
    # Precision for inla.group(z, method = "quantile", n = 13)                    340.448 2.10e+02    109.684  285.481
    # Precision for inla.group(substrate.grainsize, method = "quantile", n = 13) 1267.298 1.65e+03     97.368  769.744
    # Precision for auid                                                          398.917 3.31e+01    335.265  398.669
    # Phi for auid                                                                  0.538 7.00e-02      0.383    0.546
    # GroupRho for auid                                                             0.914 8.00e-03      0.897    0.915
    #                                                                           0.975quant     mode
    # Precision for the Gaussian observations                                      2.63e+02  255.632
    # Precision for dyri                                                           2.25e+03  661.698
    # Rho for dyri                                                                 9.32e-01    0.750
    # Precision for inla.group(t, method = "quantile", n = 13)                     3.22e+04 1784.791
    # Precision for inla.group(z, method = "quantile", n = 13)                     8.90e+02  210.519
    # Precision for inla.group(substrate.grainsize, method = "quantile", n = 13)   5.51e+03  259.117
    # Precision for auid                                                           4.65e+02  399.531
    # Phi for auid                                                                 6.53e-01    0.578
    # GroupRho for auid                                                            9.30e-01    0.915

    # Expected number of effective parameters(stdev): 921.94(40.80)
    # Number of equivalent replicates : 11.74

    # Deviance Information Criterion (DIC) ...............: -28335.04
    # Deviance Information Criterion (DIC, saturated) ....: 11769.79
    # Effective number of parameters .....................: 928.42

    # Marginal log-Likelihood:  15491.21


  }



# finished covariates ... move onto abundance index estimation



# -------------------------------------------------
# Part 8 -- Snow crab anbundance -- main mode used for production
# 9 hrs 281 configs

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



# Time used:
#     Pre = 2.67, Running = 17238, Post = 3.15, Total = 17243
# Fixed effects:
#              mean    sd 0.025quant 0.5quant 0.975quant mode kld
# (Intercept) 5.483 0.117      5.244    5.486       5.71 5.49   0

# Random effects:
#   Name	  Model
#     dyri AR1 model
#    inla.group(t, method = "quantile", n = 13) RW2 model
#    inla.group(z, method = "quantile", n = 13) RW2 model
#    inla.group(substrate.grainsize, method = "quantile", n = 13) RW2 model
#    inla.group(pca1, method = "quantile", n = 13) RW2 model
#    inla.group(pca2, method = "quantile", n = 13) RW2 model
#    auid BYM2 model

# Model hyperparameters:
#                                                                              mean     sd 0.025quant 0.5quant 0.975quant
# Precision for dyri                                                         24.644 13.295      7.198   21.987     57.758
# Rho for dyri                                                               -0.010  0.132     -0.273   -0.008      0.245
# Precision for inla.group(t, method = "quantile", n = 13)                    6.078  2.716      2.367    5.559     12.835
# Precision for inla.group(z, method = "quantile", n = 13)                    4.846  2.067      2.060    4.431     10.025
# Precision for inla.group(substrate.grainsize, method = "quantile", n = 13)  0.243  0.081      0.125    0.229      0.439
# Precision for inla.group(pca1, method = "quantile", n = 13)                 7.210  3.278      2.842    6.549     15.399
# Precision for inla.group(pca2, method = "quantile", n = 13)                20.402  9.195      8.012   18.592     43.366
# Precision for auid                                                          0.430  0.035      0.364    0.429      0.502
# Phi for auid                                                                0.688  0.056      0.571    0.691      0.790
# GroupRho for auid                                                           0.817  0.013      0.790    0.818      0.843
#                                                                              mode
# Precision for dyri                                                         16.790
# Rho for dyri                                                               -0.001
# Precision for inla.group(t, method = "quantile", n = 13)                    4.642
# Precision for inla.group(z, method = "quantile", n = 13)                    3.734
# Precision for inla.group(substrate.grainsize, method = "quantile", n = 13)  0.204
# Precision for inla.group(pca1, method = "quantile", n = 13)                 5.425
# Precision for inla.group(pca2, method = "quantile", n = 13)                15.463
# Precision for auid                                                          0.427
# Phi for auid                                                                0.696
# GroupRho for auid                                                           0.818

# Expected number of effective parameters(stdev): 1920.79(15.10)
# Number of equivalent replicates : 3.96

# Deviance Information Criterion (DIC) ...............: 47009.33
# Deviance Information Criterion (DIC, saturated) ....: 28943.24
# Effective number of parameters .....................: 1894.56

# Marginal log-Likelihood:  -23611.26
# Posterior marginals for the linear predictor and
#  the fitted values are computed


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


    sppoly = areal_units( p=p )  # to reload
    M = snowcrab_carstm( p=p, DS="carstm_inputs" )
    M$yr = M$year  # req for meanweights

    # mean weight by auidxyear
    wgts = meanweights_by_arealunit(
      set=M[M$tag=="observations",],
      AUID=as.character( sppoly$AUID ),
      yrs=p$yrs,
      fillall=TRUE,
      annual_breakdown=TRUE,
      robustify_quantiles=c(0, 0.99)  # high upper bounds are more dangerous
    )


    carstm_summary( p=p, operation="compute", carstm_model_label=p$carstm_model_label, wgts=wgts )

    RES = carstm_summary(p=p, operation="load_timeseries", carstm_model_label=p$carstm_model_label  )

    bio = carstm_summary(p=p, operation="load_spacetime_biomass", carstm_model_label=p$carstm_model_label  )
    num = carstm_summary(p=p, operation="load_spacetime_number", carstm_model_label=p$carstm_model_label  )


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
