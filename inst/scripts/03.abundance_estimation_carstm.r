

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
# (Intercept) 4.521 0.374      3.783     4.52       5.26 4.519   0

# Random effects:
#   Name	  Model
#     year_factor AR1 model
#    dyri AR1 model
#    inla.group(z, method = "quantile", n = 13) RW2 model
#    auid BYM2 model

# Model hyperparameters:
#                                                           mean    sd 0.025quant 0.5quant 0.975quant  mode
# Precision for the Gaussian observations                  0.426 0.003      0.419    0.426      0.432 0.426
# Precision for year_factor                                3.011 0.847      1.623    2.924      4.931 2.755
# Rho for year_factor                                      0.377 0.150      0.074    0.380      0.654 0.380
# Precision for dyri                                       3.127 1.127      1.382    2.982      5.743 2.683
# Rho for dyri                                             0.607 0.143      0.283    0.623      0.837 0.658
# Precision for inla.group(z, method = "quantile", n = 13) 1.130 0.337      0.603    1.086      1.917 1.002
# Precision for auid                                       0.412 0.034      0.354    0.410      0.486 0.402
# Phi for auid                                             1.000 0.000      1.000    1.000      1.000   NaN
# GroupRho for auid                                        0.722 0.025      0.667    0.724      0.766 0.731

# Expected number of effective parameters(stdev): 1541.25(23.03)
# Number of equivalent replicates : 22.78

# Deviance Information Criterion (DIC) ...............: 131164.54
# Deviance Information Criterion (DIC, saturated) ....: 36672.25
# Effective number of parameters .....................: 1543.11

# Marginal log-Likelihood:  -64636.41


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
              + f( inla.group( z, method="quantile", n=13 ), model="rw2", scale.model=TRUE, hyper=H$rw2)
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

#  Time used:
#     Pre = 2.81, Running = 1623, Post = 8.2, Total = 1634
# Fixed effects:
#              mean    sd 0.025quant 0.5quant 0.975quant  mode kld
# (Intercept) 0.144 0.034      0.076    0.145      0.202 0.146   0

# Random effects:
#   Name	  Model
#     year_factor AR1 model
#    dyri AR1 model
#    inla.group(t, method = "quantile", n = 13) RW2 model
#    inla.group(z, method = "quantile", n = 13) RW2 model
#    inla.group(substrate.grainsize, method = "quantile", n = 13) RW2 model
#    auid BYM2 model

# Model hyperparameters:
#                                                                                mean       sd 0.025quant 0.5quant
# Precision for the Gaussian observations                                    1.92e+02 2.84e+00    186.640 1.92e+02
# Precision for year_factor                                                  7.80e+02 3.94e+02    192.498 7.27e+02
# Rho for year_factor                                                        8.01e-01 1.14e-01      0.550 8.17e-01
# Precision for dyri                                                         9.03e+02 5.06e+02    247.972 8.00e+02
# Rho for dyri                                                               8.10e-02 2.35e-01     -0.352 7.20e-02
# Precision for inla.group(t, method = "quantile", n = 13)                   2.85e+04 2.94e+04   3647.288 1.99e+04
# Precision for inla.group(z, method = "quantile", n = 13)                   1.61e+02 7.49e+01     56.493 1.48e+02
# Precision for inla.group(substrate.grainsize, method = "quantile", n = 13) 2.95e+02 6.06e+02     11.110 1.30e+02
# Precision for auid                                                         1.77e+02 1.72e+01    143.759 1.77e+02
# Phi for auid                                                               9.94e-01 8.00e-03      0.974 9.97e-01
# GroupRho for auid                                                          8.68e-01 1.80e-02      0.834 8.68e-01
#                                                                            0.975quant     mode
# Precision for the Gaussian observations                                      1.98e+02  191.903
# Precision for year_factor                                                    1.67e+03  539.805
# Rho for year_factor                                                          9.67e-01    0.889
# Precision for dyri                                                           2.17e+03  594.757
# Rho for dyri                                                                 5.49e-01    0.017
# Precision for inla.group(t, method = "quantile", n = 13)                     1.06e+05 9532.957
# Precision for inla.group(z, method = "quantile", n = 13)                     3.45e+02  122.218
# Precision for inla.group(substrate.grainsize, method = "quantile", n = 13)   1.62e+03   27.315
# Precision for auid                                                           2.11e+02  176.982
# Phi for auid                                                                 1.00e+00    1.000
# GroupRho for auid                                                            9.02e-01    0.867

# Expected number of effective parameters(stdev): 769.85(33.57)
# Number of equivalent replicates : 14.06

# Deviance Information Criterion (DIC) ...............: -25408.95
# Deviance Information Criterion (DIC, saturated) ....: 11591.30
# Effective number of parameters .....................: 772.98

# Marginal log-Likelihood:  14060.02


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

  #   Time used:
  #     Pre = 4.86, Running = 1768, Post = 11.8, Total = 1785
  # Fixed effects:
  #             mean    sd 0.025quant 0.5quant 0.975quant  mode   kld
  # (Intercept) 0.013 0.032      -0.05    0.014      0.072 0.014 0.007

  # Random effects:
  #   Name	  Model
  #     year_factor AR1 model
  #   dyri AR1 model
  #   inla.group(t, method = "quantile", n = 13) RW2 model
  #   inla.group(z, method = "quantile", n = 13) RW2 model
  #   inla.group(substrate.grainsize, method = "quantile", n = 13) RW2 model
  #   auid BYM2 model

  # Model hyperparameters:
  #                                                                               mean       sd 0.025quant 0.5quant
  # Precision for the Gaussian observations                                     249.493    3.911    242.160  249.350
  # Precision for year_factor                                                  1126.871  831.853    104.778  933.654
  # Rho for year_factor                                                           0.906    0.077      0.716    0.927
  # Precision for dyri                                                          395.689  245.945     47.181  353.140
  # Rho for dyri                                                                  0.908    0.065      0.760    0.922
  # Precision for inla.group(t, method = "quantile", n = 13)                   3388.886 3692.945    500.318 2286.034
  # Precision for inla.group(z, method = "quantile", n = 13)                    342.246  169.399    116.458  308.972
  # Precision for inla.group(substrate.grainsize, method = "quantile", n = 13) 1607.268 2600.176    125.155  856.528
  # Precision for auid                                                          334.180   37.633    269.886  330.624
  # Phi for auid                                                                  0.951    0.037      0.851    0.961
  # GroupRho for auid                                                             0.759    0.037      0.673    0.764
  #                                                                           0.975quant     mode
  # Precision for the Gaussian observations                                      2.58e+02  248.886
  # Precision for year_factor                                                    3.14e+03  313.934
  # Rho for year_factor                                                          9.95e-01    0.987
  # Precision for dyri                                                           9.75e+02  156.802
  # Rho for dyri                                                                 9.94e-01    0.982
  # Precision for inla.group(t, method = "quantile", n = 13)                     1.30e+04 1188.188
  # Precision for inla.group(z, method = "quantile", n = 13)                     7.64e+02  247.987
  # Precision for inla.group(substrate.grainsize, method = "quantile", n = 13)   7.76e+03  312.351
  # Precision for auid                                                           4.17e+02  322.041
  # Phi for auid                                                                 9.92e-01    0.977
  # GroupRho for auid                                                            8.17e-01    0.779

  # Expected number of effective parameters(stdev): 859.67(27.66)
  # Number of equivalent replicates : 12.59

  # Deviance Information Criterion (DIC) ...............: -28164.77
  # Deviance Information Criterion (DIC, saturated) ....: 11648.63
  # Effective number of parameters .....................: 859.56

  # Marginal log-Likelihood:  15463.31

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

    #     Time used:
    #     Pre = 5.46, Running = 36557, Post = 14.5, Total = 36577
    # Fixed effects:
    #             mean    sd 0.025quant 0.5quant 0.975quant  mode kld
    # (Intercept) 5.56 0.204      5.166    5.558      5.971 5.556   0

    # Random effects:
    #   Name	  Model
    #     year_factor AR1 model
    #   dyri AR1 model
    #   inla.group(t, method = "quantile", n = 13) RW2 model
    #   inla.group(z, method = "quantile", n = 13) RW2 model
    #   inla.group(substrate.grainsize, method = "quantile", n = 13) RW2 model
    #   inla.group(pca1, method = "quantile", n = 13) RW2 model
    #   inla.group(pca2, method = "quantile", n = 13) RW2 model
    #   auid BYM2 model

    # Model hyperparameters:
    #                                                                             mean     sd 0.025quant 0.5quant 0.975quant
    # Precision for year_factor                                                  14.543  6.575      4.947   13.534     30.228
    # Rho for year_factor                                                         0.739  0.121      0.460    0.755      0.923
    # Precision for dyri                                                         23.107 13.229      6.935   20.089     56.893
    # Rho for dyri                                                               -0.027  0.179     -0.382   -0.022      0.312
    # Precision for inla.group(t, method = "quantile", n = 13)                    4.413  1.857      1.827    4.071      9.001
    # Precision for inla.group(z, method = "quantile", n = 13)                    5.583  2.341      2.314    5.155     11.360
    # Precision for inla.group(substrate.grainsize, method = "quantile", n = 13)  0.251  0.084      0.126    0.238      0.452
    # Precision for inla.group(pca1, method = "quantile", n = 13)                 7.630  3.286      3.051    7.031     15.734
    # Precision for inla.group(pca2, method = "quantile", n = 13)                18.666  8.216      7.176   17.198     38.801
    # Precision for auid                                                          0.408  0.041      0.332    0.407      0.494
    # Phi for auid                                                                0.764  0.049      0.660    0.766      0.852
    # GroupRho for auid                                                           0.682  0.029      0.623    0.682      0.737
    #                                                                             mode
    # Precision for year_factor                                                  11.247
    # Rho for year_factor                                                         0.798
    # Precision for dyri                                                         15.118
    # Rho for dyri                                                               -0.004
    # Precision for inla.group(t, method = "quantile", n = 13)                    3.463
    # Precision for inla.group(z, method = "quantile", n = 13)                    4.391
    # Precision for inla.group(substrate.grainsize, method = "quantile", n = 13)  0.214
    # Precision for inla.group(pca1, method = "quantile", n = 13)                 5.948
    # Precision for inla.group(pca2, method = "quantile", n = 13)                14.457
    # Precision for auid                                                          0.404
    # Phi for auid                                                                0.771
    # GroupRho for auid                                                           0.682

    # Expected number of effective parameters(stdev): 1994.79(12.63)
    # Number of equivalent replicates : 3.81

    # Deviance Information Criterion (DIC) ...............: 47097.26
    # Deviance Information Criterion (DIC, saturated) ....: 29031.17
    # Effective number of parameters .....................: 1966.99

    # Marginal log-Likelihood:  -23623.31

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
