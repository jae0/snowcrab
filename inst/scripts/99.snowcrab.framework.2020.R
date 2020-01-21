

# Snow crab --- Areal unit modelling of habitat


# -------------------------------------------------
# Initiate -- construct basic parameter list defining the main characteristics of the study
# expects bio.snowcrab/inst/scripts/03.abundance_estimation_carstm.r to have been run ..
# we are running simpler variations of that model here

  p = bio.snowcrab::snowcrab_carstm( DS="parameters", assessment.years=1999:2018 )

  # misc run params adjustments here:
  p$inla_num.threads = 2
  p$inla_blas.num.threads = 2
  p$modeldata = 'snowcrab_carstm( p=p, DS="carstm_inputs" )'




# -------------------------------------------------
# Model -- nonspatial, nontemporal -- glm
  p$carstm_model_label = "factorial_gaussian_biomass_glm"
  p$variabletomodel = "totwgt"
  p$selection = list(type = "biomass") # d efault is "number"  ... specify to override
  p$carstm_modelengine = "glm"
  p$carstm_modelcall = paste( '
    glm( formula = ', p$variabletomodel, '  ~ AUID:year -1,
         family = gaussian(link="identity"),
         data = M[ which(M$tag=="observations"), ],
    ) ' )
  M = snowcrab_carstm( p=p, DS="carstm_inputs", redo=FALSE )  # will redo if not found
  M$year = as.factor(M$year)
  M[,p$variabletomodel] = M[,p$variabletomodel] / M$data_offset  # cannot do offsets in gaussian linear model
  res = carstm_model( p=p, M=M )




# -------------------------------------------------
# Model -- nonspatial, nontemporal -- glm -- poisson
  p$carstm_model_label = "factorial_poisson_glm"
  p$variabletomodel = "totno"
  p$carstm_modelengine = "glm"
  p$carstm_modelcall = paste( '
    glm( formula = ', p$variabletomodel, '  ~ AUID:year -1 + offset( log(data_offset)),
         family = poisson(link="log"),
         data = M[ which(M$tag=="observations"), ],
    ) ' )
  M = snowcrab_carstm( p=p, DS="carstm_inputs", redo=FALSE )  # will redo if not found
  M$year = as.factor(M$year)
  res = carstm_model( p=p, M=M )




# -------------------------------------------------
  ## not working .. looks like negative values are predicted ...  these need to be filtered out

  p$carstm_model_label = "factorial_gaussian"
  p$variabletomodel = "totwgt"
  p$selection = list(type = "biomass") # d efault is "number"  ... specify to override
  p$carstm_modelengine = "inla"
  p$carstm_modelcall = paste(
    'inla( formula =', p$variabletomodel, ' ~ -1 + year_factor:auid ,
      family = "gaussian",
      data= M,
      control.compute=list(dic=TRUE, config=TRUE),
      control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
      control.predictor=list(compute=FALSE, link=1 ),
      # control.fixed = list(prec.intercept = 0.1),
      control.inla = list(  tolerance=1e-12), # increase in case values are too close to zero
    verbose=TRUE
    )'
  )
  M = snowcrab_carstm( p=p, DS="carstm_inputs", redo=FALSE )  # will redo if not found

  # remove unsampled locations in factorial methods to get sensible stats
  M$uid = paste( M$AUID, M$year, sep="." )
  withdata = unique( M$uid[which(M$tag== "observations")])
  preds = which(M$tag== "predictions")
  todrop = which( ! M$uid[ preds] %in% withdata )
  M = M[ - preds[todrop]]

  M$year_factor = factor(M$year)
  M$auid = factor( M$auid)
  M[,p$variabletomodel] = M[,p$variabletomodel] / M$data_offset  # cannot do offsets in gaussian linear model

  res = carstm_model( p=p, M=M )




# -------------------------------------------------
  p$carstm_model_label = "factorial_simple"
  p$variabletomodel = "totno"
  p$carstm_modelengine = "inla"
  p$carstm_modelcall = paste(
    'inla( formula =', p$variabletomodel,
    ' ~ -1
       + year_factor:auid
       + offset( log(data_offset)) ,
      family = "poisson",
      data= M,
      control.compute=list(dic=TRUE, config=TRUE),
      control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
      control.predictor=list(compute=FALSE, link=1 ),
      control.fixed = list(prec.intercept = 1),
      # control.fixed = H$fixed,  # priors for fixed effects, generic is ok
      control.inla = list(cmin = 0, h=1e-2, tolerance=1e-4), # increase in case values are too close to zero
      verbose=TRUE
    )'
  )
  M = snowcrab_carstm( p=p, DS="carstm_inputs", redo=FALSE )  # will redo if not found

    # remove unsampled locations in factorial methods to get sensible stats
  M$uid = paste( M$AUID, M$year, sep="." )
  withdata = unique( M$uid[which(M$tag== "observations")])
  preds = which(M$tag== "predictions")
  todrop = which( ! M$uid[ preds] %in% withdata )
  M = M[ - preds[todrop]]

  M$year_factor = factor(M$year)
  M$auid = factor( M$auid)
  M[,p$variabletomodel] = trunc( M[,p$variabletomodel] )  # poisson wants integers

  res = carstm_model( p=p, M=M )



# -------------------------------------------------
  p$carstm_model_label = "factorial_full"
  p$variabletomodel = "totno"
  p$carstm_modelengine = "inla"
  p$carstm_modelcall = paste(
    'inla( formula =', p$variabletomodel,
    ' ~ -1
       + year_factor
       + auid
       + year_factor:auid
       + offset( log(data_offset)) ,
      family = "poisson",
      data= M,
      control.compute=list(dic=TRUE, config=TRUE),
      control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
      control.predictor=list(compute=FALSE, link=1 ),
      control.fixed = list(prec.intercept = 1),
      # control.fixed = H$fixed,  # priors for fixed effects, generic is ok
      control.inla = list(cmin = 0, h=1e-2, tolerance=1e-4), # increase in case values are too close to zero
      verbose=TRUE
    )'
  )
  M = snowcrab_carstm( p=p, DS="carstm_inputs", redo=FALSE )  # will redo if not found

  # remove unsampled locations in factorial methods to get sensible stats
  M$uid = paste( M$AUID, M$year, sep="." )
  withdata = unique( M$uid[which(M$tag== "observations")])
  preds = which(M$tag== "predictions")
  todrop = which( ! M$uid[ preds] %in% withdata )
  M = M[ - preds[todrop]]

  M$year_factor = factor(M$year)
  M$auid = factor( M$auid)
  M[,p$variabletomodel] = trunc( M[,p$variabletomodel] )  # poisson wants integers


  res = carstm_model( p=p, M=M )


# -------------------------------------------------
  p$carstm_model_label = "factorial_full_covariates"
  p$variabletomodel = "totno"
  p$carstm_modelengine = "inla"
  p$carstm_modelcall = paste(
    'inla( formula =', p$variabletomodel,
    ' ~ -1
       + year_factor
       + auid
       + year_factor:auid
        + f( dyri, model="ar1", hyper=H$ar1 )
        + f( inla.group( z, method="quantile", n=13 ), model="rw2", scale.model=TRUE, hyper=H$rw2)
        + f( inla.group( substrate.grainsize, method="quantile", n=13 ), model="rw2", scale.model=TRUE, hyper=H$rw2)
        + f( inla.group( t, method="quantile", n=13 ), model="rw2", scale.model=TRUE, hyper=H$rw2)
        + f( inla.group( pca1, method="quantile", n=13 ), model="rw2", scale.model=TRUE, hyper=H$rw2)
        + f( inla.group( pca2, method="quantile", n=13 ), model="rw2", scale.model=TRUE, hyper=H$rw2)
       + offset( log(data_offset)) ,
      family = "poisson",
      data= M,
      control.compute=list(dic=TRUE, config=TRUE),
      control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
      control.predictor=list(compute=FALSE, link=1 ),
      control.fixed = list(prec.intercept = 1),
      # control.fixed = H$fixed,  # priors for fixed effects, generic is ok
      control.inla = list(cmin = 0, h=1e-2, tolerance=1e-4), # increase in case values are too close to zero
      verbose=TRUE
    )'
  )
  M = snowcrab_carstm( p=p, DS="carstm_inputs", redo=FALSE )  # will redo if not found
  # remove unsampled locations in factorial methods to get sensible stats
  M$uid = paste( M$AUID, M$year, sep="." )
  withdata = unique( M$uid[which(M$tag== "observations")])
  preds = which(M$tag== "predictions")
  todrop = which( ! M$uid[ preds] %in% withdata )
  M = M[ - preds[todrop]]

  M$year_factor = factor(M$year)
  M$auid = factor( M$auid)
  M[,p$variabletomodel] = trunc( M[,p$variabletomodel] )  # poisson wants integers


  res = carstm_model( p=p, M=M )



# -------------------------------------------------

  p$carstm_model_label = "mixed_effects_simple"
  p$variabletomodel = "totno"
  p$carstm_modelengine = "inla"
  p$carstm_modelcall = paste(
    'inla( formula =', p$variabletomodel,
      ' ~ 1
        + offset( log(data_offset))
        + year_factor
        + f( auid, model="bym2", graph=sppoly@nb, scale.model=TRUE, constr=TRUE, hyper=H$bym2),
        family = "poisson",
        data= M,
        control.compute=list(dic=TRUE, config=TRUE),
        control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
        control.predictor=list(compute=FALSE, link=1 ),
        # control.fixed = list(prec.intercept = 0.1),
        control.fixed = H$fixed,  # priors for fixed effects, generic is ok
        control.inla = list(cmin = 0, h=1e-6, tolerance=1e-12), # increase in case values are too close to zero
        # control.inla=list( strategy="laplace", cutoff=1e-6, correct=TRUE, correct.verbose=FALSE ),
        # control.inla = list( h=1e-6, tolerance=1e-12), # increase in case values are too close to zero
        # control.inla = list(h=1e-3, tolerance=1e-9, cmin=0), # restart=3), # restart a few times in case posteriors are poorly defined
        # control.mode = list( restart=TRUE, result=RES ), # restart from previous estimates
      verbose=TRUE
    )'
  )
  M = snowcrab_carstm( p=p, DS="carstm_inputs", redo=FALSE )  # will redo if not found
  M$year_factor = factor(M$year)
  M[,p$variabletomodel] = trunc( M[,p$variabletomodel] )  # poisson wants integers
  res = carstm_model( p=p, M=M )



# -------------------------------------------------
  p$carstm_model_label = "mixed_effects_static"
  p$variabletomodel = "totno"
  p$carstm_modelengine = "inla"
  p$carstm_modelcall = paste(
    'inla( formula =', p$variabletomodel,
      ' ~ 1
        + offset( log(data_offset))
        + year_factor
        + f( inla.group( z, method="quantile", n=13 ), model="rw2", scale.model=TRUE, hyper=H$rw2)
        + f( inla.group( substrate.grainsize, method="quantile", n=13 ), model="rw2", scale.model=TRUE, hyper=H$rw2)
        + f( auid, model="bym2", graph=sppoly@nb, scale.model=TRUE, constr=TRUE, hyper=H$bym2),
        family = "poisson",
        data= M,
        control.compute=list(dic=TRUE, config=TRUE),
        control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
        control.predictor=list(compute=FALSE, link=1 ),
        # control.fixed = list(prec.intercept = 0.1),
        control.fixed = H$fixed,  # priors for fixed effects, generic is ok
        control.inla = list(cmin = 0, h=1e-6, tolerance=1e-12), # increase in case values are too close to zero
        # control.inla=list( strategy="laplace", cutoff=1e-6, correct=TRUE, correct.verbose=FALSE ),
        # control.inla = list( h=1e-6, tolerance=1e-12), # increase in case values are too close to zero
        # control.inla = list(h=1e-3, tolerance=1e-9, cmin=0), # restart=3), # restart a few times in case posteriors are poorly defined
        # control.mode = list( restart=TRUE, result=RES ), # restart from previous estimates
      verbose=TRUE
    )'
  )
  M = snowcrab_carstm( p=p, DS="carstm_inputs", redo=FALSE )  # will redo if not found
  M$year_factor = factor(M$year)
  # M$auid = factor( M$auid)
  M[,p$variabletomodel] = trunc( M[,p$variabletomodel] )  # poisson wants integers
  res = carstm_model( p=p, M=M )



# -------------------------------------------------
  p$carstm_model_label = "mixed_effects_dynamic"
  p$variabletomodel = "totno"
  p$carstm_modelengine = "inla"
  p$carstm_modelcall = paste(
    'inla( formula =', p$variabletomodel,
      ' ~ 1
        + offset( log(data_offset))
        + year_factor
        + f( inla.group( z, method="quantile", n=13 ), model="rw2", scale.model=TRUE, hyper=H$rw2)
        + f( inla.group( substrate.grainsize, method="quantile", n=13 ), model="rw2", scale.model=TRUE, hyper=H$rw2)
        + f( inla.group( t, method="quantile", n=13 ), model="rw2", scale.model=TRUE, hyper=H$rw2)
        + f( inla.group( pca1, method="quantile", n=13 ), model="rw2", scale.model=TRUE, hyper=H$rw2)
        + f( inla.group( pca2, method="quantile", n=13 ), model="rw2", scale.model=TRUE, hyper=H$rw2)
        + f( auid, model="bym2", graph=sppoly@nb, scale.model=TRUE, constr=TRUE, hyper=H$bym2),
        family = "poisson",
        data= M,
        control.compute=list(dic=TRUE, config=TRUE),
        control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
        control.predictor=list(compute=FALSE, link=1 ),
        # control.fixed = list(prec.intercept = 0.1),
        control.fixed = H$fixed,  # priors for fixed effects, generic is ok
        control.inla = list(cmin = 0, h=1e-6, tolerance=1e-12), # increase in case values are too close to zero
        # control.inla=list( strategy="laplace", cutoff=1e-6, correct=TRUE, correct.verbose=FALSE ),
        # control.inla = list( h=1e-6, tolerance=1e-12), # increase in case values are too close to zero
        # control.inla = list(h=1e-3, tolerance=1e-9, cmin=0), # restart=3), # restart a few times in case posteriors are poorly defined
        # control.mode = list( restart=TRUE, result=RES ), # restart from previous estimates
      verbose=TRUE
    )'
  )
  M = snowcrab_carstm( p=p, DS="carstm_inputs", redo=FALSE )  # will redo if not found
  M$year_factor = factor(M$year)
  M[,p$variabletomodel] = trunc( M[,p$variabletomodel] )  # poisson wants integers
  res = carstm_model( p=p, M=M )


# -------------------------------------------------
  p$carstm_model_label = "separable"
  p$variabletomodel = "totno"
  p$carstm_modelengine = "inla"
  p$carstm_modelcall = paste(
    'inla( formula =', p$variabletomodel,
      ' ~ 1
        + offset( log(data_offset))
        + f( year_factor, model="ar1", hyper=H$ar1 )
        + f( inla.group( z, method="quantile", n=13 ), model="rw2", scale.model=TRUE, hyper=H$rw2)
        + f( inla.group( substrate.grainsize, method="quantile", n=13 ), model="rw2", scale.model=TRUE, hyper=H$rw2)
        + f( inla.group( t, method="quantile", n=13 ), model="rw2", scale.model=TRUE, hyper=H$rw2)
        + f( inla.group( pca1, method="quantile", n=13 ), model="rw2", scale.model=TRUE, hyper=H$rw2)
        + f( inla.group( pca2, method="quantile", n=13 ), model="rw2", scale.model=TRUE, hyper=H$rw2)
        + f( auid, model="bym2", graph=sppoly@nb, scale.model=TRUE, constr=TRUE, hyper=H$bym2),
        family = "poisson",
        data= M,
        control.compute=list(dic=TRUE, config=TRUE),
        control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
        control.predictor=list(compute=FALSE, link=1 ),
        # control.fixed = list(prec.intercept = 0.1),
        control.fixed = H$fixed,  # priors for fixed effects, generic is ok
        control.inla = list(cmin = 0, h=1e-6, tolerance=1e-12), # increase in case values are too close to zero
        # control.inla=list( strategy="laplace", cutoff=1e-6, correct=TRUE, correct.verbose=FALSE ),
        # control.inla = list( h=1e-6, tolerance=1e-12), # increase in case values are too close to zero
        # control.inla = list(h=1e-3, tolerance=1e-9, cmin=0), # restart=3), # restart a few times in case posteriors are poorly defined
        # control.mode = list( restart=TRUE, result=RES ), # restart from previous estimates
      verbose=TRUE
    )'
  )
  res = carstm_model( p=p, M=p$modeldata )



# -------------------------------------------------

  p$carstm_model_label = "separable_simple"
  p$variabletomodel = "totno"
  p$carstm_modelengine = "inla"
  p$carstm_modelcall = paste(
    'inla( formula =', p$variabletomodel,
      ' ~ 1
        + offset( log(data_offset))
        + f( year_factor, model="ar1", hyper=H$ar1 )
        + f( auid, model="bym2", graph=sppoly@nb, scale.model=TRUE, constr=TRUE, hyper=H$bym2),
        family = "poisson",
        data= M,
        control.compute=list(dic=TRUE, config=TRUE),
        control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
        control.predictor=list(compute=FALSE, link=1 ),
        # control.fixed = list(prec.intercept = 0.1),
        control.fixed = H$fixed,  # priors for fixed effects, generic is ok
        control.inla = list(cmin = 0, h=1e-6, tolerance=1e-12), # increase in case values are too close to zero
        # control.inla = list(cmin = 0 ),
        # control.inla = list( h=1e-6, tolerance=1e-12), # increase in case values are too close to zero
        # control.inla=list( strategy="laplace", cutoff=1e-6, correct=TRUE, correct.verbose=FALSE ),
        # control.inla = list( h=1e-6, tolerance=1e-12), # increase in case values are too close to zero
        # control.inla = list(h=1e-3, tolerance=1e-9, cmin=0), # restart=3), # restart a few times in case posteriors are poorly defined
        # control.mode = list( restart=TRUE, result=RES ), # restart from previous estimates
      verbose=TRUE
    )'
  )
  res = carstm_model( p=p, M=p$modeldata )





# -------------------------------------------------
  p$carstm_model_label = "nonseparable_simple" # nonseparable, basic
  p$variabletomodel = "totno"
  p$carstm_modelengine = "inla"
  p$carstm_modelcall = paste(
    'inla( formula =', p$variabletomodel,
      ' ~ 1
        + offset( log(data_offset))
        + f( auid, model="bym2", graph=sppoly@nb, group=year_factor, scale.model=TRUE, constr=TRUE, hyper=H$bym2, control.group=list(model="ar1", hyper=H$ar1_group)),
        family = "poisson",
        data= M,
        control.compute=list(dic=TRUE, config=TRUE),
        control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
        control.predictor=list(compute=FALSE, link=1 ),
        # control.fixed = list(prec.intercept = 0.1),
        control.fixed = H$fixed,  # priors for fixed effects, generic is ok
        control.inla = list(cmin = 0 ),
        # control.inla = list( h=1e-6, tolerance=1e-12), # increase in case values are too close to zero
        # control.inla=list( strategy="laplace", cutoff=1e-6, correct=TRUE, correct.verbose=FALSE ),
        # control.inla = list( h=1e-6, tolerance=1e-12), # increase in case values are too close to zero
        # control.inla = list(h=1e-3, tolerance=1e-9, cmin=0), # restart=3), # restart a few times in case posteriors are poorly defined
        # control.mode = list( restart=TRUE, result=RES ), # restart from previous estimates
      verbose=TRUE
    )'
  )
  res = carstm_model( p=p, M=p$modeldata )




# -------------------------------------------------
# Model -- nonspatial, nontemporal -- inla -- poisson
  p$carstm_model_label = "nonseparable_space-time"
  p$variabletomodel = "totno"
  p$carstm_modelengine = "inla"
  p$carstm_modelcall = paste(
    'inla( formula =', p$variabletomodel,
      ' ~ 1
        + offset( log(data_offset))
        + f( dyri, model="ar1", hyper=H$ar1 )
        + f( inla.group( t, method="quantile", n=13 ), model="rw2", scale.model=TRUE, hyper=H$rw2)
        + f( inla.group( z, method="quantile", n=13 ), model="rw2", scale.model=TRUE, hyper=H$rw2)
        + f( inla.group( substrate.grainsize, method="quantile", n=13 ), model="rw2", scale.model=TRUE, hyper=H$rw2)
        + f( inla.group( pca1, method="quantile", n=13 ), model="rw2", scale.model=TRUE, hyper=H$rw2)
        + f( inla.group( pca2, method="quantile", n=13 ), model="rw2", scale.model=TRUE, hyper=H$rw2)
        + f( auid, model="bym2", graph=sppoly@nb, group=year_factor, scale.model=TRUE, constr=TRUE, hyper=H$bym2, control.group=list(model="ar1", hyper=H$ar1_group)),
        family = "poisson",
        data= M,
        control.compute=list(dic=TRUE, config=TRUE),
        control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
        control.predictor=list(compute=FALSE, link=1 ),
        # control.fixed = list(prec.intercept = 0.1),
        control.fixed = H$fixed,  # priors for fixed effects, generic is ok
        control.inla = list(cmin = 0 ),
        # control.inla = list( h=1e-6, tolerance=1e-12), # increase in case values are too close to zero
        # control.inla=list( strategy="laplace", cutoff=1e-6, correct=TRUE, correct.verbose=FALSE ),
        # control.inla = list( h=1e-6, tolerance=1e-12), # increase in case values are too close to zero
        # control.inla = list(h=1e-3, tolerance=1e-9, cmin=0), # restart=3), # restart a few times in case posteriors are poorly defined
        # control.mode = list( restart=TRUE, result=RES ), # restart from previous estimates
      verbose=TRUE
    )'
  )
  res = carstm_model( p=p, M=p$modeldata )



# -------------------------------------------------
  p$carstm_model_label = "nonseparable_time-space"
  p$variabletomodel = "totno"
  p$carstm_modelengine = "inla"
  p$carstm_modelcall = paste(
    'inla( formula =', p$variabletomodel,
      ' ~ 1
        + offset( log(data_offset))
        + f( year_factor, model="ar1", hyper=H$ar1, group=auid, control.group=list(model="besag", graph=sppoly@nb, scale.model=TRUE ) )
        + f( dyri, model="ar1", hyper=H$ar1 )
        + f( inla.group( t, method="quantile", n=13 ), model="rw2", scale.model=TRUE, hyper=H$rw2)
        + f( inla.group( z, method="quantile", n=13 ), model="rw2", scale.model=TRUE, hyper=H$rw2)
        + f( inla.group( substrate.grainsize, method="quantile", n=13 ), model="rw2", scale.model=TRUE, hyper=H$rw2)
        + f( inla.group( pca1, method="quantile", n=13 ), model="rw2", scale.model=TRUE, hyper=H$rw2)
        + f( inla.group( pca2, method="quantile", n=13 ), model="rw2", scale.model=TRUE, hyper=H$rw2),
        family = "poisson",
        data= M,
        control.compute=list(dic=TRUE, config=TRUE),
        control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
        control.predictor=list(compute=FALSE, link=1 ),
        # control.fixed = list(prec.intercept = 0.1),
        control.fixed = H$fixed,  # priors for fixed effects, generic is ok
        control.inla = list(cmin = 0 ),
        # control.inla = list( h=1e-6, tolerance=1e-12), # increase in case values are too close to zero
        # control.inla=list( strategy="laplace", cutoff=1e-6, correct=TRUE, correct.verbose=FALSE ),
        # control.inla = list( h=1e-6, tolerance=1e-12), # increase in case values are too close to zero
        # control.inla = list(h=1e-3, tolerance=1e-9, cmin=0), # restart=3), # restart a few times in case posteriors are poorly defined
        # control.mode = list( restart=TRUE, result=RES ), # restart from previous estimates
      verbose=TRUE
    )'
  )
  res = carstm_model( p=p, M=p$modeldata )


# -------------------------------------------------
# Model -- nonspatial, nontemporal -- inla -- poisson 4.5 hrs
  p$carstm_model_label = "inla_zeroinflatedpoisson0_full"
  p$variabletomodel = "totno"
  p$carstm_modelengine = "inla"
  p$carstm_modelcall = paste(
    'inla( formula =', p$variabletomodel,
      ' ~ 1
        + offset( log(data_offset))
        + f( dyri, model="ar1", hyper=H$ar1 )
        + f( inla.group( t, method="quantile", n=13 ), model="rw2", scale.model=TRUE, hyper=H$rw2)
        + f( inla.group( z, method="quantile", n=13 ), model="rw2", scale.model=TRUE, hyper=H$rw2)
        + f( inla.group( substrate.grainsize, method="quantile", n=13 ), model="rw2", scale.model=TRUE, hyper=H$rw2)
        + f( inla.group( pca1, method="quantile", n=13 ), model="rw2", scale.model=TRUE, hyper=H$rw2)
        + f( inla.group( pca2, method="quantile", n=13 ), model="rw2", scale.model=TRUE, hyper=H$rw2)
        + f( auid, model="bym2", graph=sppoly@nb, group=year_factor, scale.model=TRUE, constr=TRUE, hyper=H$bym2, control.group=list(model="ar1", hyper=H$ar1_group)),
        family = "zeroinflatedpoisson0",
        data= M,
        control.compute=list(dic=TRUE, config=TRUE),
        control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
        control.predictor=list(compute=FALSE, link=1 ),
        # control.fixed = list(prec.intercept = 0.1),
        control.fixed = H$fixed,  # priors for fixed effects, generic is ok
        control.inla = list(cmin = 0 ),
        # control.inla = list( h=1e-6, tolerance=1e-12), # increase in case values are too close to zero
        # control.inla=list( strategy="laplace", cutoff=1e-6, correct=TRUE, correct.verbose=FALSE ),
        # control.inla = list( h=1e-6, tolerance=1e-12), # increase in case values are too close to zero
        # control.inla = list(h=1e-3, tolerance=1e-9, cmin=0), # restart=3), # restart a few times in case posteriors are poorly defined
        # control.mode = list( restart=TRUE, result=RES ), # restart from previous estimates
      verbose=TRUE
    )'
  )
  res = carstm_model( p=p, M=p$modeldata )






# -------------------------------------------------
  p$carstm_model_label = "inla_zeroinflatedpoisson1_full"  # strange
  p$variabletomodel = "totno"
  p$carstm_modelengine = "inla"
  p$carstm_modelcall = paste(
    'inla( formula =', p$variabletomodel,
      ' ~ 1
        + offset( log(data_offset))
        + f( dyri, model="ar1", hyper=H$ar1 )
        + f( inla.group( t, method="quantile", n=13 ), model="rw2", scale.model=TRUE, hyper=H$rw2)
        + f( inla.group( z, method="quantile", n=13 ), model="rw2", scale.model=TRUE, hyper=H$rw2)
        + f( inla.group( substrate.grainsize, method="quantile", n=13 ), model="rw2", scale.model=TRUE, hyper=H$rw2)
        + f( inla.group( pca1, method="quantile", n=13 ), model="rw2", scale.model=TRUE, hyper=H$rw2)
        + f( inla.group( pca2, method="quantile", n=13 ), model="rw2", scale.model=TRUE, hyper=H$rw2)
        + f( auid, model="bym2", graph=sppoly@nb, group=year_factor, scale.model=TRUE, constr=TRUE, hyper=H$bym2, control.group=list(model="ar1", hyper=H$ar1_group)),
        family = "zeroinflatedpoisson1",
        data= M,
        control.compute=list(dic=TRUE, config=TRUE),
        control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
        control.predictor=list(compute=FALSE, link=1 ),
        # control.fixed = list(prec.intercept = 0.1),
        control.fixed = H$fixed,  # priors for fixed effects, generic is ok
        control.inla = list(cmin = 0 ),
        # control.inla = list( h=1e-6, tolerance=1e-12), # increase in case values are too close to zero
        # control.inla=list( strategy="laplace", cutoff=1e-6, correct=TRUE, correct.verbose=FALSE ),
        # control.inla = list( h=1e-6, tolerance=1e-12), # increase in case values are too close to zero
        # control.inla = list(h=1e-3, tolerance=1e-9, cmin=0), # restart=3), # restart a few times in case posteriors are poorly defined
        # control.mode = list( restart=TRUE, result=RES ), # restart from previous estimates
      verbose=TRUE
    )'
  )
  res = carstm_model( p=p, M=p$modeldata )






# -------------------------------------------------
# generic calls

  res = carstm_model( p=p, DS="carstm_modelled" ) # to load currently saved res
  fit = carstm_model( p=p, DS="carstm_modelled_fit" )  # extract currently saved model fit

  if (0) {
    plot( fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE )
    plot( fit$marginals.hyperpar$"Phi for auid", type="l")  # posterior distribution of phi nonspatial dominates
    plot( fit$marginals.hyperpar$"Precision for auid", type="l")
    plot( fit$marginals.hyperpar$"Precision for setno", type="l")
  }

  # AIC(fit)
  # DIC(fit)
  s = summary(fit)
  s
  # s$dic$dic
  # s$dic$p.eff


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


  vn = paste(p$variabletomodel, "random_sample_iid", sep=".")
  if (exists(vn, res)) carstm_plot( p=p, res=res, vn=vn, time_match=list(year="1950", dyear="0.05") )

  vn = paste(p$variabletomodel, "random_auid_nonspatial", sep=".")
  if (exists(vn, res)) {
    res_dim = length( dim( res[[vn]] ) )
    if (res_dim == 1 ) time_match = NULL
    if (res_dim == 2 ) time_match = list(year="2000")
    if (res_dim == 3 ) time_match = list(year="2000", dyear="0.85" )
    carstm_plot( p=p, res=res, vn=vn, time_match=time_match )
  }

  vn = paste(p$variabletomodel, "random_auid_spatial", sep=".")
  if (exists(vn, res)) {
    res_dim = length( dim( res[[vn]] ) )
    if (res_dim == 1 ) time_match = NULL
    if (res_dim == 2 ) time_match = list(year="2000")
    if (res_dim == 3 ) time_match = list(year="2000", dyear="0.85" )
    carstm_plot( p=p, res=res, vn=vn, time_match=time_match )
  }




  # -------------------------------
  ## --- comparisons

  cc =  as.data.frame( rbind(
    c("factorial_gaussian_biomass_glm", "Factorial crossed (Gaussian) biomass                               ", "darkgreen", 26, 1, 8, "l"  ),
    c("factorial_simple",        "Factorial crossed",        "slategray",  20, 1, 4, "l"),
    c("factorial_full_covariates", "Factorial covariates",   "gray",       20, 7, 4, "l"),
    c("mixed_effects_simple",    "Mixed effects simple",     "turquoise",  21, 2, 4, "l"),
    c("mixed_effects_static",    "Mixed effects static",     "darkorange", 22, 4, 4, "l"),
    c("mixed_effects_dynamic",   "Mixed effects dynamic",    "blue",       23, 5, 4, "l"),
    c("separable",               "Separable",                "green",      24, 6, 4, "l"),
    c("separable_simple",        "Separable simple",         "darkred",    25, 7, 4, "l"),
    c("nonseparable_simple",     "Nonseparable simple",      "cyan",       26, 1, 4, "l"),
    c("nonseparable_time-space", "Nonseparable time|space",  "darkgreen",  27, 3, 4, "l"),
    c("nonseparable_space-time", "Nonseparable space|time",  "slateblue",  20, 1, 8, "l")
  ), stringsAsFactors=FALSE)

  #   c( "inla_zeroinflatedpoisson0_full", "Nonseparable overdispersed space|time")

  colnames(cc) = c("tocompare", "legend", "col", "pch", "lty", "lwd", "type" )
  cc$pch = as.numeric(cc$pch)
  cc$lty = as.numeric(cc$lty)
  cc$lwd = as.numeric(cc$lwd)

  res_ts = list()
  for ( lab in cc$tocompare ) {
    p$carstm_modelengine = "inla"
    p$variabletomodel = "totno"
    if (lab == "factorial_gaussian_biomass_glm")  {
      p$carstm_modelengine = "glm"
      p$variabletomodel = "totwgt"
    }
    res_ts[[lab]] = carstm_summary(p=p, operation="load_timeseries", carstm_model_label=lab  )
  }

  dev.new(width=11, height=7)
  i=1 ; plot( cfaall  ~ yrs, data=res_ts[[i]], lty=cc$lty[i], lwd=cc$lwd[i], col=cc$col[i], pch=cc$pch[i], type=cc$type[i], ylim=c(0,185), xlab="Year", ylab="kt")
  for (i in 2:length(res_ts)) {
    lines( cfaall ~ yrs, data=res_ts[[i]], lty=cc$lty[i], lwd=cc$lwd[i], col=cc$col[i], pch=cc$pch[i], type=cc$type[i])
  }
  legend("top", legend=cc$legend, lty=cc$lty, col=cc$col, lwd=cc$lwd, bty="n"  )










# -------------------------------------------------
# end
# -------------------------------------------------
