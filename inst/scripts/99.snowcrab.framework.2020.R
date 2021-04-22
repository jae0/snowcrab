

# Snow crab --- Areal unit modelling of habitat
# -- only AU's of interest ... not connected to a global analysis
# -- losing information outside of study area but faster
# -- results cannot be reused outside of workflow

# -------------------------------------------------
# Initiate -- construct basic parameter list defining the main characteristics of the study
# expects bio.snowcrab/inst/scripts/03.abundance_estimation_carstm.r to have been run ..
# we are running simpler variations of that model here


  # adjust based upon RAM requirements and ncores
  inla.setOption(num.threads= floor( parallel::detectCores() / 2) )
  inla.setOption(blas.num.threads= 2 )


  p = bio.snowcrab::snowcrab_parameters( project_class="carstm", assessment.years=1999:2018 )

  sppoly = areal_units( p=p, redo=TRUE )  # to reload

  # misc run params adjustments here:
  p$modeldata = 'snowcrab.db( p=p, DS="carstm_inputs" )'


# -------------------------------------------------
# Model -- nonspatial, nontemporal -- glm
  p$carstm_model_label = "factorial_gaussian_biomass_glm"
  p$variabletomodel = "totwgt"
  p$selection$type = "biomass" # d efault is "number"  ... specify to override
  p$carstm_modelengine = "glm"

  p$carstm_model_formula = as.formula( paste(
    p$variabletomodel, '  ~ AUID:year -1 ' ))

  p$carstm_model_family = gaussian(link="identity")

  M = snowcrab.db( p=p, DS="carstm_inputs", redo=FALSE )  # will redo if not found
  M$year = as.factor(M$year)
  M[,p$variabletomodel] = M[,p$variabletomodel] / M$data_offset  # cannot do offsets in gaussian linear model
  M = M[ which(M$tag=="observations"), ]

  fit = carstm_model( p=p, M=M )




# -------------------------------------------------
# Model -- nonspatial, nontemporal -- glm -- poisson
  p$carstm_model_label = "factorial_poisson_glm"
  p$variabletomodel = "totno"
  p$carstm_modelengine = "glm"

  p$carstm_model_formula = as.formula( paste(
    p$variabletomodel, ' ~ AUID:year -1 + offset( log(data_offset)) '))

  p$carstm_model_family = poisson(link="log")
  M = snowcrab.db( p=p, DS="carstm_inputs", redo=FALSE )  # will redo if not found
  M$year = as.factor(M$year)
  M = M[ which(M$tag=="observations"), ]

  fit = carstm_model( p=p, M=M )




# -------------------------------------------------
  ## not working .. looks like negative values are predicted ...

  p$carstm_model_label = "factorial_gaussian"
  p$carstm_predict_force_range = TRUE  # for factorial models this is necessary to prevent meainingless predictions
  p$variabletomodel = "totwgt"
  p$selection$type = "biomass" # d efault is "number"  ... specify to override
  p$carstm_modelengine = "inla"

  p$carstm_model_formula = as.formula( paste(
    p$variabletomodel, ' ~ year_factor:auid '))

  p$carstm_model_family = "gaussian"

  M = snowcrab.db( p=p, DS="carstm_inputs", redo=FALSE )  # will redo if not found

  # remove unsampled locations in factorial methods to get sensible stats
  M$uid = paste( M$AUID, M$year, sep="." )
  withdata = unique( M$uid[which(M$tag== "observations")])
  preds = which(M$tag== "predictions")
  todrop = which( ! M$uid[ preds] %in% withdata )
  M = M[ - preds[todrop]]

  M$year_factor = factor(M$year)
  M$auid = factor( M$auid)
  M[,p$variabletomodel] = M[,p$variabletomodel] / M$data_offset  # cannot do offsets in gaussian linear model

  fit = carstm_model( p=p, M=M )




# -------------------------------------------------
  p$carstm_model_label = "factorial_simple"
  p$carstm_predict_force_range = TRUE  # for factorial models this is necessary to prevent meainingless predictions
  p$variabletomodel = "totno"
  p$carstm_modelengine = "inla"

  p$carstm_model_formula = as.formula( paste(
    p$variabletomodel, ' ~ -1',
       ' + year_factor:auid ',
       ' + offset( log(data_offset)) ' ))

  p$carstm_model_family = "poisson"

  M = snowcrab.db( p=p, DS="carstm_inputs", redo=FALSE )  # will redo if not found

    # remove unsampled locations in factorial methods to get sensible stats
  M$uid = paste( M$AUID, M$year, sep="." )
  withdata = unique( M$uid[which(M$tag== "observations")])
  preds = which(M$tag== "predictions")
  todrop = which( ! M$uid[ preds] %in% withdata )
  M = M[ - preds[todrop]]

  M$year_factor = factor(M$year)
  M$auid = factor( M$auid)

  region.id = slot( slot(sppoly, "nb"), "region.id" )
  M$auid = match( M$AUID, region.id )

  M[,p$variabletomodel] = floor( M[,p$variabletomodel] )  # poisson wants integers

  fit = carstm_model( p=p, M=M )



# -------------------------------------------------
  p$carstm_model_label = "factorial_full"
  p$carstm_predict_force_range = TRUE  # for factorial models this is necessary to prevent meainingless predictions
  p$variabletomodel = "totno"
  p$carstm_modelengine = "inla"


  p$carstm_model_formula = as.formula( paste(
    p$variabletomodel, ' ~ -1',
       ' + year_factor',
       ' + auid',
       ' + year_factor:auid ',
       ' + offset( log(data_offset)) '
  ))

  p$carstm_model_family = "poisson"

  M = snowcrab.db( p=p, DS="carstm_inputs", redo=FALSE )  # will redo if not found

  # remove unsampled locations in factorial methods to get sensible stats
  M$uid = paste( M$AUID, M$year, sep="." )
  withdata = unique( M$uid[which(M$tag== "observations")])
  preds = which(M$tag== "predictions")
  todrop = which( ! M$uid[ preds] %in% withdata )
  M = M[ - preds[todrop]]

  M$year_factor = factor(M$year)
  M$auid = factor( M$auid)
  M[,p$variabletomodel] = floor( M[,p$variabletomodel] )  # poisson wants integers


  fit = carstm_model( p=p, M=M )


# -------------------------------------------------
  p$carstm_model_label = "factorial_full_covariates"
  p$variabletomodel = "totno"
  p$carstm_modelengine = "inla"
  p$carstm_predict_force_range = TRUE  # for factorial models this is necessary to prevent meainingless predictions

  p$carstm_model_formula = as.formula( paste(
    p$variabletomodel, ' ~ -1 ',
       ' + year_factor',
       ' + auid',
       ' + year_factor:auid',
       ' + f( dyri, model="ar1", hyper=H$ar1 )',
       ' + f( inla.group( z, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2)',
       ' + f( inla.group( substrate.grainsize, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2)',
       ' + f( inla.group( t, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2)',
       ' + f( inla.group( pca1, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2)',
       ' + f( inla.group( pca2, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2)',
       ' + offset( log(data_offset)) '
  ))
  p$carstm_model_family = "poisson"

  M = snowcrab.db( p=p, DS="carstm_inputs", redo=FALSE )  # will redo if not found
  # remove unsampled locations in factorial methods to get sensible stats
  M$uid = paste( M$AUID, M$year, sep="." )
  withdata = unique( M$uid[which(M$tag== "observations")])
  preds = which(M$tag== "predictions")
  todrop = which( ! M$uid[ preds] %in% withdata )
  M = M[ - preds[todrop]]

  M$year_factor = factor(M$year)
  M$auid = factor( M$auid)
  M[,p$variabletomodel] = floor( M[,p$variabletomodel] )  # poisson wants integers

  fit = carstm_model( p=p, M=M )



# -------------------------------------------------

# -------------------------------------------------
  p$carstm_model_label = "covariates_only"
  p$variabletomodel = "totno"
  p$carstm_modelengine = "inla"

  p$carstm_model_formula = as.formula( paste(
    p$variabletomodel, ' ~ -1 ',
      ' + f( inla.group( z, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
      ' + f( inla.group( substrate.grainsize, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
      ' + f( inla.group( t, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
      ' + f( inla.group( pca1, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
      ' + f( inla.group( pca2, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
      ' + offset( log(data_offset))  '
  ))

  p$carstm_model_family = "poisson"

  M = snowcrab.db( p=p, DS="carstm_inputs", redo=FALSE )  # will redo if not found
  # remove unsampled locations in factorial methods to get sensible stats
  M$uid = paste( M$AUID, M$year, sep="." )
  withdata = unique( M$uid[which(M$tag== "observations")])
  preds = which(M$tag== "predictions")
  todrop = which( ! M$uid[ preds] %in% withdata )
  M = M[ - preds[todrop]]

  M$year_factor = factor(M$year)
  M$auid = factor( M$auid)
  M[,p$variabletomodel] = floor( M[,p$variabletomodel] )  # poisson wants integers


  fit = carstm_model( p=p, M=M )



# -------------------------------------------------


  p$carstm_model_label = "mixed_effects_simple"
  p$variabletomodel = "totno"
  p$carstm_modelengine = "inla"


  p$carstm_model_formula = as.formula( paste(
    p$variabletomodel, ' ~ 1 ',
      ' + offset( log(data_offset)) ',
      ' + year_factor',
      ' + f( auid, model="bym2", graph=slot(sppoly, "nb"), scale.model=TRUE, constr=TRUE, hyper=H$bym2)'
  ))

  p$carstm_model_family = "poisson"

  M = snowcrab.db( p=p, DS="carstm_inputs", redo=FALSE )  # will redo if not found
  M$year_factor = factor(M$year)
  M[,p$variabletomodel] = floor( M[,p$variabletomodel] )  # poisson wants integers
  fit = carstm_model( p=p, M=M )



# -------------------------------------------------
  p$carstm_model_label = "mixed_effects_static"
  p$variabletomodel = "totno"
  p$carstm_modelengine = "inla"

  p$carstm_model_formula = as.formula( paste(
    p$variabletomodel, ' ~ 1 ',
      ' + offset( log(data_offset)) ',
      ' + year_factor',
      ' + f( inla.group( z, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2)' ,
      ' + f( inla.group( substrate.grainsize, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
      ' + f( auid, model="bym2", graph=slot(sppoly, "nb"), scale.model=TRUE, constr=TRUE, hyper=H$bym2) '
  ))
  p$carstm_model_family = "poisson"

  M = snowcrab.db( p=p, DS="carstm_inputs", redo=FALSE )  # will redo if not found
  M$year_factor = factor(M$year)
  # M$auid = factor( M$auid)
  M[,p$variabletomodel] = floor( M[,p$variabletomodel] )  # poisson wants integers
  fit = carstm_model( p=p, M=M )



# -------------------------------------------------
  p$carstm_model_label = "mixed_effects_dynamic"
  p$variabletomodel = "totno"
  p$carstm_modelengine = "inla"

  p$carstm_model_formula = as.formula( paste(
    p$variabletomodel, ' ~ 1 ',
      ' + offset( log(data_offset)) ',
      ' + year_factor',
      ' + f( inla.group( z, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2)',
      ' + f( inla.group( substrate.grainsize, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2)',
      ' + f( inla.group( t, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2)',
      ' + f( inla.group( pca1, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2)',
      ' + f( inla.group( pca2, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2)',
      ' + f( auid, model="bym2", graph=slot(sppoly, "nb"), scale.model=TRUE, constr=TRUE, hyper=H$bym2)'
  ))

  p$carstm_model_family = "poisson"

  M = snowcrab.db( p=p, DS="carstm_inputs", redo=FALSE )  # will redo if not found
  M$year_factor = factor(M$year)
  M[,p$variabletomodel] = floor( M[,p$variabletomodel] )  # poisson wants integers
  fit = carstm_model( p=p, M=M )


# -------------------------------------------------
  p$carstm_model_label = "separable"
  p$variabletomodel = "totno"
  p$carstm_modelengine = "inla"
  p$carstm_model_formula = as.formula( paste(
    p$variabletomodel, ' ~ 1 ',
      ' + offset( log(data_offset)) ',
      ' + f( year_factor, model="ar1", hyper=H$ar1 )',
      ' + f( inla.group( z, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2)',
      ' + f( inla.group( substrate.grainsize, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2)',
      ' + f( inla.group( t, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2)',
      ' + f( inla.group( pca1, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2)',
      ' + f( inla.group( pca2, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2)',
      ' + f( auid, model="bym2", graph=slot(sppoly, "nb"), scale.model=TRUE, constr=TRUE, hyper=H$bym2)'
  ))

  p$carstm_model_family = "poisson"
  fit = carstm_model( p=p, M=p$modeldata )



# -------------------------------------------------

  p$carstm_model_label = "separable_simple"
  p$variabletomodel = "totno"
  p$carstm_modelengine = "inla"
  p$carstm_model_formula = as.formula( paste(
    p$variabletomodel, ' ~ 1 ',
      ' + offset( log(data_offset)) ',
      ' + f( year_factor, model="ar1", hyper=H$ar1 ) ',
      ' + f( auid, model="bym2", graph=slot(sppoly, "nb"), scale.model=TRUE, constr=TRUE, hyper=H$bym2) '
  ))

  p$carstm_model_family = "poisson"
  fit = carstm_model( p=p, M=p$modeldata )



# -------------------------------------------------
  p$carstm_model_label = "nonseparable_simple" # nonseparable, basic
  p$variabletomodel = "totno"
  p$carstm_modelengine = "inla"
  p$carstm_model_formula = as.formula( paste(
    p$variabletomodel, ' ~ 1 ',
      ' + offset( log(data_offset)) ',
      ' + f( auid, model="bym2", graph=slot(sppoly, "nb"), group=year_factor, scale.model=TRUE, constr=TRUE, hyper=H$bym2, control.group=list(model="ar1", hyper=H$ar1_group)) '
  ))

  p$carstm_model_family = "poisson"
  fit = carstm_model( p=p, M=p$modeldata )



# -------------------------------------------------
# Model -- nonspatial, nontemporal -- inla -- poisson
  p$carstm_model_label = "nonseparable_space-time"
  p$variabletomodel = "totno"
  p$carstm_modelengine = "inla"
  p$carstm_model_formula = as.formula( paste(
    p$variabletomodel, ' ~ 1 ',
      ' + offset( log(data_offset)) ',
      ' + f( dyri, model="ar1", hyper=H$ar1 ) ',
      ' + f( inla.group( t, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
      ' + f( inla.group( z, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
      ' + f( inla.group( substrate.grainsize, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
      ' + f( inla.group( pca1, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
      ' + f( inla.group( pca2, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
      ' + f( auid, model="bym2", graph=slot(sppoly, "nb"), group=year_factor, scale.model=TRUE, constr=TRUE, hyper=H$bym2, control.group=list(model="ar1", hyper=H$ar1_group)) '
  ))

  p$carstm_model_family = "poisson"
  fit = carstm_model( p=p, M=p$modeldata )



# -------------------------------------------------
# Model -- nonspatial, nontemporal -- inla -- poisson
  p$carstm_model_label = "nonseparable_space-time_no_pca"
  p$variabletomodel = "totno"
  p$carstm_modelengine = "inla"
  p$carstm_model_formula = as.formula( paste(
    p$variabletomodel, ' ~ 1 ',
      ' + offset( log(data_offset)) ',
      ' + f( dyri, model="ar1", hyper=H$ar1 ) ',
      ' + f( inla.group( t, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
      ' + f( inla.group( z, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
      ' + f( inla.group( substrate.grainsize, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
      ' + f( auid, model="bym2", graph=slot(sppoly, "nb"), group=year_factor, scale.model=TRUE, constr=TRUE, hyper=H$bym2, control.group=list(model="ar1", hyper=H$ar1_group)) '
  ))

  p$carstm_model_family = "poisson"
  fit = carstm_model( p=p, M=p$modeldata )


# -------------------------------------------------
  p$carstm_model_label = "nonseparable_time-space"
  p$variabletomodel = "totno"
  p$carstm_modelengine = "inla"
  p$carstm_model_formula = as.formula( paste(
    p$variabletomodel, ' ~ 1 ',
      ' + offset( log(data_offset)) ',
      ' + f( year_factor, model="ar1", hyper=H$ar1, group=auid, control.group=list(model="besag", graph=slot(sppoly, "nb"), scale.model=TRUE ) ) ',
      ' + f( dyri, model="ar1", hyper=H$ar1 ) ',
      ' + f( inla.group( t, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
      ' + f( inla.group( z, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
      ' + f( inla.group( substrate.grainsize, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
      ' + f( inla.group( pca1, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
      ' + f( inla.group( pca2, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2) '
  ))

  p$carstm_model_family = "poisson"
  fit = carstm_model( p=p, M=p$modeldata )


# -------------------------------------------------
# Model -- nonspatial, nontemporal -- inla -- poisson 4.5 hrs
  p$carstm_model_label = "inla_zeroinflatedpoisson0_full"
  p$variabletomodel = "totno"
  p$carstm_modelengine = "inla"
  p$carstm_model_formula = as.formula( paste(
    p$variabletomodel, ' ~ 1 ',
      ' + offset( log(data_offset)) ',
      ' + f( dyri, model="ar1", hyper=H$ar1 ) ',
      ' + f( inla.group( t, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
      ' + f( inla.group( z, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
      ' + f( inla.group( substrate.grainsize, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
      ' + f( inla.group( pca1, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
      ' + f( inla.group( pca2, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
      ' + f( auid, model="bym2", graph=slot(sppoly, "nb"), group=year_factor, scale.model=TRUE, constr=TRUE, hyper=H$bym2, control.group=list(model="ar1", hyper=H$ar1_group)) '
  ))
  p$carstm_model_family =  "zeroinflatedpoisson0"
    fit = carstm_model( p=p, M=p$modeldata )



# -------------------------------------------------
  p$carstm_model_label = "inla_zeroinflatedpoisson1_full"  # strange
  p$variabletomodel = "totno"
  p$carstm_modelengine = "inla"
  p$carstm_model_formula = as.formula( paste(
    p$variabletomodel, ' ~ 1 ',
      ' + offset( log(data_offset)) ',
      ' + f( dyri, model="ar1", hyper=H$ar1 ) ',
      ' + f( inla.group( t, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
      ' + f( inla.group( z, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
      ' + f( inla.group( substrate.grainsize, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
      ' + f( inla.group( pca1, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
      ' + f( inla.group( pca2, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
      ' + f( auid, model="bym2", graph=slot(sppoly, "nb"), group=year_factor, scale.model=TRUE, constr=TRUE, hyper=H$bym2, control.group=list(model="ar1", hyper=H$ar1_group)) '
  ))

  p$carstm_model_family =  "zeroinflatedpoisson1"
  fit = carstm_model( p=p, M=p$modeldata )






# -------------------------------------------------
# generic calls

  if (0) {
    m = get("inla.models", INLA:::inla.get.inlaEnv())
    m$latent$rw2$min.diff = NULL
    assign("inla.models", m, INLA:::inla.get.inlaEnv())
  }


# extract results and examine
  fit =  carstm_model( p=p, DS="carstm_modelled_fit" )  # extract currently saved model fit
  s = summary(fit)
  s
  # s$dic$dic
  # s$dic$p.eff
  # AIC(fit)
  # DIC(fit)
  -sum(log(fit$cpo$cpo), na.rm=TRUE)




  # for mapping
  res = carstm_model( p=p, carstm_model_label=p$carstm_model_label, DS="carstm_modelled_summary" )


  if (0) {
    plot( fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE )
    plot( fit$marginals.hyperpar$"Phi for auid", type="l")  # posterior distribution of phi nonspatial dominates
    plot( fit$marginals.hyperpar$"Precision for auid", type="l")
    plot( fit$marginals.hyperpar$"Precision for setno", type="l")
  }


  sppoly = areal_units( p=p )  # to reload
  # M = snowcrab.db( p=p, DS="carstm_inputs" )
  # M$yr = M$year  # req for meanweights

  # # mean weight by auidxyear
  # wgts = meanweights_by_arealunit(
  #   set=M[M$tag=="observations",],
  #   AUID=as.character( sppoly$AUID ),
  #   yrs=p$yrs,
  #   fillall=TRUE,
  #   annual_breakdown=TRUE,
  #   robustify_quantiles=c(0, 0.99)  # high upper bounds are more dangerous
  # )


  # aggregate time series


  weight_year = meanweights_by_arealunit_modelled( p=p, redo=TRUE )  ## needed in carstm_output_compute

  RES = snowcrab.db(p=p, DS="carstm_output_compute"   )

  RES = snowcrab.db(p=p, DS="carstm_output_timeseries"  )
  bio = snowcrab.db(p=p, DS="carstm_output_spacetime_biomass"  )
  num = snowcrab.db(p=p, DS="carstm_output_spacetime_number"  )

  # plots with mean and 95% CI

    plot( cfaall_median ~ yrs, data=RES, lty="solid", lwd=4, pch=20, col="slateblue", type="b", ylab="Biomass (kt)", xlab="")
    lines( cfaall_lb ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )
    lines( cfaall_ub ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )

    plot( cfasouth_median ~ yrs, data=RES, lty="solid", lwd=4, pch=20, col="slateblue", type="b", ylab="Biomass (kt)", xlab="")
    lines( cfasouth_lb ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )
    lines( cfasouth_ub ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )

    plot( cfanorth_median ~ yrs, data=RES, lty="solid", lwd=4, pch=20, col="slateblue", type="b", ylab="Biomass (kt)", xlab="")
    lines( cfanorth_lb ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )
    lines( cfanorth_ub ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )

    plot( cfa4x_median ~ yrs, data=RES, lty="solid", lwd=4, pch=20, col="slateblue", type="b", ylab="Biomass (kt)", xlab="")
    lines( cfa4x_lb ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )
    lines( cfa4x_ub ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )

  p$boundingbox = list( xlim=p$corners$lon, ylim=p$corners$lat) # bounding box for plots using spplot
  p$coastLayout = aegis.coastline::coastline_layout(p=p)
  p$mypalette=RColorBrewer::brewer.pal(9, "YlOrRd")
  sppoly = areal_units( p=p )  # to reload

  # map it ..mean density

  yrc = as.character( 2000 )

  vn = "pred"
  dev.new()
  sppoly[,vn] = bio[,yrc]
  brks = interval_break(X= sppoly[[vn]], n=length(p$mypalette), style="quantile")
  spplot( sppoly, vn, col.regions=p$mypalette, main=vn, at=brks, sp.layout=p$coastLayout, col="transparent" )


  time_match = NULL
  time_match = list(year=yrc)
  time_match = list(year=yrc, dyear="0.85" )
  
  vn = paste(p$variabletomodel, "random_sample_iid", sep=".")
  carstm_map( res=res, vn=vn, time_match=time_match, 
    # at=seq(-2, 10, by=2),          
    sp.layout = p$coastLayout, 
    col.regions = p$mypalette, 
    main=paste( vn, paste0(time_match, collapse="-") )  )


  vn = paste(p$variabletomodel, "random_auid_nonspatial", sep=".")
  carstm_map( res=res, vn=vn,  time_match=time_match, 
    # at=seq(-2, 10, by=2),          
    sp.layout = p$coastLayout, 
    col.regions = p$mypalette, 
    main=paste( vn, paste0(time_match, collapse="-") )  )

  
  vn = paste(p$variabletomodel, "random_auid_spatial", sep=".")
  carstm_map( res=res, vn=vn,  time_match=time_match, 
    # at=seq(-2, 10, by=2),          
    sp.layout = p$coastLayout, 
    col.regions = p$mypalette, 
    main=paste( vn, paste0(time_match, collapse="-") )  )

  



  # -------------------------------
  ## --- comparisons

  cc =  as.data.frame( rbind(
    c("factorial_gaussian_biomass_glm", "Factorial crossed (Gaussian) biomass                                               ", "darkgreen", 26, "solid", 8, "l"  ),
    c("factorial_simple",        "Factorial crossed",        "slategray",  20, "solid", 4, "l"),
    c("factorial_full_covariates", "Factorial covariates",   "darkred",       20, "solid", 4, "l"),
    c("mixed_effects_simple",    "Mixed effects simple",     "turquoise",  21, "dotdash", 4, "l"),
    c("mixed_effects_static",    "Mixed effects static",     "darkorange", 22, "dotdash", 4, "l"),
    c("mixed_effects_dynamic",   "Mixed effects dynamic",    "blue",       23, "dotdash", 4, "l"),
    c("separable",               "Separable",                "lightgreen",      24, "dotted", 6, "l"),
    c("separable_simple",        "Separable simple",         "orange",    25, "dotted", 6, "l"),
    c("nonseparable_simple",     "Nonseparable simple",      "cyan",       26, "dashed", 6, "l"),
    c("nonseparable_space-time_no_pca",     "Nonseparable space|time no PCA",  "darkgray",  24, "dashed", 6, "l"),
    c("nonseparable_time-space", "Nonseparable time|space",  "darkgreen",  27, "dashed", 6, "l"),
    c("nonseparable_space-time", "Nonseparable space|time",  "slateblue",  20, "dashed", 6, "l")
  ), stringsAsFactors=FALSE)


  # cc =  as.data.frame( rbind(
  #   c("factorial_gaussian_biomass_glm", "Factorial crossed (Gaussian) biomass                               ", "darkgreen", 26, 1, 8, "l"  ),
  #   c("mixed_effects_dynamic",   "Mixed effects dynamic",    "blue",       23, 5, 4, "l"),
  #   c("covariates_only",     "Ecosystem covariates only", "darkorange", 26, 1, 4, "l"),
  #   c("nonseparable_simple",     "Nonseparable simple",      "cyan",       26, 1, 4, "l"),
  #   c("nonseparable_space-time_no_pca",     "Nonseparable space|time no PCA",      "green",      24, 6, 4, "l"),
  #   c("nonseparable_space-time", "Nonseparable space|time",  "slateblue",  20, 1, 8, "l")
  # ), stringsAsFactors=FALSE)



#    c( "inla_zeroinflatedpoisson0_full", "Nonseparable overdispersed space|time" ,  "orange",  20, 1, 8, "l")

  colnames(cc) = c("tocompare", "legend", "col", "pch", "lty", "lwd", "type" )
  cc$pch = as.numeric(cc$pch)
  # cc$lty = as.character(cc$lty)
  cc$lwd = as.numeric(cc$lwd)

  res_ts = list()
  for ( lab in cc$tocompare ) {
    p$carstm_modelengine = "inla"
    p$variabletomodel = "totno"
    if (lab == "factorial_gaussian_biomass_glm")  {
      p$carstm_modelengine = "glm"
      p$variabletomodel = "totwgt"
    }
    p$carstm_model_label=lab
    res_ts[[lab]] = snowcrab.db(p=p, DS="carstm_output_timeseries"  )
  }

  dev.new(width=11, height=7)
  i=1 ; plot( cfaall_median  ~ yrs, data=res_ts[[i]], lty=cc$lty[i], lwd=cc$lwd[i], col=cc$col[i], pch=cc$pch[i], type=cc$type[i], ylim=c(0,185), xlab="Year", ylab="kt")
  for (i in 2:length(res_ts)) {
      lines( cfaall_median ~ yrs, data=res_ts[[i]], lty=cc$lty[i], lwd=cc$lwd[i], col=cc$col[i], pch=cc$pch[i], type=cc$type[i])
  }
  legend("top", legend=cc$legend, lty=cc$lty, col=cc$col, lwd=cc$lwd, bty="n"  )




#############


# Presence-absence analysis
# .. these require alternate data inputs for each class of snow crab and so data needs to be created for each



  p = bio.snowcrab::snowcrab_parameters( project_class="carstm", assessment.years=1999:2018 )

  p$modeldata = 'snowcrab.db( p=p, DS="carstm_inputs" )'



# -------------------------------------------------
# Model -- nonspatial, nontemporal -- inla -- binomial -- immature all
  p$selection$type = "presence_absence"
  p$selection$biologicals=list(
    spec_bio=bio.taxonomy::taxonomy.recode( from="spec", to="parsimonious", tolookup=p$groundfish_species_code ),
    # sex=1, # female
    mat=0, # do not use maturity status in groundfish data as it is suspect ..
    len= c( 1, 50 )/10, #  mm -> cm ; aegis_db in cm
    ranged_data="len"
  )
  p$selection$survey=list(
    data.source = c("snowcrab"),
    yr = p$assessment.years,      # time frame for comparison specified above
    settype = 1, # same as geartype in groundfish_survey_db
    polygon_enforce=TRUE,  # make sure mis-classified stations or incorrectly entered positions get filtered out
    strata_toremove = NULL #,  # emphasize that all data enters analysis initially ..
    # ranged_data = c("dyear")  # not used .. just to show how to use range_data
  )
  p$variabletomodel = "pa"
  p$carstm_model_label = "nonseparable_space-time_pa_immature_small"
  p$carstm_modelengine = "inla"
  p$carstm_model_formula = as.formula( paste(
    p$variabletomodel, ' ~ 1 ',
      ' + f( dyri, model="ar1", hyper=H$ar1 ) ',
      ' + f( inla.group( t, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
      ' + f( inla.group( z, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
      ' + f( inla.group( substrate.grainsize, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
      ' + f( inla.group( pca1, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
      ' + f( inla.group( pca2, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
      ' + f( auid, model="bym2", graph=slot(sppoly, "nb"), group=year_factor, scale.model=TRUE, constr=TRUE, hyper=H$bym2, control.group=list(model="ar1", hyper=H$ar1_group)) '
  ))

  p$carstm_model_family = "binomial"  # "nbinomial", "betabinomial", "zeroinflatedbinomial0" , "zeroinflatednbinomial0"
  p$carstm_model_inla_control_familiy = list(control.link=list(model='logit'))

  M = snowcrab.db( p=p, DS="carstm_inputs", redo=TRUE )  # will redo if not found
  M = NULL; gc()
  fit = carstm_model( p=p, M=p$modeldata )





# -------------------------------------------------
# Model -- nonspatial, nontemporal -- inla -- binomial -- fishable component
  p$selection = list()
  p$selection$type = "presence_absence"
  p$selection$biologicals=list(
    spec_bio=bio.taxonomy::taxonomy.recode( from="spec", to="parsimonious", tolookup=p$groundfish_species_code ),
    sex=0, # male
    mat=1, # do not use maturity status in groundfish data as it is suspect ..
    len= c( 95, 200 )/10, #  mm -> cm ; aegis_db in cm
    ranged_data="len"
  )
  p$selection$survey=list(
    data.source = c("snowcrab", "groundfish", "logbook"),
    yr = p$assessment.years,      # time frame for comparison specified above
    settype = 1, # same as geartype in groundfish_survey_db
    polygon_enforce=TRUE,  # make sure mis-classified stations or incorrectly entered positions get filtered out
    strata_toremove = NULL #,  # emphasize that all data enters analysis initially ..
    # ranged_data = c("dyear")  # not used .. just to show how to use range_data
  )
  p$variabletomodel = "pa"
  p$carstm_model_label = "nonseparable_space-time_pa_fishable"
  p$carstm_modelengine = "inla"
  p$carstm_model_formula = as.formula( paste(
    p$variabletomodel, ' ~ 1 ',
      ' + f( dyri, model="ar1", hyper=H$ar1 ) ',
      ' + f( inla.group( t, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
      ' + f( inla.group( z, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
      ' + f( inla.group( substrate.grainsize, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
      ' + f( inla.group( pca1, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
      ' + f( inla.group( pca2, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
      ' + f( auid, model="bym2", graph=slot(sppoly, "nb"), group=year_factor, scale.model=TRUE, constr=TRUE, hyper=H$bym2, control.group=list(model="ar1", hyper=H$ar1_group)) '
  ))

  p$carstm_model_family = "binomial"  # "nbinomial", "betabinomial", "zeroinflatedbinomial0" , "zeroinflatednbinomial0"
  p$carstm_model_inla_control_familiy = list(control.link=list(model='logit'))

  M = snowcrab.db( p=p, DS="carstm_inputs", redo=TRUE )  # will redo if not found
  M = NULL; gc()
  fit = carstm_model( p=p, M=p$modeldata )




# -------------------------------------------------
# Model -- nonspatial, nontemporal -- inla -- binomial -- mature females
  p$selection$type = "presence_absence"
  p$selection$biologicals=list(
    spec_bio=bio.taxonomy::taxonomy.recode( from="spec", to="parsimonious", tolookup=p$groundfish_species_code ),
    sex=1, # female
    mat=1, # do not use maturity status in groundfish data as it is suspect ..
    len= c( 30, 96 )/10, #  mm -> cm ; aegis_db in cm
    ranged_data="len"
  )
  p$selection$survey=list(
    data.source = c("snowcrab"),
    yr = p$assessment.years,      # time frame for comparison specified above
    settype = 1, # same as geartype in groundfish_survey_db
    polygon_enforce=TRUE,  # make sure mis-classified stations or incorrectly entered positions get filtered out
    strata_toremove = NULL #,  # emphasize that all data enters analysis initially ..
    # ranged_data = c("dyear")  # not used .. just to show how to use range_data
  )
  p$variabletomodel = "pa"
  p$carstm_model_label = "nonseparable_space-time_pa_mature_females"
  p$carstm_modelengine = "inla"
  p$carstm_model_formula = as.formula( paste(
    p$variabletomodel, ' ~ 1 ',
      ' + f( dyri, model="ar1", hyper=H$ar1 ) ',
      ' + f( inla.group( t, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
      ' + f( inla.group( z, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
      ' + f( inla.group( substrate.grainsize, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
      ' + f( inla.group( pca1, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
      ' + f( inla.group( pca2, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
      ' + f( auid, model="bym2", graph=slot(sppoly, "nb"), group=year_factor, scale.model=TRUE, constr=TRUE, hyper=H$bym2, control.group=list(model="ar1", hyper=H$ar1_group)) '
  ) )

  p$carstm_model_family = "binomial"  # "nbinomial", "betabinomial", "zeroinflatedbinomial0" , "zeroinflatednbinomial0"
  p$carstm_model_inla_control_familiy = list(control.link=list(model='logit'))

  M = snowcrab.db( p=p, DS="carstm_inputs", redo=TRUE )  # will redo if not found
  M = NULL; gc()
  fit = carstm_model( p=p, M=p$modeldata )





# -------------------------------------------------
# Model -- nonspatial, nontemporal -- inla -- binomial -- immature all
  p$selection$type = "presence_absence"
  p$selection$biologicals=list(
    spec_bio=bio.taxonomy::taxonomy.recode( from="spec", to="parsimonious", tolookup=p$groundfish_species_code ),
    # sex=1, # female
    mat=0, # do not use maturity status in groundfish data as it is suspect ..
    len= c( 1, 200 )/10, #  mm -> cm ; aegis_db in cm
    ranged_data="len"
  )
  p$selection$survey=list(
    data.source = c("snowcrab"),
    yr = p$assessment.years,      # time frame for comparison specified above
    settype = 1, # same as geartype in groundfish_survey_db
    polygon_enforce=TRUE,  # make sure mis-classified stations or incorrectly entered positions get filtered out
    strata_toremove = NULL #,  # emphasize that all data enters analysis initially ..
    # ranged_data = c("dyear")  # not used .. just to show how to use range_data
  )
  p$variabletomodel = "pa"
  p$carstm_model_label = "nonseparable_space-time_pa_immature"
  p$carstm_modelengine = "inla"
  p$carstm_model_formula = as.formula( paste(
    p$variabletomodel, ' ~ 1 ',
      ' + f( dyri, model="ar1", hyper=H$ar1 ) ',
      ' + f( inla.group( t, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
      ' + f( inla.group( z, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
      ' + f( inla.group( substrate.grainsize, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
      ' + f( inla.group( pca1, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
      ' + f( inla.group( pca2, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
      ' + f( auid, model="bym2", graph=slot(sppoly, "nb"), group=year_factor, scale.model=TRUE, constr=TRUE, hyper=H$bym2, control.group=list(model="ar1", hyper=H$ar1_group)) '
  ) )

  p$carstm_model_family = "binomial"  # "nbinomial", "betabinomial", "zeroinflatedbinomial0" , "zeroinflatednbinomial0"
  p$carstm_model_inla_control_familiy = list(control.link=list(model='logit'))

  M = snowcrab.db( p=p, DS="carstm_inputs", redo=TRUE )  # will redo if not found
  M = NULL; gc()
  fit = carstm_model( p=p, M=p$modeldata )




# -------------------------------------------------
# Model -- nonspatial, nontemporal -- inla -- negative binomial -- fishable component
  p$selection$type = "presence_absence"
  p$selection$biologicals=list(
    spec_bio=bio.taxonomy::taxonomy.recode( from="spec", to="parsimonious", tolookup=p$groundfish_species_code ),
    sex=0, # male
    mat=1, # do not use maturity status in groundfish data as it is suspect ..
    len= c( 95, 200 )/10, #  mm -> cm ; aegis_db in cm
    ranged_data="len"
  )
  p$selection$survey=list(
    data.source = c("snowcrab", "groundfish", "logbook"),
    yr = p$assessment.years,      # time frame for comparison specified above
    settype = 1, # same as geartype in groundfish_survey_db
    polygon_enforce=TRUE,  # make sure mis-classified stations or incorrectly entered positions get filtered out
    strata_toremove = NULL #,  # emphasize that all data enters analysis initially ..
    # ranged_data = c("dyear")  # not used .. just to show how to use range_data
  )
  p$variabletomodel = "pa"
  p$carstm_model_label = "nonseparable_space-time_pa_fishable_nbinomial"
  p$carstm_modelengine = "inla"
  p$carstm_model_formula = as.formula( paste(
    p$variabletomodel, ' ~ 1 ',
        + f( dyri, model="ar1", hyper=H$ar1 )
        + f( inla.group( t, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2)
        + f( inla.group( z, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2)
        + f( inla.group( substrate.grainsize, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2)
        + f( inla.group( pca1, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2)
        + f( inla.group( pca2, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2)
        + f( auid, model="bym2", graph=slot(sppoly, "nb"), group=year_factor, scale.model=TRUE, constr=TRUE, hyper=H$bym2, control.group=list(model="ar1", hyper=H$ar1_group)),
  ) )

  p$carstm_model_family = "nbinomial"  # "nbinomial", "betabinomial", "zeroinflatedbinomial0" , "zeroinflatednbinomial0"
  p$carstm_model_inla_control_familiy = NULL

  M = snowcrab.db( p=p, DS="carstm_inputs", redo=TRUE )  # will redo if not found
  M = NULL; gc()
  fit = carstm_model( p=p, M=p$modeldata )





# -------------------------------------------------
# Model -- nonspatial, nontemporal -- inla -- negative binomial -- fishable component
  p$selection$type = "presence_absence"
  p$selection$biologicals=list(
    spec_bio=bio.taxonomy::taxonomy.recode( from="spec", to="parsimonious", tolookup=p$groundfish_species_code ),
    sex=0, # male
    mat=1, # do not use maturity status in groundfish data as it is suspect ..
    len= c( 95, 200 )/10, #  mm -> cm ; aegis_db in cm
    ranged_data="len"
  )
  p$selection$survey=list(
    data.source = c("snowcrab", "groundfish", "logbook"),
    yr = p$assessment.years,      # time frame for comparison specified above
    settype = 1, # same as geartype in groundfish_survey_db
    polygon_enforce=TRUE,  # make sure mis-classified stations or incorrectly entered positions get filtered out
    strata_toremove = NULL #,  # emphasize that all data enters analysis initially ..
    # ranged_data = c("dyear")  # not used .. just to show how to use range_data
  )
  p$variabletomodel = "pa"
  p$carstm_model_label = "nonseparable_space-time_pa_fishable_zeroinflatedbinomial0"
  p$carstm_modelengine = "inla"
  p$carstm_model_formula = as.formula( paste(
    p$variabletomodel, ' ~ 1 ',
      ' + f( dyri, model="ar1", hyper=H$ar1 ) ',
      ' + f( inla.group( t, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
      ' + f( inla.group( z, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
      ' + f( inla.group( substrate.grainsize, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
      ' + f( inla.group( pca1, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
      ' + f( inla.group( pca2, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
      ' + f( auid, model="bym2", graph=slot(sppoly, "nb"), group=year_factor, scale.model=TRUE, constr=TRUE, hyper=H$bym2, control.group=list(model="ar1", hyper=H$ar1_group))  '
  ) )

  p$carstm_model_family  = "zeroinflatedbinomial0", #  "binomial",  # "nbinomial", "betabinomial", "zeroinflatedbinomial0" , "zeroinflatednbinomial0"
  p$carstm_model_inla_control_familiy = NULL

  M = snowcrab.db( p=p, DS="carstm_inputs", redo=TRUE )  # will redo if not found
  M = NULL; gc()
  fit = carstm_model( p=p, M=p$modeldata )





# -------------------------------------------------
# Model -- nonspatial, nontemporal -- inla -- negative binomial -- fishable component
  p$selection$type = "presence_absence"
  p$selection$biologicals=list(
    spec_bio=bio.taxonomy::taxonomy.recode( from="spec", to="parsimonious", tolookup=p$groundfish_species_code ),
    sex=0, # male
    mat=1, # do not use maturity status in groundfish data as it is suspect ..
    len= c( 95, 200 )/10, #  mm -> cm ; aegis_db in cm
    ranged_data="len"
  )
  p$selection$survey=list(
    data.source = c("snowcrab", "groundfish", "logbook"),
    yr = p$assessment.years,      # time frame for comparison specified above
    settype = 1, # same as geartype in groundfish_survey_db
    polygon_enforce=TRUE,  # make sure mis-classified stations or incorrectly entered positions get filtered out
    strata_toremove = NULL #,  # emphasize that all data enters analysis initially ..
    # ranged_data = c("dyear")  # not used .. just to show how to use range_data
  )
  p$variabletomodel = "pa"
  p$carstm_model_label = "nonseparable_space-time_pa_fishable_zeroinflatedbinomial1"
  p$carstm_modelengine = "inla"
  p$carstm_model_formula = as.formula( paste(
    p$variabletomodel, ' ~ 1 ',
      ' + f( dyri, model="ar1", hyper=H$ar1 ) ',
      ' + f( inla.group( t, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
      ' + f( inla.group( z, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
      ' + f( inla.group( substrate.grainsize, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
      ' + f( inla.group( pca1, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
      ' + f( inla.group( pca2, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
      ' + f( auid, model="bym2", graph=slot(sppoly, "nb"), group=year_factor, scale.model=TRUE, constr=TRUE, hyper=H$bym2, control.group=list(model="ar1", hyper=H$ar1_group)) '
  ) )

  p$carstm_model_family  = "zeroinflatedbinomial1", #  "binomial",  # "nbinomial", "betabinomial", "zeroinflatedbinomial0" , "zeroinflatednbinomial0"
  p$carstm_model_inla_control_familiy = NULL



  M = snowcrab.db( p=p, DS="carstm_inputs", redo=TRUE )  # will redo if not found
  M = NULL; gc()
  fit = carstm_model( p=p, M=p$modeldata )



# -------------------------------------------------

# -------------------------------------------------



# -------------------------------------------------
# generic calls


# extract results and examine
  fit =  carstm_model( p=p, DS="carstm_modelled_fit" )  # extract currently saved model fit

  res = carstm_model( p=p, carstm_model_label=p$carstm_model_label, DS="carstm_modelled_summary" )


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
  #
  -sum(log(fit$cpo$cpo), na.rm=TRUE)



  # surface area weighted average
  RES = snowcrab.db(p=p, DS="carstm_output_compute"  )
  RES = snowcrab.db(p=p, DS="carstm_output_timeseries"  )

  pa = snowcrab.db(p=p, DS="carstm_output_spacetime_pa"  )

# plots with 95% PI

    plot( cfaall ~ yrs, data=RES, lty="solid", lwd=4, pch=20, col="slateblue", type="b", ylab="Prob of observing snow crab", xlab="", ylim=c(0.2,0.4))
    lines( cfaall_lb ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )
    lines( cfaall_ub ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )

    plot( cfasouth ~ yrs, data=RES, lty="solid", lwd=4, pch=20, col="slateblue", type="b", ylab="Prob of observing snow crab", xlab="", ylim=c(0,1))
    lines( cfasouth_lb ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )
    lines( cfasouth_ub ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )

    plot( cfanorth ~ yrs, data=RES, lty="solid", lwd=4, pch=20, col="slateblue", type="b", ylab="Prob of observing snow crab", xlab="", ylim=c(0,1))
    lines( cfanorth_lb ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )
    lines( cfanorth_ub ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )

    plot( cfa4x ~ yrs, data=RES, lty="solid", lwd=4, pch=20, col="slateblue", type="b", ylab="Prob of observing snow crab", xlab="", ylim=c(0,1))
    lines( cfa4x_lb ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )
    lines( cfa4x_ub ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )


  p$boundingbox = list( xlim=p$corners$lon, ylim=p$corners$lat) # bounding box for plots using spplot
  p$coastLayout = aegis.coastline::coastline_layout(p=p)
  p$mypalette=RColorBrewer::brewer.pal(9, "YlOrRd")
  sppoly = areal_units( p=p )  # to reload

  # map it ..mean density
  vn = "pred"
  sppoly[,vn] = pa[,"2017"]
  brks = interval_break(X= sppoly[[vn]], n=length(p$mypalette), style="quantile")
  spplot( sppoly, vn, col.regions=p$mypalette, main=vn, at=brks, sp.layout=p$coastLayout, col="transparent", sub= p$carstm_model_label )





  time_match = NULL
  time_match = list(year=yrc)
  time_match = list(year=yrc, dyear="0.85" )
  
  vn = paste(p$variabletomodel, "random_sample_iid", sep=".")
  carstm_map( res=res, vn=vn, time_match=time_match, 
    # at=seq(-2, 10, by=2),          
    sp.layout = p$coastLayout, 
    col.regions = p$mypalette, 
    main=paste( vn, paste0(time_match, collapse="-") )  )


  vn = paste(p$variabletomodel, "random_auid_nonspatial", sep=".")
  carstm_map( res=res, vn=vn, time_match=time_match, 
    # at=seq(-2, 10, by=2),          
    sp.layout = p$coastLayout, 
    col.regions = p$mypalette, 
    main=paste( vn, paste0(time_match, collapse="-") )  )

  
  vn = paste(p$variabletomodel, "random_auid_spatial", sep=".")
  carstm_map( res=res, vn=vn, time_match=time_match, 
    # at=seq(-2, 10, by=2),          
    sp.layout = p$coastLayout, 
    col.regions = p$mypalette, 
    main=paste( vn, paste0(time_match, collapse="-") )  )





  # -------------------------------
  ## --- comparisons ... not yet chosen/completed .. just a copy from abundance ..

  if (0) {

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
        p$carstm_model_label=lab
        res_ts[[lab]] = snowcrab.db(p=p, DS="carstm_output_timeseries"  )
      }

      dev.new(width=11, height=7)
      i=1 ; plot( cfaall  ~ yrs, data=res_ts[[i]], lty=cc$lty[i], lwd=cc$lwd[i], col=cc$col[i], pch=cc$pch[i], type=cc$type[i], ylim=c(0,185), xlab="Year", ylab="kt")
      for (i in 2:length(res_ts)) {
        lines( cfaall ~ yrs, data=res_ts[[i]], lty=cc$lty[i], lwd=cc$lwd[i], col=cc$col[i], pch=cc$pch[i], type=cc$type[i])
      }
      legend("top", legend=cc$legend, lty=cc$lty, col=cc$col, lwd=cc$lwd, bty="n"  )

  }







## Plot maps of residuals of numbers per set obs vs pred


year.assessment = 2018
#To add a title to any carstm_map, please see below example
#carstm_map( res=res, vn=vn, main=list(label="my plot title", cex=2) )


# -------------------------------------------------
# Part 1 -- construct basic parameter list defining the main characteristics of the study
# require(aegis)

 p = bio.snowcrab::snowcrab_parameters( project_class="carstm", assessment.years=1999:year.assessment )


  plot.dir=paste(p$modeldir,"prediction.plots", year.assessment, sep="/" )



# extract results and examine



  fit =  carstm_model( p=p, DS="carstm_modelled_fit" )  # extract currently saved model fit
  summary(fit)

  res = carstm_model( p=p, DS="carstm_modelled_summary"  )



  # prediction surface
  crs_lonlat = st_crs(projection_proj4string("lonlat_wgs84"))

  sppoly = areal_units( p=p )  # will redo if not found
  sppoly = st_transform(sppoly, crs=crs_lonlat )


  # do this immediately to reduce storage for sppoly (before adding other variables)
  M = snowcrab.db( p=p, DS="biological_data" )  # will redo if not found .. not used here but used for data matching/lookup in other aegis projects that use bathymetry
  M$tiyr=lubridate::decimal_date(M$timestamp)

  # M$totno = M$totno_adjusted / M$cf_set_no   # convert density to counts
  # M$totwgt = M$totwgt_adjusted / M$cf_set_mass # convert density to total wgt

  # M$data_offset = 1 / M$cf_set_no  ## offset only used in poisson model


  # reduce size
  M = M[ which( M$lon > p$corners$lon[1] & M$lon < p$corners$lon[2]  & M$lat > p$corners$lat[1] & M$lat < p$corners$lat[2] ), ]
  # levelplot(z.mean~plon+plat, data=M, aspect="iso")

  M$AUID = st_points_in_polygons(
    pts = st_as_sf( M, coords=c("lon","lat"), crs=crs_lonlat ),
    polys = sppoly[, "AUID"],
    varname="AUID"
  )

  M = M[!is.na(M$AUID),]

  names(M)[which(names(M)=="yr") ] = "year"
  # M = M[ which(M$year %in% p$yrs), ]
  # M$tiyr = lubridate::decimal_date ( M$timestamp )
  # M$dyear = M$tiyr - M$year

  MM = res$M

  obsMM = MM[MM$tag=="observations",]
  plocs = MM[MM$tag=="predictions",]

  obs = M
  if (p$variabletomodel=="totno") {
    obs$density = obs$totno / obs$data_offset
    rr = as.data.frame.table(res$totno.predicted)
  }
  if (p$variabletomodel=="totwgt") {
    obs$density = obs$totwgt / obs$data_offset
    rr = as.data.frame.table(res$totwgt.predicted)
  }

  rr$AUID = as.character( rr$AUID)


region.id = slot( slot(sppoly, "nb"), "region.id" )
rr$auid = match( rr$AUID, region.id )

rr$AUID_character = rr$AUID
rr$AUID = rr$auid  -- temp fix to get numbers to match bym geometry


  rr$year = as.numeric( as.character( rr$year) )

  obs = merge( obs, rr, by=c("year", "AUID"), all.x=TRUE, all.y=FALSE )
  obs$Freq[ !is.finite(obs$Freq) ] = 0
  obs$resid =  obs$Freq - obs$density
  obs$resid_per_set = obs$resid * obs$data_offset
  obs$yr = obs$year


  vn = "resid"
  vn = "resid_per_set"
  #er = range( obs[,vn], na.rm=T) * c(0.95, 1.05)
  er = c(-100, 100)

  resol = p$pres

  B = bathymetry_db(p=p, DS="baseline")  # 1 km (p$pres )

  for ( y in  2000:2018 ) {
      ii = which( obs$yr==y & is.finite(obs[,vn] ))
      if ( length(ii) > 3 ) {
      dir.create( file.path( p$project.outputdir, "residuals", p$carstm_model_label), recursive=TRUE, showWarnings =FALSE)
      fn = file.path( p$project.outputdir, "residuals", p$carstm_model_label, paste( "residuals", y, "png", sep=".") )
      png( filename=fn, width=3072, height=2304, pointsize=40, res=300 )
       lp = map_simple( toplot=obs[ ii, c("plon","plat", vn) ], plotarea=B, resol=1, theta=15, filterdistances=7.5, vn=vn, annot=paste("Residuals", y), er=er )
       print(lp)
      dev.off()
      print(fn)
    }
  }




## ---------
#### final estimation of biomass via fishery models and associated figures and tables:

#Pick whichever year reference below is correct (most often year.assessment...-1)
if (!exists("year.assessment")) {
   year.assessment=lubridate::year(Sys.Date())
   year.assessment=lubridate::year(Sys.Date()) - 1
}

year.assessment = 2018
p = bio.snowcrab::load.environment(
  year.assessment=year.assessment,
  assessment_years = 2000:year.assessment,
  vars.tomodel="R0.mass"
)



 
  p$fishery_model = fishery_model( DS = "logistic_parameters", p=p, tag=p$areal_units_type )
  p$fishery_model$stancode = stan_initialize( stan_code=fishery_model( p=p, DS="stan_surplus_production" ) )

  # override : 
  p$fishery_model$outdir = file.path(project.datadirectory('bio.snowcrab'), "assessments", p$year.assessment, "nonseparable_simple" )

  # force use of specific carstm results
  p$carstm_model_label = "nonseparable_simple"
  p$fishery_model$standata = bio.snowcrab::fishery_model( DS="data_aggregated_timeseries", p=p, assessment_years=p$assessment_years )

  str( p$fishery_model)

  res = fishery_model( 
    DS="logistic_model", 
    p=p, 
    tag=p$areal_units_type,
    # from here down are params for cmdstanr::sample()
    data=p$fishery_model$standata, 
    iter_warmup = 14000,
    iter_sampling = 5000,
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






# -------------------------------------------------
# end
# -------------------------------------------------



