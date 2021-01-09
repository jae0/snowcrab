

# snow crab using alt AUs
# using simplified work flow (AU's used consistently across all covariates)
# -- only AU's of interest ... not connected to a global analysis
# -- losing information outside of study area but faster
# -- results cannot be reused outside of workflow
#
# --- note this only predicts to time slices where it is required
# -------------------------------------------------
# set up paramerters

# -------------------------------------------------
# 1. set up polygon parameters

  yrs = 1999:2020


  areal_units_type = "tesselation"
  areal_units_resolution_km = 1
  for ( areal_units_constraint_nmin in c( 0, 1, 3, 5, 10, 15, 20, 25, 30, 40, 50  ) ) {
      p = bio.snowcrab::snowcrab_carstm(
        DS="parameters",
        assessment.years=yrs,
        modeldir = project.datadirectory("bio.snowcrab", "modelled", "testing" ),  ## <--- important: specify save location
        carstm_model_label = paste( "testing", areal_units_type, areal_units_resolution_km, areal_units_constraint_nmin, sep="_" ),
        inputdata_spatial_discretization_planar_km = 1,
        boundingbox = list( xlim = c(-70.5, -56.5), ylim=c(39.5, 47.5)), # bounding box for plots using spplot
        areal_units_proj4string_planar_km = projection_proj4string("utm20"), # set up default map projection
        areal_units_constraint = "snowcrab",
        areal_units_constraint_nmin = areal_units_constraint_nmin,
        areal_units_type= areal_units_type,
        areal_units_resolution_km = areal_units_resolution_km,
        sa_threshold_km2 = 5,
        inla_num.threads = 4,
        inla_blas.num.threads = 4
      )
    sppoly = areal_units( p=p, redo=TRUE )  # to create
  }



  areal_units_type = "lattice"
  for ( areal_units_resolution_km in c( 5, 10, 15, 20, 25, 30, 40, 50, 60, 80, 100 ) ) {
    for ( areal_units_constraint_nmin in c( 0, 3, 5, 10, 15  ) ) {
        p = bio.snowcrab::snowcrab_carstm(
          DS="parameters",
          assessment.years=yrs,
          modeldir = project.datadirectory("bio.snowcrab", "modelled", "testing" ),  ## <--- important: specify save location
          carstm_model_label = paste( "testing", areal_units_type, areal_units_resolution_km, areal_units_constraint_nmin, sep="_" ),
          inputdata_spatial_discretization_planar_km = 1,
          boundingbox = list( xlim = c(-70.5, -56.5), ylim=c(39.5, 47.5)), # bounding box for plots using spplot
          areal_units_proj4string_planar_km = projection_proj4string("utm20"), # set up default map projection
          areal_units_constraint = "snowcrab",
          areal_units_constraint_nmin = areal_units_constraint_nmin,
          areal_units_type= areal_units_type,
          areal_units_resolution_km = areal_units_resolution_km,
          sa_threshold_km2 = areal_units_resolution_km/2,
          inla_num.threads = 4,
          inla_blas.num.threads = 4
        )
      sppoly = areal_units( p=p, redo=TRUE )  # to create
    }
  }



if (0) {
  W.nb = poly2nb(sppoly, row.names=sppoly$internal_id, queen=TRUE)  # slow .. ~1hr?
  plot(W.nb, st_geometry(sppoly))
}


  #  areal_units_type = "lattice"
  areal_units_type = "tesselation"

  # 5 - 10 works well .. mean and variances stabilize
  # 3 is too low for ts analysis,...  esp for temperature
  # areal_units_constraint_nmin  = floor(length(yrs) / 3) # = 6
  # areal_units_constraint_nmin = 3
  # areal_units_constraint_nmin = length(yrs)
  # areal_units_constraint_nmin = 4
  # areal_units_constraint_nmin = 10
  # areal_units_constraint_nmin = 15
  # areal_units_constraint_nmin = 30
  areal_units_constraint_nmin = 10

  #   areal_units_resolution_km = 2
  #   areal_units_resolution_km = 5
  #   areal_units_resolution_km = 10 # too disconnected (at nmin=3)
  #   areal_units_resolution_km = 15 # well connected at nmin=3,
  #  areal_units_resolution_km =  20
  #   areal_units_resolution_km =  25  # default snow crab
  #   areal_units_resolution_km =  30
  #   areal_units_resolution_km =  40
  areal_units_resolution_km = 1

  p = bio.snowcrab::snowcrab_carstm(
    DS="parameters",
    assessment.years=1999:2019,
    modeldir = project.datadirectory("bio.snowcrab", "modelled", "testing" ),  ## <--- important: specify save location
    carstm_model_label = paste( "testing", areal_units_type, areal_units_resolution_km, areal_units_constraint_nmin, sep="_" ),
    inputdata_spatial_discretization_planar_km = 1,
    boundingbox = list( xlim = c(-70.5, -56.5), ylim=c(39.5, 47.5)), # bounding box for plots using spplot
    areal_units_proj4string_planar_km = projection_proj4string("utm20"), # set up default map projection
    areal_units_constraint = "snowcrab",
    areal_units_constraint_nmin = areal_units_constraint_nmin,
    areal_units_type= areal_units_type,
    areal_units_resolution_km = areal_units_resolution_km,
    sa_threshold_km2 = 5,
    inla_num.threads = 1,
    inla_blas.num.threads = 1
  )


  (p$carstm_model_label ) # required to use for abundance estimation
  # [1] carstm_model_label="testing_snowcrab_polygons_tesselation_1_3"

  if (0) coastLayout = aegis.coastline::coastline_layout( p=p, redo=TRUE )

  p = c(p, aegis.coastline::coastline_layout( p=p ) )
  # p$mypalette = RColorBrewer::brewer.pal(9, "YlOrRd")
  p$mypalette = RColorBrewer::brewer.pal(9, "Spectral")[9:1]


  # p$areal_units_proj4string_planar_km = projection_proj4string("omerc_nova_scotia")   # oblique mercator, centred on Scotian Shelf rotated by 325 degrees



  # create polygons
  if (0) {
    # ensure if polys exist and create if required
    for (au in c("cfanorth", "cfasouth", "cfa4x", "cfaall" )) {
      plot(polygon_managementareas( species="snowcrab", au))
    }
  }

  sppoly = areal_units( p=p )


  if (0) {
    # verify that params are correct
    sppoly = areal_units( p=p, redo=TRUE )  # to create

    plot(sppoly[,"AUID"], col="orange")
    plot( slot(sppoly, "nb"), coords=st_centroid(st_geometry( as(sppoly, "sf")) ),  col="green", add=T )

    dev.new(); spplot( sppoly, "au_sa_km2", main="SA_km2", sp.layout=p$coastLayout )
    dev.new(); spplot( sppoly, "au_sa_km2", main="AUID", sp.layout=p$coastLayout,  col.regions=RColorBrewer::brewer.pal(8, "Accent") )
    length(sppoly)
    summary(sppoly$au_sa_km2)

  }




# -------------------------------------------------
# Part 3 -- create covariate field for bathymetry
# bathymetry -- ensure the data assimilation in bathymetry is first completed :: 01.bathymetry_data.R

  pB = aegis.bathymetry::bathymetry_parameters( p=parameters_reset(p), project_class="carstm"  )
  # M = bathymetry_db( p=pB, DS="aggregated_data" , redo=TRUE )  # not used if raw data are being used
  M = NULL; gc()
  M = bathymetry_db( p=pB, DS="carstm_inputs", redo=TRUE  ) # will redo if not found
  M = NULL; gc()

  fit = carstm_model( p=pB, M='bathymetry_db( p=p, DS="carstm_inputs" )', DS="redo"  ) # run model and obtain predictions

  if (0) {
    # to use a saved instance
    fit = carstm_model( p=pB, DS="carstm_modelled_fit" )  # extract currently saved model fit
    summary(fit)

    res = carstm_summary( p=pB )  # load summary

    # maps of some of the results
    vn = paste(pB$variabletomodel, "predicted", sep=".")
    zplot = carstm_plot( p=pB, res=res, vn=vn )
    print(zplot)

     #to save map of predicted

    if (0) {
      fn=paste("z.predicted", year.assesment,"pdf", sep="." )
      pdf(zplot, file=paste(plot.dir, fn, sep="/"))
    }

    vn = paste(pB$variabletomodel, "predicted_se", sep=".")
    zplot=carstm_plot( p=pB, res=res, vn=vn )
    print(zplot)

    vn = paste(pB$variabletomodel, "random_auid_nonspatial", sep=".")
    carstm_plot( p=pB, res=res, vn=vn )

    vn = paste(pB$variabletomodel, "random_auid_spatial", sep=".")
    carstm_plot( p=pB, res=res, vn=vn )

    plot(fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE, single=TRUE )

  }


# -------------------------------------------------
# Part 4 -- create covariate field for  substrate
# ensure the data assimilation in substrate is first completed :: 01.substrate_data.R

  pS = substrate_parameters( p=parameters_reset(p), project_class="carstm"  )
  M = substrate_db( p=pS, DS="aggregated_data", redo=TRUE )  # used for data matching/lookup in other aegis projects that use substrate
  M = substrate_db( p=pS, DS="carstm_inputs", redo=TRUE )  # will redo if not found
  M = NULL; gc()
  fit = carstm_model( p=pS, M='substrate_db( p=p, DS="carstm_inputs")', DS="redo" )  # run model and obtain predictions

  if (0) {
    # to use a saved instance
    fit = carstm_model( p=pS, DS="carstm_modelled_fit" )  # extract currently saved model fit
    summary(fit)

    res = carstm_summary( p=pS )  # load summary

    vn = paste(pS$variabletomodel, "predicted", sep=".")
    carstm_plot( p=pS, res=res, vn=vn ) # maps of some of the results

    vn = paste(pS$variabletomodel, "predicted_se", sep=".")
    carstm_plot( p=pS, res=res, vn=vn )

    vn = paste(pS$variabletomodel, "random_auid_nonspatial", sep=".")
    carstm_plot( p=pS, res=res, vn=vn )

    vn = paste(pS$variabletomodel, "random_auid_spatial", sep=".")
    carstm_plot( p=pS, res=res, vn=vn )

    plot(fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE, single=TRUE )
  }




# -------------------------------------------------
# Part 5 -- create covariate field for temperature
# ensure the data assimilation in temperature is first completed :: 01.temperature_data.R

  pT = temperature_parameters( p=parameters_reset(p), project_class="carstm"  )
  M = temperature_db( p=pT, DS="aggregated_data", redo=TRUE )  #  used for data matching/lookup in other aegis projects that use temperature
  M = temperature_db( p=pT, DS="carstm_inputs", redo=TRUE )  # will redo if not found
  M = NULL; gc()
  fit = carstm_model( p=pT, M='temperature_db( p=p, DS="carstm_inputs")', DS="redo"  ) # run model and obtain predictions

  if (0) {
    # to use a saved instance
    fit = carstm_model( p=pT, DS="carstm_modelled_fit" )  # extract currently saved model fit
    summary(fit)

    res = carstm_summary( p=pT )

    vn = paste(pT$variabletomodel, "predicted", sep=".")
    carstm_plot( p=pT, res=res, vn=vn, time_match=list(year="2000", dyear="0.65" ) )       # maps of some of the results -- note maps of other dyears do not exists as they are not predicted upon

    vn = paste(pT$variabletomodel, "predicted_se", sep=".")
    carstm_plot( p=pT, res=res, vn=vn, time_match=list(year="2000", dyear="0.65" ) )       # maps of some of the results

    vn = paste(pT$variabletomodel, "random_auid_nonspatial", sep=".")
    carstm_plot( p=pT, res=res, vn=vn, time_match=list(year="2000"  ) )       # maps of some of the results
    vn = paste(pT$variabletomodel, "random_auid_spatial", sep=".")
    carstm_plot( p=pT, res=res, vn=vn, time_match=list(year="2000"  ) )       # maps of some of the results


    plot(fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE, single=TRUE )


    # recent=as.character((year.assessment-6): year.assessment)
    # vn = paste(pT$variabletomodel, "predicted", sep=".")

    # for (x in recent){
    #   fn=paste(x,"t",  "pdf", sep=".")
    #   outfile=paste(plot.dir, fn, sep="/")
    #   each.plot=   carstm_plot( p=pT, res=res, vn=vn, time_match=list(year=x, dyear="0.65" ) )
    #   pdf(outfile)
    #   print(each.plot)
    #   dev.off()
    # }

  }


# -------------------------------------------------
# Part 6 -- create covariate field for species composition 1
# ensure that survey data is assimilated : bio.snowcrab::01snowcb_data.R, aegis.survey::01.surveys.data.R , etc.

  pPC1 = speciescomposition_parameters( p=parameters_reset(p), project_class="carstm", variabletomodel="pca1" )
  M = speciescomposition_db( p=pPC1, DS="carstm_inputs",  redo=TRUE )  # will redo if not found
  M = NULL; gc()
  fit = carstm_model( p=pPC1, M='speciescomposition_db( p=p, DS="carstm_inputs" )', DS="redo"   ) # run model and obtain predictions

  if (0) {
      # to use a saved instance
      fit = carstm_model( p=pPC1, DS="carstm_modelled_fit" )  # extract currently saved model fit
      summary(fit)

      res = carstm_summary( p=pPC1  )

      plot( fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE )

      vn = paste(pPC1$variabletomodel, "predicted", sep=".")
      carstm_plot( p=pPC1, res=res, vn=vn, time_match=list(year="2017" ), dyear="0.65" ) # maps of some of the results
      vn = paste(pPC1$variabletomodel, "predicted_se", sep=".")
      carstm_plot( p=pPC1, res=res, vn=vn, time_match=list(year="2017" ), dyear="0.65" ) # maps of some of the results

      vn = paste(pPC1$variabletomodel, "random_auid_nonspatial", sep=".")
      carstm_plot( p=pPC1, res=res, vn=vn, time_match=list(year="2017" ), dyear="0.65" )       # maps of some of the results , dyear="0.65"
      vn = paste(pPC1$variabletomodel, "random_auid_spatial", sep=".")
      carstm_plot( p=pPC1, res=res, vn=vn, time_match=list(year="2017" ), dyear="0.65" )       # maps of some of the results , dyear="0.65"
  }



# -------------------------------------------------
# Part 7 -- create covariate field for species composition 2
# ensure that survey data is assimilated : bio.snowcrab::01snowcb_data.R, aegis.survey::01.surveys.data.R ,

  pPC2 = speciescomposition_parameters( p=parameters_reset(p), project_class="carstm", variabletomodel="pca2")
  M = speciescomposition_db( p=pPC2, DS="carstm_inputs", redo=TRUE )  # will redo if not found
  M = NULL; gc()
  fit = carstm_model( p=pPC2, M='speciescomposition_db( p=p, DS="carstm_inputs" )', DS="redo"  ) # run model and obtain predictions
  if (0) {

    # to use a saved instance
    fit = carstm_model( p=pPC2, DS="carstm_modelled_fit" )  # extract currently saved model fit
    summary(fit)

    res = carstm_summary( p=pPC2  )

    plot( fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE )

    vn = paste(pPC2$variabletomodel, "predicted", sep=".")
    carstm_plot( p=pPC2, res=res, vn=vn, time_match=list(year="2017" ), dyear="0.65" )       # maps of some of the results
    vn = paste(pPC2$variabletomodel, "predicted_se", sep=".")
    carstm_plot( p=pPC2, res=res, vn=vn, time_match=list(year="2017" ), dyear="0.65" ) # maps of some of the results

    vn = paste(pPC2$variabletomodel, "random_auid_nonspatial", sep=".")
    carstm_plot( p=pPC2, res=res, vn=vn, time_match=list(year="2017" ), dyear="0.65" )       # maps of some of the results , dyear="0.65"
    vn = paste(pPC2$variabletomodel, "random_auid_spatial", sep=".")
    carstm_plot( p=pPC2, res=res, vn=vn, time_match=list(year="2017" ), dyear="0.65" )       # maps of some of the results , dyear="0.65"

  }

# finished covariates ... move onto abundance index estimation


# -------------------------------------------------
# Part 8 -- Snow crab anbundance -- main mode used for production


if (0) {
  # use alternate model -- zero-inflated1

          p$carstm_model_label = "zeroinflated"  #unique to this project ... to permit alt model forms/variations within the same overall carstm

          p$carstm_model_call = paste(
          'inla( formula =', p$variabletomodel,
          ' ~ 1
            + offset( log(data_offset))
            + f( dyri, model="ar1", hyper=H$ar1 )
            + f( inla.group( t, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2)
            + f( inla.group( z, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2)
            + f( inla.group( substrate.grainsize, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2)
            + f( inla.group( pca1, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2)
            + f( inla.group( pca2, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2)
            + f( auid, model="bym2", graph=slot(sppoly, "nb"), group=year_factor, scale.model=TRUE, constr=TRUE, hyper=H$bym2, control.group=list(model="ar1", hyper=H$ar1_group)),
            family = "zeroinflatedpoisson1",
            data= M,
            control.compute=list(cpo=TRUE, waic=TRUE, dic=TRUE, config=TRUE),
            control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
            control.predictor=list(compute=FALSE, link=1 ),
            #control.fixed = list(prec.intercept = 0.1),
            control.fixed = H$fixed,  # priors for fixed effects, generic is ok
            control.inla = list(cmin = 0, h=1e-4, tolerance=1e-9, strategy="adaptive", optimise.strategy="smart"),
            # control.inla = list( h=1e-6, tolerance=1e-12), # increase in case values are too close to zero
            #control.inla=list( strategy="laplace", cutoff=1e-6, correct=TRUE, correct.verbose=FALSE ),
            # control.inla = list( h=1e-6, tolerance=1e-12), # increase in case values are too close to zero
            # control.inla = list(h=1e-3, tolerance=1e-9, cmin=0), # restart=3), # restart a few times in case posteriors are poorly defined
            # control.mode = list( restart=TRUE, result=RES ), # restart from previous estimates
          verbose=TRUE
          )'
        )


}


  M = snowcrab_carstm( p=p, DS="carstm_inputs", redo=TRUE )  # will redo if not found
  fit = carstm_model( p=p, M='snowcrab_carstm( p=p, DS="carstm_inputs" )' ) # 151 configs and long optim .. 19 hrs
  RES = snowcrab_carstm(p=p, DS="carstm_output_compute" )

  if (0) {

      fit =  carstm_model( p=p, DS="carstm_modelled_fit" )  # extract currently saved model fit
      summary(fit)

      res = carstm_summary( p=p )

      vn = paste(p$variabletomodel, "predicted", sep=".")
      carstm_plot( p=p, res=res, vn=vn, time_match=list(year="2019" ) )     # maps of some of the results

      plot(fit)
      plot(fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE, single=TRUE )
      s = summary(fit)
      s$dic$dic
      s$dic$p.eff

      # maps of some of the results
      vn = paste(p$variabletomodel, "random_sample_iid", sep=".")
      if (exists(vn, res)) carstm_plot( p=p, res=res, vn=vn, time_match=list(year="2019", dyear="0.65") )

      vn = paste(p$variabletomodel, "random_auid_nonspatial", sep=".")
      if (exists(vn, res)) {
        res_dim = dim( res[[vn]] )
        if (res_dim == 1 ) time_match = NULL
        if (res_dim == 2 ) time_match = list(year="2000")
        if (res_dim == 3 ) time_match = list(year="2000", dyear="0.65" )
        carstm_plot( p=p, res=res, vn=vn, time_match=time_match )
      }

      vn = paste(p$variabletomodel, "random_auid_spatial", sep=".")
      if (exists(vn, res)) {
        res_dim = dim( res[[vn]] )
        if (res_dim == 1 ) time_match = NULL
        if (res_dim == 2 ) time_match = list(year="2000")
        if (res_dim == 3 ) time_match = list(year="2000", dyear="0.65" )
        carstm_plot( p=p, res=res, vn=vn, time_match=time_match )
      }

      RES = snowcrab_carstm(p=p, DS="carstm_output_timeseries" )
      bio = snowcrab_carstm(p=p, DS="carstm_output_spacetime_biomass" )
      num = snowcrab_carstm(p=p, DS="carstm_output_spacetime_number" )


      plot( cfaall ~ yrs, data=RES, lty="solid", lwd=4, pch=20, col="slateblue", type="b", ylab="Biomass (kt)", xlab="")
      lines( cfaall_lb ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )
      lines( cfaall_ub ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )

      plot( cfasouth ~ yrs, data=RES, lty="solid", lwd=4, pch=20, col="slateblue", type="b", ylab="Biomass (kt)", xlab="")
      lines( cfasouth_lb ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )
      lines( cfasouth_ub ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )

      plot( cfanorth ~ yrs, data=RES, lty="solid", lwd=4, pch=20, col="slateblue", type="b", ylab="Biomass (kt)", xlab="")
      lines( cfanorth_lb ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )
      lines( cfanorth_ub ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )

      plot( cfa4x ~ yrs, data=RES, lty="solid", lwd=4, pch=20, col="slateblue", type="b", ylab="Biomass (kt)", xlab="")
      lines( cfa4x_lb ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )
      lines( cfa4x_ub ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )


      p$boundingbox = list( xlim=p$corners$lon, ylim=p$corners$lat) # bounding box for plots using spplot

      p$coastLayout = aegis.coastline::coastline_layout(p=p)
      # p$mypalette=RColorBrewer::brewer.pal(9, "YlOrRd")

      sppoly = areal_units( p=p )  # to reload

      # map it ..mean density
      vn = "pred"
      slot(sppoly, "data")[,vn] = bio[,"2019"]
      brks = interval_break(X= sppoly[[vn]], n=length(p$mypalette), style="quantile")
      spplot( sppoly, vn, col.regions=p$mypalette, main=vn, at=brks, sp.layout=p$coastLayout, col="transparent" )


      plot( fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE )
      plot( fit$marginals.hyperpar$"Phi for auid", type="l")  # posterior distribution of phi nonspatial dominates
      plot( fit$marginals.hyperpar$"Precision for auid", type="l")
      plot( fit$marginals.hyperpar$"Precision for setno", type="l")
  }





# ----

# estimate numbers for "recruitment"

  p$selection$type = "number"

  p$selection$biologicals = list(
    spec_bio=bio.taxonomy::taxonomy.recode( from="spec", to="parsimonious", tolookup=p$groundfish_species_code ),
    sex=0, # male
    # mat=1, # do not use maturity status in groundfish data as it is suspect ..
    len= c( 80, 94.99 )/10, #  mm -> cm ; aegis_db in cm
    ranged_data="len"
  )


  M = snowcrab_carstm( p=p, DS="carstm_inputs", redo=TRUE )  # will redo if not found
  fit = carstm_model( p=p, M='snowcrab_carstm( p=p, DS="carstm_inputs" )' ) # 151 configs and long optim .. 19 hrs
  fit =  carstm_model( p=p, DS="carstm_modelled_fit" )  # extract currently saved model fit
  summary(fit)
  res = carstm_summary( p=p )
  RES = snowcrab_carstm(p=p, DS="carstm_output_compute" )



# for the abundance run :


