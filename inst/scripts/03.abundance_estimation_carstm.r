

# Snow crab --- Areal unit modelling of habitat  -- no reliance upon stmv fields


# -------------------------------------------------
# Part 1 -- construct basic parameter list defining the main characteristics of the study
# require(aegis)

  year.assessment = 2020

  p = bio.snowcrab::snowcrab_parameters( 
    project_class="carstm", 
    assessment.years=2000:year.assessment, 
    areal_units_type="lattice",
    selection = list(type = "number")
 )

# ------------------------------------------------
# Part 2 -- spatiotemporal statistical model 

  if ( spataiotemporal_model ) {

    if (0) {
      # polygon structure:: create if not yet made
      # adjust based upon RAM requirements and ncores
      
      p$areal_units_constraint_nmin = 3
      p$areal_units_constraint_ntarget = 10

      for (au in c("cfanorth", "cfasouth", "cfa4x", "cfaall" )) plot(polygon_managementareas( species="snowcrab", au))
      xydata = snowcrab.db( p=p, DS="areal_units_input", redo=TRUE )

      sppoly = areal_units( p=p, redo=TRUE )  # create constrained polygons with neighbourhood as an attribute
      MS = NULL
    }

    sppoly = areal_units( p=p )  # to reload
    plot( sppoly[, "au_sa_km2"]  )


    # -------------------------------------------------
    M = snowcrab.db( p=p, DS="carstm_inputs", redo=TRUE )  # will redo if not found
    M = NULL; gc()

    fit = carstm_model( p=p, M='snowcrab.db( p=p, DS="carstm_inputs" )' ) # 151 configs and long optim .. 19 hrs
    # fit = carstm_model( p=p, DS="carstm_modelled_fit")

      # extract results
      if (0) {
        # very large files .. slow 
        fit = carstm_model( p=p, DS="carstm_modelled_fit" )  # extract currently saved model fit
        plot(fit)
        plot(fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE )
      }


    res = carstm_model( p=p, DS="carstm_modelled_summary"  ) # to load currently saved results
    res$summary$dic$dic
    res$summary$dic$p.eff
    res$dyear


    plot_crs = p$aegis_proj4string_planar_km
    coastline=aegis.coastline::coastline_db( DS="eastcoast_gadm", project_to=plot_crs )
    isobaths=aegis.bathymetry::isobath_db( depths=c(50, 100, 200, 400 ), project_to=plot_crs )
    managementlines = aegis.polygons::area_lines.db( DS="cfa.regions", returntype="sf", project_to=plot_crs )
  
    time_match = list( year=as.character(2020)  )
    carstm_map(  res=res, 
        vn=paste(p$variabletomodel, "predicted", sep="."), 
        time_match=time_match, 
        coastline=coastline,
        managementlines=managementlines,
        isobaths=isobaths,
        main=paste("Predicted abundance", paste0(time_match, collapse="-") )  
    )
      

    # map all :
    vn = paste(p$variabletomodel, "predicted", sep=".")

    outputdir = file.path( p$modeldir, p$carstm_model_label, "predicted.numerical.densitites" )
    if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )

    brks = pretty(  quantile(res[[vn]], probs=c(0,0.975))  )

    for (y in res$year ){

        time_match = list( year=as.character(y)  )
        fn_root = paste("Predicted_abundance", paste0(time_match, collapse="-"), sep="_")
        fn = file.path( outputdir, paste(fn_root, "png", sep=".") )

          carstm_map(  res=res, vn=vn, time_match=time_match, 
            breaks =brks,
            coastline=coastline,
            isobaths=isobaths,
            managementlines=managementlines,
            main=paste("Predicted numerial abundance", paste0(time_match, collapse="-") ),
            outfilename=fn
          )  

    }
    


    snowcrab.db(p=p, DS="carstm_output_compute" )
    
    RES = snowcrab.db(p=p, DS="carstm_output_timeseries" )

    bio = snowcrab.db(p=p, DS="carstm_output_spacetime_biomass" )
    num = snowcrab.db(p=p, DS="carstm_output_spacetime_number" )

    outputdir = file.path( p$modeldir, p$carstm_model_label, "aggregated_biomass_timeseries" )

    if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )


    png( filename=file.path( outputdir, "cfa_all.png"), width=3072, height=2304, pointsize=12, res=300 )
      plot( cfaall ~ yrs, data=RES, lty="solid", lwd=4, pch=20, col="slateblue", type="b", ylab="Biomass index (kt)", xlab="")
      lines( cfaall_lb ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )
      lines( cfaall_ub ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )
    dev.off()


    png( filename=file.path( outputdir, "cfa_south.png"), width=3072, height=2304, pointsize=12, res=300 )
      plot( cfasouth ~ yrs, data=RES, lty="solid", lwd=4, pch=20, col="slateblue", type="b", ylab="Biomass index (kt)", xlab="")
      lines( cfasouth_lb ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )
      lines( cfasouth_ub ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )
    dev.off()


    png( filename=file.path( outputdir, "cfa_north.png"), width=3072, height=2304, pointsize=12, res=300 )
      plot( cfanorth ~ yrs, data=RES, lty="solid", lwd=4, pch=20, col="slateblue", type="b", ylab="Biomass index (kt)", xlab="")
      lines( cfanorth_lb ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )
      lines( cfanorth_ub ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )
    dev.off()


    png( filename=file.path( outputdir, "cfa_4x.png"), width=3072, height=2304, pointsize=12, res=300 )
      plot( cfa4x ~ yrs, data=RES, lty="solid", lwd=4, pch=20, col="slateblue", type="b", ylab="Biomass index (kt)", xlab="")
      lines( cfa4x_lb ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )
      lines( cfa4x_ub ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )
    dev.off()





    # map it ..mean density

    sppoly = areal_units( p=p )  # to reload

    plot_crs = p$aegis_proj4string_planar_km
    coastline=aegis.coastline::coastline_db( DS="eastcoast_gadm", project_to=plot_crs )
    isobaths=aegis.bathymetry::isobath_db( depths=c(50, 100, 200, 400, 800), project_to=plot_crs  )
    managementlines = aegis.polygons::area_lines.db( DS="cfa.regions", returntype="sf", project_to=plot_crs )
  
    vn = paste("biomass", "predicted", sep=".")

    outputdir = file.path( p$modeldir, p$carstm_model_label, "predicted.biomass.densitites" )

    if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )

    brks = pretty(  quantile( bio[], probs=c(0,0.975) )* 10^6 )

    for (i in 1:length(p$yrs) ){
      y = as.character( p$yrs[i] )
      sppoly[,vn] = bio[,y] * 10^6
      fn = file.path( outputdir , paste( "biomass", y, "png", sep=".") )

        carstm_map(  sppoly=sppoly, vn=vn,    
          breaks=brks, 
          coastline=coastline,
          isobaths=isobaths,
          managementlines=managementlines,
          main=paste("Predicted biomass density", y ),  
          outfilename=fn
        )
    }

    plot( fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE )
    plot( fit$marginals.hyperpar$"Phi for auid", type="l")  # posterior distribution of phi nonspatial dominates
    plot( fit$marginals.hyperpar$"Precision for auid", type="l")
    plot( fit$marginals.hyperpar$"Precision for setno", type="l")



    (res$summary)


    Fixed effects:
                mean    sd 0.025quant 0.5quant 0.975quant  mode kld
    (Intercept) 6.602 0.102      6.403    6.602      6.801 6.602   0

    Random effects:
      Name	  Model
        dyri AR1 model
      inla.group(t, method = "quantile", n = 11) RW2 model
      inla.group(z, method = "quantile", n = 11) RW2 model
      inla.group(substrate.grainsize, method = "quantile", n = 11) RW2 model
      inla.group(pca1, method = "quantile", n = 11) RW2 model
      inla.group(pca2, method = "quantile", n = 11) RW2 model
      auid BYM2 model

    Model hyperparameters:
                                                                                mean    sd 0.025quant 0.5quant 0.975quant
    Precision for dyri                                                         60.314 0.379     59.581   60.320     60.927
    Rho for dyri                                                                0.763 0.001      0.762    0.763      0.764
    Precision for inla.group(t, method = "quantile", n = 11)                   54.500 0.251     54.118   54.477     54.967
    Precision for inla.group(z, method = "quantile", n = 11)                   54.221 0.198     53.941   54.189     54.651
    Precision for inla.group(substrate.grainsize, method = "quantile", n = 11) 59.579 0.237     59.204   59.565     59.998
    Precision for inla.group(pca1, method = "quantile", n = 11)                55.907 0.503     55.023   55.885     56.787
    Precision for inla.group(pca2, method = "quantile", n = 11)                50.613 0.264     50.211   50.587     51.108
    Precision for auid                                                         32.051 0.028     31.992   32.052     32.102
    Phi for auid                                                                0.060 0.000      0.060    0.060      0.060
    GroupRho for auid                                                           0.850 0.000      0.850    0.850      0.852
                                                                                mode
    Precision for dyri                                                         60.114
    Rho for dyri                                                                0.762
    Precision for inla.group(t, method = "quantile", n = 11)                   54.147
    Precision for inla.group(z, method = "quantile", n = 11)                   53.983
    Precision for inla.group(substrate.grainsize, method = "quantile", n = 11) 59.252
    Precision for inla.group(pca1, method = "quantile", n = 11)                55.481
    Precision for inla.group(pca2, method = "quantile", n = 11)                50.261
    Precision for auid                                                         32.058
    Phi for auid                                                                0.060
    GroupRho for auid                                                           0.850

    Expected number of effective parameters(stdev): 440.55(1.45)
    Number of equivalent replicates : 15.80 

    Deviance Information Criterion (DIC) ...............: 43123.14
    Deviance Information Criterion (DIC, saturated) ....: 27013.84
    Effective number of parameters .....................: 374.76

    Watanabe-Akaike information criterion (WAIC) ...: 44938.34
    Effective number of parameters .................: 1957.84

    Marginal log-Likelihood:  -21151.77 
    Posterior marginals for the linear predictor and
    the fitted values are computed




  }


  ########## 


  if (fishery_model) {

      p$fishery_model = fishery_model( DS = "logistic_parameters", p=p, tag=p$areal_units_type )
      p$fishery_model$stancode = stan_initialize( stan_code=fishery_model( p=p, DS="stan_surplus_production" ) )

      str( p$fishery_model)

      res = fishery_model( 
        DS="logistic_model", 
        p=p, 
        tag=p$areal_units_type,
        # from here down are params for cmdstanr::sample()
        data=p$fishery_model$standata, 
        iter_warmup = 12000,
        iter_sampling = 5000,
        seed = 123,
        chains = 3,
        parallel_chains = 3,  # The maximum number of MCMC chains to run in parallel.
        max_treedepth = 18,
        adapt_delta = 0.99,
        refresh = 1000
      )
 
      # res = fishery_model( p=p, DS="samples", tag=p$areal_units_type )  # to get samples
  
      if (0) {
        fit = fishery_model( p=p, DS="fit", tag=p$areal_units_type )  # to get samples
      
        print( fit, max_rows=30 )
        # fit$summary("K", "r", "q")
        
        fit$cmdstan_diagnose()
        fit$cmdstan_summary()
      }


      color_scheme_set("gray")
      mcmc_dens(res$mcmc, regex_pars="K",  facet_args = list(nrow = 3, labeller = ggplot2::label_parsed ) ) + facet_text(size = 14 )   
      # mcmc_hist( fit$draws("K"))


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
  }


# end
