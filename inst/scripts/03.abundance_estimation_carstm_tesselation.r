

# Snow crab --- Areal unit modelling of habitat  -- no reliance upon stmv fields


# -------------------------------------------------
# Part 1 -- construct basic parameter list defining the main characteristics of the study
# require(aegis)

  year.assessment = 2020

  p = bio.snowcrab::snowcrab_parameters( 
    project_class="carstm", 
    assessment.years=2000:year.assessment, 
    areal_units_type="tesselation",
    carstm_model_label = "tesselation",   # default is the name of areal_units_type  
    selection = list(type = "number")
 )
  

# ------------------------------------------------
# Part 2 -- spatiotemporal statistical model 

  if ( spatiotemporal_model ) {

      if (0) {
      # polygon structure:: create if not yet made
      # adjust based upon RAM requirements and ncores
        # create if not yet made
        for (au in c("cfanorth", "cfasouth", "cfa4x", "cfaall" )) plot(polygon_managementareas( species="snowcrab", au))
        xydata = snowcrab.db( p=p, DS="areal_units_input", redo=TRUE )

        sppoly = areal_units( p=p, hull_alpha=16, redo=TRUE, verbose=TRUE )  # create constrained polygons with neighbourhood as an attribute
        plot( sppoly[, "npts"]  )
        
        MS = NULL

      }

      sppoly = areal_units( p=p )  # to reload

    
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
      isobaths=aegis.bathymetry::isobath_db( depths=c(50, 100, 200, 400, 800), project_to=plot_crs )
      managementlines = aegis.polygons::area_lines.db( DS="cfa.regions", returntype="sf", project_to=plot_crs )
    
      time_match = list( year=as.character(2020)  )
      carstm_map(  res=res, 
          vn=paste(p$variabletomodel, "predicted", sep="."), 
          time_match=time_match, 
          coastline=coastline,
          managementlines=managementlines,
          isobaths=isobaths,
          main=paste("Predicted numerical abundance", paste0(time_match, collapse="-") )  
      )
        

      # map all :
      vn = paste(p$variabletomodel, "predicted", sep=".")

      outputdir = file.path( p$modeldir, p$carstm_model_label, "predicted.numerical.densitites" )

      if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )

      brks = pretty(  quantile(res[[vn]], probs=c(0,0.975))  )

      for (y in res$year ){

          time_match = list( year=as.character(y)  )
          fn_root = paste("Predicted_numerical_abundance", paste0(time_match, collapse="-"), sep="_")
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


      ( fn = file.path( outputdir, "cfa_all.png") )
      png( filename=fn, width=3072, height=2304, pointsize=12, res=300 )
        plot( cfaall ~ yrs, data=RES, lty="solid", lwd=4, pch=20, col="slateblue", type="b", ylab="Biomass index (kt)", xlab="")
        lines( cfaall_lb ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )
        lines( cfaall_ub ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )
      dev.off()

      
      ( fn = file.path( outputdir, "cfa_south.png") )
      png( filename=fn, width=3072, height=2304, pointsize=12, res=300 )
        plot( cfasouth ~ yrs, data=RES, lty="solid", lwd=4, pch=20, col="slateblue", type="b", ylab="Biomass index (kt)", xlab="")
        lines( cfasouth_lb ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )
        lines( cfasouth_ub ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )
      dev.off()

      ( fn = file.path( outputdir, "cfa_north.png") )
      png( filename=fn, width=3072, height=2304, pointsize=12, res=300 )
        plot( cfanorth ~ yrs, data=RES, lty="solid", lwd=4, pch=20, col="slateblue", type="b", ylab="Biomass index (kt)", xlab="")
        lines( cfanorth_lb ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )
        lines( cfanorth_ub ~ yrs, data=RES, lty="dotted", lwd=2, col="slategray" )
      dev.off()

      ( fn = file.path( outputdir, "cfa_4x.png") ) 
      png( filename=fn, width=3072, height=2304, pointsize=12, res=300 )
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

Number of fn-calls= 784 with 8.496 sec/fn-call
Parallel efficiency for 'Preparations' and 'Approx inference':
    Accumulated CPU-time is equivalent to 5.93 threads running at 100%
    Efficiency using 12 threads = 49.42%

Saving carstm fit (this can be slow): /home/jae/bio.data/bio.snowcrab/modelled/tesselation/bio.snowcrab|snowcrab|tesselation|1|snowcrab|10|3|snowcrab_managementareas.rdata|carstm_modelled_fit|totno|inla|1000|30.rdata
Computing summaries ...
carstm summary saved as: /home/jae/bio.data/bio.snowcrab/modelled/tesselation/bio.snowcrab|snowcrab|tesselation|1|snowcrab|10|3|snowcrab_managementareas.rdata|carstm_modelled_summary|totno|inla|1000|30.rdata

Call:
   c("inla(formula = p$carstm_model_formula, family = p$carstm_model_family, ", " data = M, verbose = TRUE, 
   control.compute = list(dic = TRUE, ", " waic = TRUE, cpo = FALSE, config = TRUE), control.predictor = 
   list(compute = TRUE, ", " link = 1), control.family = p$options.control.family, ", " control.inla = 
   p$options.control.inla[[civ]], control.results = list(return.marginals.random = TRUE, ", " 
   return.marginals.predictor = TRUE))") 
Time used:
    Pre = 4.36, Running = 6664, Post = 5.65, Total = 6674 
Fixed effects:
             mean    sd 0.025quant 0.5quant 0.975quant  mode kld
(Intercept) 7.035 0.106      6.827    7.035      7.244 7.035   0

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
Precision for dyri                                                         53.187 0.109     53.038   53.167     53.440
Rho for dyri                                                                0.754 0.000      0.753    0.754      0.755
Precision for inla.group(t, method = "quantile", n = 11)                   58.617 0.184     58.300   58.613     58.926
Precision for inla.group(z, method = "quantile", n = 11)                   55.026 0.389     54.505   54.952     55.948
Precision for inla.group(substrate.grainsize, method = "quantile", n = 11) 58.669 0.213     58.307   58.665     59.024
Precision for inla.group(pca1, method = "quantile", n = 11)                55.032 0.393     54.505   54.957     55.964
Precision for inla.group(pca2, method = "quantile", n = 11)                56.259 0.251     55.894   56.224     56.772
Precision for auid                                                         58.480 0.002     58.475   58.479     58.985
Phi for auid                                                                0.048 0.000      0.048    0.048      0.049
GroupRho for auid                                                           0.785 0.000      0.785    0.785      0.785
                                                                             mode
Precision for dyri                                                         53.066
Rho for dyri                                                                0.753
Precision for inla.group(t, method = "quantile", n = 11)                   58.417
Precision for inla.group(z, method = "quantile", n = 11)                   54.624
Precision for inla.group(substrate.grainsize, method = "quantile", n = 11) 58.418
Precision for inla.group(pca1, method = "quantile", n = 11)                54.624
Precision for inla.group(pca2, method = "quantile", n = 11)                55.943
Precision for auid                                                         58.480
Phi for auid                                                                0.048
GroupRho for auid                                                           0.785

Expected number of effective parameters(stdev): 608.49(1.62)
Number of equivalent replicates : 11.88 

Deviance Information Criterion (DIC) ...............: 56306.57
Deviance Information Criterion (DIC, saturated) ....: 38406.06
Effective number of parameters .....................: -1193.28

Watanabe-Akaike information criterion (WAIC) ...: 65141.77
Effective number of parameters .................: 5098.54

Marginal log-Likelihood:  -31355.14 
Posterior marginals for the linear predictor and
 the fitted values are computed




  }


  ########## 


  if (fishery_model) {

      p$fishery_model = fishery_model( DS = "logistic_parameters", p=p, tag=p$areal_units_type )
      
      p$fishery_model$stancode_compiled = rstan::stan_model( model_code=p$fishery_model$stancode )

      if (0) {
        # testing:
        # later:::ensureInitialized()  # solve mode error
        str( p$fishery_model$standata )
        p$fishery_model$standata$Kmu = c(4, 40, 1)
        p$fishery_model$standata$rmu = c(1, 1, 1)
        p$fishery_model$standata$qmu = c(1, 1, 1)
        p$fishery_model$standata$Ksd = c(0.5, 0.5, 0.5) * p$fishery_model$standata$Kmu  # c( 2, 20, 0.5)
        p$fishery_model$standata$rsd = c(0.5, 0.5, 0.5) * p$fishery_model$standata$rmu  # rep( 0.3, 3)
        p$fishery_model$standata$qsd = c(0.5, 0.5, 0.5) * p$fishery_model$standata$qmu  # rep( 0.3, 3)

        # to recompile
      }


      str( p$fishery_model)

      res = fishery_model( p=p, DS="logistic_model", tag=p$areal_units_type,
        chains=3, cores=3, iter=20000, warmup=12000, refresh=1000,
        control=list(adapt_delta=0.99, max_treedepth=18)
      )
      # res = fishery_model( p=p, DS="logistic_samples", tag=p$areal_units_type )  # to get samples



      p$fishery_model$outdir = file.path( p$modeldir, p$carstm_model_label, "fishery_model_results" )


      # frequency density of key parameters
      figure.mcmc( "K", res=res, fn=file.path(p$fishery_model$outdir, "K.density.png" ) )
      figure.mcmc( "r", res=res, fn=file.path(p$fishery_model$outdir, "r.density.png" ) )
      figure.mcmc( "q", res=res, fn=file.path(p$fishery_model$outdir, "q.density.png" ) ,xrange=c(0,5.5))
      figure.mcmc( "FMSY", res=res, fn=file.path(p$fishery_model$outdir, "FMSY.density.png" ) )
      # figure.mcmc( "bosd", res=res, fn=file.path(p$fishery_model$outdir, "bosd.density.png" ) )
      # figure.mcmc( "bpsd", res=res, fn=file.path(p$fishery_model$outdir, "bpsd.density.png" ) )

      # timeseries
      figure.mcmc( type="timeseries", vname="biomass", res=res, fn=file.path(p$fishery_model$outdir, "biomass.timeseries.png" ), save.plot=T )
      figure.mcmc( type="timeseries", vname="fishingmortality", res=res, fn=file.path(p$fishery_model$outdir, "fishingmortality.timeseries.png" ) )

      # Summary table of mean values for inclusion in document
      biomass.summary.table(x)

      # Harvest control rules
      figure.mcmc( type="hcr", vname="default", res=res, fn=file.path(p$fishery_model$outdir, "hcr.default.png" ), save.plot=T  )
      figure.mcmc( type="hcr", vname="simple", res=res, fn=file.path(p$fishery_model$outdir, "hcr.simple.png" ) )

      # diagnostics
      # figure.mcmc( type="diagnostic.errors", res=res, fn=file.path(p$fishery_model$outdir, "diagnostic.errors.png" ) )
      # figure.mcmc( type="diagnostic.phase", res=res, fn=file.path(p$fishery_model$outdir, "diagnostic.phase.png" ) )


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
